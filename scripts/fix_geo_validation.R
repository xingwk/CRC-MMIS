#!/usr/bin/env Rscript
# External validation fix: GSE39582 (and diagnose GSE33113)
library(survival); library(dplyr); library(readr)

BASE <- "C:/Users/xingw/WorkBuddy/2026-06-04-22-17-17"
GEO <- file.path(BASE, "data", "geo")
RESULTS <- file.path(BASE, "results", "revisions_r2")
dir.create(RESULTS, showWarnings=FALSE, recursive=TRUE)

# --- Load TCGA-trained coefficients ---
selected <- read_csv(file.path(BASE, "results/v31_run01_seed42/selected_genes.csv"),
                     show_col_types=FALSE)
coefs_run1 <- setNames(selected$coefficient, selected$gene)
selected2 <- read_csv(file.path(BASE, "results/v31_run02_seed123/selected_genes.csv"),
                      show_col_types=FALSE)
coefs_run2 <- setNames(selected2$coefficient, selected2$gene)

# --- Load GEO gene-symbol matrices ---
load_geo_matrix <- function(cohort) {
  f <- file.path(GEO, paste0(cohort, "_series_matrix.txt.gz"))
  mat <- read.table(gzfile(f), header=TRUE, sep="\t", row.names=1, 
                    check.names=FALSE, comment.char="!", fill=TRUE)
  cat(sprintf("  %s: %d probes x %d samples\n", cohort, nrow(mat), ncol(mat)))
  mat
}

# --- Load probe→symbol mapping ---
probe_map <- read_csv(file.path(GEO, "GPL570_probe_to_symbol.csv"), show_col_types=FALSE)
cat(sprintf("Loaded probe map: %d probes\n", nrow(probe_map)))

# --- Convert probe matrix to gene matrix ---
probe_to_gene <- function(pmat, pmap) {
  symbols <- pmap$gene_symbol[match(rownames(pmat), pmap$probe_id)]
  valid <- !is.na(symbols) & symbols != ""
  pmat <- pmat[valid, ]; gene_sym <- symbols[valid]
  # Collapse by mean
  gene_mat <- aggregate(pmat, by=list(gene_sym), FUN=mean, na.rm=TRUE)
  rownames(gene_mat) <- gene_mat$Group.1
  gene_mat$Group.1 <- NULL
  as.matrix(gene_mat)
}

# --- Process a single cohort with a given coefficient set ---
process_cohort <- function(cohort, coefs, run_label) {
  cat(sprintf("\n=== %s (Run: %s) ===\n", cohort, run_label))
  
  # Load survival
  surv <- read_csv(file.path(GEO, paste0(cohort, "_survival.csv")), show_col_types=FALSE)
  cat(sprintf("  Survival data: %d rows, cols: %s\n", nrow(surv),
              paste(names(surv), collapse=", ")))
  
  # Detect time/event columns
  time_col <- names(surv)[grepl("os_time|rfs_time|time_to|time\\b", names(surv), ignore.case=TRUE)][1]
  event_col <- names(surv)[grepl("os_event|rfs_event|event|status", names(surv), ignore.case=TRUE)][1]
  
  if (is.na(time_col)) {
    cat("  ⚠️ No time column found!\n")
    return(NULL)
  }
  cat(sprintf("  Time col: %s, Event col: %s\n", time_col, event_col))
  
  # Count events  
  if (!is.na(event_col)) {
    n_events <- sum(surv[[event_col]] == 1, na.rm=TRUE)
    n_total <- sum(!is.na(surv[[time_col]]) & surv[[time_col]] > 0)
  } else {
    # Check if time can be used as recurrence endpoint (event = any non-NA time)
    n_events <- sum(!is.na(surv[[time_col]]) & surv[[time_col]] < 2000, na.rm=TRUE)
    n_total <- sum(!is.na(surv[[time_col]]), na.rm=TRUE)
  }
  cat(sprintf("  N=%d, events=%d\n", n_total, n_events))
  
  if (n_events < 5) {
    cat(sprintf("  ⚠️ Insufficient events (%d < 5) — skipping\n", n_events))
    return(list(cohort=cohort, run=run_label, n=n_total, genes=length(coefs),
                HR=NA, p=NA, CI_low=NA, CI_high=NA, events=n_events, status="insufficient_events"))
  }
  
  # Load and convert expression
  pmat <- load_geo_matrix(cohort)
  gmat <- probe_to_gene(pmat, probe_map)
  cat(sprintf("  Gene matrix: %d genes x %d samples\n", nrow(gmat), ncol(gmat)))
  
  # Match genes
  common_genes <- intersect(names(coefs), rownames(gmat))
  cat(sprintf("  Matching genes: %d/%d\n", length(common_genes), length(coefs)))
  if (length(common_genes) < 4) {
    cat("  ⚠️ Too few matching genes\n")
    return(NULL)
  }
  
  # Compute risk score
  risk <- rep(0, ncol(gmat))
  for (g in common_genes) {
    risk <- risk + coefs[g] * gmat[g, ]
  }
  risk <- scale(risk)[,1]
  
  # Match samples
  sample_col <- names(surv)[1]  # first column is sample ID
  surv_samples <- surv[[sample_col]]
  expr_samples <- colnames(gmat)
  
  # Find intersect
  common_samples <- intersect(surv_samples, expr_samples)
  cat(sprintf("  Sample overlap: %d\n", length(common_samples)))
  
  if (length(common_samples) < 20) {
    cat("  ⚠️ Too few overlapping samples\n")
    return(NULL)
  }
  
  # Build analysis DF
  surv_idx <- match(common_samples, surv_samples)
  expr_idx <- match(common_samples, expr_samples)
  
  time_vals <- as.numeric(surv[[time_col]][surv_idx])
  event_vals <- if (is.na(event_col)) rep(1, length(common_samples)) else as.numeric(surv[[event_col]][surv_idx])
  risk_vals <- risk[expr_idx]
  
  # Remove NA
  valid <- !is.na(time_vals) & time_vals > 0 & !is.na(event_vals) & !is.na(risk_vals)
  cat(sprintf("  Valid samples: %d, events: %.0f\n", sum(valid), sum(event_vals[valid])))
  
  if (sum(valid) < 20 || sum(event_vals[valid]) < 5) {
    cat("  ⚠️ Not enough valid samples/events\n")
    return(NULL)
  }
  
  # Cox
  fit <- coxph(Surv(time_vals[valid], event_vals[valid]) ~ risk_vals[valid])
  s <- summary(fit)
  hr <- s$conf.int[1, 1]
  hr_ci <- s$conf.int[1, 3:4]
  p <- s$waldtest["pvalue"]
  c_idx <- s$concordance["C"]
  
  cat(sprintf("  HR=%.2f [%.2f-%.2f], p=%.4f, C=%.3f\n", hr, hr_ci[1], hr_ci[2], p, c_idx))
  
  list(cohort=cohort, run=run_label, n=sum(valid), genes=length(common_genes),
       HR=hr, p=p, CI_low=hr_ci[1], CI_high=hr_ci[2], C_index=c_idx,
       events=sum(event_vals[valid]), status="ok")
}

# --- Run all cohorts ---
all_results <- list()

# Run 1 (seed=42) — 10-gene signature
for (cohort in c("GSE17536", "GSE14333", "GSE38832", "GSE33113", "GSE39582")) {
  r <- process_cohort(cohort, coefs_run1, "Run1_seed42")
  if (!is.null(r)) all_results[[length(all_results) + 1]] <- r
}

# Run 2 (seed=123) — 14-gene signature (different gene count from different analysis)
for (cohort in c("GSE17536", "GSE14333", "GSE38832", "GSE33113", "GSE39582")) {
  r <- process_cohort(cohort, coefs_run2, "Run2_seed123")
  if (!is.null(r)) all_results[[length(all_results) + 1]] <- r
}

# --- Save results ---
results_df <- bind_rows(lapply(all_results, as_tibble))
write_csv(results_df, file.path(RESULTS, "geo_validation_fixed.csv"))
cat(sprintf("\n=== RESULTS SAVED to %s ===\n", file.path(RESULTS, "geo_validation_fixed.csv")))

# Print summary table
cat("\n=== SUMMARY ===\n")
for (i in 1:nrow(results_df)) {
  r <- results_df[i,]
  cat(sprintf("%-10s %-12s N=%-4d Ev=%d HR=%.2f [%.2f-%.2f] p=%.4f %s\n",
              r$cohort, r$run, r$n, r$events, r$HR, r$CI_low, r$CI_high, r$p, r$status))
}

# --- Recompute meta-analysis ---
library(meta)
cat("\n=== META-ANALYSIS (Random Effects) ===\n")
valid <- results_df %>% filter(!is.na(HR), !is.na(CI_low), status == "ok")
n_valid <- nrow(valid)
cat(sprintf("Valid studies: %d\n", n_valid))
if (n_valid >= 2) {
  meta_res <- metagen(TE = log(valid$HR), seTE = (log(valid$CI_high) - log(valid$CI_low))/3.92,
                       studlab = paste(valid$cohort, valid$run, sep="_"),
                       sm = "HR", common = FALSE, random = TRUE)
  cat(sprintf("RE pooled HR: %.3f [%.3f-%.3f], I2=%.1f%%\n",
              exp(meta_res$TE.random), exp(meta_res$lower.random), exp(meta_res$upper.random),
              meta_res$I2))
}
