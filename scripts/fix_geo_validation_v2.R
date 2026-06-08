#!/usr/bin/env Rscript
# External validation fix v2 — correct OS endpoint selection
library(survival); library(dplyr); library(readr)

BASE <- "C:/Users/xingw/WorkBuddy/2026-06-04-22-17-17"
GEO <- file.path(BASE, "data", "geo")
RESULTS <- file.path(BASE, "results", "revisions_r2")
dir.create(RESULTS, showWarnings=FALSE, recursive=TRUE)

# --- Survival endpoint mapping (manually verified per cohort) ---
ENDPOINTS <- list(
  GSE17536 = list(time = "os_time", event = "os_event", type = "OS"),
  GSE14333 = list(time = "dfs_time", event = "dfs_event", type = "DFS*"),
  GSE38832 = list(time = "dss_time", event = "dss_event", type = "DSS*"),
  GSE33113 = list(time = "time_to_recurrence", event = NA, type = "Recurrence*"),
  GSE39582 = list(time = "os_time", event = "os_event", type = "OS")
)

cat(sprintf("Endpoint mapping:\n"))
for (cohort in names(ENDPOINTS)) {
  ep <- ENDPOINTS[[cohort]]
  cat(sprintf("  %-12s: %s (%s)\n", cohort, ep$time, ep$type))
}

# --- Load coefficients ---
selected <- read_csv(file.path(BASE, "results/v31_run01_seed42/selected_genes.csv"),
                     show_col_types=FALSE)
coefs_r1 <- setNames(selected$coefficient, selected$gene)
selected2 <- read_csv(file.path(BASE, "results/v31_run02_seed123/selected_genes.csv"),
                      show_col_types=FALSE)
coefs_r2 <- setNames(selected2$coefficient, selected2$gene)

# --- Load probe→symbol ---
probe_map <- read_csv(file.path(GEO, "GPL570_probe_to_symbol.csv"), show_col_types=FALSE)

# --- Functions ---
load_geo_matrix <- function(cohort) {
  f <- file.path(GEO, paste0(cohort, "_series_matrix.txt.gz"))
  read.table(gzfile(f), header=TRUE, sep="\t", row.names=1, 
             check.names=FALSE, comment.char="!", fill=TRUE)
}

probe_to_gene <- function(pmat, pmap) {
  symbols <- pmap$gene_symbol[match(rownames(pmat), pmap$probe_id)]
  valid <- !is.na(symbols) & symbols != ""
  pmat <- pmat[valid, ]; gene_sym <- symbols[valid]
  gmat <- aggregate(pmat, by=list(gene_sym), FUN=mean, na.rm=TRUE)
  rownames(gmat) <- gmat$Group.1; gmat$Group.1 <- NULL
  as.matrix(gmat)
}

# --- Process cohort ---
process_cohort <- function(cohort, coefs, run_label) {
  ep <- ENDPOINTS[[cohort]]
  cat(sprintf("\n=== %s (%s, endpoint=%s) ===\n", cohort, run_label, ep$type))
  
  # Load survival
  surv <- read_csv(file.path(GEO, paste0(cohort, "_survival.csv")), show_col_types=FALSE)
  sample_col <- names(surv)[1]
  
  # Get time and event
  time_vals <- as.numeric(surv[[ep$time]])
  event_vals <- if (is.na(ep$event)) rep(1, nrow(surv)) else as.numeric(surv[[ep$event]])
  surv_samples <- surv[[sample_col]]
  
  n_total <- sum(!is.na(time_vals) & time_vals > 0)
  n_events <- sum(event_vals == 1, na.rm=TRUE)
  cat(sprintf("  Data: %d samples, %d events\n", n_total, n_events))
  
  # Load expression
  pmat <- load_geo_matrix(cohort)
  gmat <- probe_to_gene(pmat, probe_map)
  
  common_genes <- intersect(names(coefs), rownames(gmat))
  if (length(common_genes) < 4) return(NULL)
  
  # Risk score
  risk <- rep(0, ncol(gmat))
  for (g in common_genes) risk <- risk + coefs[g] * gmat[g, ]
  risk <- scale(risk)[,1]
  
  # Match
  common_s <- intersect(surv_samples, colnames(gmat))
  si <- match(common_s, surv_samples)
  ei <- match(common_s, colnames(gmat))
  
  t2 <- time_vals[si]; e2 <- event_vals[si]; r2 <- risk[ei]
  valid <- !is.na(t2) & t2 > 0 & !is.na(e2) & !is.na(r2)
  
  if (sum(valid) < 20) return(NULL)
  
  # Cox
  fit <- coxph(Surv(t2[valid], e2[valid]) ~ r2[valid])
  s <- summary(fit)
  hr <- s$conf.int[1,1]; hr_ci <- s$conf.int[1,3:4]
  p <- s$waldtest["pvalue"]; c_idx <- s$concordance["C"]
  n_ev <- sum(e2[valid])
  
  cat(sprintf("  HR=%.2f [%.2f-%.2f], p=%.4f, C=%.3f, events=%d\n",
              hr, hr_ci[1], hr_ci[2], p, c_idx, n_ev))
  
  list(cohort=cohort, run=run_label, n=sum(valid), genes=length(common_genes),
       HR=hr, p=p, CI_low=hr_ci[1], CI_high=hr_ci[2], C_index=c_idx,
       events=n_ev, endpoint=ep$type, status="ok")
}

# --- Run ---
all_results <- list()
for (cohort in names(ENDPOINTS)) {
  r <- process_cohort(cohort, coefs_r1, "Run1_seed42")
  if (!is.null(r)) all_results[[length(all_results)+1]] <- r
}
for (cohort in names(ENDPOINTS)) {
  r <- process_cohort(cohort, coefs_r2, "Run2_seed123")
  if (!is.null(r)) all_results[[length(all_results)+1]] <- r
}

# --- Save ---
results_df <- bind_rows(lapply(all_results, as_tibble))
write_csv(results_df, file.path(RESULTS, "geo_validation_fixed_v2.csv"))

# --- Meta-analysis (OS-only: GSE17536 + GSE39582) ---
cat("\n================================\n")
cat("OS-only cohorts: GSE17536 + GSE39582\n")
os_data <- results_df %>% filter(endpoint == "OS", !is.na(HR))
library(meta)
if (nrow(os_data) >= 2) {
  meta_res <- metagen(TE=log(os_data$HR), 
                       seTE=(log(os_data$CI_high)-log(os_data$CI_low))/3.92,
                       studlab=paste(os_data$cohort, os_data$run, sep="_"),
                       sm="HR", common=FALSE, random=TRUE)
  cat(sprintf("RE pooled HR: %.3f [%.3f-%.3f], I2=%.1f%%\n",
      exp(meta_res$TE.random), exp(meta_res$lower.random), exp(meta_res$upper.random), meta_res$I2))
}

cat("\n================================\n")
cat("All 5 cohorts (mixed endpoints):\n")
all_data <- results_df %>% filter(!is.na(HR))
meta_all <- metagen(TE=log(all_data$HR),
                     seTE=(log(all_data$CI_high)-log(all_data$CI_low))/3.92,
                     studlab=paste(all_data$cohort, all_data$run, sep="_"),
                     sm="HR", common=FALSE, random=TRUE)
cat(sprintf("RE pooled HR: %.3f [%.3f-%.3f], I2=%.1f%%\n",
    exp(meta_all$TE.random), exp(meta_all$lower.random), exp(meta_all$upper.random), meta_all$I2))

cat("\nDone!\n")
