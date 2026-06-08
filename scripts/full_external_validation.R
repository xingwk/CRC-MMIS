#!/usr/bin/env Rscript
# CRC-MMIS Full External Validation — R survival::coxph
sink("/c/Users/xingw/WorkBuddy/2026-06-04-22-17-17/results/revisions_r2/validation_R.log")
library(survival); library(dplyr); library(readr)

BASE <- "C:/Users/xingw/WorkBuddy/2026-06-04-22-17-17"
GEO <- file.path(BASE, "data", "geo")
OUT <- file.path(BASE, "results", "revisions_r2")
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

# ──────────────────────────
# Load coefficients (both seeds)
# ──────────────────────────
coefs_r1 <- setNames(
  read_csv(file.path(BASE, "results/v31_run01_seed42/selected_genes.csv"), show_col_types=FALSE)$coefficient,
  read_csv(file.path(BASE, "results/v31_run01_seed42/selected_genes.csv"), show_col_types=FALSE)$gene)
coefs_r2 <- setNames(
  read_csv(file.path(BASE, "results/v31_run02_seed123/selected_genes.csv"), show_col_types=FALSE)$coefficient,
  read_csv(file.path(BASE, "results/v31_run02_seed123/selected_genes.csv"), show_col_types=FALSE)$gene)

cat(sprintf("Run1 (seed=42): %d genes\n", length(coefs_r1)))
cat(sprintf("Run2 (seed=123): %d genes\n", length(coefs_r2)))

# ──────────────────────────
# Load probe→symbol mapping
# ──────────────────────────
probes <- read_csv(file.path(GEO, "GPL570_probe_to_symbol.csv"), show_col_types=FALSE)
pmap <- setNames(probes$gene_symbol, probes$probe_id)

# ──────────────────────────
# Endpoint definitions
# ──────────────────────────
ep_list <- list(
  GSE17536 = list(time="os_time", event="os_event", type="OS"),
  GSE14333 = list(time="dfs_time", event="dfs_event", type="DFS"),
  GSE38832 = list(time="dss_time", event="dss_event", type="DSS"),
  GSE33113 = list(time="time_to_recurrence", event=NA, type="Recurrence"),
  GSE39582 = list(time="os_time", event="os_event", type="OS")
)

# ──────────────────────────
# Helper: load GEO matrix and convert to gene symbols
# ──────────────────────────
load_geo <- function(cohort) {
  f <- file.path(GEO, paste0(cohort, "_series_matrix.txt.gz"))
  cat(sprintf("  Loading %s...", cohort))
  raw <- read.table(gzfile(f), header=TRUE, sep="\t", row.names=1,
                    check.names=FALSE, comment.char="!", fill=TRUE)
  # Map probes to genes
  syms <- pmap[rownames(raw)]
  valid <- !is.na(syms) & syms != ""
  raw <- raw[valid, ]
  syms <- syms[valid]
  # Collapse by mean
  gmat <- as.data.frame(raw)
  gmat$gene <- syms
  gmat <- aggregate(. ~ gene, data=gmat, FUN=mean, na.rm=TRUE)
  rownames(gmat) <- gmat$gene; gmat$gene <- NULL
  cat(sprintf(" %d genes × %d samples\n", nrow(gmat), ncol(gmat)))
  as.matrix(gmat)
}

# ──────────────────────────
# Process one cohort × seed combination
# ──────────────────────────
process <- function(cohort, coefs, seed_label) {
  ep <- ep_list[[cohort]]
  cat(sprintf("\n  %s / %s (%s)...", cohort, seed_label, ep$type))
  
  # Load survival
  surv <- read_csv(file.path(GEO, paste0(cohort, "_survival.csv")), show_col_types=FALSE)
  sids <- surv[[1]]  # first column = sample ID
  
  t_vals <- as.numeric(surv[[ep$time]])
  e_vals <- if (is.na(ep$event)) rep(1, nrow(surv)) else as.numeric(surv[[ep$event]])
  
  # Load expression
  gmat <- load_geo(cohort)
  
  # Match genes
  common <- intersect(names(coefs), rownames(gmat))
  cat(sprintf(" %d/%d genes", length(common), length(coefs)))
  if (length(common) < 4) return(NULL)
  
  # Risk score
  risk <- rep(0, ncol(gmat))
  for (g in common) risk <- risk + coefs[g] * gmat[g, ]
  risk <- scale(risk)[, 1]
  
  # Match samples
  expr_ids <- colnames(gmat)
  common_s <- intersect(sids, expr_ids)
  si <- match(common_s, sids)
  ei <- match(common_s, expr_ids)
  
  t2 <- t_vals[si]; e2 <- e_vals[si]; r2 <- risk[ei]
  ok <- !is.na(t2) & t2 > 0 & !is.na(e2) & !is.na(r2)
  n_ev <- sum(e2[ok])
  cat(sprintf(" N=%d E=%d", sum(ok), n_ev))
  
  if (sum(ok) < 20 || n_ev < 3) {
    cat(" → skipped\n")
    return(NULL)
  }
  
  # Cox
  fit <- coxph(Surv(t2[ok], e2[ok]) ~ r2[ok])
  s <- summary(fit)
  hr <- s$conf.int[1, 1]
  hr_lo <- s$conf.int[1, 3]
  hr_hi <- s$conf.int[1, 4]
  pval <- s$waldtest["pvalue"]
  c_idx <- s$concordance["C"]
  cat(sprintf(" HR=%.2f [%.2f-%.2f] p=%.4f C=%.3f\n", hr, hr_lo, hr_hi, pval, c_idx))
  
  data.frame(cohort=cohort, run=seed_label, n=sum(ok), genes=length(common),
             HR=hr, HR_low=hr_lo, HR_high=hr_hi, p_value=pval,
             C_index=c_idx, events=n_ev, endpoint=ep$type, stringsAsFactors=FALSE)
}

# ──────────────────────────
# RUN ALL
# ──────────────────────────
cat("\n==================== FULL VALIDATION ====================\n")
results <- list()
for (cohort in names(ep_list)) {
  r <- process(cohort, coefs_r1, "Run1")
  if (!is.null(r)) results[[length(results)+1]] <- r
}
for (cohort in names(ep_list)) {
  r <- process(cohort, coefs_r2, "Run2")
  if (!is.null(r)) results[[length(results)+1]] <- r
}

all_r <- bind_rows(results)
write_csv(all_r, file.path(OUT, "geo_validation_full_R.csv"))
cat(sprintf("\nSaved: %d results\n", nrow(all_r)))

# ──────────────────────────
# SUMMARY TABLE
# ──────────────────────────
cat("\n==================== SUMMARY ====================\n")
cat(sprintf("%-12s %-6s %5s %4s %8s %-15s %7s %s\n",
            "Cohort","Run","N","Ev","HR","95% CI","p","Endpoint"))
cat(paste(rep("-",75), collapse=""),"\n")
for (i in 1:nrow(all_r)) {
  r <- all_r[i,]
  sig <- if(r$p_value < 0.05) " *" else ""
  cat(sprintf("%-12s %-6s %5d %4d %6.2f [%5.2f-%5.2f] %6.4f %s%s\n",
              r$cohort, r$run, r$n, r$events, r$HR, r$HR_low, r$HR_high, r$p_value, r$endpoint, sig))
}

# ──────────────────────────
# META-ANALYSIS (manual DerSimonian-Laird)
# ──────────────────────────
DSL_meta <- function(hr_vec, se_vec) {
  ok <- !is.na(hr_vec) & !is.na(se_vec)
  log_hr <- log(hr_vec[ok])
  se <- se_vec[ok]
  w <- 1 / se^2
  fe_pooled <- sum(w * log_hr) / sum(w)
  Q <- sum(w * (log_hr - fe_pooled)^2)
  df <- length(log_hr) - 1
  C <- sum(w) - sum(w^2) / sum(w)
  tau2 <- max(0, (Q - df) / C)
  w_re <- 1 / (se^2 + tau2)
  re_pooled <- sum(w_re * log_hr) / sum(w_re)
  re_se <- sqrt(1 / sum(w_re))
  I2 <- max(0, (Q - df) / Q * 100)
  if (is.nan(I2)) I2 <- 0
  list(HR = exp(re_pooled), HR_low = exp(re_pooled - 1.96*re_se), 
       HR_high = exp(re_pooled + 1.96*re_se), I2=I2, Q=Q, df=df, tau2=tau2, n=length(log_hr))
}

cat("\n==================== META ====================\n")
# OS-only
os <- all_r %>% filter(endpoint == "OS", !is.na(HR))
if (nrow(os) >= 2) {
  m_os <- DSL_meta(os$HR, (log(os$HR_high)-log(os$HR_low))/3.92)
  cat(sprintf("\n--- OS-only (%d studies: %s) ---\n", m_os$n, paste(os$cohort, os$run, sep="/", collapse=", ")))
  cat(sprintf("RE pooled HR: %.3f [%.3f-%.3f]\n", m_os$HR, m_os$HR_low, m_os$HR_high))
  cat(sprintf("I2 = %.1f%%, Q = %.1f, df = %d\n", m_os$I2, m_os$Q, m_os$df))
}

# All cohorts
all_d <- all_r %>% filter(!is.na(HR))
if (nrow(all_d) >= 2) {
  m_all <- DSL_meta(all_d$HR, (log(all_d$HR_high)-log(all_d$HR_low))/3.92)
  cat(sprintf("\n--- All %d cohorts (mixed endpoints) ---\n", m_all$n))
  cat(sprintf("RE pooled HR: %.3f [%.3f-%.3f]\n", m_all$HR, m_all$HR_low, m_all$HR_high))
  cat(sprintf("I2 = %.1f%%, Q = %.1f, df = %d\n", m_all$I2, m_all$Q, m_all$df))
}

# ──────────────────────────
# Generate forest plot CSV (matching old format)
# ──────────────────────────
fcsv <- all_d %>%
  mutate(logHR = log(HR), se = (log(HR_high) - log(HR_low))/3.92) %>%
  select(cohort, run, n, genes=endpoint, HR, p_value, logHR, se, HR_low, HR_high, events, endpoint)
write_csv(fcsv, file.path(OUT, "forest_plot_corrected.csv"))

# Also generate the simplified forest_plot_meta.csv format
simple <- all_d %>%
  mutate(cohort_display = paste0(cohort, " (", endpoint, ")"),
         p_display = round(p_value, 3)) %>%
  select(cohort = cohort_display, run, n, genes = events, HR, p = p_display,
         HR_low, HR_high)
write_csv(simple, file.path(OUT, "forest_plot_meta_updated.csv"))

cat("\n=== DONE ===\n")
sink()
cat(sprintf("Full results: %s\n", file.path(OUT, "geo_validation_full_R.csv")))
cat(sprintf("Forest plot CSV: %s\n", file.path(OUT, "forest_plot_meta_updated.csv")))
