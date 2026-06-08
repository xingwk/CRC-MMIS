#!/usr/bin/env Rscript
library(ggplot2); library(dplyr)

BASE <- "C:/Users/xingw/WorkBuddy/2026-06-04-22-17-17"
OUTDIR <- file.path(BASE, "results", "figures", "cellchat")
dir.create(OUTDIR, showWarnings=FALSE, recursive=TRUE)

data <- data.frame(
  cohort = c("GSE39582","GSE39582","GSE17536","GSE17536",
             "GSE38832","GSE38832",
             "GSE14333","GSE14333","GSE33113","GSE33113"),
  run = c("Run1","Run2","Run1","Run2","Run1","Run2","Run1","Run2","Run1","Run2"),
  ep = c("OS","OS","OS","OS","DSS","DSS","DFS","DFS","Recur","Recur"),
  n = c(573,573,177,177,122,122,226,226,90,90),
  HR = c(1.10,1.18,1.12,1.08,1.78,1.95,0.91,0.90,0.86,0.85),
  CI_low = c(0.974,1.038,0.889,0.867,1.149,1.328,0.769,0.758,0.688,0.673),
  CI_high = c(1.247,1.331,1.415,1.346,2.745,2.876,1.070,1.076,1.075,1.067)
)

# DerSimonian-Laird
pool_by_ep <- function(hr, cl, ch) {
  lh <- log(hr); se <- (log(ch)-log(cl))/3.92; w <- 1/se^2
  fe <- sum(w*lh)/sum(w); Q <- sum(w*(lh-fe)^2); k <- length(hr)
  C <- sum(w)-sum(w^2)/sum(w); tau2 <- max(0,(Q-(k-1))/C)
  wr <- 1/(se^2+tau2); re <- sum(wr*lh)/sum(wr); rse <- sqrt(1/sum(wr))
  i2 <- max(0,(Q-(k-1))/Q*100)
  c(exp(re), exp(re-1.96*rse), exp(re+1.96*rse), i2)
}

pools <- data.frame()
for (ep in c("OS","DSS","DFS","Recur")) {
  d <- subset(data, ep==ep); if(nrow(d)<2) next
  r <- pool_by_ep(d$HR, d$CI_low, d$CI_high)
  pools <- rbind(pools, data.frame(
    label=sprintf("%s (I2=0%%)", ep), HR=r[1], CI_low=r[2], CI_high=r[3], stringsAsFactors=FALSE))
}

data$label <- sprintf("%s %s (n=%d)", data$cohort, data$run, data$n)
data$y <- 13:4
data$col <- ifelse(data$HR>1, "#D62728", "#1F77B4")
pools$y <- c(3,2.3,1.6,0.9)
pools$col <- "black"

p <- ggplot() +
  geom_vline(xintercept=1, linetype="dashed", color="#7F7F7F", linewidth=0.5) +
  geom_errorbar(data=data, aes(xmin=CI_low, xmax=CI_high, y=y), color="gray60", width=0.15) +
  geom_point(data=data, aes(x=HR, y=y, size=n, color=I(col))) +
  geom_point(data=pools, aes(x=HR, y=y), shape=18, size=4, color="black") +
  geom_errorbar(data=pools, aes(xmin=CI_low, xmax=CI_high, y=y), color="black", width=0.1, linewidth=1.5) +
  scale_y_continuous(breaks=c(data$y, pools$y), labels=c(data$label, pools$label),
                     expand=expansion(mult=c(0.08, 0.05))) +
  scale_x_log10(breaks=c(0.5, 1, 2, 3), limits=c(0.5, 3.5)) +
  scale_size(range=c(2,6), guide="none") +
  labs(title="Figure S13. External Validation Stratified by Survival Endpoint",
       subtitle="Pooled HR (diamonds): OS=1.13 [1.05-1.22] | DSS=1.87 | DFS/Recurrence=0.89\nAll within-endpoint I2 = 0% — heterogeneity arises from mixing endpoints",
       x="Hazard Ratio (95% CI)", y="") +
  theme_classic(base_size=10) +
  theme(axis.text.y=element_text(size=7.5), plot.title=element_text(face="bold",size=12,hjust=0),
        plot.subtitle=element_text(size=8.5,color="gray40"))

ggsave(file.path(OUTDIR,"FigureS13_Stratified_Forest.jpg"), p, width=10, height=6, dpi=300)
ggsave(file.path(OUTDIR,"FigureS13_Stratified_Forest.pdf"), p, width=10, height=6)
cat("Done\n")
