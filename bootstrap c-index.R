Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
library(magrittr)
library(tibble)
library(stringr)
library("glmnet")
library("survival")
library(ggplot2)
library(survminer) 
library(timeROC)

load("./3.lasso_cox.rdata")


index = c('CCR7','MS4A1','IRF4','ADAMDEC1','IGHV3-30','AIM2','SLAMF7','GAPT','IGKC','UBD','ADAM28','MMP7','ICOS','IL7R','CCL19','VTCN1','MMP9','SEMA3C','PROM1','CFTR','STMN2')
index %in% rownames(Tpm_TLS_ornot) %>% table()

choose_df = Tpm_TLS_ornot[index,]
choose_df = as.data.frame(t(choose_df))
choose_df = choose_df[match(Clinical_TLS_ornot$sample,rownames(choose_df)),]
choose_df %>% dim()
Clinical_TLS_ornot %>% dim()


set.seed(123)          # 可重复
B <- 1000              # Bootstrap 次数
n <- nrow(mysurv)        # 样本数

# 记录每个基因被选中的次数
gene_freq <- rep(0, ncol(t(myexpr)))
names(gene_freq) <- colnames(t(myexpr))
expr = t(myexpr)
# 用于 Cox 的 y 对象

for (b in 1:B) {
  
  # 1) 有放回抽样
  idx <- sample(1:n, n, replace = TRUE)
  x_b <- expr[idx, , drop = FALSE]
  y_b <- mysurv$group[idx]
  
  # 2) 交叉验证选 lambda
  cvfit_b <- tryCatch(
    cv.glmnet(x_b,  y_b, 
              #10倍交叉验证，非必须限定条件，这篇文献有，其他文献大多没提
              family='cox', alpha = 1
    ),
    error = function(e) NULL
  )
  if (is.null(cvfit_b)) next
  
  # 3) 提取非零系数基因
  coef_b <- coef(cvfit_b, s = cvfit_b$lambda.min)
  selected_b <- rownames(coef_b)[as.numeric(coef_b) != 0]
  
  # 4) 累计入选次数
  gene_freq[selected_b] <- gene_freq[selected_b] + 1
  
  if (b %% 100 == 0) cat("Bootstrap:", b, "\n")
}

# 入选频率
gene_freq_pct <- gene_freq / B
gene_freq_pct <- sort(gene_freq_pct, decreasing = TRUE)

# 保留入选频率 > 50% 的基因
stable_genes <- names(gene_freq_pct[gene_freq_pct > 0.5])
print(stable_genes)


{
  gene_freq_pct_df <- data.frame(
    gene  = names(gene_freq_pct),
    freq  = as.numeric(gene_freq_pct),
    stringsAsFactors = FALSE
  )
  
  # 按频率从高到低排序，固定基因顺序
  gene_freq_pct_df$gene <- factor(gene_freq_pct_df$gene, levels = gene_freq_pct_df$gene[order(gene_freq_pct_df$freq, decreasing = TRUE)])
  
  # 标记是否 > 50%
  gene_freq_pct_df$selected <- ifelse(gene_freq_pct_df$freq > 0.5, "Stable (>50%)", "Unstable (≤50%)")
  gene_freq_pct_df$selected <- factor(gene_freq_pct_df$selected, levels = c("Stable (>50%)", "Unstable (≤50%)"))
  
  p <- ggplot(gene_freq_pct_df, aes(x = gene, y = freq, fill = selected)) +
    geom_col(width = 0.7, color = "grey20", linewidth = 0.3) +
    geom_hline(yintercept = 0.5, linetype = "dashed",
               color = "red", linewidth = 0.7) +
    geom_text(aes(label = sprintf("%.2f", freq)),
              vjust = -0.4, size = 3, color = "grey20") +
    scale_fill_manual(values = c("Stable (>50%)"   = "#2C7FB8",
                                 "Unstable (≤50%)" = "#BDBDBD")) +
    scale_y_continuous(limits = c(0, 1.08),
                       breaks = seq(0, 1, 0.2),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(
      x = NULL,
      y = "Bootstrap selection frequency",
      fill = NULL,
      title = "Bootstrap stability of LASSO-Cox selected genes",
      subtitle = "1000 bootstrap resamples; dashed line = 50% cutoff"
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.text.x        = element_text(angle = 45, hjust = 1,
                                        color = "black", face = "italic"),
      axis.text.y        = element_text(color = "black"),
      legend.position    = c(0.78, 0.85),
      legend.background  = element_rect(fill = "white", color = "grey70"),
      plot.title         = element_text(face = "bold", hjust = 0.5),
      plot.subtitle      = element_text(hjust = 0.5, color = "grey30")
    )
  
  p
}

# -------------------------- 生存分析
# TCGA 队列
{
  tcga_info = signature %>% cbind(.,Clinical_TLS_ornot)
  risk_score = tcga_info$.
  cutoff <- median(risk_score)
  risk_group <- ifelse(risk_score > cutoff, "High", "Low")
  risk_group <- factor(risk_group, levels = c("Low", "High"))
  
  surv_df <- data.frame(
    time       = tcga_info$OS.time,
    status     = tcga_info$OS,
    risk_score = risk_score,
    risk_group = risk_group
  )
  
  fit_km <- survfit(Surv(time, status) ~ risk_group, data = surv_df)
  
  ggsurvplot(
    fit_km,
    data = surv_df,
    pval = TRUE,             # log-rank p
    risk.table = TRUE,
    conf.int = TRUE,
    xlab = "Time",
    ylab = "Survival probability",
    legend.title = "Risk",
    palette = c("#2E9FDF", "#E7B800")
  )
  
  sd <- survdiff(Surv(time, status) ~ risk_group, data = surv_df)
  pval <- 1 - pchisq(sd$chisq, df = length(sd$n) - 1)
  cat("Log-rank p =", pval, "\n")
  
  cindex_train <- concordance(Surv(time, status) ~ risk_score, data = surv_df)
  cat("C-index (training) =", cindex_train$concordance, "\n")
}
# ICGC队列
{
  # 验证集：expr_val, surv_val（样本顺序一致）
  risk_score_val <- for_xt$V1
  names(risk_score_val) <- rownames(for_xt)
  
  surv_val_df <- data.frame(
    time       = for_xt$SUR,
    status     = for_xt$STATUS,
    risk_score = risk_score_val,
    risk_group = factor(ifelse(risk_score_val > cutoff, "High", "Low"),
                        levels = c("Low", "High"))
  )
  
  # KM
  fit_km_val <- survfit(Surv(time, status) ~ risk_group, data = surv_val_df)
  ggsurvplot(fit_km_val, data = surv_val_df, pval = TRUE, risk.table = TRUE)
  
  # log-rank
  sd_val <- survdiff(Surv(time, status) ~ risk_group, data = surv_val_df)
  pval_val <- 1 - pchisq(sd_val$chisq, df = length(sd_val$n) - 1)
  
  # C-index
  cindex_val <- concordance(Surv(time, status) ~ I(-risk_score),
                            data = surv_val_df)$concordance
  cat("Validation C-index =", cindex_val, "\n")
}

# ---------------------- c-index绘图
{
  library(survival)
  
  # 计算训练集 C-index
  c_train <- concordance(Surv(time, status) ~ I(-risk_score),
                         data = surv_df)
  c_train_val <- c_train$concordance
  c_train_se  <- sqrt(c_train$var)
  c_train_lo  <- c_train_val - 1.96 * c_train_se
  c_train_hi  <- c_train_val + 1.96 * c_train_se
  
  # 计算验证集 C-index
  c_val <- concordance(Surv(time, status) ~ I(-risk_score),
                       data = surv_val_df)
  c_val_val <- c_val$concordance
  c_val_se  <- sqrt(c_val$var)
  c_val_lo  <- c_val_val - 1.96 * c_val_se
  c_val_hi  <- c_val_val + 1.96 * c_val_se
  
  cindex_df <- data.frame(
    Cohort  = c("TCGA-LIHC (Training)", "ICGC-LIRI-JP (Validation)"),
    Cindex  = c(c_train_val, c_val_val),
    Lower   = c(c_train_lo, c_val_lo),
    Upper   = c(c_train_hi, c_val_hi)
  )
  
  ggplot(cindex_df, aes(x = Cohort, y = Cindex, fill = Cohort)) +
    geom_col(width = 0.6, color = "grey20") +
    geom_errorbar(aes(ymin = Lower, ymax = Upper),
                  width = 0.15, color = "grey20") +
    geom_text(aes(label = sprintf("%.3f", Cindex)),
              vjust = -0.6, size = 4) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
    scale_fill_manual(values = c("#2C7FB8", "#F03B20"), guide = "none") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(x = NULL, y = "C-index (95% CI)",
         title = "Model discrimination in training and validation cohorts") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          axis.text = element_text(color = "black"))
}

#时间依赖性auc
{
  roc_res <- timeROC(
    T      = surv$time,
    delta  = surv$status,
    marker = risk_score,
    cause  = 1,
    times  = c(12, 24, 36, 60),
    iid    = TRUE
  )
  
  # 提取 AUC
  auc_df <- data.frame(
    Time = c(12, 24, 36, 60),
    AUC  = roc_res$AUC
  )
  print(auc_df)
  
  # 画 ROC 曲线
  plot(roc_res, time = 12, col = "#E41A1C", title = FALSE, lwd = 2)
  plot(roc_res, time = 24, col = "#377EB8", add = TRUE, lwd = 2)
  plot(roc_res, time = 36, col = "#4DAF4A", add = TRUE, lwd = 2)
  plot(roc_res, time = 60, col = "#984EA3", add = TRUE, lwd = 2)
  legend("bottomright",
         legend = sprintf("%d-month AUC = %.3f",
                          c(12, 24, 36, 60), roc_res$AUC),
         col = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"),
         lwd = 2, bty = "n")
}
