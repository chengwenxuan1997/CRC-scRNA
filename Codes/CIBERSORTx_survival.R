# path <- "i:/genomicdata/jiaoda/CIBERSORTx/workflow";setwd(path)
# sig.mat <- read.table("Step1/Step1_SignatureMatrix.txt", sep = "\t", row.names = 1, header = T)
# mix.mat <- read.table("Step2/bulk.txt", sep = "\t", row.names = 1, header = T)
# 
# source("Step2/CIBERSORT.R")
# library(tibble)
# library(preprocessCore)
# library(e1071)
# Fraction <- CIBERSORT(sig_matrix = sig.mat,
#                     mixture_file = mix.mat,
#                     perm = 10,
#                     absolute = F,
#                     QN = T)

path <- "I:/genomicdata/jiaoda/CIBERSORTx";setwd(path)
data <- readRDS("I:/genomicdata/jiaoda/pseudotime/data.rds")
CD8T_info <- readRDS("I:/genomicdata/jiaoda/pseudotime/CD8T_info.rds")
all(rownames(CD8T_info)==colnames(data))
data <- data[, CD8T_info$celltypes4 != "Undetermined"]
data <- round(data)
label <- t(CD8T_info$celltypes4[CD8T_info$celltypes4 != "Undetermined"])
# data <- data[, CD8T_info$celltypes4 %in% c("CD103+ TRM", "ITGB2+ TRM", "non-TRM")]
# label <- t(CD8T_info$celltypes4[CD8T_info$celltypes4 %in% c("CD103+ TRM", "ITGB2+ TRM", "non-TRM")])
data <- rbind(label, as.matrix(data))
write.table(data, file = "data.txt", sep = "\t", quote = F, row.names = T, col.names = F)

library(Seurat)
CD8T0826 <- readRDS("I:/genomicdata/jiaoda/pseudotime/CD8T0826.rds")
CD8T_info <- CD8T0826@meta.data
data <- GetAssayData(CD8T0826, slot = "data")
data <- data[, CD8T_info$celltypes4 != "Undetermined"]
# data <- exp(data)-1
data <- 2^data-1
# data <- exp(data)
label <- t(CD8T_info$celltypes4[CD8T_info$celltypes4 != "Undetermined"])
rm(CD8T0826)
data <- rbind(label, as.matrix(round(data, digits = 3)))
write.table(data, file = "data2.txt", sep = "\t", quote = F, row.names = T, col.names = F)


library(survival)
library(survminer)

# Fraction <- read.table("Step2/CIBERSORTx_Job4_Results.txt", sep = "\t", header = T, row.names = 1)
Fraction <- read.table("Step2/CIBERSORTx_Job6_Adjusted.txt", sep = "\t", header = T, row.names = 1)
# Fraction <- read.table("Step2/FractionWithBatchEffectRemoved.txt", sep = "\t", header = T, row.names = 1)
Fraction <- read.table("Step2/tcga_gep_corrected.txt", sep = "\t", header = T, row.names = 1)

num = 698
plotdata <- data.frame(
  "celltypes" = c(rep("ITGB2+ TRM", num), rep("non TRM", num), rep("CD103+ TRM", num), rep("CD8+ Exhausted", num), rep("CD8+ Naive", num), rep("CD8+ Cycling", num)),
  "fraction" = c(Fraction[,1], Fraction[,2], Fraction[,3], Fraction[,4], Fraction[,5], Fraction[,6])
)
plotdata$celltypes <- factor(x = plotdata$celltypes,
                             levels = c("ITGB2+ TRM", "non TRM", "CD103+ TRM", "CD8+ Exhausted", "CD8+ Naive", "CD8+ Cycling"))
ggplot(data = plotdata, aes(x = celltypes, y = fraction, fill = celltypes))+
  geom_violin(alpha = 0.4, position = position_dodge(width = .75), 
              size = 0.8, color="white") +
  geom_point(shape = 21, size=1.2, # 点的性状和大小
             position = position_jitterdodge(), # 让点散开
             aes(color=celltypes), alpha = 0.6) +
  geom_boxplot(notch = TRUE, outlier.size = -1, 
               color="black", lwd = 0.6, alpha = 0.7) +
  theme_bw()

# Fraction <- read.table("Step2/FractionWithBatchEffectRemoved.txt", sep = "\t", header = T, row.names = 1)
# GSE39582.surv <- readRDS("i:/genomicdata/jiaoda/CIBERSORTx/workflow/Step2/GSE39582.surv.rds")
surv <- readRDS("i:/genomicdata/jiaoda/CIBERSORTx/workflow/Step2/tcga.surv.rds")
rownames(Fraction) <- gsub("\\.","-", rownames(Fraction))
all(rownames(surv) %in% rownames(Fraction))
Fraction <- Fraction[rownames(surv),]

tmp <- cbind(Fraction[,c(1:6)], surv)
rfs.cox.res <- lapply(colnames(tmp)[1:6], function(x){
  coxph(Surv(RFS.time,RFS) ~ tmp[[x]], data = tmp)
})
rfs.cox.summary <- lapply(rfs.cox.res, summary)
unlist(lapply(rfs.cox.summary, function(x){x$logtest["pvalue"]}))

pfi.cox.res <- lapply(colnames(tmp)[1:6], function(x){
  coxph(Surv(PFI.time,PFI) ~ tmp[[x]], data = tmp)
})
pfi.cox.summary <- lapply(pfi.cox.res, summary)
unlist(lapply(pfi.cox.summary, function(x){x$logtest["pvalue"]}))

os.cox.res <- lapply(colnames(tmp)[1:6], function(x){
  coxph(Surv(OS.time,OS) ~ tmp[[x]], data = tmp)
})
os.cox.summary <- lapply(os.cox.res, summary)
unlist(lapply(os.cox.summary, function(x){x$logtest["pvalue"]}))

plotdata <- data.frame(
  row.names = rownames(Fraction),
  "score" = Fraction[,3],
  "Group" = ifelse(Fraction[,3] > quantile(Fraction[,3], 0.5),"High","Low"),
  surv[rownames(Fraction),]
)

summary(coxph(Surv(OS.time, OS) ~ plotdata$Group, data = plotdata))$logtest[["pvalue"]]

fit <- survfit(Surv(OS.time, OS) ~ Group,
                data = plotdata,
                type = "kaplan-meier",
                error     = "greenwood",
                conf.type = "plain",
                na.action = na.exclude)
ggsurvplot(fit = fit,
           data = plotdata,
           conf.int          = FALSE,
           risk.table        = TRUE,
           risk.table.col    = "strata",
           palette           = "jco",
           size              = 1,
           legend.title      = "",
           pval              = TRUE,
           surv.median.line  = "hv",
           xlab              = "Time (Days)",
           ylab              = "Survival probability (%)",
           risk.table.y.text = FALSE)



##########
a<-lapply(seq(0.20,0.8,0.01), function(x){
  plotdata <- data.frame(
    row.names = NULL,
    # "score" = data,
    # "group" = ifelse(data$TRM > quantile(data$TRM, probs = x), "high", "low"),
    "group" = ifelse(data$ITGB2 > quantile(data$ITGB2, probs = x), "high", "low"),
    data
  )
  
  summary(coxph(Surv(OS.time, OS) ~ group, data = plotdata))$logtest[["pvalue"]]
  # summary(coxph(Surv(RFS.time, RFS) ~ group, data = plotdata))$logtest[["pvalue"]]
  # summary(coxph(Surv(PFI.time, PFI) ~ group, data = plotdata))$logtest[["pvalue"]]
})
plot(x = seq(0.20,0.8,0.01), y = a)
