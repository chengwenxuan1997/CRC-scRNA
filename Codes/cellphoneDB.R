# export h5ad for cellphone
path <- "~/genomicdata/jiaoda/CellphoneDB/"
library(SeuratDisk)
library(Seurat)
library(magrittr)
library(reshape2)
# use function Prep10X, postCellphoneDB and ReadCellphoneDB

expr <- readRDS("~/genomicdata/jiaoda/expr0617.rds")
cell_info <- readRDS("~/genomicdata/jiaoda/cell_info.rds")
gene_inference <- read.table("~/genomicdata/jiaoda/raw/genes.tsv", sep = "\t", header = T)
B <- table(subset(cell_info,type == "B")["celltypes"])
B <- names(B[B>0])
epi <- table(subset(cell_info,type == "epithelic")["celltypes"])
epi <- names(epi[epi>0])
mye <- table(subset(cell_info,type == "myeloid")["celltypes"])
mye <- names(mye[mye>0])
stro <- table(subset(cell_info,type == "stromal")["celltypes"])
stro <- names(stro[stro>0])
CD8T <- table(subset(cell_info,division == "CD8+ T cells")["celltypes"])
CD8T <- names(CD8T[CD8T>0])
panT <- table(subset(cell_info,type == "panT")["celltypes"])
panT <- names(panT[panT>0])

sample_info <- subset(cell_info, type %in% c("myeloid","panT"))
sample <- expr[,rownames(sample_info)]
write.table(sample,paste0(path,"count.txt"),sep = "\t")
sample_info <- sample_info["celltypes"]
write.table(sample_info,paste0(path,"info.txt"))

expr <- read.table("~/genomicdata/learning/cellphoneDB/test_counts.txt", header = T, row.names = 1)
seu <- CreateSeuratObject(expr)
seu <- NormalizeData(seu)
SaveH5Seurat(seu,"~/genomicdata/learning/cellphoneDB/test_counts")
Convert("~/genomicdata/learning/cellphoneDB/test_counts.h5seurat",dest = "h5ad")
expr <- GetAssayData(seu, slot = "count")
expr <- as.matrix(expr)
expr <- cbind(rownames(expr),expr)
colnames(expr)[1] <- "Gene"
rownames(expr) <- NULL
write.table(expr,"~/genomicdata/learning/cellphoneDB/test1.txt",sep = "\t",quote = F,row.names = F)

info <- data.frame("Cell" = rownames(sample_info),
                   "cell_type" = sample_info$celltypes)
write.table(info, paste0(path,"info.txt"), sep = "\t", quote = F, row.names = F)
seu <- CreateSeuratObject(sample)%>%NormalizeData()
norm.count <- seu@assays$RNA@data
rownames(norm.count)<-gene_inference$gene_ids[match(rownames(norm.count),gene_inference$gene)]

dir.create(paste0(path,"test1"))
features <- rownames(norm.count)
write.table(features,paste0(path,"test1/features.tsv"),sep = "\t", quote = F, row.names = F, col.names = F)
barcodes <- colnames(norm.count)
write.table(barcodes,paste0(path,"test1/barcodes.tsv"),sep = "\t", quote = F, row.names = F, col.names = F)
matrix <- as(norm.count,"dgTMatrix")
rownames(matrix) <- NULL;colnames(matrix) <- NULL
Matrix::writeMM(matrix,paste0(path,"test1/matrix.mtx"))

# preparation for cellphoneDB
Prep10X <- function(project,cells){
  if (dir.exists(paste0(path,project))) {
    unlink(paste0(path,project),recursive = T, force = T)
    dir.create(paste0(path,project))
  }else{dir.create(paste0(path,project))}
  
  sample_info <- subset(cell_info, celltypes %in% cells)
  if (strsplit(project,"_")[[1]][3]=="NC") sample_info <- subset(sample_info, group == "NC")
  if (strsplit(project,"_")[[1]][3]=="TM") sample_info <- subset(sample_info, group != "NC")
  info <- data.frame("Cell" = rownames(sample_info),
                     "cell_type" = sample_info$celltypes)
  sample <- expr[,rownames(sample_info)]
  seu <- CreateSeuratObject(sample)%>%NormalizeData()
  norm.count <- seu@assays$RNA@data
  rownames(norm.count)<-gene_inference$gene_ids[match(rownames(norm.count),gene_inference$gene)]
  features <- rownames(norm.count)
  barcodes <- colnames(norm.count)
  matrix <- as(norm.count,"dgTMatrix")
  
  print(dim(matrix))
  print("create metadata")
  write.table(info, paste0(path,project,"/info.txt"), sep = "\t", quote = F, row.names = F)
  dir.create(paste0(path,project,"/data"))
  print("create features.tsv")
  write.table(features,paste0(path,project,"/data/features.tsv"),sep = "\t", quote = F, row.names = F, col.names = F)
  print("create barcodes.tsv")
  write.table(barcodes,paste0(path,project,"/data/barcodes.tsv"),sep = "\t", quote = F, row.names = F, col.names = F)
  print("create matrix")
  Matrix::writeMM(matrix,paste0(path,project,"/data/matrix.mtx"))
  print(paste0("cd ~/genomicdata/jiaoda/CellphoneDB/",project))
  print(paste0("cellphonedb method statistical_analysis info.txt data"))
}

Prep10X("epi_mye",c(epi,mye))
Prep10X("epi_panT",c(epi,panT))
Prep10X("stro_mye",c(stro,mye))
Prep10X("stro_panT",c(stro,panT))
Prep10X("stro_B",c(stro,B))
Prep10X("B_CD8T_NC",c(B,CD8T))
Prep10X("B_CD8T_TM",c(B,CD8T))
Prep10X("B_mye_NC",c(B,mye))
Prep10X("B_mye_TM",c(B,mye))
deconvoluted <- read.table(paste0(path,"epi_mye/out/deconvoluted.txt"), sep = "\t", header = T)
means <- read.table(paste0(path,"epi_mye/out/means.txt"), sep = "\t", header = T)
pvalues <- read.table(paste0(path,"epi_mye/out/pvalues.txt"), sep = "\t", header = T)
sig <- read.table(paste0(path,"epi_mye/out/significant_means.txt"), sep = "\t", header = T)

sig_melt <- melt(sig, id.vars = colnames(sig)[1:12])
pvalues_melt <- melt(pvalues, id.vars = colnames(pvalues)[1:11])
sig_melt$id <- paste0(sig_melt$id_cp_interaction,sig_melt$variable)
pvalues_melt$id <- paste0(pvalues_melt$id_cp_interaction,pvalues_melt$variable)
sum(duplicated(pvalues_melt$id))
sig_melt$pvalue <- pvalues_melt$value[match(sig_melt$id,pvalues_melt$id)]
sig_melt <- subset(sig_melt, !is.na(value))
a<- subset(sig_melt, gene_a != "" & gene_b != "")
sig_melt$gene_a <- gene_inference$gene[match(sig_melt$gene_a,gene_inference$gene_ids)]
sig_melt$gene_b <- gene_inference$gene[match(sig_melt$gene_b,gene_inference$gene_ids)]
comb <- as.data.frame(t(combn(c(epi,mye),2)))
inter <- data.frame("source" = c(epi,mye,comb$V1,comb$V2), "target" = c(epi,mye,comb$V2,comb$V1))
inter$interaction <- paste0(inter$source,".",inter$target)
inter$interaction <- gsub("\\+|\\-|\\ ",".",inter$interaction)
sig_melt <- cbind(sig_melt,inter[match(sig_melt$variable,inter$interaction),])
res <- data.frame(
  "id_cp_interaction" = sig_melt$id_cp_interaction, "interacting_pair" = sig_melt$interacting_pair,
  "partner_a" = sig_melt$partner_a, "partner_b" = sig_melt$partner_b,
  "gene_a" = sig_melt$gene_a, "gene_b" = sig_melt$gene_b,
  "secreted" = sig_melt$secreted, "receptor_a" = sig_melt$receptor_a, "receptor_b" = sig_melt$receptor_b,
  "annotation_strategy" = sig_melt$annotation_strategy, "is_integrin" = sig_melt$is_integrin,
  "rank" = sig_melt$rank, "source" = sig_melt$source, "target" = sig_melt$target, "means" = sig_melt$value, "pvalue" = sig_melt$pvalue
)
result <- list(
  "LR" = subset(res, source %in% epi & target %in% mye),
  "RL" = subset(res, source %in% mye & target %in% epi)
)
openxlsx::write.xlsx(x = result, file = "~/genomicdata/jiaoda/CellphoneDB/epi_mye/epi_mye.xlsx")

postCellphoneDB <- function(project, group1 = NULL, group2 = NULL){
  dir <- paste0(path,project)
  pvalues <- read.table(paste0(dir,"/out/pvalues.txt"), sep = "\t", header = T)
  sig <- read.table(paste0(dir,"/out/significant_means.txt"), sep = "\t", header = T)
  
  sig_melt <- melt(sig, id.vars = colnames(sig)[1:12])
  pvalues_melt <- melt(pvalues, id.vars = colnames(pvalues)[1:11])
  sig_melt$id <- paste0(sig_melt$id_cp_interaction,sig_melt$variable)
  pvalues_melt$id <- paste0(pvalues_melt$id_cp_interaction,pvalues_melt$variable)
  sig_melt$pvalue <- pvalues_melt$value[match(sig_melt$id,pvalues_melt$id)]
  sig_melt <- subset(sig_melt, !is.na(value))
  sig_melt$gene_a <- gene_inference$gene[match(sig_melt$gene_a,gene_inference$gene_ids)]
  sig_melt$gene_b <- gene_inference$gene[match(sig_melt$gene_b,gene_inference$gene_ids)]
  
  if (is.null(group1)) group1 <- strsplit(project,"_")[[1]][1]%>%get()
  if (is.null(group2)) group2 <- strsplit(project,"_")[[1]][2]%>%get()
  comb <- as.data.frame(t(combn(c(group1,group2),2)))
  inter <- data.frame("source" = c(group1,group2,comb$V1,comb$V2), "target" = c(group1,group2,comb$V2,comb$V1))
  inter$interaction <- paste0(inter$source,".",inter$target)
  inter$interaction <- gsub("\\+|\\-|\\ ",".",inter$interaction)
  sig_melt <- cbind(sig_melt,inter[match(sig_melt$variable,inter$interaction),])
  res <- data.frame(
    "id_cp_interaction" = sig_melt$id_cp_interaction, "interacting_pair" = sig_melt$interacting_pair,
    "partner_a" = sig_melt$partner_a, "partner_b" = sig_melt$partner_b,
    "gene_a" = sig_melt$gene_a, "gene_b" = sig_melt$gene_b,
    "secreted" = sig_melt$secreted, "receptor_a" = sig_melt$receptor_a, "receptor_b" = sig_melt$receptor_b,
    "annotation_strategy" = sig_melt$annotation_strategy, "is_integrin" = sig_melt$is_integrin,
    "rank" = sig_melt$rank, "source" = sig_melt$source, "target" = sig_melt$target, "means" = sig_melt$value, "pvalue" = sig_melt$pvalue
  )
  result <- list(
    "LR" = subset(res, source %in% group1 & target %in% group2),
    "RL" = subset(res, source %in% group2 & target %in% group1)
  )
  openxlsx::write.xlsx(x = result, file = paste0(dir,"/",project,".xlsx"))
}

postCellphoneDB("epi_panT")
postCellphoneDB("epi_mye")
postCellphoneDB("stro_mye")
postCellphoneDB("stro_panT")
postCellphoneDB("stro_B")
postCellphoneDB("B_CD8T_NC",group1 = "IgA+ plasma cells")
postCellphoneDB("B_CD8T_TM",group1 = "IgA+ plasma cells")
postCellphoneDB("B_mye_NC",group1 = "IgA+ plasma cells")
postCellphoneDB("B_mye_TM",group1 = "IgA+ plasma cells")


# graph
ReadCellphoneDB <- function(project){
  dir <- paste0(path,project,"/out/")
  means <- read.table(paste0(dir,"means.txt"), sep = "\t", header = T)
  pvalues <- read.table(paste0(dir,"pvalues.txt"), sep = "\t", header = T)
  group1 <- strsplit(project,"_")[[1]][1]%>%get()
  group2 <- strsplit(project,"_")[[1]][2]%>%get()
  comb <- as.data.frame(t(combn(c(group1,group2),2)))
  inter <- data.frame("source" = c(group1,group2,comb$V1,comb$V2), "target" = c(group1,group2,comb$V2,comb$V1))
  inter$interaction <- paste0(inter$source,".",inter$target)
  inter$interaction <- gsub("\\+|\\-|\\ ",".",inter$interaction)
  
  data <- list("means" = means, "pvalues" = pvalues)
  data <- lapply(data, function(x){
    x = melt(x, id.vars = colnames(x)[1:11])
    x$gene_a <- gene_inference$gene[match(x$gene_a,gene_inference$gene_ids)]
    x$gene_b <- gene_inference$gene[match(x$gene_b,gene_inference$gene_ids)]
    x$id <- paste0(x$id_cp_interaction,x$variable)
    x <- cbind(x,inter[match(x$variable,inter$interaction),])
    return(x)
  })
  means.melt <- data[["means"]]; pvalues.melt <- data[["pvalues"]]
  if (all(means.melt$id==pvalues.melt$id)){
    res <- data.frame(
      "id_cp_interaction" = means.melt$id_cp_interaction, "interacting_pair" = means.melt$interacting_pair,
      "partner_a" = means.melt$partner_a, "partner_b" = means.melt$partner_b,
      "gene_a" = means.melt$gene_a, "gene_b" = means.melt$gene_b,
      "secreted" = means.melt$secreted, "receptor_a" = means.melt$receptor_a, "receptor_b" = means.melt$receptor_b,
      "annotation_strategy" = means.melt$annotation_strategy, "is_integrin" = means.melt$is_integrin,
      "source" = means.melt$source, "target" = means.melt$target, "means" = means.melt$value, "pvalue" = pvalues.melt$value
    )
  }else{res = "fail to read"}
  return(res)
}


stro_B <- ReadCellphoneDB("stro_B")
stro_mye <- ReadCellphoneDB("stro_mye")
stro_panT <- ReadCellphoneDB("stro_panT")

data <- rbind(stro_mye,stro_B,stro_panT)
# data <- subset(data, gene_a == "CXCL12" & gene_b == "CXCR4" & source %in% stro & !target %in% stro)
data <- subset(data, gene_a == "CCL5" & gene_b == "CCR5" & source %in% stro & !target %in% stro)
data <- data[,c("source","target","means","pvalue")]
data$P <- cut(data$pvalue,c(1,0.05,0.01,0),include.lowest = T)
data$P <- factor(data$P,levels = c("(0.05,1]","(0.01,0.05]","[0,0.01]"))
data$source <- factor(data$source,
                      levels = rev(c("Enteric glial cells", "Pericytes", "Vascular SMCs",
                                     "EPCs", "EC-2", "EC-1", "Tip-like ECs", "Stalk-like ECs",
                                     "Myofibroblasts", "Stromal-3", "Stromal-2", "Stromal-1", "PDGRFA+ fibroblasts")
                                 ))
data$target <- factor(data$target,
                      levels = c(panT, B, mye))
data$means[data$means>2]=2
# data <- subset(data, target != "CD8+ Naïve")
data <- subset(data, target %in% mye)

# grDevices::cairo_pdf(paste0(path,"CXCL12-CXCR4.pdf"), width = 876/100, height = 618/100)
pdf(paste0(path,"CCL5-CCR5.pdf"), width = 348/100, height = 601/100)
ggplot(data, aes(x = target, y = source, col = means, size = P))+
  geom_point()+
  scale_colour_gradientn(colors = colorRampPalette(c("darkblue","yellow","red"))(99), na.value = "white") +
  scale_y_discrete(limits = rev)+
  theme_linedraw() + theme(panel.grid.major = element_blank())+
  theme(axis.text.y = element_text(color = "black", size = 10), axis.title = element_blank(),
        axis.text.x = element_text(color = "black", size = 10, angle = 90, hjust = 1, vjust = 0.5))+
  geom_vline(xintercept = c(10.5,14.5))+
  ggtitle("CCL5-CCR5")
invisible(dev.off())

ggsave(paste0(path,"CXCL12-CXCR4.pdf"), plot, width = 876/100, height = 618/100)
