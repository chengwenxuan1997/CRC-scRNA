library(Seurat)
library(scMetabolism)
library(ComplexHeatmap)

# computation -------------------------------------------------------------


path <- "i:/genomicdata/jiaoda/scMetabolism";setwd(path)

count <- readRDS("~/data/cwx/jiaoda/raw/countCD8T.rds")
KEGG<-sc.metabolism(countexp = as.matrix(count),
                    method = "VISION", 
                    imputation = F, 
                    ncores = 10, 
                    metabolism.type = "KEGG")
REACTOME<-sc.metabolism(countexp = as.matrix(count),
                        method = "VISION", 
                        imputation = F, 
                        ncores = 10, 
                        metabolism.type = "REACTOME")

CD8Tmeta <- readRDS("I:/genomicdata/jiaoda/scMetabolism/CD8Tmeta.rds")
CD8Tres <- readRDS("I:/genomicdata/jiaoda/scMetabolism/CD8Tres.rds")
CD8Tmeta <- subset(CD8Tmeta, celltypes3 != "Undetermined")
CD8Tmeta$celltype <- CD8Tmeta$celltypes3
cell_info <- CD8Tmeta  

avg.scMetabolism <- function(res, cell_info){
  celltypes <- names(table(cell_info$celltype))[table(cell_info$celltype)>0]
  
  CtMat <- lapply(celltypes, function(x){
    as.numeric(cell_info$celltype == x)
  })
  CtMat <- do.call(cbind, CtMat)
  rownames(CtMat) = rownames(cell_info)
  colnames(CtMat) = celltypes
  rownames(CtMat) = gsub("-", "\\.", rownames(CtMat))
  if (!all(rownames(CtMat)==colnames(res$KEGG)) ){print("NAME ERROR")}
  if (!all(table(cell_info$celltype)[table(cell_info$celltype)>0] == colSums(CtMat))){print("NUM ERROR")}
  res.avg <- lapply(res, function(x){
    sum = as.matrix(x) %*% CtMat
    avg = sum
    for (i in 1:length(celltypes)){
      avg[,i] = sum[,i]/colSums(CtMat)[i]
    }
    rownames(avg) = rownames(x)
    colnames(avg) = celltypes
    return(avg)
  })
  return(res.avg)
}


# plot --------------------------------------------------------------------


pathway <- c("Retinol metabolism", "Histidine metabolism", "Primary bile acid biosynthesis", "Glycolysis / Gluconeogenesis",
             "Citrate cycle (TCA cycle)", "Vitamin B2 riboflavin metabolism", "Nicotinate metabolism",
             "Biotin transport and metabolism", "Metabolism of porphyrins", "Arachidonic acid metabolism")
# celltypes <- c("CD8+ Naïve", "CD8+ Effecor-1", "CD8+ Effecor-2", "CD8+ Exhausted", "CD8+ Cycling", "CD8+ NKT")
celltypes <- c("CD8+ Effecor-2", "CD8+ Exhausted", "CD8+ Effecor-1", "CD8+ NKT")# , "CD8+ Naïve", "CD8+ Cycling"
rownames(plotdata)[grep("Sialic", rownames(plotdata))]
res.avg <- avg.scMetabolism(res = CD8Tres, cell_info = CD8Tmeta)
plotdata <- rbind(res.avg$KEGG, res.avg$REACTOME)
all(pathway%in%rownames(plotdata))
all(celltypes%in%colnames(plotdata))
plotdata <- plotdata[pathway, celltypes]
range(plotdata)
plot(density(plotdata))
plotdata <- t(scale(t(plotdata)))
plotdata <- plotdata[!is.na(plotdata[,1]),]
plot(density(plotdata))
plotdata[plotdata<(-2)]=-2
plotdata[plotdata>2]=2
scales::show_col(colors)
pdf("CD8TscMetabolism0112.pdf", width = 600/100, height = 575/100)
pheatmap(mat = plotdata, 
         cluster_rows = F, cluster_cols = F,
         border_color = "white", 
         color = colorRampPalette(colors = c("#77A4D0","white","red"))(100),
         cellwidth = 24, cellheight = 12,
         angle_col = "45",
         fontsize_col = 12, fontsize_row = 10,
         fontsize = 15,
         gaps_row = 5)
invisible(dev.off())
