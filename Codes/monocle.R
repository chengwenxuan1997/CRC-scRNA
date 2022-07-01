
# monocle 2021/12/21 ------------------------------------------------------

## "CD8+ Naïve", "GZMK+ Effector(1)", "CD103+ Trm"

library(Seurat)
library(monocle)
library(cowplot)

work.path <- "i:/genomicdata/jiaoda/pseudotime/1217";setwd(work.path)
raw.path <- "i:/genomicdata/jiaoda/raw"

# load SeuratObject, find markers
load(file.path(raw.path, "seurat_CD8T加naive(1).Rdata"))
CD8T@meta.data$celltypes4<-CD8T@meta.data$seurat_clusters
Idents(CD8T) <- "celltypes4" 
CD8T<- RenameIdents(CD8T, '4'='CD8+ Naïve','8'='CD8+ Naïve',
                    '10'='GZMK+ Effector(1)','1'='GZMK+ Effector(1)','0'='CD103+ Trm','9'='CD103+ Trm','2'='CD103+ Trm','6'='CD103+ Trm',
                    '5'='CD8+ Exhausted','11'='CD8+ Cycling','3'='GZMK+ Effector(2)','7'='Undetermined')
CD8T$celltypes4<-Idents(CD8T)
# CD8T <- subset(CD8T, celltypes4 %in% c("CD8+ Naïve", "GZMK+ Effector(1)", "CD103+ Trm", "CD8+ Exhausted"))
# CD8T <- subset(CD8T, celltypes4 %in% c("CD8+ Naïve", "GZMK+ Effector(1)", "CD103+ Trm"))
# CD8T <- subset(CD8T, celltypes4 %in% c("CD8+ Naïve", "GZMK+ Effector(1)", "CD103+ Trm", "GZMK+ Effector(2)"))
CD8T <- subset(CD8T, celltypes4 %in% c("CD8+ Naïve", "CD103+ Trm", "CD8+ Exhausted"))

seu <- CD8T;rm(CD8T)
markers <- FindAllMarkers(seu)
markers <- subset(markers, p_val_adj < 0.05)


# create monocle object
count <- readRDS(file.path(raw.path, "countCD8T.rds"))
colnames(count) <- gsub("\\.", "-", colnames(count));sum(colnames(seu)%in%colnames(count))
expr <- count[, colnames(seu)]
Barcodes <- seu@meta.data[,-c(1:3)]
Genes <- data.frame(row.names = rownames(seu), rownames(seu));colnames(Genes)="gene_short_name"
pd <- new("AnnotatedDataFrame", data = Barcodes)
fd <- new("AnnotatedDataFrame", data = Genes)
rm(count);rm(seu)

cds <- newCellDataSet(expr, phenoData = pd, featureData = fd)
cds <- detectGenes(cds)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
expressed_genes <- rownames(subset(fData(cds), num_cells_expressed >= 10))

# dim reduction
DisTable <- dispersionTable(cds)
DisTable <- subset(DisTable, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)
# ordering_genes <- unique(union(DisTable$gene_id, markers$gene))
ordering_genes <- unique(markers$gene)
cds <- setOrderingFilter(cds, ordering_genes)
plot_ordering_genes(cds)
cds <- reduceDimension(cds, max_components = 2, method = "DDRTree")

# find root state
plot_cell_trajectory(cds, color_by = "celltypes4") + facet_wrap(~celltypes4, ncol = 3)+theme(legend.position = "right") +
  # ggtitle("union of genes from Seurat::FindAllMarker and monocle::dispersionTable")
  ggtitle("genes from Seurat::FindAllMarker")
cds <- orderCells(cds)
cds <- orderCells(cds, root_state = 3)
plot_cell_trajectory(cds, color_by = "Pseudotime")
plot_cell_trajectory(cds, color_by = "State") + facet_wrap(~State, ncol = 3)

MaxMin <- function(vec, minmum = quantile(vec, 0.05), maximum = quantile(vec, 0.95)){
  vec[vec < minmum] = minmum
  vec[vec > maximum] = maximum
  (vec-minmum)/(maximum-minmum)
}
cds$Pseudotime[cds$State == 3] = MaxMin(cds$Pseudotime[cds$State == 3]) * 5
cds$Pseudotime[cds$State == 1] = MaxMin(cds$Pseudotime[cds$State == 1]) * 5 + 5
cds$Pseudotime[cds$State == 2] = MaxMin(cds$Pseudotime[cds$State == 2]) * 5 + 5
  
# saveRDS(cds, "01064celltypetmp.rds")
saveRDS(cds, "01063celltype.rds")

# plot for pseudotime
pdf(file.path(work.path, "CD8T_pseudotime_20220106.pdf"), width = 485/100, height = 575/100)
plot_grid(plot_cell_trajectory(cds, color_by = "celltypes4", cell_size = 0.5)+theme(legend.position = "right")+ 
            scale_color_manual(breaks = c("CD8+ Naïve", "GZMK+ Effector(1)", "CD103+ Trm", "CD8+ Exhausted"), 
                               values=c("#E41A1C", "#377EB8", "#4DAF4A", "#CB2EE3")),
          plot_cell_trajectory(cds, color_by = "Pseudotime", cell_size = 0.5)+theme(legend.position = "right"),
          ncol = 1)
invisible(dev.off())

# heatmap

BEAM.res <- BEAM(cds = cds, branch_point = 2)
BEAM.res <- BEAM.res[order(BEAM.res$qval),]
# plot.gene <- rownames(subset(BEAM.res, qval < 1e-4 & use_for_ordering == TRUE), )
plot.gene <- head(BEAM.res$gene_short_name, 50)
plot_genes_branched_heatmap(cds[plot.gene,],
                            branch_point = 2,
                            num_clusters = 4,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)
# test <- branchTest(cds = cds, branch_point = 2)

plot_cell_trajectory(cds, color_by = "State") + facet_wrap(~State, ncol = 3)
plot_genes_branched_pseudotime(cds[c("CD2", "NEAT1"),], branch_point = 2, color_by = "State")
plot_cell_trajectory(cds, color_by = "celltypes4")

# CD2 state5 high; NEAT1 state234 high


BEAM.res <- BEAM(cds = cds, branch_point = 1)
BEAM.res <- BEAM.res[order(BEAM.res$qval),]
# plot.gene <- rownames(subset(BEAM.res, qval < 1e-4 & use_for_ordering == TRUE), )
plot.gene <- head(BEAM.res$gene_short_name, 50)
plot.gene <- c("ITGAE", "IL2RB", "GZMB", "IFNG", "GZMK", "PDCD1", "CTLA4", "TIGIT")
pdf("pseudotime0118.pdf", width = 920/100, height = 890/100)
plot_genes_branched_heatmap(cds[plot.gene,],
                            branch_point = 1,
                            num_clusters = 5,
                            cores = 1,
                            use_gene_short_name = T,
                            branch_labels = c("state2(up)", "state1(left)"),
                            show_rownames = T)
invisible(dev.off())
plot_cell_trajectory(cds, color_by = "State") + facet_wrap(~State, ncol = 3)
plot_genes_branched_pseudotime(cds[c("CREM", "ATP5F1E"),], branch_point = 1, color_by = "State")
plot_cell_trajectory(cds, color_by = "celltypes4")

# ATP5F1E state 1; CREM state 2
plot.gene <- c("GZMB", "IFNG")
plot_genes_branched_pseudotime(cds[c("GZMB", "IFNG"),], branch_point = 1, color_by = "State")
plot_genes_branched_pseudotime(cds[c("GZMB", "IFNG"),], branch_point = 1, color_by = "celltypes4")
plot_genes_branched_pseudotime(cds[c("GZMB", "IFNG"),], branch_point = 1, color_by = "group", cell_size = 2)
plot_grid(plot_cell_trajectory(cds, color_by = "celltypes4")+ 
            scale_color_manual(breaks = c("CD8+ Naïve", "GZMK+ Effector(1)", "CD103+ Trm", "CD8+ Exhausted"), 
                               values=c("#E41A1C", "#377EB8", "#4DAF4A", "#CB2EE3")), 
          plot_cell_trajectory(cds, color_by = "State"))
plot_cell_trajectory(cds, color_by = "State") + facet_wrap(~State, ncol = 3)
plot_genes_branched_heatmap(cds[plot.gene,],
                            branch_point = 1,
                            cluster_rows = F, 
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)
plot_cell_trajectory(cds, markers = c("GZMB", "IFNG"), use_color_gradient = TRUE)

# plot(x = NULL, y = NULL, xlim = c(0, 5), ylim = c(0, 5))
# plot(function(x) 0.8*x, 0, 5, col = "red", lwd = 4, add = T)
# plot(function(x) 1*x, 0, 4, col = "green", lwd = 4, add = T)
# plot(function(x) 1*x, 4, 5, col = "green", lwd = 4, add = T, lty = 2)
# abline(a = 4, b = 0, lty =2)
# text(x = 3, y = 3.5, label = "CD8 exhaust", col = "dark green")
# text(x = 3, y = 2, label = "CD103", col = "red")
# text(x = 3, y = 4.5, label = "CD8 exhaust prediction", col = "dark green")
