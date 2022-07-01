library(Seurat)
library(infercnv)
library(ggplot2)
library(ggpubr)
library(circlize)


# Reference file preparation ----------------------------------------------

path <- "~/genomicdata/jiaoda/inferCNV/"
expr <- readRDS("~/genomicdata/jiaoda/expr0617.rds")
auc.seu <- readRDS("~/genomicdata/jiaoda/SCENIC/auc.seu.rds")
cell_info<-FetchData(auc.seu, vars = c("group","patient","type","division","celltypes"))
cell_info$group[cell_info$group=="SSLD"] <- "SSL"
cell_info$GP[cell_info$group=="NC"]="NC"
cell_info$GP[cell_info$group!="NC"]="TM"
cell_info$ident <- paste0(cell_info$patient,"_",cell_info$celltypes)
expr <- expr[,rownames(cell_info)]
if (!file.exists(paste0(path,"gencode_v21_gen_pos.complete.txt"))){
  download.file("https://data.broadinstitute.org/Trinity/CTAT/cnv/gencode_v21_gen_pos.complete.txt",
                paste0(path,"gencode_v21_gen_pos.complete.txt"))
}
if (!file.exists(paste0(path,"gencode_v19_gene_pos.txt"))){
  download.file("https://data.broadinstitute.org/Trinity/CTAT/cnv/gencode_v19_gene_pos.txt",
                paste0("gencode_v19_gene_pos.txt"))
}
if (!file.exists(paste0(path,"gencode.v38.annotation.gtf.gz"))){
  download.file("http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz",
                paste0(path,"gencode.v38.annotation.gtf.gz"))
}
gene<-read.table("~/genomicdata/jiaoda/genes.tsv", sep = "\t", header = T)
rownames(expr)<-gene$gene_ids[match(rownames(expr),gene$gene)]
gene_order_file21 <- read.table(paste0(path,"gencode_v21_gen_pos.complete.txt"), 
                              sep = "\t", header = F,
                              row.names = 1,
                              check.names = F)
gene_order_file19 <- read.table(paste0(path,"gencode_v19_gene_pos.txt"), 
                                sep = "\t", header = F,
                                row.names = 1,
                                check.names = F)

##########GRCh38
# download.file("https://raw.githubusercontent.com/broadinstitute/infercnv/master/scripts/gtf_to_position_file.py","gtf_to_position_file.py")
gene_order_file38 <- read.table(paste0(path,"gencode.v38.txt"),
                                sep = "\t", row.names = 1,
                                header = F, check.names = F)
gene38<-strsplit(rownames(gene_order_file38),"\\.")
gene38<-do.call(rbind,gene38)[,1]
gene_order_file <- gene_order_file38[!duplicated(gene38),]
rownames(gene_order_file) <-gene38[!duplicated(gene38)]
# compare<-data.frame("GRCH21"=gene$gene_ids[match(rownames(expr),gene$gene)],
#                     "GRCH38"=gene$gene_ids[match(rownames(expr),gene$gene)])
# download.file("http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.chr_patch_hapl_scaff.annotation.gtf.gz",
#               paste0(path,"gencode.v38.chr_patch_hapl_scaff.annotation.gtf.gz"))



# Run inferCNV ------------------------------------------------------------


# epi_info <- cell_info[grep("Epi-1|Epi-2|Epi-3",cell_info$ident),]
# epi_info <- cell_info[cell_info$celltypes%in%c("Epi-1","Epi-2","Epi-3"),]
epi_info <- cell_info[cell_info$celltypes%in%c("Epi-1","Epi-2","Epi-3","Mature colonocytes"),]
# epi_info <- cell_info[grep("Epi-2",cell_info$ident),]
# epi_info <- cell_info[cell_info$type=="epithelic",]
# a<-epi_info$ident[substr(epi_info$ident,1,2)!="NC"]
# a<-do.call(rbind,strsplit(a,"_"))
# epi_info$ident[substr(epi_info$ident,1,2)!="NC"] <- a[,1]
# epi_info$ident <- gsub("^NC[0-9]_","NC_",epi_info$ident)
epi_info$ident <- paste0(epi_info$GP,"_",epi_info$celltypes)

epi_expr <- expr[,rownames(epi_info)]

epi_info1 <- as.data.frame(as.character(epi_info$ident))
rownames(epi_info1) <- rownames(epi_info)
epi_cnv <- CreateInfercnvObject(raw_counts_matrix = epi_expr,
                                     annotations_file = epi_info1,
                                     gene_order_file = gene_order_file,
                                     ref_group_names = names(table(epi_info1))[substr(names(table(epi_info1)),0,2)=="NC"])
                                     # ref_group_names = NULL)
epi_cnv <- infercnv::run(infercnv_obj = epi_cnv,
                              cutoff = 0.1,
                              out_dir = paste0(path,"pan_epi_by_celltype_nosub"),
                              cluster_by_groups = T,
                              denoise = T,
                              # noise_filter = 0.2,
                              # sd_amplifier = 2,
                              HMM = T,
                              # analysis_mode = "subclusters",
                              analysis_mode = "sample",
                              num_threads = 20)


# Do plot -----------------------------------------------------------------



final_infercnv_obj = readRDS(paste0(path,"epioutput1/run.final.infercnv_obj"))
compare_obj = readRDS(paste0(path, "epioutput1/19_HMM_pred.repr_intensitiesHMMi6.hmm_mode-samples.Pnorm_0.5.infercnv_obj"))
new_gene_order = data.frame()
for (chr_name in c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22")) {
  new_gene_order = rbind(new_gene_order, infercnv_obj@gene_order[which(infercnv_obj@gene_order[["chr"]] == chr_name) , , drop=FALSE])
}
names(new_gene_order) <- c("chr", "start", "stop")
copy_infercnv_obj@gene_order = new_gene_order
copy_infercnv_obj@expr.data = infercnv_obj@expr.data[rownames(new_gene_order), , drop=FALSE]



## Boxplot -----------------------------------------------------------------

# cnv <- readRDS(paste0(path,"pan_epi/run.final.infercnv_obj"))
cnv <- readRDS(paste0(path,"pan_epi_by_celltype/run.final.infercnv_obj"))
cnv_state <- readRDS(paste0(path,"pan_epi/17_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.infercnv_obj"))
data <- cnv@expr.data-1
data <- apply(data,2,function(x){sum(x^2)})
sum(names(data)==rownames(epi_info))
plotdata <- data.frame("CNV_score" = data,
                       "celltype" = epi_info$celltypes,
                       "group" = as.factor(as.character(epi_info$group)),
                       "CNV_signal" = data/dim(cnv@expr.data)[1])
compare_type <- list(c("Epi-2","Mature colonocytes"))
compare_group <- list(c("NC","SSL"))
# pdf(paste0(path,"pan_epi_by_celltype/boxplot.pdf"))
ggplot(data = plotdata, mapping = aes(x = celltype, y = CNV_score, fill = celltype))+
  geom_boxplot(width=0.4)+
  theme(legend.position = "top")+
  theme_classic()+
  scale_color_manual(values = c("#0A6BAD","#33A02C","#E31A1C","#B1DE89"),aesthetics = "fill")+
  stat_compare_means(comparisons = compare_type)

ggplot(data = plotdata, mapping = aes(x = group, y = CNV_signal, fill = group))+
  geom_boxplot(width=0.4)+
  theme(legend.position = "top")+
  theme_classic()+
  stat_compare_means(comparisons = compare_group)
  # scale_color_manual(values = c("#0A6BAD","#33A02C","#E31A1C","#B1DE89"),aesthetics = "fill")
# invisible(dev.off())


## boxviolin plot -----------------------------------------------------------------

cell.col <- c("Epi-2" = "#0A6BAD",
              "Epi-1" = "#33A02C",
              "Epi-3" = "#E10C0D",
              "Mature colonocytes" = "#B2DF8A")
dat <- read.csv("~/genomicdata/jiaoda/inferCNV/inferCNV_plotdata.csv",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
dat$celltype <- factor(dat$celltype,levels = c("Epi-2","Epi-1","Epi-3","Mature colonocytes"))
gg<-ggplot(data = dat,aes(x = celltype, #分组列名
                      y = log2(CNV_score + 1) , #连续变量列名
                      fill = celltype))+ #按分组填充颜色
  scale_fill_manual(values = cell.col) + #用自定义颜色填充
  geom_violin(alpha = 0.4, position = position_dodge(width = .75), 
              size = 0.8, color="white") +
  
  geom_point(shape = 21, size=1.2, # 点的性状和大小
             position = position_jitterdodge(), # 让点散开
             aes(color=celltype), alpha = 0.6) +
  scale_color_manual(values = cell.col) + #用自定义颜色填充
  geom_boxplot(notch = TRUE, outlier.size = -1, 
               color="black", lwd = 0.6, alpha = 0.7) +
  theme_bw() + 
  ylab("log2(CNV score + 1)") +
  xlab("") + ylim(0,5.5) +
  theme(axis.text.x = element_text(hjust = 1, size = 10, angle = 45,color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.position = "top",
        legend.title = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10)) +
  stat_compare_means(method = "kruskal.test", label.y = 5.3, label.x = 1.5)
gg+guides(fill=guide_legend(nrow=2,byrow=TRUE))
ggsave("~/genomicdata/jiaoda/inferCNV/cnv score of epi cells.pdf",width = 3,height = 5)

# Heatmap -----------------------------------------------------------------
library(ComplexHeatmap)

cnv <- readRDS(paste0(path,"pan_epi_by_celltype/19_HMM_pred.repr_intensitiesHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.infercnv_obj"))
gene_annotation <- cnv@gene_order$chr
names(gene_annotation) <- rownames(cnv@gene_order)
genecol <- colorRampPalette(RColorBrewer::brewer.pal(12,"Set3"))(22)
names(genecol) <- names(table(gene_annotation))
anncol <- HeatmapAnnotation(gene = gene_annotation,
                            col = list(gene = genecol))
ha <- HeatmapAnnotation(gene = anno_block(gp = gpar(fill = colorRampPalette(RColorBrewer::brewer.pal(12,"Set3"))(22)),
                                          labels = paste0("chr",seq(1:22))))

ref <- c("NC_Mature colonocytes","NC_Epi-1","NC_Epi-3","NC_Epi-2")
ref_plotdata <- lapply(ref, function(x){
  expr = cnv@expr.data[,cnv@reference_grouped_cell_indices[[x]]]
  x = expr[,order.dendrogram(as.dendrogram(cnv@tumor_subclusters$hc[[x]]))]
})
ref_plotdata <- do.call(cbind,ref_plotdata)
# ref_annotation <- epi_info1[colnames(ref_plotdata),];names(ref_annotation) = colnames(ref_annotation)
ref_annotation <- factor(epi_info1[colnames(ref_plotdata),], levels = unique(epi_info1[colnames(ref_plotdata),]))
ref_annrow <- rowAnnotation(celltype = ref_annotation,
                          col = list(
                            celltype = c(
                              "NC_Epi-2" = "#0A6BAD",
                              "NC_Epi-1" = "#33A02C",
                              "NC_Epi-3" = "#E10C0D",
                              "NC_Mature colonocytes" = "#B2DF8A"
                            )
                          ))
tp <- Heatmap(t(ref_plotdata),
              width = ncol(t(ref_plotdata))*unit(0.04, "mm"),
              height = nrow(t(ref_plotdata))*unit(0.03, "mm"),
              col = colorRamp2(c(0,1,3),c("darkblue","white","darkred")),
              row_split = ref_annotation, row_gap = unit(0, "mm"), row_title = "NC",
              column_split = gene_annotation, column_gap = unit(0, "mm"), column_title = NULL,
              border = T, 
              cluster_rows = F, cluster_columns = F,
              show_row_names = F, show_column_names = F,
              left_annotation = ref_annrow
              # bottom_annotation =  ha
)

obs <- c("TM_Mature colonocytes","TM_Epi-1","TM_Epi-3","TM_Epi-2")
obs_plotdata <- lapply(obs, function(x){
  expr = cnv@expr.data[,cnv@observation_grouped_cell_indices[[x]]]
  x = expr[,order.dendrogram(as.dendrogram(cnv@tumor_subclusters$hc[[x]]))]
})
obs_plotdata <- do.call(cbind,obs_plotdata)
# obs_annotation <- epi_info1[colnames(obs_plotdata),];names(obs_annotation) = colnames(obs_annotation)
obs_annotation <- factor(epi_info1[colnames(obs_plotdata),], levels = unique(epi_info1[colnames(obs_plotdata),]))
obs_annrow <- rowAnnotation(celltype = obs_annotation,
                            col = list(
                              celltype = c(
                                "TM_Epi-2" = "#0A6BAD",
                                "TM_Epi-1" = "#33A02C",
                                "TM_Epi-3" = "#E10C0D",
                                "TM_Mature colonocytes" = "#B2DF8A"
                              )
                            ))
bp <- Heatmap(t(obs_plotdata),
              width = ncol(t(obs_plotdata))*unit(0.04, "mm"),
              height = nrow(t(obs_plotdata))*unit(0.03, "mm"),
              col = colorRamp2(c(0,1,3),c("darkblue","white","darkred")),
              row_split = obs_annotation, row_gap = unit(0, "mm"), row_title = "TM",
              column_split = gene_annotation, column_gap = unit(0, "mm"), column_title = NULL,
              border = T, 
              cluster_rows = F, cluster_columns = F,
              show_row_names = F, show_column_names = F,
              left_annotation = obs_annrow,
              top_annotation =  ha
)
pdf(paste0(path,"test2.pdf"), width = 38*0.39, height = 26*0.39)
draw(tp%v%bp, column_title = "inferCNV")
invisible(dev.off())

