
# cellcycle ---------------------------------------------------------------

path <- "i:/genomicdata/jiaoda/cellcycle";setwd(path)
rawData.path <- "i:/genomicdata/jiaoda/raw"
seurat <- CellCycleScoring(seurat, 
                           s.features = cc.genes$s.genes, 
                           g2m.features = cc.genes$g2m.genes, 
                           set.ident = TRUE)
panT_info <- FetchData(seurat, vars = c("celltypes","group","Phase"))

load(file.path(rawData.path, "seurat_CD8T加naive(1).Rdata"))
CD8T@meta.data$celltypes4<-CD8T@meta.data$seurat_clusters
Idents(CD8T) <- "celltypes4" 
CD8T<- RenameIdents(CD8T, '4'='CD8+ Naïve','8'='CD8+ Naïve',
                    '10'='GZMK+ Effector(1)','1'='GZMK+ Effector(1)','0'='CD103+ Trm','9'='CD103+ Trm','2'='CD103+ Trm','6'='CD103+ Trm',
                    '5'='CD8+ Exhausted','11'='CD8+ Cycling','3'='GZMK+ Effector(2)','7'='Undetermined')
CD8T$celltypes4<-Idents(CD8T)
CD8T <- CellCycleScoring(CD8T, 
                         s.features = cc.genes$s.genes, 
                         g2m.features = cc.genes$g2m.genes, 
                         set.ident = TRUE)

CD8T_info <- FetchData(CD8T, vars = c("celltypes4","group","Phase"))
colnames(CD8T_info) <- c("celltypes","group","Phase")
panT_info <- panT_info[!rownames(panT_info)%in%rownames(CD8T_info),]

data <- rbind(CD8T_info, panT_info)
data$group <- as.character(data$group)
data$group[data$group != "NC"] <- "TM"
celltype <- names(table(data$celltypes)[table(data$celltypes)>0])
celltype <- c("CD8+ Naïve", "GZMK+ Effector(1)", "CD103+ Trm", "CD8+ Exhausted", "CD8+ Cycling", "GZMK+ Effector(2)")

plotlist <- lapply(celltype, function(x){
  plotdata <- subset(data, celltypes == x)
  plotdata <- as.data.frame(table(plotdata$group, plotdata$Phase))
  colnames(plotdata) <- c("Group", "Phase", "Num")
  plotdata$Phase <- factor(as.character(plotdata$Phase), levels = c("G1", "S", "G2M"))
  plotdata$title <- x
  plotdata$sum[plotdata$Group == "TM"] <- sum(plotdata$Num[plotdata$Group == "TM"])
  plotdata$sum[plotdata$Group == "NC"] <- sum(plotdata$Num[plotdata$Group == "NC"])
  plotdata$pct <- round(plotdata$Num/plotdata$sum*100, 1)
  plotdata$color[plotdata$Phase %in% c("G1","G2M")] <- "white"
  plotdata$color[plotdata$Phase %in% "S"] <- "black"
  
  ggplot(data = plotdata, aes(x = Group, y = pct, fill = Phase, label = paste0(pct, "%"), color = color)) +
    geom_bar(position = "stack", stat = "identity", alpha = 0.8, color = NA)  +
    scale_fill_manual(values=c("G1" = "#A9030D", "S" = "#EBEBED", "G2M" = "#003893")) +
    scale_color_manual(values = c("black" = "black", "white" = "white")) +
    geom_text(size = 3.5, position = position_stack(vjust = 0.5)) +
    facet_grid(~title) +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
          panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5, colour = "red", vjust = -10),
          strip.background = element_rect(fill="grey", color = "black", size = 1),
          strip.text = element_text(size = 9, colour="black"),
          axis.title = element_blank())
})
legend <- get_legend(plotlist[[1]] + theme(legend.position="right"))
plotlist <- lapply(plotlist, function(x){x + theme(legend.position = "none")})
plotlist_new <- lapply(plotlist, function(x){x + theme(axis.text.y = element_blank())})
plotlist_new[c(1)] = plotlist[c(1)]
raw_plot <- plot_grid(plotlist = plotlist_new, ncol = 6, nrow = 1, axis = l, rel_widths = c(3.5,3,3))
plot <- plot_grid(raw_plot, legend, nrow = 1, rel_widths = c(8,1))
grDevices::cairo_pdf(paste0(path, "/cellcycle1213.pdf"), width = 826/100, height = 231/100)
plot
invisible(dev.off())
