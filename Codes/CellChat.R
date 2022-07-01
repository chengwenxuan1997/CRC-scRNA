# Curated CellChat Script

########################################
#############standard pipeline
########################################
cellchat <- createCellChat(object = sub_expr, meta = cell_ident, group.by = "ident")
CellChatDB <-CellChatDB.human
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
# future::plan("multiprocess", workers = 4) #unstable in rstudio
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat)
# cellchat <- filterCommunication(cellchat, min.cells = 10)
# cellchat <- computeCommunProbPathway(cellchat)

# circle plot
df.net <- subsetCommunication(cellchat)
subdata1 <- subset(df.net, source %in% group1 & target %in% group2)
subdata2 <- subset(df.net, source %in% group2 & target %in% group1)
write.csv(subdata1, paste0(path, "output/stro_CD8T.csv"))
write.csv(subdata2, paste0(path, "output/CD8T_stro.csv"))
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
mat <- cellchat@net$weight
mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
mat2["IgA+ plasma cells", ] <- mat["IgA+ plasma cells", ]
mat2[,"IgA+ plasma cells"] <- mat[,"IgA+ plasma cells"]
netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = "IgA+ plasma cells")
# hierachy plot
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = "wnt",  vertex.receiver = vertex.receiver)
# bubble plot
netVisual_bubble(cellchat, sources.use = 7, targets.use = c(1:6), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = c(1:6), targets.use = 7, remove.isolate = FALSE)
# chord diagram
netVisual_chord_gene(cellchat, sources.use = 7, targets.use = c(1:6,8:11), lab.cex = 0.5,legend.pos.y = 30)
# plotGeneExpression(cellchat, signaling = "CXCL")
# system analysis
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2


##############################
#####project
##############################

path <- "i:/genomicdata/jiaoda/cellchat/";setwd(path)

library(CellChat)
library(patchwork)
library(openxlsx)
library(Seurat)
library(cowplot)

runCellChat <- function(expr,cell_ident){
  cellchat <- createCellChat(object = expr, meta = cell_ident, group.by = "ident")
  CellChatDB <-CellChatDB.human
  CellChatDB.use <- CellChatDB
  cellchat@DB <- CellChatDB.use
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.human)
  cellchat <- computeCommunProb(cellchat)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  return(cellchat)
}
pairCellChat <- function(group1, group2){
  sub_info <- subset(cell_info, celltypes %in% c(group1, group2))
  sub_expr <- expr[,rownames(sub_info)]
  
  NC_info <- subset(sub_info, group == "NC")
  NC_expr <- expr[,rownames(NC_info)]
  NC_ident <- data.frame("ident" = as.character(NC_info$celltypes),
                         row.names = rownames(NC_info))
  NC_cellchat <- runCellChat(NC_expr, NC_ident)
  
  TM_info <- subset(sub_info, group != "NC")
  TM_expr <- expr[,rownames(TM_info)]
  TM_ident <- data.frame("ident" = as.character(TM_info$celltypes),
                         row.names = rownames(TM_info))
  TM_cellchat <- runCellChat(TM_expr, TM_ident)
  result <- list("NC" = NC_cellchat, "TM" = TM_cellchat)
  return(result)
}
allCellChat <- function(group1, group2){
  sub_info <- subset(cell_info, celltypes %in% c(group1, group2))
  sub_expr <- expr[,rownames(sub_info)]
  
  ALL_info <- sub_info
  ALL_expr <- expr[,rownames(ALL_info)]
  ALL_ident <- data.frame("ident" = as.character(ALL_info$celltypes),
                          row.names = rownames(ALL_info))
  ALL_cellchat <- runCellChat(ALL_expr, ALL_ident)
  return(ALL_cellchat)
}

expr <- readRDS("i:/genomicdata/jiaoda/raw/expr0617.rds")
seu <- CreateSeuratObject(expr)
seu <- NormalizeData(seu)
expr <- seu@assays$RNA@data
cell_info <- readRDS("i:/genomicdata/jiaoda/raw/cell_info.rds")
load("i:/genomicdata/jiaoda/raw/CD8T_metadata.Rdata")
CD8T_metadata <- subset(CD8T_metadata, celltypes5 != "Undetermined")
rownames(CD8T_metadata) <- gsub("-", "\\.", rownames(CD8T_metadata))
expr <- expr[,rownames(cell_info)]
all(rownames(cell_info)%in%colnames(expr))
all(rownames(CD8T_metadata)%in%rownames(cell_info))

cell_info <- data.frame(row.names = rownames(cell_info), lapply(cell_info, as.character), stringsAsFactors = F)
cell_info$celltypes[match(rownames(CD8T_metadata), rownames(cell_info))] <- CD8T_metadata$celltypes5
cell_info$celltypes <- gsub("TRM", "Trm", cell_info$celltypes)

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

B_mye <- pairCellChat(B,mye)
B_CD8T <- pairCellChat(B,CD8T)
# stro_mye <- pairCellChat(stro,mye)
# stro_CD8T <- pairCellChat(stro,CD8T)
# mye_panT <- pairCellChat(mye,panT)
# B_panT <- pairCellChat(B,panT)
# epi_mye <- pairCellChat(epi,mye)
# epi_panT <- pairCellChat(epi,panT)
# epi_mye_all <- allCellChat(epi,mye)
# epi_panT_all <- allCellChat(epi,panT)
stro_mye_all <- allCellChat(stro,mye)
stro_panT_all <- allCellChat(stro,panT)
stro_B_all <- allCellChat(stro,B)
saveRDS(B_mye,paste0(path,"B_mye.rds"))
saveRDS(B_CD8T,paste0(path,"B_CD8T.rds"))
# saveRDS(stro_mye,paste0(path,"stro_mye.rds"))
# saveRDS(stro_CD8T,paste0(path,"stro_CD8T.rds"))
# saveRDS(mye_panT,paste0(path,"mye_panT.rds"))
# saveRDS(B_panT,paste0(path,"B_panT.rds"))
# saveRDS(epi_mye,paste0(path,"epi_mye.rds"))
# saveRDS(epi_panT,paste0(path,"epi_panT.rds"))
saveRDS(stro_mye_all,paste0(path,"stro_mye_all.rds"))
saveRDS(stro_panT_all,paste0(path,"stro_panT_all.rds"))
saveRDS(stro_B_all,paste0(path,"stro_B_all.rds"))

############################
########post process(get data frame which describe significant interaction)
############################
postpro <- function(cellchat,path){
  name <- as.character(substitute(cellchat))
  group1 <- strsplit(name, "_")[[1]][1] %>% get()
  group2 <- strsplit(name, "_")[[1]][2] %>% get()
  df.net <- lapply(cellchat, subsetCommunication)
  result <- list()
  result[["NC_LR"]] <- subset(df.net$NC, source %in% group1 & target %in% group2)
  result[["NC_RL"]] <- subset(df.net$NC, source %in% group2 & target %in% group1)
  result[["TM_LR"]] <- subset(df.net$TM, source %in% group1 & target %in% group2)
  result[["TM_RL"]] <- subset(df.net$TM, source %in% group2 & target %in% group1)
  write.xlsx(result,paste0(path,substitute(cellchat),".xlsx"))
  return(result)
}
Allpostpro <- function(cellchat,path){
  name <- as.character(substitute(cellchat))
  group1 <- strsplit(name, "_")[[1]][1] %>% get()
  group2 <- strsplit(name, "_")[[1]][2] %>% get()
  df.net <- subsetCommunication(cellchat)
  result <- list()
  result[["LR"]] <- subset(df.net, source %in% group1 & target %in% group2)
  result[["RL"]] <- subset(df.net, source %in% group2 & target %in% group1)
  write.xlsx(result,paste0(path,substitute(cellchat),".xlsx"))
  return(result)
}

############################
########dot plot
############################

############### 
GetInteractionScore <- function(interlist, group1, group2, group){
  # interaction <- unique(as.character(df.net$interaction_name))
  # gene <- lapply(interaction, FUN = function(x){unlist(strsplit(x,"_"))})[[1]]
  # source_expr <- expr[gene[1],cell_info$celltypes == group1 & cell_info$group == "NC"]
  # target_expr <- expr[gene[-1],cell_info$celltypes == group2[1] & cell_info$group == "NC"]
  # value <-  (mean(source_expr)+mean(target_expr)*(length(gene)-1))/length(gene)
  # score <- sapply(group2, FUN = function(x){
  #   source_expr <- expr[gene[1],cell_info$celltypes == group1 & cell_info$group == "NC"]
  #   target_expr <- expr[gene[-1],cell_info$celltypes == x & cell_info$group == "NC"]
  #   value <- mean(c(as.vector(source_expr),as.vector(target_expr)))
  #   value <- (mean(source_expr)+mean(target_expr)*(length(gene)-1))/length(gene)
  # })
  gene <- lapply(interlist, FUN = function(x){unlist(strsplit(x,"_"))})
  score <- lapply(gene, FUN = function(g){
    sapply(group1, FUN = function(x){
      sapply(group2, FUN = function(y){
        source_expr <- expr[g[1],cell_info$celltypes == x & cell_info$group %in% group]
        target_expr <- expr[g[-1],cell_info$celltypes == y & cell_info$group %in% group]
        # value <- mean(c(as.vector(source_expr),as.vector(target_expr)))
        value <- (mean(source_expr)+mean(target_expr)*(length(g)-1))/length(g)
      })
    })
  })
  score <- do.call(cbind,score)
  colnames(score) =  interlist
  return(score)
}
turn <- function(mtx, name){
  record <- lapply(rownames(mtx), FUN = function(x){
    lapply(colnames(mtx), FUN = function(y){
      c(x,y,mtx[x,y])
    })
  })
  record <- lapply(record, FUN = function(x){do.call(rbind,x)})
  record <- do.call(rbind, record)
  record <- as.data.frame(record)
  colnames(record) <- name
  return(record)
}
FetchPvalue <- function(data,array){
  index <- dimnames(array)
  Pvalue <- array[match("IgA+ plasma cells", index[[1]]),
                  match(data$target, index[[2]]),
                  match(data$interaction_name1, index[[3]])]
  return(Pvalue)
}

## figure1
interlist <- c("MIF_CD74_CXCR4","MDK_ITGA4_ITGB1","MDK_NCL",
               "GAS6_MERTK","CLEC2B_KLRB1","CD69_KLRB1",
               "ICAM2_ITGAL_ITGB2","NAMPT_ITGA5_ITGB1",
               "SEMA4D_PLXNB2","SEMA4A_PLXNB2")

DataToPlot <- function(group,  division, xorder = NULL){
  score <- GetInteractionScore(interlist, "IgA+ plasma cells", c(CD8T,mye), group)
  if (is.null(xorder)){xorder <- rownames(score)}
  yorder <- interlist
  plotdata <- turn(score, c("target","interaction_name","score"))
  plotdata$interaction_name1 <- plotdata$interaction_name
  plotdata$interaction_name1[plotdata$interaction_name1=="CD69_KLRB1"] = "CLEC2C_KLRB1"
  plotdata$target <- factor(plotdata$target, levels = xorder)
  plotdata$interaction_name <- factor(plotdata$interaction_name, levels = yorder)
  plotdata$score <- as.numeric(plotdata$score)
  plotdata$Pvalue1 <- sapply(rownames(plotdata), FUN = function(x){
    FetchPvalue(plotdata[x,],B_mye[[division]]@net$pval)
  })
  plotdata$Pvalue2 <- sapply(rownames(plotdata), FUN = function(x){
    FetchPvalue(plotdata[x,],B_CD8T[[division]]@net$pval)
  })
  plotdata$Pvalue <- paste0(plotdata$Pvalue1,plotdata$Pvalue2)
  plotdata$Pvalue <- as.numeric(gsub("NA","",plotdata$Pvalue))
  plotdata$P <- cut(plotdata$Pvalue,c(1,0.05,0.01,0),include.lowest = T)
  plotdata$P <- factor(plotdata$P,levels = c("(0.05,1]","(0.01,0.05]","[0,0.01]"))
  plotdata$score[plotdata$score>5]=5
  plotdata <- plotdata[!is.na(plotdata$P),]
}
plotorder <- c("CD8+ Naïve", "CD103+ TRM", "ITGB2+ TRM", "non-TRM CD8+ Activated", "CD8+ Exhausted", "CD8+ Cycling",
               mye)
TMdata <- DataToPlot(c("HP","SSL","SSLD","TSA"), "TM", plotorder)
NCdata <- DataToPlot("NC", "NC", plotorder)
TMdata$target <- plyr::mapvalues(TMdata$target,
                                 from = c("CD103+ TRM", "ITGB2+ TRM", "non-TRM CD8+ Activated"),
                                 to = c("CD103+ Trm", "ITGB2+ Trm", "non-Trm CD8+ Activated"))
NCdata$target <- plyr::mapvalues(NCdata$target,
                                 from = c("CD103+ TRM", "ITGB2+ TRM", "non-TRM CD8+ Activated"),
                                 to = c("CD103+ Trm", "ITGB2+ Trm", "non-Trm CD8+ Activated"))

p1 <- ggplot(NCdata, aes(x = target, y = interaction_name, col = score, size = P)) +
  ggtitle("normal") +
  geom_point()+
  scale_colour_gradientn(colors = colorRampPalette(c("darkblue","yellow","red"))(101), limits = c(0, 5)) +
  scale_y_discrete(limits = rev)+
  theme_linedraw() + theme(panel.grid.major = element_blank())+
  theme(axis.text.y = element_text(color = "#70B7E6"), axis.title = element_blank(),
        axis.text.x = element_text(color = "#70B7E6", angle = 90, hjust = 1))+
  geom_vline(xintercept = 6.5)

p2 <- ggplot(TMdata, aes(x = target, y = interaction_name, col = score, size = P)) +
  ggtitle("serrated lesions") +
  geom_point()+
  scale_colour_gradientn(colors = colorRampPalette(c("darkblue","yellow","red"))(99), na.value = "white", limits = c(0, 5)) +
  scale_y_discrete(limits = rev)+
  theme_linedraw() + theme(panel.grid.major = element_blank())+
  theme(axis.text.y = element_text(color = "#70B7E6"), axis.title = element_blank(),
        axis.text.x = element_text(color = "#70B7E6", angle = 90, hjust = 1))+
  geom_vline(xintercept = 6.5)

legend <- get_legend(p1)

plot_grid(p1 + theme(legend.position = "none"), 
          legend, 
          p2 + theme(legend.position = "none"), 
          ncol = 3, rel_widths = c(4,1,4))
ggsave("lgA.pdf", width = 876/100, height = 580/100)

#########################
##### dotplot for one pair in multiple celltypes
#########################
ArrayToDataFrame <- function(array,pair = NULL){
  source = dimnames(array)[[1]]
  target = dimnames(array)[[2]]
  if (is.null(pair)) pair = dimnames(array)[[3]]
  data = lapply(source,function(i){
    lapply(target,function(j){
      t(as.matrix(sapply(pair,function(k){
        c(i, j, k, as.numeric(array[i,j,k]))#source, target, pair, value
      })))
    })
  })
  data = lapply(data,function(x) do.call(rbind,x))
  data = do.call(rbind,data)
  data = as.data.frame(data)
  rownames(data) = NULL
  colnames(data) = c("source", "target", "pair", "value")
  return(data)
}
ExtractCellChat <- function(cellchat,pair){
  prob <- cellchat@net$prob
  prob <- ArrayToDataFrame(prob,pair)
  pval <- cellchat@net$pval
  pval <- ArrayToDataFrame(pval,pair)
  if (all(prob[,1:3]==pval[,1:3])){
    res = data.frame(
      "source" = prob$source,
      "target" = prob$target,
      "pair" = prob$pair,
      "prob" = prob$value,
      "pval" = pval$value
    )
  }else{print("index not match")}
}
FetchInteractionScoreForDataFrame <- function(data){
  data = split(data,f = paste0(data$source, data$target, data$pair))
  data = lapply(data, function(x){
    source = as.character(x[1]); target = as.character(x[2])
    ligand = strsplit(as.character(x[3]), "_")[[1]][1]
    receptor = strsplit(as.character(x[3]), "_")[[1]][2]
    source.expr = expr[ligand,cell_info$celltypes%in%source]
    target.expr = expr[receptor,cell_info$celltypes%in%target]
    score = (mean(source.expr)+mean(target.expr))/2
    x["score"] = score
    return(x)
  })
  data = as.data.frame(do.call(rbind,data))
  rownames(data) = NULL
  colnames(data) = c("source", "target", "pair", "prob", "pval", "score")
  return(data)
}

project <- list("stro_panT_all","stro_B_all","stro_mye_all")
cellchatdata <- lapply(project, function(x){
  cellchat = readRDS(paste0(path,x,".rds"))
  data = ExtractCellChat(cellchat,"CXCL12_CXCR4")
  # data = ExtractCellChat(cellchat,"CCL5_CCR1")
  data = subset(data, source%in%stro & target%in%c(panT,B,mye))
})
cellchatdata <- do.call(rbind,cellchatdata)
cellchatdata <- FetchInteractionScoreForDataFrame(cellchatdata)

plotdata <- cellchatdata
plotdata$source <- factor(plotdata$source,
                          levels = c("PDGRFA+ fibroblasts", "Stromal-1", "Stromal-2",
                                     "Stromal-3", "Myofibroblasts", "Stalk-like ECs", 
                                     "Tip-like ECs", "EC-1", "EC-2", "EPCs", "Vascular SMCs",
                                     "Pericytes", "Enteric glial cells"))
# plotdata$target <- factor(plotdata$target, levels = c(panT, B, "FCN1+ inflammatory macrophages","Anti-inflammatory macrophages","Conventional DCs","Plasmacytoid DCs"))
plotdata$P <- cut(as.numeric(plotdata$pval),c(1,0.05,0.01,0),include.lowest = T)
plotdata$P <- factor(plotdata$P,levels = c("(0.05,1]","(0.01,0.05]","[0,0.01]"))
plotdata$score[plotdata$score>5]=5

CXCL12 = plotdata
CXCL12 = subset(CXCL12, !target=="CD8+ Naïve")
CXCL12$target <- factor(CXCL12$target, 
                        levels = c("CD8+ Naïve", "CD103+ Trm", "ITGB2+ Trm", "non-Trm CD8+ Activated", "CD8+ Exhausted", "CD8+ Cycling", "CD4+ T cells",
                                   "Naïve T cells", "Regulatory T cells","γδ T cells", "NK cells",
                                   B, c("FCN1+ inflammatory macrophages","Anti-inflammatory macrophages","Conventional DCs","Plasmacytoid DCs"))) 
grDevices::cairo_pdf(paste0(path,"CXCL12-CXCR40926.pdf"), width = 798/100, height = 490/100)
ggplot(CXCL12, aes(x = target, y = source, col = score, size = P))+
  geom_point()+
  scale_colour_gradientn(colors = colorRampPalette(c("darkblue","yellow","red"))(99), na.value = "white") +
  scale_y_discrete(limits = rev)+
  theme_linedraw() + theme(panel.grid.major = element_blank())+
  theme(axis.text.y = element_text(color = "black", size = 10), axis.title = element_blank(),
        axis.text.x = element_text(color = "black", size = 10, angle = 45, hjust = 1, vjust = 1))+
  geom_vline(xintercept = c(10.5,14.5))+
  ggtitle("CXCL12-CXCR4")
invisible(dev.off())

# CCL5=plotdata
# CCL5 = subset(CCL5, target%in%mye)
# CCL5$target <- factor(CCL5$target, 
#                       levels = c(panT, B, c("FCN1+ inflammatory macrophages","Anti-inflammatory macrophages","Conventional DCs","Plasmacytoid DCs"))) 
# grDevices::cairo_pdf(paste0(path,"CCL5-CCR1.pdf"), width = 379/100, height = 622/100)
# ggplot(CCL5, aes(x = target, y = source, col = score, size = P))+
#   geom_point()+
#   scale_colour_gradientn(colors = colorRampPalette(c("darkblue","yellow","red"))(99), na.value = "white") +
#   scale_y_discrete(limits = rev)+
#   theme_linedraw() + theme(panel.grid.major = element_blank())+
#   theme(axis.text.y = element_text(color = "black", size = 10), axis.title = element_blank(),
#         axis.text.x = element_text(color = "black", size = 10, angle = 45, hjust = 1, vjust = 1))+
#   ggtitle("CCL5-CCR1")
# invisible(dev.off())