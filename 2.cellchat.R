library(CellChat)
library(patchwork)
source('scRNA_scripts/lib.R')
setwd("./24.5")

sce.all@meta.data$group = "NA"
library(readr)
group <- read_csv("../sample.csv")


for(i in 1:nrow(group)){
  sce.all@meta.data[which(sce.all@meta.data$orig.ident == group$sample[i]),
                    'group'] <- group$state[i]
}
Idents(sce.all)=sce.all$group
future::plan("multisession", workers = 12)

data.dir <- './comparison'
setwd(data.dir)

sce.all.Pre_R = sce.all[,sce.all@meta.data$group %in% c('Pre_R')]
sce.all.Pre_NR = sce.all[,sce.all@meta.data$group %in% c('Pre_NR')]
sce.all.Post_NR = sce.all[,sce.all@meta.data$group %in% c('Post_NR')]
sce.all.Post_R = sce.all[,sce.all@meta.data$group %in% c('Post_R')]


sub_list <- list(sce.all.Pre_NR,sce.all.Post_NR,sce.all.Pre_R,sce.all.Post_R)
names(sub_list) <- c('sce.all.Pre_NR','sce.all.Post_NR','sce.all.Pre_R','sce.all.Post_R')

cellchat_list <- lapply(sub_list, function(sce) {
  cellchat <- createCellChat(sce, group.by='celltype')
  cellchat
  levels( cellchat@idents )
  
  #设置受配体互作数据库####
  ## A list includes the ligand-receptor interactions
  CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
  showDatabaseCategory(CellChatDB)
  dplyr::glimpse(CellChatDB$interaction)
  
  unique(CellChatDB$interaction$pathway_name)
  
  # CellChatDB.use <- CellChatDB
  ## key : the name of the variable in CellChatDB interaction_input
  ## CellChatDB.ss <- subsetDB(CellChatDB, search = "Secreted Signaling", key='annotation')
  
  # set the used database in the cellchatect
  cellchat@DB <- CellChatDB
  
  #预处理表达数据####
  #数据处理过程包括： - 鉴定一个细胞类型中的过表达的配体或受体 - 鉴定过表达的配受体互作
  #基于PPI网络平滑基因表达（可选：适用于处理低深度单细胞测序处理的dropout）
  # subset表达数据，提取仅在互作数据库中的基因，减少下游分析数据量
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
  # future::plan("multiprocess", workers = 4)
  # 识别在单个细胞类型中过表达配/受体
  cellchat <- identifyOverExpressedGenes(cellchat)
  
  # 识别过表达互作对
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  # 平滑表达值（目的是消除dropout影响，可选不用）
  # We also provide a function to project gene expression data onto protein-protein interaction (PPI) network. Specifically, a diffusion process is used to smooth genes’ expression values based on their neighbors’ defined in a high-confidence experimentally validated protein-protein network. This function is useful when analyzing single-cell data with shallow sequencing depth because the projection reduces the dropout effects of signaling genes, in particular for possible zero expression of subunits of ligands/receptors. 
  cellchat <- projectData(cellchat, PPI.human)
  
  #计算和推断细胞间通讯网络####
  # 互作可能性计算
  cellchat <- computeCommunProb(cellchat, type = "triMean")    
  
  # 过滤表达细胞比例低的互作对
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  #推断通路水平的互作网络####
  cellchat <- computeCommunProbPathway(cellchat )
  cellchat@netP
  #汇总及展示细胞通讯网络整体结果####
  # 整合通讯网络结果
  cellchat <- aggregateNet(cellchat)
  groupSize <- as.numeric(table(cellchat@idents))
  # 计算网络中心性
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
  selectK(cellchat, pattern = "outgoing")
  #识别信号模式
  #分泌信号类型分为：outgoing和incoming
  #信号模式有助于探索多个细胞类型或信号通路的协同作用
  nPatterns = 2
  cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
  return(cellchat)
  })

library(ggalluvial)
library(NMF)
names(cellchat_list) <- c('Pre_NR','Post_NR','Pre_R','Post_R')

saveRDS(cellchat_list,file = "cellchat_list.RDS")

cellchat_list_r <- list(cellchat_list[[3]],cellchat_list[[4]])
names(cellchat_list_r) <- c('Pre_R','Post_R')

cellchat_R <- mergeCellChat(cellchat_list, 
                            add.names = c('Pre_R','Post_R'))
cellchat <- mergeCellChat(cellchat_list, 
                          add.names = c('Pre_NR','Post_NR','Pre_R','Post_R'))
cellchat
#Compare the total number of interactions and interaction strength
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

#Circle plot showing differential number of interactions or interaction strength among different cell populations across two datasets

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2

weight.max <- getMaxWeight(cellchat_list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(cellchat_list)) {
  netVisual_circle(cellchat_list[[i]]@net$count, weight.scale = T, label.edge= F, 
                   edge.weight.max = weight.max[2], edge.width.max = 12, 
                   title.name = paste0("Number of interactions - ", names(cellchat_list)[i]))
}

mat <- cellchat_list[[1]]@net$weight
groupSize <- as.numeric(table(cellchat_list[[1]]@idents))
par(mfrow = c(3,5), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), label.edge= F, title.name = rownames(mat)[i])
}
mat <- cellchat_list[[2]]@net$weight
groupSize <- as.numeric(table(cellchat_list[[2]]@idents))
par(mfrow = c(3,5), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), label.edge= F, title.name = rownames(mat)[i])
}
mat <- cellchat_list[[3]]@net$weight
groupSize <- as.numeric(table(cellchat_list[[3]]@idents))
par(mfrow = c(3,5), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), label.edge= F, title.name = rownames(mat)[i])
}
mat <- cellchat_list[[4]]@net$weight
groupSize <- as.numeric(table(cellchat_list[[4]]@idents))
par(mfrow = c(3,5), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), label.edge= F, title.name = rownames(mat)[i])
}

#Compare the major sources and targets in a 2D space
num.link <- sapply(cellchat_list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(cellchat_list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(cellchat_list[[i]], 
                                               title = names(cellchat_list)[i], 
                                               weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)

# Identify the signaling changes of specific cell populations
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Cancer stem cell", 
                                            signaling.exclude = "MIF")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "SELL+ stem T", 
                                            signaling.exclude = c("MIF"))

patchwork::wrap_plots(plots = list(gg1,gg2))



gg1 <- netAnalysis_signalingChanges_scatter(cellchat_R, idents.use = "Cancer stem cell", 
                                            signaling.exclude = "MIF")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat_R, idents.use = "SELL+ stem T", 
                                            signaling.exclude = c("MIF"))

patchwork::wrap_plots(plots = list(gg1,gg2))

#Identify signaling groups based on their functional similarity
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)

reticulate::py_install(packages = 'umap-learn')

gg1 <- rankNet(cellchat_R, mode = "comparison", measure = "weight", 
               sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat_R, mode = "comparison", measure = "weight", 
               sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)

gg1 + gg2

library(ComplexHeatmap)

i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(cellchat_list_r[[i]]@netP$pathways, cellchat_list_r[[i+1]]@netP$pathways)
pathway.union <- pathway.union[1:20]
ht1 = netAnalysis_signalingRole_heatmap(cellchat_list_r[[i]], pattern = "outgoing", 
                                        signaling = pathway.union, title = names(cellchat_list_r)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(cellchat_list_r[[i+1]], pattern = "outgoing", 
                                        signaling = pathway.union, title = names(cellchat_list_r)[i+1], width = 5, height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(cellchat_list_r[[i]], pattern = "incoming", 
                                        signaling = pathway.union, title = names(cellchat_list_r)[i], width = 5, height = 6, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(cellchat_list_r[[i+1]], pattern = "incoming", 
                                        signaling = pathway.union, title = names(cellchat_list_r)[i+1], width = 5, height = 6, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(cellchat_list_r[[i]], pattern = "all", signaling = pathway.union, title = names(cellchat_list_r)[i], width = 5, height = 6, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(cellchat_list_r[[i+1]], pattern = "all", signaling = pathway.union, title = names(cellchat_list_r)[i+1], width = 5, height = 6, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

netVisual_bubble(cellchat_R, sources.use = 4, targets.use = c(5:11),  
                 comparison = c(1, 2), angle.x = 45)
# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "Post_R"
# define a char name used for storing the results of differential expression analysis
features.name = paste0(pos.dataset, ".merged")

# perform differential expression analysis 
# Of note, compared to CellChat version < v2, CellChat v2 now performs an ultra-fast Wilcoxon test using the presto package, which gives smaller values of logFC. Thus we here set a smaller value of thresh.fc compared to the original one (thresh.fc = 0.1). Users can also provide a vector and dataframe of customized DEGs by modifying the cellchat@var.features$LS.merged and cellchat@var.features$LS.merged.info. 

cellchat_R <- identifyOverExpressedGenes(cellchat_R, group.dataset = "datasets", 
                                         pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.05,thresh.p = 0.05) 
#> Use the joint cell labels from the merged CellChat object

# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat_R, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat_R, net = net, datasets = "Post_R",ligand.logFC = 0.05, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated receptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat_R, net = net, datasets = "Pre_R",ligand.logFC = -0.05, receptor.logFC = NULL)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat_R)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat_R)


pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat_R, pairLR.use = pairLR.use.up, sources.use = 3, 
                        targets.use = 9, comparison = c(1, 2),  angle.x = 90, 
                        remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(cellchat_list_r)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat_R, pairLR.use = pairLR.use.down, sources.use = 3, 
                        targets.use = 9, comparison = c(1, 2),  angle.x = 90, 
                        remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(cellchat_list_r)[2]))
#> Comparing communications on a merged object
gg1 + gg2

gg1 <- netVisual_bubble(cellchat_R, sources.use = c(3), targets.use = c(9),  
                        comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in Post_R", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat_R, sources.use = c(3), targets.use = c(9),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in Post_R", angle.x = 45, remove.isolate = T)
gg1 + gg2

# Chord diagram
par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(cellchat_list_r[[2]], sources.use = 4, targets.use = c(5:11), 
                     slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, 
                     title.name = paste0("Up-regulated signaling in ", names(cellchat_list_r)[2]))
netVisual_chord_gene(cellchat_list_r[[1]], sources.use = 4, targets.use = c(5:11), 
                     slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, 
                     title.name = paste0("Down-regulated signaling in ", names(cellchat_list_r)[2]))

computeEnrichmentScore(net.up, species = 'human', variable.both = TRUE)

par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(cellchat_list_r)) {
  netVisual_chord_gene(cellchat_list_r[[i]], sources.use = c(3), targets.use = c(9),  
                       title.name = paste0("Signaling received by SELL+ stem T - ", 
                                           names(cellchat_list_r)[i]), legend.pos.x = 10)
}


for (i in 1:length(cellchat_list_r)) {
  netVisual_chord_gene(cellchat_list_r[[i]], sources.use = 3, targets.use = c(9), 
                       lab.cex = 0.5, title.name = paste0("Signaling from CSC - ", 
                                                          names(cellchat_list_r)[i]))
}


par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(cellchat_list_r)) {
  netVisual_chord_gene(cellchat_list_r[[i]], sources.use = 3, targets.use = c(9),slot.name = "netP", title.name = paste0("Signaling pathways sending from CSC - ", names(cellchat_list_r)[i]), legend.pos.x = 10)
}

cellchat_R@meta$datasets = factor(cellchat_R@meta$datasets, levels = c("Pre_R", "Post_R")) # set factor level
plotGeneExpression(cellchat_R, signaling = "CXCL", split.by = "datasets", colors.ggplot = T, type = "violin")

save(cellchat_R, file = "cellchat_merged_R.RData")



