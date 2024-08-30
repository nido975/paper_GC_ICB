dir <- "D:/bioinfo/GSE189926/gc"
samples=list.files( dir ,pattern = 'gz')
samples <- sapply(samples, function(x) strsplit(x, "_")[[1]][1])

do.call(rbind,lapply(seu_list, dim))
names(seu_list) =  samples
gc_list <- seu_list[c(1:6,9:13,16:22)]
samples <- samples[c(1:6,9:13,16:22)]
sce.all=merge(x=gc_list[[1]],
              y=gc_list[ -1 ],
              add.cell.ids = samples  ) 

names(sce.all@assays$RNA@layers)
sce.all[["RNA"]]$counts 
# Alternate accessor function with the same result
LayerData(sce.all, assay = "RNA", layer = "counts")
sce.all <- JoinLayers(sce.all)
dim(sce.all[["RNA"]]$counts )

###### step2: QC质控 ######
#去除双胞

sep <- as.SingleCellExperiment(sce.all)
library(BiocParallel)
sep <- scDblFinder(sep, samples="orig.ident", 
                   BPPARAM = SnowParam(workers = 12, type = "SOCK"))
table(sep$scDblFinder.class)
sce.all <- as.Seurat(sep,counts = "counts",data = NULL)
sce.all[["RNA5"]] <- as(object = sce.all[["RNA"]], Class = "Assay5")
names(sce.all)
DefaultAssay(sce.all)='RNA5'
sce.all[['RNA']]=NULL  

sce.all <- sce.all[,sce.all$scDblFinder.class %in% "singlet"]

library(Seurat)
library(decontX)
library(ggplot2)
library(cowplot)
library(beepr)

decontx_result <- decontX(sce.all@assays[["RNA5"]]@layers[["counts"]])
sce.all$Contamination <- decontx_result$contamination

sce.all.filt <- subset(sce.all, Contamination < 0.2)
 
dir.create("./1-subQC")
setwd("./1-subQC")
# 如果过滤的太狠，就需要去修改这个过滤代码
source('../scRNA_scripts/qc.R')
sce.all.filt = basic_qc(sce.all.filt)
print(dim(sce.all.filt))
setwd('..')
getwd()

saveRDS(sce.all.filt,file = "sce.all.filt.rds")

mito_genes=rownames(sce.all.filt)[grep("^MT-", rownames(sce.all.filt),ignore.case = T)] 
print(mito_genes) #可能是13个线粒体基因
ribo_genes=rownames(sce.all.filt)[grep("^Rp[sl]", rownames(sce.all.filt),ignore.case = T)]
print(ribo_genes)
# 合并线粒体基因和核糖体基因
remove_genes = c(mito_genes, ribo_genes)
# 假设你的Seurat对象是seurat_obj
# 计算所有基因的平均表达量
avg_expression <- Matrix::rowMeans(sce.all.filt@assays[["RNA5"]]@layers[["counts"]]) 
names(avg_expression) <- rownames(sce.all.filt@assays[["RNA5"]])
# 找到表达量最高的前20个基因
top_genes <- names(sort(avg_expression, decreasing = TRUE))[1:20]

# 使用VlnPlot函数绘制这些基因的表达量
VlnPlot(sce.all.filt, features = top_genes)

ig_genes <- grep("^IG", rownames(sce.all.filt@assays[["RNA5"]]), value = TRUE)
remove_genes = c(mito_genes, ribo_genes,ig_genes)
# 去除基因
sce.all.filt = subset(sce.all.filt, features = setdiff(rownames(sce.all.filt), 
                                                       remove_genes))


sce.all.filt <- NormalizeData(sce.all.filt, 
                           normalization.method = "LogNormalize",
                           scale.factor = 1e4) 
sce.all.filt <- FindVariableFeatures(sce.all.filt)
sce.all.filt <- ScaleData(sce.all.filt)
sce.all.filt <- RunPCA(sce.all.filt, features = VariableFeatures(object = sce.all.filt))

top10 <- head(VariableFeatures(sce.all.filt), 10)
plot1 <- VariableFeaturePlot(sce.all.filt)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 
plot2
sce.all.filt <- RunHarmony(sce.all.filt, "orig.ident")
names(sce.all.filt@reductions)

sce.all.filt <- FindNeighbors(sce.all.filt, dims = 1:30)

sce.all.filt <- RunUMAP(sce.all.filt, dims = 1:30, reduction = "harmony",
                        reduction.name = "UMAP_dim30", reduction.key = "UMAP_dim30_")

#设置不同的分辨率，观察分群效果(选择哪一个？)
for (res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5,0.8,1)) {
  sce.all.filt=FindClusters(sce.all.filt, #graph.name = "CCA_snn", 
                             resolution = res, algorithm = 1)
}
p2_tree=clustree(sce.all.filt@meta.data, prefix = "RNA5_snn_res.")
p2_tree
options(future.seed = TRUE)


DotPlot(sce.all.filt, features = marker)+ 
  coord_flip() + 
  RotatedAxis()


DimPlot(sce.all.filt,reduction = "UMAP_dim30",label=T,group.by = "RNA5_snn_res.1") 

# ggsave(filename='umap-by-orig.ident-after-harmony',plot = p)
markers <- FindAllMarkers(sce.all.filt, only.pos = TRUE, 
                          min.pct = 0.5, logfc.threshold = 0.25) 
write.csv(markers, "markers.csv")

top20 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.csv(top20, file = "top20.markers.csv")


#细胞周期分析
str(cc.genes) # Seurat 内置的细胞周期 markers

s_genes = cc.genes$s.genes # S 期 markers, 43个

g2m_genes = cc.genes$g2m.genes # G2M 期 markers, 54个

sce.all.filt = CellCycleScoring(sce.all.filt, s.features = s_genes, g2m.features = g2m_genes)
DimPlot(sce.all.filt, reduction = "pca", group.by = "Phase")


# install.packages("devtools")
devtools::install_github("campbio/decontX")



marker <- read.table("../MARKER.txt")
marker <- marker$V1

DotPlot(sce.all.filt, features = marker)+ 
  coord_flip() + 
  RotatedAxis()



sce.int <- subset(sce.all.filt, 
                      subset = RNA5_snn_res.1 != c(0,1,2,7,11,12,14,17,19,21,23,27,28,35,40))

sce.int[["percent.mt"]] <- PercentageFeatureSet(sce.int, pattern = "^MT-")
sce.int <- NormalizeData(sce.int)
sce.int <- FindVariableFeatures(sce.int)
sce.int <- ScaleData(sce.int)
sce.int <- RunPCA(sce.int, features = VariableFeatures(object = sce.int))
sce.int <- RunHarmony(sce.int, "orig.ident")
sce.int <- RunUMAP(sce.int, dims = 1:30, reduction = "harmony",
                   reduction.name = "UMAP_dim30", reduction.key = "UMAP_dim30_")
sce.int <- FindNeighbors(sce.int, dims = 1:30)
for (res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5,0.8,1)) {
  sce.int=FindClusters(sce.int, #graph.name = "CCA_snn", 
                            resolution = res, algorithm = 1)
}
p2_tree=clustree(sce.int@meta.data, prefix = "RNA5_snn_res.")
p2_tree

DotPlot(sce.int, features = marker)+ 
  coord_flip() + 
  RotatedAxis()

DimPlot(sce.int,reduction = "UMAP_dim30",label=T,group.by = "RNA5_snn_res.1")





#########手动给各个单细胞亚群命名##########


sce.all.int <- sce.all.int8
source('scRNA_scripts/lib.R')
dir.create('check-by-0.1')
setwd('check-by-0.1')
sel.clust = "RNA5_snn_res.0.1"
sce.all.int <- SetIdent(sce.all.int, value = sel.clust)
table(sce.all.int@active.ident) 
source('../scRNA_scripts/check-all-markers.R')
setwd('../') 
getwd()
DimPlot(sce.all.int,reduction = "UMAP_dim30",label=T,group.by = "RNA5_snn_res.0.1" ) 
DotPlot(sce.all.int, features = marker,group.by = "RNA5_snn_res.0.1")+ 
  coord_flip() + 
  RotatedAxis()

markers <- FindAllMarkers(sce.all.int, only.pos = TRUE, 
                          min.pct = 0.5, logfc.threshold = 0.25) 
write.csv(markers, "markers.csv")

top20 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.csv(top20, file = "top20.markers.csv")


library(readxl)
celltype <- read_excel("anno.xlsx")


head(celltype)
celltype
table(celltype$celltype)
sce.all.int@meta.data$celltype = "NA"

for(i in 1:nrow(celltype)){
  sce.all.int@meta.data[which(sce.all.int@meta.data$RNA5_snn_res.0.5 == celltype$cluster[i]),'celltype'] <- celltype$celltype[i]}
Idents(sce.all.int)=sce.all.int$celltype

DimPlot(sce.all.int,group.by = "celltype",label = T,raster=F)

ggsave('umap_by_celltype.pdf',width=10,height=8)

sel.clust = "celltype"
sce.all.int <- SetIdent(sce.all.int, value = sel.clust)
table(sce.all.int@active.ident) 

dir.create('check-by-celltype')
setwd('check-by-celltype')
source('../scRNA_scripts/check-all-markers.R')
setwd('../') 
getwd()
phe=sce.all.int@meta.data
save(phe,file = 'phe.Rdata')
pdf('celltype-vs-orig.ident.pdf',width = 10)

gplots::balloonplot(table(sce.all.int$celltype,sce.all.int$orig.ident))
dev.off()

saveRDS(sce.all.int,file = "sce.anno.rds")

sce.all@meta.data[["combined"]] <- ""
sce.all@meta.data[["combined"]] <- paste(sce.all$group, sce.all$celltype, sep = "_")
Idents(sce.all)=sce.all$combined

markers <- FindAllMarkers(sce.all, only.pos = TRUE, 
                          min.pct = 0.5, logfc.threshold = 0.25)
top20 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.csv(top20, file = "top20.markers.celltype.group.csv")
