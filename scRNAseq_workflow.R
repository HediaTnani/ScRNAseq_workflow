# Example sensitive resistant from public data 
library(Seurat)
library(tidyverse)
library(magrittr)
library(future)
library(stringr)
library(harmony)
#options(future.globals.maxSize = 891289600)
plan("multisession")
options(future.globals.maxSize = 8000 * 1024 ^ 5)
dirname <- "~/Documents/scRNAseq/test_data"
in.list <- list.files(path=dirname, recursive = F, full.names = T) %>% sort()
in.list
samples.ls <- in.list %>% furrr::future_map( ~ .x)
names(samples.ls) <- sapply(in.list, basename) %>% str_sub(.,1,9)  %>% purrr::set_names() 
seurat.ls <- samples.ls %>% furrr::future_map( ~ Read10X(.x)) %>% furrr::future_map( ~ CreateSeuratObject(.x,min.features = 100), project = names(samples.ls))
for (i in seq_along(names(samples.ls))){
  seurat.ls[[i]]@meta.data$orig.ident <- factor(rep(names(samples.ls)[i], nrow(seurat.ls[[i]]@meta.data)))
  Idents(seurat.ls[[i]]) <- seurat.ls[[i]]@meta.data$orig.ident
}
## Plotting QC Metrics
for (i in names(samples.ls)){
  seurat.ls[[i]] <- PercentageFeatureSet(object = seurat.ls[[i]], pattern = "^MT-", col.name = "percent.MT")
  seurat.ls[[i]] <- PercentageFeatureSet(object = seurat.ls[[i]], pattern = "^RP[SL]", col.name = "percent.Ribo")
  feats <- c("nFeature_RNA", "nCount_RNA", "percent.MT", "percent.Ribo")
  p1 <- VlnPlot(seurat.ls[[i]], group.by = "orig.ident", features = feats, pt.size = 0, ncol = 4) +
    NoLegend()
  p2<- RidgePlot(seurat.ls[[i]], features=feats, ncol = 2)
  plot(p1)
}
# Filtering step using mean and SD
for (i in names(samples.ls)){
  count.max <- round(mean(seurat.ls[[i]]@meta.data$nCount_RNA) + 2 * sd(seurat.ls[[i]]@meta.data$nCount_RNA), digits = -2)
  count.min <- round(mean(seurat.ls[[i]]@meta.data$nCount_RNA) - 2 * sd(seurat.ls[[i]]@meta.data$nCount_RNA), digits = -2)
  feat.max <- round(mean(seurat.ls[[i]]@meta.data$nFeature_RNA) + 2 * sd(seurat.ls[[i]]@meta.data$nFeature_RNA), digits = -2)
  feat.min <- round(mean(seurat.ls[[i]]@meta.data$nFeature_RNA) - 2 * sd(seurat.ls[[i]]@meta.data$nFeature_RNA), digits = -2)
  mt.min <- round(mean(seurat.ls[[i]]@meta.data$percent.MT) - 2 * sd(seurat.ls[[i]]@meta.data$percent.MT))
  mt.max <- round(mean(seurat.ls[[i]]@meta.data$percent.MT) + 2 * sd(seurat.ls[[i]]@meta.data$percent.MT))
  
  ## Set minimum parameters to 0 if negative value
  if (count.min < 0){
    count.min <- 0
  } else {
    count.min <- count.min
  }
  
  if (feat.min < 0){
    feat.min <- 0
  } else {
    feat.min <- feat.min
  }
  
  ## Filter cells
  seurat.ls[[i]] <- subset(seurat.ls[[i]], subset = nFeature_RNA > feat.min & nFeature_RNA < feat.max & nCount_RNA < count.max & nCount_RNA > count.min & percent.MT < 20)
}
# QC after filtering
for (i in names(samples.ls)){
  feats <- c("nFeature_RNA", "nCount_RNA", "percent.MT", "percent.Ribo")
  p1 <- VlnPlot(seurat.ls[[i]], group.by = "orig.ident", features = feats, pt.size = 0, ncol = 4) +
    NoLegend()
  p2 <- RidgePlot(seurat.ls[[i]], features=feats, ncol = 2)
  plot(p1)
  plot(p2)
}

# Normalization and merge
# SCTransform
seurat.ls <- lapply(X = seurat.ls, FUN = function(x) {
  x <- SCTransform(object = x , conserve.memory = TRUE)
})
features <- SelectIntegrationFeatures(object.list = seurat.ls, nfeatures = 3000)
# Merge
combined <- merge(x = seurat.ls[[1]],
                  y = seurat.ls[2:length(seurat.ls)],
                  add.cell.ids = names(samples.ls)) 
# QC after merging
feats <- c("nFeature_RNA", "nCount_RNA", "percent.MT", "percent.Ribo")
p1 <-  VlnPlot(combined, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 4) +
  NoLegend()
p2 <- RidgePlot(combined, features=feats, ncol = 2)
plot(p1)
plot(p2)
# Dim reduction and finding clusters
options(future.globals.onReference = "ignore")
VariableFeatures(combined) <- features
combined <- RunPCA(object = combined, assay = "SCT", npcs = 30)
combined<- RunHarmony(object = combined,
                      assay.use = "SCT",
                      reduction = "pca",
                      dims.use = 1:30,
                      group.by.vars = "orig.ident",
                      plot_convergence = TRUE)
combined <- RunUMAP(object = combined, assay = "SCT", reduction = "harmony", dims = 1:30)
options(future.globals.onReference = "ignore")
combined <- FindNeighbors(object = combined, assay = "SCT", reduction = "harmony", dims = 1:30)
combined <- FindClusters(object = combined, resolution = 0.5)
combined <- RunTSNE(combined, dims = 1:30)
# Dim reduction visualization
plot_grid(ncol = 2,
          p1 <- DimPlot(combined, reduction = "umap", group.by = "orig.ident")+NoAxes()+ggtitle("UMAP Grouped"),
          p2 <-DimPlot(combined, reduction = "umap", split.by = "orig.ident")+NoAxes()+ggtitle("UMAP splitted"),
          p3 <- DimPlot(combined, reduction = "umap", group.by = "seurat_clusters", label = TRUE,repel = TRUE)+NoAxes()+ggtitle("Clusters"))
plot(p1+p2+p3)
# Featureplot con quality metrics
p5 <- FeaturePlot(combined, features= c("nFeature_RNA", "nCount_RNA", "percent.MT", "percent.Ribo"))
plot(p5)
# Conserved markers
library(metap)
library(multtest)
DefaultAssay(combined) <- "RNA"
get_conserved <- function(cluster){
  FindConservedMarkers(combined,
                       ident.1 = cluster,
                       grouping.var = "orig.ident",
                       only.pos = TRUE) %>%
    rownames_to_column(var = "gene") %>%
    cbind(cluster_id = cluster, .)
}
n = length(levels(combined$seurat_clusters)) -1 
conserved_markers <- map_dfr(0:n, get_conserved)
# GO per cluster
require(clusterProfiler)
require(org.Hs.eg.db) 
n = length(levels(combined$seurat_clusters)) -1 
for (cl in 0:n){
  clustern <-filter(conserved_markers, cluster_id == cl)
  genes = clustern$gene
  gene.df <- bitr(genes, fromType = "SYMBOL",
                  toType = c("ENSEMBL", "ENTREZID"),
                  OrgDb = org.Hs.eg.db)
  enriched <- enrichGO(gene         = gene.df$ENSEMBL,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "ALL",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05) 
}
