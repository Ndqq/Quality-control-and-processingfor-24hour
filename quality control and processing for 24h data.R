
Sys.setenv(Language="en")
library(ggplot2)
library(R.utils)
library(scDblFinder)
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(Matrix)
library(SoupX)
library(harmony)
library(SeuratDisk)
library(DoubletFinder)
library(scds)
library(grid)
library(SingleCellExperiment)
library(viridis)
library(slingshot)
library(scater)
library(RColorBrewer)
library(dplyr)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(ggpubr)   
library(patchwork)

###counts after soupx
H24a <- H24a_postSoupX_20250417
H24b <- H24b_postSoupX_20250417
C24a <- C24a_postSoupX_20250417
C24b <- C24b_postSoupX_20250417

###SCDlbfinder_H24a
sce_H24a <- as.SingleCellExperiment(H24a)
sce_H24a <- scDblFinder(sce_H24a)
H24a@meta.data$doublet_score <- sce_H24a$scDblFinder.score
H24a@meta.data$doublet_class <- sce_H24a$scDblFinder.class
H24a_singlet <- subset(H24a, subset = doublet_class == "singlet")

H24a <- NormalizeData(H24a)
H24a <- FindVariableFeatures(H24a)
H24a <- ScaleData(H24a)
H24a <- RunPCA(H24a)
H24a <- RunUMAP(H24a, dims=1:35)
H24a <- FindNeighbors(H24a, dims = 1:35)
H24a <- FindClusters(H24a, resolution = 0.9) 
DimPlot(H24a, group.by = "doublet_class", cols = c("singlet" = "#21908CFF", "doublet" = "#FDE725FF")) +
  ggtitle("Doublet Classification_H24a (scDblFinder)") +
  theme_minimal()

H24a@meta.data$seurat_clusters <- Idents(H24a)
df_cluster_doublet <- H24a@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(mean_score = mean(doublet_score, na.rm = TRUE),
            doublet_rate = mean(doublet_class == "doublet"))
ggplot(df_cluster_doublet, aes(x = seurat_clusters, y = mean_score, fill = doublet_rate)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "#56B1F7", high = "#132B43", name = "Doublet Rate")+
  labs(title = "H24a_Average Doublet Score per Cluster",
       x = "Cluster", y = "Mean Doublet Score") +
  theme_minimal()

####we find cluster5 has highest percent of doublet
cluster5_doublet <- subset(H24a, idents = "5")

DefaultAssay(cluster5_doublet) <- "RNA"

FeaturePlot(cluster5_doublet, 
            features = c("mpx", "mpeg1.1", "krt17", "elavl3", "sox10"), 
            reduction = "umap", cols = c("lightgrey", "red"), pt.size = 1.2) +
  ggtitle("Potential Lineage Mix in Cluster 5")
###cluster5 is composed by the neutrophil and macrophage
###not very sure whether is the real expression of neutrophil subtype
###but I decided to clear all the doublet according to the scdblfinder


###SCDlbfinder_H24b
sce_H24b <- as.SingleCellExperiment(H24b)
sce_H24b <- scDblFinder(sce_H24b)
H24b@meta.data$doublet_score <- sce_H24b$scDblFinder.score
H24b@meta.data$doublet_class <- sce_H24b$scDblFinder.class
H24b_singlet <- subset(H24b, subset = doublet_class == "singlet")

H24b <- NormalizeData(H24b)
H24b <- FindVariableFeatures(H24b)
H24b <- ScaleData(H24b)
H24b <- RunPCA(H24b)
H24b <- RunUMAP(H24b, dims=1:35)
H24b <- FindNeighbors(H24b, dims = 1:35)
H24b <- FindClusters(H24b, resolution = 0.9) 
DimPlot(H24b, group.by = "doublet_class", cols = c("singlet" = "#21908CFF", "doublet" = "#FDE725FF")) +
  ggtitle("Doublet Classification_H24b (scDblFinder)") +
  theme_minimal()

H24b@meta.data$seurat_clusters <- Idents(H24b)
df_cluster_doublet <- H24b@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(mean_score = mean(doublet_score, na.rm = TRUE),
            doublet_rate = mean(doublet_class == "doublet"))
ggplot(df_cluster_doublet, aes(x = seurat_clusters, y = mean_score, fill = doublet_rate)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "#56B1F7", high = "#132B43", name = "Doublet Rate")+
  labs(title = "H24b_Average Doublet Score per Cluster",
       x = "Cluster", y = "Mean Doublet Score") +
  theme_minimal()

###we find cluster13 has highest percent of doublet
cluster13_doublet <- subset(H24b, idents = "13")

DefaultAssay(cluster13_doublet) <- "RNA"

FeaturePlot(cluster13_doublet, 
            features = c("mpx", "mpeg1.1", "krt17", "elavl3", "sox10"), 
            reduction = "umap", cols = c("lightgrey", "red"), pt.size = 1.2) +
  ggtitle("Potential Lineage Mix in Cluster 13")
###cluster15 is composed by the neutrophil and macrophage
###not very sure whether is the real expression of neutrophil subtype
###but I decided to clear all the doublet according to the scdblfinder

###SCDblfinder_C24a
sce_C24a <- as.SingleCellExperiment(C24a)
sce_C24a <- scDblFinder(sce_C24a)
C24a@meta.data$doublet_score <- sce_C24a$scDblFinder.score
C24a@meta.data$doublet_class <- sce_C24a$scDblFinder.class
C24a_singlet <- subset(C24a, subset = doublet_class == "singlet")
C24a <- NormalizeData(C24a)
C24a <- FindVariableFeatures(C24a)
C24a <- ScaleData(C24a)
C24a <- RunPCA(C24a)
C24a <- RunUMAP(C24a, dims=1:35)
C24a <- FindNeighbors(C24a, dims = 1:35)
C24a <- FindClusters(C24a, resolution = 0.9) 
DimPlot(C24a, group.by = "doublet_class", cols = c("singlet" = "#21908CFF", "doublet" = "#FDE725FF")) +
  ggtitle("Doublet Classification_C24a (scDblFinder)") +
  theme_minimal()

C24a@meta.data$seurat_clusters <- Idents(C24a)
df_cluster_doublet <- C24a@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(mean_score = mean(doublet_score, na.rm = TRUE),
            doublet_rate = mean(doublet_class == "doublet"))
ggplot(df_cluster_doublet, aes(x = seurat_clusters, y = mean_score, fill = doublet_rate)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "#56B1F7", high = "#132B43", name = "Doublet Rate")+
  labs(title = "C24a_Average Doublet Score per Cluster",
       x = "Cluster", y = "Mean Doublet Score") +
  theme_minimal()

###we find cluster15 has highest percent of doublet
cluster15_doublet <- subset(C24a, idents = "15")
DefaultAssay(cluster15_doublet) <- "RNA"
FeaturePlot(cluster15_doublet, 
            features = c("mpx", "mpeg1.1", "krt17", "elavl3", "sox10"), 
            reduction = "umap", cols = c("lightgrey", "red"), pt.size = 1.2) +
  ggtitle("Potential Lineage Mix in Cluster 15")
###cluster15 is composed by the neutrophil and macrophage
###not very sure whether is the real expression of neutrophil subtype
###but I decided to clear all the doublet according to the scdblfinder

###SCDblfinder_C24b
sce_C24b <- as.SingleCellExperiment(C24b)
sce_C24b <- scDblFinder(sce_C24b)
C24b@meta.data$doublet_score <- sce_C24b$scDblFinder.score
C24b@meta.data$doublet_class <- sce_C24b$scDblFinder.class
C24b_singlet <- subset(C24b, subset = doublet_class == "singlet")

C24b <- NormalizeData(C24b)
C24b <- FindVariableFeatures(C24b)
C24b <- ScaleData(C24b)
C24b <- RunPCA(C24b)
C24b <- RunUMAP(C24b, dims=1:35)
C24b <- FindNeighbors(C24b, dims = 1:35)
C24b <- FindClusters(C24b, resolution = 0.9) 
DimPlot(C24b, group.by = "doublet_class", cols = c("singlet" = "#21908CFF", "doublet" = "#FDE725FF")) +
  ggtitle("Doublet Classification_C24b (scDblFinder)") +
  theme_minimal()
C24b@meta.data$seurat_clusters <- Idents(C24b)
df_cluster_doublet <- C24b@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(mean_score = mean(doublet_score, na.rm = TRUE),
            doublet_rate = mean(doublet_class == "doublet"))
ggplot(df_cluster_doublet, aes(x = seurat_clusters, y = mean_score, fill = doublet_rate)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "#56B1F7", high = "#132B43", name = "Doublet Rate")+
  labs(title = "C24b_Average Doublet Score per Cluster",
       x = "Cluster", y = "Mean Doublet Score") +
  theme_minimal()

###merged 4 samples for next qc
###This is the cell cycle gene from the thesis of Abigail
s_genes_zf <- c("ccnd1", "ccnd2a", "ccnd2", "ccnd3", "ccne1", "ccne2", "meis1b",
                "ccna1", "ccna2", "cdc6", "mcm2", "mcm3", "mcm4", "mcm5", "mcm6", 
                "mcm7", "mcm10", "orc1", "orc6")

g2m_genes_zf <- c("ccna1", "ccna2", "foxm1", "skp2", "cdk1", "ccnb1", "ccnb2", "ccnb3", 
                  "mki67", "top2a", "oip5", "nusap1", "kif20a", "kif23", "racgap1", 
                  "sgo1", "anln")

###qc for H24a
H24a <- H24a_singlet
H24a <- NormalizeData(H24a)
H24a <- FindVariableFeatures(H24a)
H24a <- ScaleData(H24a)
H24a <- RunPCA(H24a)
H24a <- RunUMAP(H24a, dims = 1:35)
H24a <- FindNeighbors(H24a, dims = 1:35)
H24a <- FindClusters(H24a, resolution = 0.7)
DimPlot(H24a, group.by = "seurat_clusters", label = TRUE) + ggtitle("H24a singlet - Cluster UMAP")
Idents(H24a) <- "seurat_clusters"

mt_genes <- grep("^mt-", rownames(H24a), value = TRUE, ignore.case = TRUE)
H24a[["percent.mt"]] <- PercentageFeatureSet(H24a, features = mt_genes)

VlnPlot(H24a, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        pt.size = 0.1, ncol = 3 )

summary(H24a$nFeature_RNA)
summary(H24a$nCount_RNA)
summary(H24a$percent.mt)

H24a <- subset(H24a,
               subset = nFeature_RNA > 800 & 
                 nFeature_RNA < 8000 &
                 nCount_RNA > 800 &
                 nCount_RNA < 100000 & 
                 percent.mt < 12)
DefaultAssay(H24a) <- "RNA"
H24a <- NormalizeData(H24a)
H24a <- CellCycleScoring(H24a, s.features = s_genes_zf, g2m.features = g2m_genes_zf)
FeaturePlot(H24a, features = c("krt5", "mpeg1.1", "lyz","sox2","fgfrl1a"), 
            reduction = "umap", label = TRUE)


###qc for H24b
H24b <- H24b_singlet
H24b <- NormalizeData(H24b)
H24b <- FindVariableFeatures(H24b)
H24b <- ScaleData(H24b)
H24b <- RunPCA(H24b)
H24b <- RunUMAP(H24b, dims = 1:35)
H24b <- FindNeighbors(H24b, dims = 1:35)
H24b <- FindClusters(H24b, resolution = 0.7)
DimPlot(H24b, group.by = "seurat_clusters", label = TRUE) + ggtitle("H24b Clean - Cluster UMAP")
Idents(H24b) <- "seurat_clusters"


mt_genes <- grep("^mt-", rownames(H24b), value = TRUE, ignore.case = TRUE)
H24b[["percent.mt"]] <- PercentageFeatureSet(H24b, features = mt_genes)

VlnPlot(H24b, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        pt.size = 0.1, ncol = 3)

summary(H24b$nFeature_RNA)
summary(H24b$nCount_RNA)
summary(H24b$percent.mt)

H24b <- subset(H24b,
               subset = nFeature_RNA > 800 & 
                 nFeature_RNA < 9500 & 
                 nCount_RNA > 800 & 
                 nCount_RNA < 160000 & 
                 percent.mt < 12)
DefaultAssay(H24b) <- "RNA"
H24b <- NormalizeData(H24b)
H24b <- CellCycleScoring(H24b, s.features = s_genes_zf, g2m.features = g2m_genes_zf)
FeaturePlot(H24b, features = c("krt5", "mpeg1.1", "lyz","sox2","fgfrl1a"), 
            reduction = "umap", label = TRUE)


###qc for C24a
C24a <- C24a_singlet
C24a <- NormalizeData(C24a)
C24a <- FindVariableFeatures(C24a)
C24a <- ScaleData(C24a)
C24a <- RunPCA(C24a)
C24a <- RunUMAP(C24a, dims = 1:30)
C24a <- FindNeighbors(C24a, dims = 1:30)
C24a <- FindClusters(C24a, resolution = 0.7)
DimPlot(C24a, group.by = "seurat_clusters", label = TRUE) + ggtitle("C24a Clean - Cluster UMAP")
Idents(H24a) <- "seurat_clusters"


mt_genes <- grep("^mt-", rownames(C24a), value = TRUE, ignore.case = TRUE)
C24a[["percent.mt"]] <- PercentageFeatureSet(C24a, features = mt_genes)

VlnPlot(C24a, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        pt.size = 0.1, ncol = 3 )

summary(C24a$nFeature_RNA)
summary(C24a$nCount_RNA)
summary(C24a$percent.mt)

C24a <- subset(C24a,
               subset = nFeature_RNA > 800 & 
                 nFeature_RNA < 5200 & 
                 nCount_RNA > 800 &
                 nCount_RNA < 50000 & 
                 percent.mt < 12)
DefaultAssay(C24a) <- "RNA"
C24a <- NormalizeData(C24a)
C24a <- CellCycleScoring(C24a, s.features = s_genes_zf, g2m.features = g2m_genes_zf)
FeaturePlot(C24a, features = c("krt5", "mpeg1.1", "lyz","sox2","fgfrl1a"), 
            reduction = "umap", label = TRUE)



###qc for C24b

C24b <- C24b_singlet
C24b <- NormalizeData(C24b)
C24b <- FindVariableFeatures(C24b)
C24b <- ScaleData(C24b)
C24b <- RunPCA(C24b)
C24b <- RunUMAP(C24b, dims = 1:35)
C24b <- FindNeighbors(C24b, dims = 1:35)
C24b <- FindClusters(C24b, resolution = 0.7)
DimPlot(C24b, group.by = "seurat_clusters", label = TRUE) + ggtitle("C24b Clean - Cluster UMAP")

mt_genes <- grep("^mt-", rownames(C24b), value = TRUE, ignore.case = TRUE)
C24b[["percent.mt"]] <- PercentageFeatureSet(C24b, features = mt_genes)

VlnPlot(C24b, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        pt.size = 0.1, ncol = 3 )

summary(C24b$nFeature_RNA)
summary(C24b$nCount_RNA)
summary(C24b$percent.mt)

C24b <- subset(C24b,
               subset = nFeature_RNA > 800 & 
                 nFeature_RNA < 9000 & 
                 nCount_RNA > 800 & 
                 nCount_RNA < 100000 & 
                 percent.mt < 12)
DefaultAssay(C24b) <- "RNA"
C24b <- NormalizeData(C24b)
C24b <- CellCycleScoring(C24b, s.features = s_genes_zf, g2m.features = g2m_genes_zf)
FeaturePlot(C24b, features = c("krt5", "mpeg1.1", "lyz","sox2","fgfrl1a"), 
            reduction = "umap", label = TRUE)

setwd("D:/project2_cGAS-STING/scRNA_seq/merge_8_24_36_48h")
saveRDS(C24a, file = "C24a_postqc.rds")
saveRDS(H24a, file = "H24a_postqc.rds")
saveRDS(C24b, file = "C24b_postqc.rds")
saveRDS(H24b, file = "H24b_postqc.rds")

###merge 4 samples for next batch effect clear 

H24a$sample <- "H24a"
H24b$sample <- "H24b"
C24a$sample <- "C24a"
C24b$sample <- "C24b"

H24a$batch <- "a"
H24b$batch <- "b"
C24a$batch <- "a"
C24b$batch <- "b"

H24a$genotype <- "HRAS"
H24b$genotype <- "HRAS"
C24a$genotype <- "CAAX"
C24b$genotype <- "CAAX"

###merged
seurat_merged <- merge(H24a, y = list(H24b, C24a, C24b), 
                       add.cell.ids = c("H24a", "H24b", "C24a", "C24b"), 
                       project = "24_merged")

###SCT and cell cycle, percent mt regression
seurat_merged <- SCTransform(seurat_merged, 
                             vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"), 
                             verbose = FALSE)

###bacth effect
seurat_merged <- RunPCA(seurat_merged, verbose = FALSE)
seurat_merged <- RunHarmony(seurat_merged, group.by.vars = "batch")
ElbowPlot(seurat_merged, ndims = 50)
seurat_merged <- RunUMAP(seurat_merged, reduction = "harmony", dims = 1:35)
seurat_merged <- FindNeighbors(seurat_merged, reduction = "harmony", dims = 1:35)
seurat_merged <- FindClusters(seurat_merged, resolution = 0.7)

DimPlot(seurat_merged, group.by = "sample") +
  ggtitle("Sample origin after Harmony_24h")
DimPlot(seurat_merged, group.by = "Phase") +
  ggtitle("Cell Cycle Phase_24h")
DimPlot(seurat_merged, group.by = "batch") +
  ggtitle("Batch distribution (A vs B)")
DimPlot(seurat_merged, reduction = "umap", label = TRUE) +
  ggtitle("Clusters after Harmony integration_24h")
DimPlot(seurat_merged, group.by = "genotype") +
  ggtitle("24h_Genotype: HRAS vs CAAX")

###annotation
FeaturePlot(seurat_merged, features = c("krt5", "mpeg1.1", "lyz","sox2","fgfrl1a","angptl7"), 
            reduction = "umap", label = TRUE)
VlnPlot(seurat_merged, features = c("agr2", "muc5.2", "fev", "gcm2","fgfrl1a"), group.by = "seurat_clusters")
seurat_merged$celltype_anno <- plyr::mapvalues(
  x = seurat_merged$seurat_clusters,
  from = names(cluster.ids),
  to = unname(cluster.ids)
)

DimPlot(seurat_merged, group.by = "celltype_anno", label = TRUE, repel = TRUE) +
  ggtitle("UMAP with Annotated Cell Types_24h") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  )

###about the cluster which is not annotated
target_clusters <- c(7, 9, 10, 11, 13, 18, 20, 21, 22)
Idents(seurat_merged) <- "seurat_clusters"
seurat_merged <- PrepSCTFindMarkers(seurat_merged)

markers_target_list <- list()
for (cl in target_clusters) {
  markers_target_list[[as.character(cl)]] <- FindMarkers(
    object = seurat_merged,
    ident.1 = cl,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25
  ) %>% 
    tibble::rownames_to_column(var = "gene")  
  markers_target_list[[as.character(cl)]]$cluster <- cl
}


markers_target <- bind_rows(markers_target_list)

top10_markers_target <- markers_target %>%
  group_by(cluster) %>%
  top_n(10, avg_log2FC)
View(top10_markers_target)

VlnPlot(seurat_merged, features = c("ifng1"))

###subset krt
krt_subset <- subset(seurat_merged, idents= c("3","5","12","16"))
krt_subset$celltype <- "Krt"
krt_subset <- SCTransform(krt_subset, 
                          vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"), 
                          verbose = FALSE)
krt_subset <- RunPCA(krt_subset, verbose = FALSE)
krt_subset <- RunUMAP(krt_subset, dims = 1:35) 
krt_subset <- FindNeighbors(krt_subset, dims = 1:35)
krt_subset <- FindClusters(krt_subset, resolution =0.7) 
setwd("D:/project2_cGAS-STING/senescence")
saveRDS(krt_subset, file = "krt_subset.rds")
###UMAP
DimPlot(krt_subset, label = TRUE, group.by = "seurat_clusters") +
  ggtitle("Krt Reclustering by SCT")
DimPlot(krt_subset, label = FALSE, group.by = "genotype") +
  ggtitle("krt Reclustering by SCT")
###TSNE
krt_subset <- RunTSNE(krt_subset, dims = 1:35, perplexity = 15)
DimPlot(krt_subset, reduction = "tsne", label = TRUE, group.by = "seurat_clusters") +
  ggtitle("Krt Reclustering - tSNE")
DimPlot(krt_subset, reduction = "tsne", label = TRUE, group.by = "genotype") +
  ggtitle("Krt Reclustering - tSNE")
