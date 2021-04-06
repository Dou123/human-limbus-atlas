# 1) Doublet removal and quality control
# A total of 6 samples/libraries were used in this study: 13y.L, 53y.L, 53y.R, 64y.L, 64y.R, 95y.L; Each library was performed on DoubletFinder to remove doublets, independently. Here we take 53y.R as an example: 
library(Seurat)
library(dplyr)
library(DoubletFinder)

Limbus.53y.R <- Read10X(".Limbus.53y.R/outs/filtered_feature_bc_matrix")
Limbus.53y.R <- CreateSeuratObject(Limbus.53y.Ry.R, project = "Limbus.53y.R", min.cells = 5, min.features = 200)
Limbus.53y.R <- RenameCells(Limbus.53y.R, add.cell.id = " Limbus.53y.R")
Limbus.53y.R[["percent.mt"]] <- PercentageFeatureSet(Limbus.53y.R,  pattern = "^MT-")
VlnPlot(Limbus.53y.R, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.01)
Limbus.53y.R <- subset(Limbus.53y.R, subset = nFeature_RNA > 500 & nFeature_RNA < 7000 & percent.mt < 10)
Limbus.53y.R <- NormalizeData(Limbus.53y.R,normalization.method = "LogNormalize", scale.factor = 10000)
Limbus.53y.R <- FindVariableFeatures(Limbus.53y.R, selection.method = "vst", nfeatures = 4000)
Limbus.53y.R <- ScaleData(Limbus.53y.R, vars.to.regress = "percent.mt")
Limbus.53y.R <- RunPCA(Limbus.53y.R, features = VariableFeatures(Limbus.53y.R))
ElbowPlot(Limbus.53y.R)
DimHeatmap(Limbus.53y.R, dims = 1:20, cells = 500, balanced = TRUE)
Limbus.53y.R <- FindNeighbors(Limbus.53y.R, dims = 1:13)
Limbus.53y.R <- FindClusters(Limbus.53y.R, resolution = 0.3)
Limbus.53y.R <- RunTSNE(Limbus.53y.R, dims = 1:13)
DimPlot(Limbus.53y.R, label = T) + NoLegend()

# Doublet removal
sweep.res.list_Limbus.53y.R <- paramSweep_v3(Limbus.53y.R, PCs = 1:13, sct = FALSE)
sweep.stats_Limbus.53y.R <- summarizeSweep(sweep.res.list_Limbus.53y.R, GT = FALSE)
bcmvn_Limbus.53y.R <- find.pK(sweep.stats_Limbus.53y.R)
pK_value <- as.numeric(as.character(bcmvn_Limbus.53y.R$pK[bcmvn_Limbus.53y.R$BCmetric == max(bcmvn_Limbus.53y.R$BCmetric)]))
annotations <- Limbus.53y.R@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.05*dim(Limbus.53y.R)[2])
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_", pN_value,"_", pK_value,'_', nExp_poi)
Limbus.53y.R <- doubletFinder_v3(Limbus.53y.R, PCs = 1:13, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
Limbus.53y.R <- doubletFinder_v3(Limbus.53y.R, PCs = 1:13, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value, sct = FALSE)
Limbus.53y.R[["doublet"]] <- Limbus.53y.R@meta.data[,ncol(Limbus.53y.R@meta.data)]
Limbus.53y.R <- subset(Limbus.53y.R,doublet=="Singlet")

# Diccociation-associated cluster
FeaturePlot(Limbus.53y.R, features = c("FOSB","HSPA1A","DCN","FOS"), pt.size = 1, min.cutoff = "q9", cols = c("#EDEDED", "red"))
Limbus.53y.R <- subset(Limbus.53y.R, idents="6",invert=T)


# 2) CCA integration 
limbus.53 <- Seurat:::merge(Limbus.53y.L, Limbus.53y.R) # The same batch
limbus.64 <- Seurat:::merge(Limbus.64y.L, Limbus.64y.R) # The same batch
limbus.13 <- FindVariableFeatures(limbus.13, nfeatures = 4000)
limbus.53 <- FindVariableFeatures(limbus.53, nfeatures = 4000)
limbus.64 <- FindVariableFeatures(limbus.64, nfeatures = 4000)
limbus.95 <- FindVariableFeatures(limbus.95, nfeatures = 4000)
human.anchors <- FindIntegrationAnchors(object.list = list(limbus.13, limbus.53, limbus.64, limbus.95), dims = 1:20, anchor.features = 4000)
human.combined <- IntegrateData(anchorset = human.anchors, dims = 1:20)
Limbus <- ScaleData(Limbus, vars.to.regress = "percent.mt")
Limbus <- RunPCA(Limbus, features = VariableFeatures(object = Limbus))
DimPlot(Limbus, reduction = "pca")
ElbowPlot(Limbus)
DimHeatmap(Limbus, dims = 1:20, cells = 500, balanced = TRUE)
Limbus <- FindNeighbors(Limbus, dims = 1:16)
Limbus <- FindClusters(Limbus, resolution =0.3)
Limbus <- RunTSNE(Limbus, dims = 1:16, seed.use = 5)
DimPlot(Limbus, reduction = "tsne", label = TRUE, pt.size = 0.5, label.size = 10) + NoLegend()
Top.markers <- FindAllMarkers(Limbus, logfc.threshold = 0.25, only.pos = TRUE, min.pct=0.25)


# 3) Cell type annotation
Limbus.name <- Limbus
new.cluster.ids <- c("LEpC", "LEpC", "StC", "StC", "LEpC", "VEnC", "LEpC", "PeC", "ImC", "MeC", "StC", "LEnC", "ImC", "ScC")
names(new.cluster.ids) <- levels(Limbus.name)
Limbus.name <- RenameIdents(Limbus.name, new.cluster.ids)
DimPlot(Limbus.name, reduction = "tsne", label = TRUE, pt.size = 0.5, label.size = 10) + NoLegend()
