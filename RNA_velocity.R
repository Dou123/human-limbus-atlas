# For epithelium subclusters
library(velocyto.R) 
library(Seurat)
library(SeuratWrappers)
library(scales)
seurat_ob <- readRDS("limbus.Epi.subcluster.rds")
loom_list <- list("./Limbus.13y.loom", "./Limbus.53y.L.loom", "./Limbus.53y.R.loom", "./Limbus.64y.L.loom", "./Limbus.64y.R.loom", "./Limbus.95y.R.loom")
ldat_seurat_list = lapply(loom_list, function(loomx){
  ldat = ReadVelocity(loomx)
  ldat_seurat = as.Seurat(ldat)})
ldat_seurat = Seurat:::merge.Seurat(x = ldat_seurat_list[[1]],
                                    y = ldat_seurat_list[2:length(ldat_seurat_list)],
                                    merge.data = T)
seurat_ob[["spliced"]] = ldat_seurat[["spliced"]]
seurat_ob[["unspliced"]] = ldat_seurat[["unspliced"]]
seurat_ob <- RunVelocity(seurat_ob, deltaT = 1, kCells = 25, fit.quantile = 0.02, verbose=T)
DimPlot(seurat_ob,label=T)
hue_pal()(6)
all_colors = c("#F8766D", "#00BA38", "#619CFF")
ident.colors <- all_colors[1:length(unique(Idents(seurat_ob)))]
names(ident.colors) <- levels(seurat_ob)
cell.colors <- ident.colors[Idents(seurat_ob)]
names(cell.colors) <- colnames(seurat_ob)
emb = Embeddings(seurat_ob, reduction = "tsne")
show.velocity.on.embedding.cor(emb, 
                               vel = Tool(object = seurat_ob, slot = "RunVelocity"),
                               n = 200, 
                               scale = "sqrt", 
                               min.arrow.size = 0,
                               cell.colors = ac(x = cell.colors, alpha = 0.5),
                               cex = 0.8, 
                               arrow.scale = 10, 
                               show.grid.flow = TRUE,
                               min.grid.cell.mass = 0.5, 
                               grid.n = 40, 
                               arrow.lwd = 1.5,
                               do.par = F,
                               n.cores = 24, 
                               cell.border.alpha = 0)
