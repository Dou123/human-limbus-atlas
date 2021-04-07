library(Seurat)
library(dplyr)
library(Matrix)
library(monocle)

## Pseudotime analysis
seurat_ob <- readRDS("limbua.Epi.rds")

data <- as(as.matrix(seurat_ob@assays$RNA@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = seurat_ob@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)  

cds <- newCellDataSet(data,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size());
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds = detectGenes(cds,min_expr = 1)
pData(cds)$Total_mRNA = Matrix::colSums(exprs(cds))

clustering_DEGs = differentialGeneTest(cds,fullModelFormulaStr = paste0("~", "seurat_clusters"))   
ordering_genes <- row.names (subset(clustering_DEGs, qval < 0.01))
cds = setOrderingFilter(cds,ordering_genes = ordering_genes)

cds = reduceDimension(cds,max_components = 2,verbose = T,check_duplicates = F,method = 'DDRTree')
cds = orderCells(cds,reverse = F)

## Visualization

plot_cell_trajectory(cds, color_by = "Pseudotime", cell_size = 1) + scale_colour_viridis_c(option = "inferno")
plot_cell_trajectory(cds, color_by = "celltype")

seurat_ob@meta.data$Pseudotime = ""
seurat_ob@meta.data[rownames(cds@phenoData),]$Pseudotime = cds$Pseudotime
seurat_ob@meta.data$Pseudotime = as.numeric(seurat_ob@meta.data$Pseudotime)
FeaturePlot(seurat_ob, features= "Pseudotime", reduction = "tsne", pt.size= 1) + scale_colour_viridis_c(option = "inferno")

plot_cell_trajectory(cds, markers="TP63", cell_size=0.5, use_color_gradient=T) + 
  scale_color_gradient2(low="#EDEDED",mid="yellow",high="red") 


## Heatmap

to_be_tested <- row.names(subset(fData(cds),gene_short_name %in% c("TP63", "MYC", "PARP1")))
cds_subset <- cds[to_be_tested,]
plot_genes_in_pseudotime(cds_subset,color_by="seurat_clusters")  

genes=as.factor(subset(cds@featureData@data, use_for_ordering == TRUE)$gene_short_name)
to_be_tested <- row.names(subset(fData(cds),gene_short_name %in% levels(genes)))
cds_sub <- cds[to_be_tested,]
heatmap = plot_pseudotime_heatmap(cds_sub, cluster_rows = T, num_clusters = 3, show_rownames = F, return_heatmap = T )

## TF heatmap
genes=read.table("TFlist.txt")
to_be_tested <- row.names(subset(fData(cds),gene_short_name %in% levels(genes[,1])))
cds_sub <- cds[to_be_tested,]
heatmap = plot_pseudotime_heatmap(cds_sub, cluster_rows = T, num_clusters = 3, show_rownames = T, return_heatmap = T )
