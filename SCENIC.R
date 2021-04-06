library(dplyr)
library(stringr)
library(tidyr)
library(ggrepel)
library(cowplot)
library(Seurat)
library(tidyverse)
library(svMisc)
library(Biobase)

RemoveOutlier <- function(
  metric,
  nmads = 5,
  type = c("both", "lower", "higher"),
  subset = NULL,
  min_diff = NA
) {if (!is.null(subset)) {
    submetric <- metric[subset]
    if (length(submetric) == 0L) {
      warning("no observations remaining after subsetting")
    }
  } else {
    submetric <- metric
  }
  cur.med <- median(submetric, na.rm = TRUE)
  cur.mad <- mad(submetric, center = cur.med, na.rm = TRUE)
  diff.val <- max(min_diff, nmads * cur.mad, na.rm = TRUE)
  upper.limit <- cur.med + diff.val
  lower.limit <- cur.med - diff.val
  type <- match.arg(type)
  if (type == "lower") {
    upper.limit <- Inf
  } else if (type == "higher") {
    lower.limit <- -Inf
  }
  kx = metric < lower.limit | upper.limit < metric
  val = c(lower=lower.limit, higher=upper.limit)
  if ( log ){
    val <- 10^val
  }
  attr(kx, "thresholds") <- val
  return( kx )
}

regulons = read.csv("regulons.csv")[-c(1:2),]
colnames(regulons) =  c("TF","MotifID","AUC","Annotation","Context","MotifSimilarityQvalue","NES","OrthologousIdentity","RankAtMax","TargetGenes")
regulons_res = regulons %>% group_by(TF) %>% do(data.frame(target = paste(.$TargetGenes, collapse = "_"))) %>% ungroup()
regulons_res = as.data.frame(regulons_res )
for(i in 1:nrow(regulons_res)){
  target_str = str_extract_all(regulons_res$target[i], "\\'(\\w+)\\'" )
  target_genes = gsub("'", "", unique(target_str[[1]]))
  regulons_res[i,3] = paste(target_genes[1:length(target_genes)], collapse = ", ")
  regulons_res[i,4] = length(target_genes)
}
colnames(regulons_res) = c("TF", "Target_raw", "Target_genes", "Target_num")
regulons_tmp = gsub("^","(", regulons_res$Target_num)
regulons_res$Target_num = gsub("$","g)", regulons_tmp)
regulons_res$name <- paste(regulons_res$TF, regulons_res$Target_num, sep=" ")
regulons_num = regulons_res$name
auc = read.csv("auc_mtx.csv", header=1, row.names = 1)
colnames(auc) = regulons_num
regulonAUC_mat = t(auc)
seurat_ob = readRDS("******************************.rds")
DefaultAssay(seurat_ob) <- "RNA"
metadata = as.matrix(Idents(seurat_ob))
colnames(metadata) = "subtype"
seurat_ave = AverageExpression(seurat_ob, assays = "RNA")$RNA
colnames(seurat_ave) = paste0(colnames(seurat_ave), "_ave")
regulonActivity_byCellType <- sapply(split(rownames(metadata), as.data.frame(metadata)[["subtype"]]), function(cells){
  rowmean = apply(regulonAUC_mat[,cells], 1, function(x){
    mean(x[!RemoveOutlier(x, nmads = 3, type = "higher")])
  })
})
if( length(unique(metadata[,"subtype"])) ==  2){
  regulonActivity_byCellType_processed <- t(scale(t(regulonActivity_byCellType), center = T, scale=F))
  regulonActivity_byCellType_processed <- na.omit(regulonActivity_byCellType_processed)
}else{
  regulonActivity_byCellType_processed <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
  regulonActivity_byCellType_processed <- na.omit(regulonActivity_byCellType_processed)
}
regulonActivity_byCellType_processed <- regulonActivity_byCellType_processed[rowMax(regulonActivity_byCellType_processed) > 0,]
regulonActivity_byCellType_processed <- regulonActivity_byCellType_processed[order(regulonActivity_byCellType_processed[,1], decreasing = T),]
regulonActivity_byCellType_processed <- regulonActivity_byCellType_processed[,c(2,1)]
activiy_heatmap = pheatmap::pheatmap( regulonActivity_byCellType_processed,
                                      cluster_rows = F,
                                      color = colorRampPalette(c("navy", "white", "firebrick3"))(20),
                                      angle_col = 45,
                                      treeheight_col = 10,
                                      treeheight_row = 0,
                                      border_color = NA,
                                      cellheight = 13, 
                                      cellwidth = 13)
