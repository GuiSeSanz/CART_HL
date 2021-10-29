# library(GSVA)
library(Seurat)
library(parallel)

# get the signatures
gene_sets <- list()
gene_sets[['High_CD4']] <- as.character(read.delim(file="./Data/signature/BatchK_CD4_SIGGenes_BASAL_BOTH_LowvsHigh.tsv",sep="\t",header=T)$sigGenes_symbol)
gene_sets[['High_CD8']]<- as.character(read.delim(file="./Data/signature/BatchK_CD8_SIGGenes_BASAL_BOTH_LowvsHigh.tsv",sep="\t",header=T)$sigGenes_symbol)

signatures_path <- '/home/sevastopol/data/gserranos/CART_HL/Data/signature/OtherSignatures'
for (sgn in list.files(signatures_path)){
    sgn_name <- stringr::str_extract(sgn, '(?<=_)[\\w]+')
    gene_sets[[sgn_name]] <- as.character(read.delim(paste0(signatures_path, '/', sgn))[,1])
}

GeneSets_HL <- readRDS('./Data/GeneSets_HL.rds')
for (sgn in names(GeneSets_HL)){
    gene_sets[[sgn]] <-  as.character(GeneSets_HL[[sgn]][complete.cases(GeneSets_HL[[sgn]]),])
}

rds_test <- readRDS('/home/sevastopol/data/mcallejac/JuanRo_SimiC/data/CART_HIGHLOW/Scores_Improved_Apr/HighLowCod_ctrl_integrated_seurat_cd4cd8_clusters.rds')
normData_ScaledNorm <- as.data.frame(rds_test@assays$RNA@scale.data)
normData_Norm <- as.data.frame(as.matrix(rds_test@assays$RNA@data))
coords_SC <- as.data.frame(rds_test@reductions$umap@cell.embeddings)

norm_matrix <- readRDS('./Data/Norm_counts_zscaled_rlog_HL.rds')
CD8_BULK <- read.delim(file="./Data/signature/BatchK_CD8_output_BASAL_BOTH_LowvsHigh.tsv",sep="\t",header=T)
CD4_BULK <- read.delim(file="./Data/signature/BatchK_CD4_output_BASAL_BOTH_LowvsHigh.tsv",sep="\t",header=T)


rownames(CD8_BULK) <- CD8_BULK$GeneID
CD8_BULK <- CD8_BULK[, grepl('^D', colnames(CD8_BULK))]

rownames(CD4_BULK) <- CD4_BULK$GeneID
CD4_BULK <- CD4_BULK[, grepl('^D', colnames(CD4_BULK))]


gsva_ScaledNorm <- GSVA::gsva(as.matrix(normData_ScaledNorm), gene_sets, verbose=TRUE)
saveRDS(gsva_ScaledNorm, './Data/Signatures_gsva_ScaledNorm.rds')
gsva_Norm       <- GSVA::gsva(as.matrix(normData_Norm), gene_sets, verbose=TRUE)
saveRDS(gsva_Norm, './Data/Signatures_gsva_Norm.rds')


gsva_CD8Bulk <- GSVA::gsva(as.matrix(CD8_BULK), gene_sets, verbose=TRUE)
saveRDS(gsva_CD8Bulk, './Data/Signatures_gsva_CD8Bulk.rds')

gsva_CD4_BULK <- GSVA::gsva(as.matrix(CD4_BULK), gene_sets, verbose=TRUE)
saveRDS(gsva_CD4_BULK, './Data/Signatures_gsva_CD4_BULK.rds')

# 
# GSVA scores are calculated non-parametrically using a KS-like random walk statistic and a negative value for a particular sample and gene set, means that the gene set has a lower expression than the same gene set with 
# a positive value in a different sample, or than another different gene set with a positive value. Whether you should expect positive or negative values for a particular gene set depends on the expression levels of the 
# genes that form that gene set with respect to the expression levels of the genes outside that gene set. Think of a ranking of genes by decreasing expression levels in a particular sample, so that the top of the ranking 
# contains the genes with higher expression levels and the bottom of the ranking contains the genes with lower expression levels. Now imagine a gene set whose genes are mostly located at the bottom of that ranking. That 
# gene set is likely to get a negative KS-like random walk statistic (what we simply call a 'GSVA score'). Here below, you should see the image of Figure 1 of the paper, step three of the sketched GSVA algorithm shows a 
# toy KS-like random walk for a gene set with 3 genes whose expression values are at the top of the ranking. This will lead to a positive GSVA score, but if you imagine the gene set at the bottom of the ranking, the two 
# random walk step CDFs (inside and outside the gene set) would be inverted and the GSVA score would be negative."
# 
