#!/usr/bin/env Rscript
#===============================================================================
#' Author: Maria E. Calleja
#' Date: 2021/MAY
#===============================================================================
####Following my annotations to sync with sevastopol and follow Guille's instructions for
## dev SimiC

#rsync -av --progress jrrodriguez@hpclogin.unav.es:/home/jrrodriguez/SC_HighLow/temp_mar/Integrated/Scores_Improved_Apr /home/sevastopol/data/mcallejac/JuanRo_SimiC/data/CART_HIGHLOW/
#rsync -av --progress /home/mecc/clustermecc/JuanRo_SimiC/src/MECC_Test_launch_SimiC_CART_HIGHLOW_sevastopol.py sevastopol@159.237.148.112:~/data/mcallejac/JuanRo_SimiC/src/
#rsync -av --progress /home/mecc/clustermecc/JuanRo_SimiC/src/MECC_SimicData_CART_HIGHLOW.R sevastopol@159.237.148.112:~/data/mcallejac/JuanRo_SimiC/src/

###NEW SIMIC route, change cross validation to TRUE and then proceed as usual. 
#@/home/sevastopol/data/gserranos/SimiC/Dev_SimiC/simicLASSO_git/code/simiclasso"

##/home/mecc/clustermecc/JuanRo_SimiC/src/MECC_Test_launch_SimiC_CART_HIGHLOW_sevastopol.py
#~/data/mcallejac/JuanRo_SimiC/src$ 
#rsync -av --progress ./MECC_Test_launch_SimiC_CART_HIGHLOW_sevastopol.py /home/sevastopol/data/gserranos/SimiC/Dev_SimiC/simicLASSO_git/code/simiclasso/
#python /home/sevastopol/data/gserranos/SimiC/Dev_SimiC/simicLASSO_git/code/simiclasso/MECC_Test_launch_SimiC_CART_HIGHLOW_sevastopol.py >> ./LaunchSimic_CART_HIGHLOW_outerr.log 2>&1




CD4_signature_genes <-read.delim(file="/home/sevastopol/data/gserranos/CART_HL/Data/signature/BatchK_CD4_SIGGenes_BASAL_BOTH_LowvsHigh.tsv",sep="\t",header=T)
CD8_signature_genes <-read.delim(file="/home/sevastopol/data/gserranos/CART_HL/Data/signature/BatchK_CD8_SIGGenes_BASAL_BOTH_LowvsHigh.tsv",sep="\t",header=T)

CD4_signature <- read.delim(file="/home/sevastopol/data/gserranos/CART_HL/Data/signature/BatchK_CD4_output_BASAL_BOTH_LowvsHigh.tsv",sep="\t",header=T)
CD8_signature <- read.delim(file="/home/sevastopol/data/gserranos/CART_HL/Data/signature/BatchK_CD8_output_BASAL_BOTH_LowvsHigh.tsv",sep="\t",header=T)


apply_signature <- function(dataset, signature, genes, sample_name, cargene = 'CAR-pCCL-BCMA'){
    results <- data.frame(cell_id = colnames(dataset))
    signature <- as.data.frame(signature)
    signature <- signature[signature$GeneID %in% genes$sigGenes_symbol, c('GeneID', 'log2FoldChange', 'padj')]
    # we keep the genes present on the signature list
    # bulkNormalized$gene_id <- rownames(bulkNormalized)
    CAR_gene_expr <- as.data.frame(t(dataset[grepl(cargene, rownames(dataset)), ]))
    CAR_gene_expr$cell_id <- rownames(CAR_gene_expr)
    dataset <- dataset[rownames(dataset) %in% genes$sigGenes_symbol,]

    logplusone <- function(x) {log(x + 0.5)}
    l <- as.data.frame(apply(dataset, 2, logplusone))
    # head(l[,1:5])
    zscore <- as.data.frame(scale(l))

    # dim(zscore)
    # head(zscore[,1:5])
    zscore$GeneID <- rownames(zscore)
    ann_markers <- merge(zscore, signature, by='GeneID')
    DT.high <- ann_markers[, !colnames(ann_markers) %in% c('GeneID', 'log2FoldChange', 'padj')] * ann_markers[, 'log2FoldChange']
    DT.low <- ann_markers[, !colnames(ann_markers) %in% c('GeneID', 'log2FoldChange', 'padj')] * -ann_markers[, 'log2FoldChange']
    DT.high <- cbind(DT.high, ann_markers[, c('GeneID', 'log2FoldChange', 'padj')])
    DT.low  <- cbind(DT.low , ann_markers[, c('GeneID', 'log2FoldChange', 'padj')])
    DT.cd4.high <- setNames(as.data.frame(colSums(DT.high[, !colnames(DT.high) %in% c('GeneID', 'log2FoldChange', 'padj')], )), c('High_pondered'))
    DT.cd4.high$Patient_ID <- rownames(DT.cd4.high)
    DT.cd4.low  <- setNames(as.data.frame(colSums(DT.low[, !colnames(DT.high) %in% c('GeneID', 'log2FoldChange', 'padj')], )), c('Low_pondered'))
    DT.cd4.low$Patient_ID <- rownames(DT.cd4.low)


    all_signatures <- merge(DT.cd4.high, DT.cd4.low, by='Patient_ID')
    # stoped at line 202....
    set.seed(123)
    pdf('./Plots/Hist_signatures.pdf')
    cowplot::plot_grid(
        ggplot2::ggplot(all_signatures , ggplot2::aes(x=High_pondered)) + ggplot2::geom_histogram (bins=50, fill = '#7F7F7F', alpha = 0.6)+ geom_vline(xintercept = 0, 
                color = "red") + theme_classic(), 
        ggplot2::ggplot(all_signatures , ggplot2::aes(x=Low_pondered)) + ggplot2::geom_histogram (bins=50, fill = '#DBBE78', alpha = 0.6)+ geom_vline(xintercept = 0, 
                color = "red", size=1)+ theme_classic(), 
        ggplot2::ggplot(all_signatures , ggplot2::aes(x=scale(High_pondered))) + ggplot2::geom_histogram (bins=50, fill = '#7F7F7F', alpha = 0.6)+ geom_vline(xintercept = 0, 
                color = "red", size=1)+ theme_classic(), 
        ggplot2::ggplot(all_signatures , ggplot2::aes(x=scale(Low_pondered))) + ggplot2::geom_histogram (bins=50, fill = '#DBBE78', alpha = 0.6)+ geom_vline(xintercept = 0, 
                color = "red", size=1)+ theme_classic(), 
    ncol=2)
    dev.off()

    #highness
    car_exp <- scale(all_signatures$High_pondered)
    # car_exp <- all_signatures$High_pondered
    car_p33 <- quantile(car_exp[car_exp>0], probs = c(0.25))
    car_p66 <- quantile(car_exp[car_exp>0], probs = c(0.75))
    # quantile(car_exp[car_exp>0], probs = c(0.99))
    car_exp.hl <- as.data.frame(car_exp)
    car_exp.hl$CAR_level <- ifelse((car_exp>0 & car_exp<=car_p33), 1 , 
                            ifelse((car_exp>car_p33 & car_exp<=car_p66) ,2, ifelse( (car_exp > car_p66), 3, 0 )))
    car_exp.hl$CAR_High_level_COD <- ifelse((car_exp>0 & car_exp<=car_p33), 'Low' , 
                            ifelse((car_exp>car_p33 & car_exp<=car_p66) ,'Med_to_high', ifelse( (car_exp > car_p66), 'High', 'Negative' )))

    car_exp.hl$Patient_ID <- all_signatures$Patient_ID
    results$CAR_HIGH_SCORE_LEVEL_HL <- car_exp.hl$CAR_level
    results$CAR_HIGH_SCORE_COD <- car_exp.hl$CAR_High_level_COD
    results$CAR_HIGH_SCORE_pondered_scaled <- car_exp.hl$V1
    # lowness
    car_exp <- scale(all_signatures$Low_pondered)
    car_p33 <- quantile(car_exp[car_exp>0], probs = c(0.25))
    car_p66 <- quantile(car_exp[car_exp>0], probs = c(0.75))
    # quantile(car_exp[car_exp>0], probs = c(0.99))
    car_exp.hl <- as.data.frame(car_exp)
    car_exp.hl$CAR_level <- ifelse((car_exp>0 & car_exp<=car_p33), 1 , 
                            ifelse((car_exp>car_p33 & car_exp<=car_p66) ,2, ifelse( (car_exp > car_p66), 3, 0 )))
    car_exp.hl$CAR_High_level_COD <- ifelse((car_exp>0 & car_exp<=car_p33), 'Low' , 
                            ifelse((car_exp>car_p33 & car_exp<=car_p66) ,'Med_to_high', ifelse( (car_exp > car_p66), 'High', 'Negative' )))

    car_exp.hl$Patient_ID <- all_signatures$Patient_ID



 

    results$CAR_LOW_SCORE_LEVEL_HL <- car_exp.hl$CAR_level
    results$CAR_LOW_SCORE_COD <- car_exp.hl$CAR_High_level_COD
    results$CAR_LOW_SCORE_pondered_scaled <- car_exp.hl$V1

    results$Overall_score <- (results$CAR_HIGH_SCORE_LEVEL_HL>0)*(results$CAR_HIGH_SCORE_LEVEL_HL +3) + results$CAR_LOW_SCORE_LEVEL_HL
    results$High_pondered_bin <- ifelse(all_signatures$High_pondered >0, 'High', 'Low')
    results$High_pondered <- all_signatures$High_pondered

    results <- merge(results, CAR_gene_expr, by = 'cell_id')
    results$Sample <- sample_name
    return(results)
}




library(reticulate)
#reticulate::use_python("/usr/bin/python3")
py_discover_config("magic")
library(Rmagic)
library(Seurat)
library(cowplot)
library(dplyr)
library(future)
library(viridis)
library(ggplot2)

options(stringsAsFactors = FALSE)
#===============================================================================
## GLOBAL VARIABLES.

#PARENT_DIR<-"/home/mcallejac"

#===============================================================================
## FUNCTIONS.
##
#===============================================================================
## MAIN.
################################################################################

hemato_data <- readRDS("./Data/HighLowCod_ctrl_integrated_seurat_cd4cd8_clusters.rds")

cell_populations <- c('ctrl')
#cell_populations <- c('stim')

Idents(hemato_data) <- hemato_data$sample
hemato_data <- subset(hemato_data, idents=cell_populations)

labels_2_keep <- c(0:1,3,4:6,8,9,10,12,15:17,19,21)
cells_2_keep_1 <- names(hemato_data$integrated_snn_res.0.8 [hemato_data$integrated_snn_res.0.8 %in% labels_2_keep])

# calculate the signature
exprMatrix_cd4 <-read.csv(file="/home/sevastopol/data/gserranos/CART_HL/Data/signature/ForGuille_exprMatrix_cd4.csv", row.names = "X")
colnames(exprMatrix_cd4) <- sub('\\.', '-', colnames(exprMatrix_cd4))
exprMatrix_cd8 <-read.csv(file="/home/sevastopol/data/gserranos/CART_HL/Data/signature/ForGuille_exprMatrix_cd8.csv", row.names = "X")
colnames(exprMatrix_cd8) <- sub('\\.', '-', colnames(exprMatrix_cd8))

results_CD4_SC_FP <- apply_signature(exprMatrix_cd4, CD4_signature, CD4_signature_genes, 'CD4', 'CAR-pCCL-BCMA')
results_CD8_SC_FP <- apply_signature(exprMatrix_cd8, CD8_signature, CD8_signature_genes, 'CD8', 'CAR-pCCL-BCMA')
cell_scores <- rbind(results_CD4_SC_FP[, c('Sample', 'High_pondered', 'cell_id')], results_CD8_SC_FP[, c('Sample', 'High_pondered', 'cell_id')])
cell_scores$BinScore <- ifelse(cell_scores$High_pondered >0 , 'High', 'Low')
rownames(cell_scores) <- cell_scores$cell_id

# labels_2_keep <- c("High", "Low")
# cells_2_keep_2 <- names(hemato_data$SCORE_COD_HIGH_LOW [hemato_data$SCORE_COD_HIGH_LOW %in% labels_2_keep]) 
### Aqui tendrias para elegir cualquiera de los tres: SCORE_CONCENS_COD_HIGH_LOW ,  SCORE_COD_HIGH_LOW, SCORE_COD_HIGH_LOW_ENHANCED
### Habría que decidir si además se incoporan los "Else", que en principio no son ni muy high, pero tampoco muy low. 
### He dejado solo high y low

# cells_2_keep <- intersect(cells_2_keep_1, cells_2_keep_2)
cells_2_keep <- cells_2_keep_1
# length(cells_2_keep)

### data from seurat ####
hemato_data_raw <- as.data.frame(hemato_data@assays$RNA@counts)
# remove unexpressed genes
hemato_data_raw <- hemato_data_raw[, colnames(hemato_data_raw) %in% cells_2_keep]
unexpresed_genes <- names(which(rowSums(abs(hemato_data_raw))<1e-6))
hemato_data_raw <- hemato_data_raw[ !rownames(hemato_data_raw) %in% unexpresed_genes, ]
dim(hemato_data_raw)

# remove undesired genes
#mt_genes <- names(which(str_detect(rownames(hemato_data_raw), "^MT-")))
mt_genes <- grep("^MT-", rownames(hemato_data_raw), value = TRUE)
rps_genes <- grep("^RP[SL]", rownames(hemato_data_raw), value = TRUE)
genes_2_discard <- union(mt_genes, rps_genes)
hemato_data_raw <- hemato_data_raw[ !rownames(hemato_data_raw) %in% genes_2_discard, ]

saveRDS(hemato_data_raw, "/home/sevastopol/data/gserranos/CART_HL/SimiC/Data/CART_HIGHLOW_INT_raw.rds")

# cell assignments
cell_populations <- c("Low", "High")
cell_group_idx <- as.data.frame(cell_scores$BinScore )
rownames(cell_group_idx) <- cell_scores$cell_id
assignment_2_post <- setNames(cell_group_idx[rownames(cell_group_idx) %in% cells_2_keep,, drop=FALSE], 'Score')
cell_group_idx <- as.character(cell_group_idx[rownames(cell_group_idx) %in% cells_2_keep,])
assignment <- as.character(seq(0,length(cell_populations)-1))
names(assignment) <- cell_populations
cluster_assignments <- assignment[cell_group_idx]
assignment_2_post$BinScore <- ifelse(assignment_2_post$Score  == 'Low', 0,1)
saveRDS(assignment_2_post, './Data/Assignments_high_low.rds')

### RUN MAGIC #### #MAGIC works on a cells x genes matrix, seurat gives a genes x cells matrix
hemato_data_raw<-Rmagic::library.size.normalize(t(hemato_data_raw))
hemato_data_raw <- sqrt(hemato_data_raw)

#reticulate::virtualenv_remove('r-reticulate')
data_MAGIC <- magic(hemato_data_raw,genes='all_genes') 

data_MAGIC_df <- as.data.frame(data_MAGIC)

data_MAGIC_df <- data_MAGIC_df

analysis_root_dir <- "/home/sevastopol/data/gserranos/CART_HL/SimiC/Data/"
#file_idx <- 'simic_val_CAR_INTEGRATED_filtered'
st=format(Sys.time(), "_%Y-%m-%d_%H_%M")
file_idx <- 'CART_HIGHLOW_INT'
DF_f <- paste0(analysis_root_dir,file_idx, st,".DF.pickle") #a genes x cells matrix
TFs_f <- paste0(analysis_root_dir,file_idx, st,".TFs.pickle") #name of genes to be use as drivers
cluster_assignment_f <- paste0(analysis_root_dir,file_idx, st,".clustAssign.txt") #phenotype of interest of the cells



### SELECT TFs and TARGETS ####
MAX_NUM_TFs=300
MAX_NUM_TARGETS=3000

TFs_list <- py_load_object("/home/sevastopol/data/mcallejac/JuanRo_SimiC/human_TF_TFchkpnt.pickle")
TFs_list <- c(TFs_list,"CAR-pCCL-BCMA")
TFs <- intersect(colnames(data_MAGIC_df),TFs_list)

MAD_TFs <- order(apply(data_MAGIC_df[,TFs],2,mad), decreasing = TRUE)
top_MAD_tfs <- na.omit(TFs[MAD_TFs[1:MAX_NUM_TFs]])

target_genes <- setdiff(colnames(data_MAGIC_df),TFs_list)

MAD_targets <- order(apply(data_MAGIC_df[,target_genes],2,mad), decreasing = TRUE)
top_MAD_targets <- na.omit(target_genes[MAD_targets[1:MAX_NUM_TARGETS]])

input_data <- data_MAGIC_df[,c(top_MAD_tfs,top_MAD_targets)]
TFs <- top_MAD_tfs

sum(is.na(TFs))
### WRITE DATA TO DISK ####

write(cluster_assignments,file = cluster_assignment_f)

sum(is.nan(as.matrix(input_data)))
sum(is.infinite(as.matrix(input_data)))

reticulate::py_save_object(as.data.frame(input_data), filename = DF_f)
reticulate::py_save_object(TFs,filename = TFs_f)

################################################################################
## Saving and etcs 
################################################################################
save.image (file =file.path(PARENT_DIR,"JuanRo_SimiC/RSession/220621_CART_HIGHLOW.RData"))
savehistory(file =file.path(PARENT_DIR,"JuanRo_SimiC/RSession/220621_CART_HIGHLOW.Rhistory"))
sink(file =file.path(PARENT_DIR,"JuanRo_SimiC/RSession/220621_CART_HIGHLOW.txt"))
toLatex(sessionInfo())
sink(NULL)
