
##########################
##########################
##########################

metadata <- read.table('./Data/METADATA_JRR.csv', sep=',' , header=TRUE)
bulkNormalized <- as.data.frame(readRDS('./Data/Normalized_Counts_Bulk.rds'))
# exprMatrix_cd4<-read.csv(file="./Data/signature/ForGuille_exprMatrix_cd4.csv")
# exprMatrix_cd8<-read.csv(file="./Data/signature/ForGuille_exprMatrix_cd8.csv")
# cd4_metadata <- read.csv(file="./Data/signature/ForGuille_metadata_cd4.csv", row.names = "X")
# cd8_metadata <- read.csv(file="./Data/signature/ForGuille_metadata_cd8.csv", row.names = "X")

CD4_signature_genes <-read.delim(file="./Data/signature/BatchK_CD4_SIGGenes_BASAL_BOTH_LowvsHigh.tsv",sep="\t",header=T)
CD8_signature_genes <-read.delim(file="./Data/signature/BatchK_CD8_SIGGenes_BASAL_BOTH_LowvsHigh.tsv",sep="\t",header=T)

CD4_signature<-read.delim(file="./Data/signature/BatchK_CD4_output_BASAL_BOTH_LowvsHigh.tsv",sep="\t",header=T)
CD8_signature<-read.delim(file="./Data/signature/BatchK_CD8_output_BASAL_BOTH_LowvsHigh.tsv",sep="\t",header=T)


apply_signature <- function(dataset, signature, genes){
    res.df <- as.data.frame(signature)
    res2.cd4 <- res.df[res.df$GeneID %in% genes$sigGenes_symbol, c('GeneID', 'log2FoldChange', 'padj')]
    # we keep the genes present on the signature list
    # bulkNormalized$gene_id <- rownames(bulkNormalized)
    dataset <- dataset[rownames(dataset) %in% genes$sigGenes_symbol,]

    logplusone <- function(x) {log(x + 0.5)}
    l <- as.data.frame(apply(dataset, 2, logplusone))
    # head(l[,1:5])
    zscore <- as.data.frame(scale(l))

    # dim(zscore)
    # head(zscore[,1:5])
    zscore$GeneID <- rownames(zscore)
    ann_markers <- merge(zscore, res2.cd4, by='GeneID')
    DT.high <- ann_markers[, !colnames(ann_markers) %in% c('GeneID', 'log2FoldChange', 'padj')] * ann_markers[, 'log2FoldChange']
    DT.low <- ann_markers[, !colnames(ann_markers) %in% c('GeneID', 'log2FoldChange', 'padj')] * -ann_markers[, 'log2FoldChange']
    DT.high <- cbind(DT.high, ann_markers[, c('GeneID', 'log2FoldChange', 'padj')])
    DT.low  <- cbind(DT.low , ann_markers[, c('GeneID', 'log2FoldChange', 'padj')])
    DT.cd4.high <- setNames(as.data.frame(colSums(DT.high[, !colnames(DT.high) %in% c('GeneID', 'log2FoldChange', 'padj')], )), c('High_pondered'))
    DT.cd4.high$Patient_ID <- rownames(DT.cd4.high)
    DT.cd4.low  <- setNames(as.data.frame(colSums(DT.low[, !colnames(DT.high) %in% c('GeneID', 'log2FoldChange', 'padj')], )), c('Low_pondered'))
    DT.cd4.low$Patient_ID <- rownames(DT.cd4.low)

    # pdf('./Plots/Hist_signatures.pdf')
    # ggplot2::ggplot(DT.cd4.high , ggplot2::aes(x=High_pondered)) + ggplot2::geom_histogram (bins=50)
    # ggplot2::ggplot(DT.cd4.low , ggplot2::aes(x=Low_pondered)) + ggplot2::geom_histogram (bins=50)
    # dev.off()

    all_signatures <- merge(DT.cd4.high, DT.cd4.low, by='Patient_ID')
    # stoped at line 202....
    set.seed(123)
    #highness
    car_exp <- scale(all_signatures$High_pondered)
    car_p33 <- quantile(car_exp[car_exp>0], probs = c(0.25))
    car_p66 <- quantile(car_exp[car_exp>0], probs = c(0.75))
    # quantile(car_exp[car_exp>0], probs = c(0.99))
    car_exp.hl <- as.data.frame(car_exp)
    car_exp.hl$CAR_level <- ifelse((car_exp>0 & car_exp<=car_p33), 1 , 
                            ifelse((car_exp>car_p33 & car_exp<=car_p66) ,2, ifelse( (car_exp > car_p66), 3, 0 )))
    car_exp.hl$CAR_High_level_COD <- ifelse((car_exp>0 & car_exp<=car_p33), 'Low' , 
                            ifelse((car_exp>car_p33 & car_exp<=car_p66) ,'Med_to_high', ifelse( (car_exp > car_p66), 'High', 'Negative' )))

    car_exp.hl$Patient_ID <- all_signatures$Patient_ID

    metadata$CAR_HIGH_SCORE_LEVEL_HL <- car_exp.hl$CAR_level
    metadata$CAR_HIGH_SCORE_COD <- car_exp.hl$CAR_High_level_COD
    metadata$CAR_HIGH_SCORE_pondered_scaled <- car_exp.hl$V1
    # lowness
    car_exp <- scale(all_signatures$Low_pondered)
    # car_exp <- all_signatures$Low_pondered
    car_p33 <- quantile(car_exp[car_exp>0], probs = c(0.25))
    car_p66 <- quantile(car_exp[car_exp>0], probs = c(0.75))
    # quantile(car_exp[car_exp>0], probs = c(0.99))
    car_exp.hl <- as.data.frame(car_exp)
    car_exp.hl$CAR_level <- ifelse((car_exp>0 & car_exp<=car_p33), 1 , 
                            ifelse((car_exp>car_p33 & car_exp<=car_p66) ,2, ifelse( (car_exp > car_p66), 3, 0 )))
    car_exp.hl$CAR_High_level_COD <- ifelse((car_exp>0 & car_exp<=car_p33), 'Low' , 
                            ifelse((car_exp>car_p33 & car_exp<=car_p66) ,'Med_to_high', ifelse( (car_exp > car_p66), 'High', 'Negative' )))

    car_exp.hl$Patient_ID <- all_signatures$Patient_ID
    metadata$CAR_LOW_SCORE_LEVEL_HL <- car_exp.hl$CAR_level
    metadata$CAR_LOW_SCORE_COD <- car_exp.hl$CAR_High_level_COD
    metadata$CAR_LOW_SCORE_pondered_scaled <- car_exp.hl$V1
    metadata$Overall_score <- (metadata$CAR_HIGH_SCORE_LEVEL_HL>0)*(metadata$CAR_HIGH_SCORE_LEVEL_HL +3) + metadata$CAR_LOW_SCORE_LEVEL_HL
    return(metadata)
}



test_CD4 <- apply_signature(bulkNormalized, CD4_signature, CD4_signature_genes)
test_CD8 <- apply_signature(bulkNormalized, CD8_signature, CD8_signature_genes)

test_CD8$new_OS <- as.character(test_CD8$OS)
test_CD8$new_OS[test_CD8$new_OS == 'PR'] <- 'NR'
test_CD8$new_OS[test_CD8$new_OS == 'PRTD'] <- 'CR'
test_CD4$new_OS <- as.character(test_CD4$OS)
test_CD4$new_OS[test_CD4$new_OS == 'PR'] <- 'NR'
test_CD4$new_OS[test_CD4$new_OS == 'PRTD'] <- 'CR'

plot_HM <- function(dataset, title, gaps_col=NULL){
    tmp <- pheatmap::pheatmap(dataset, cluster_cols=FALSE, cluster_rows=FALSE, main = title, display_numbers = TRUE, cellheigh = 20, cellwidth=20, silent=TRUE, legend=FALSE, gaps_col = gaps_col)[[4]]
    return(tmp)
}



pdf('./Plots/Signature.pdf', 5, 7)
test_CD4$bin_overallScore <- ifelse(test_CD4$Overall_score > 3, 'High', 'Low')
test_CD8$bin_overallScore <- ifelse(test_CD8$Overall_score > 3, 'High', 'Low')
test_CD4$bin_overallScore <- factor(test_CD4$bin_overallScore, levels=c('Low', 'High'))
test_CD8$bin_overallScore <- factor(test_CD8$bin_overallScore, levels=c('Low', 'High'))

epivalcd4 <- round(Epi::twoby2(table(test_CD4$new_OS, test_CD4$bin_overallScore))$p.value[[2]], 3)
epivalcd8 <- round(Epi::twoby2(table(test_CD8$new_OS, test_CD8$bin_overallScore))$p.value[[2]], 3)

cowplot::plot_grid(
plot_HM(table(test_CD4$new_OS, test_CD4$Overall_score), 'CD4 signature', 3),
plot_HM(table(test_CD8$new_OS, test_CD8$Overall_score), 'CD8 signature', 3),
plot_HM(round(prop.table(as.table(table(test_CD4$new_OS, test_CD4$bin_overallScore)), 1)*100, 3), paste0('CD4 signature\nRelative Risk Pval:', epivalcd4)), 
plot_HM(round(prop.table(as.table(table(test_CD8$new_OS, test_CD8$bin_overallScore)), 1)*100, 3), paste0('CD8 signature\nRelative Risk Pval:', epivalcd8)), 
# plot_HM(table(test_CD4$new_OS, test_CD4$bin_overallScore), paste0('CD4 signature\nRelative Risk Pval:', epivalcd4)),
# plot_HM(table(test_CD8$new_OS, test_CD8$bin_overallScore), paste0('CD8 signature\nRelative Risk Pval:', epivalcd8)), 
ncol =2
)

dev.off()
