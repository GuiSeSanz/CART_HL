
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


apply_signature <- function(dataset, signature, genes, annotation){
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

    annotation$CAR_HIGH_SCORE_LEVEL_HL <- car_exp.hl$CAR_level
    annotation$CAR_HIGH_SCORE_COD <- car_exp.hl$CAR_High_level_COD
    annotation$CAR_HIGH_SCORE_pondered_scaled <- car_exp.hl$V1
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
    annotation$CAR_LOW_SCORE_LEVEL_HL <- car_exp.hl$CAR_level
    annotation$CAR_LOW_SCORE_COD <- car_exp.hl$CAR_High_level_COD
    annotation$CAR_LOW_SCORE_pondered_scaled <- car_exp.hl$V1
    annotation$Overall_score <- (annotation$CAR_HIGH_SCORE_LEVEL_HL>0)*(annotation$CAR_HIGH_SCORE_LEVEL_HL +3) + annotation$CAR_LOW_SCORE_LEVEL_HL
    annotation$High_pondered <- all_signatures$High_pondered
    return(annotation)
}



test_CD4 <- apply_signature(bulkNormalized, CD4_signature, CD4_signature_genes, metadata)
test_CD8 <- apply_signature(bulkNormalized, CD8_signature, CD8_signature_genes, metadata)

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




pdf('./Plots/Signature_HighScore.pdf', 5, 7)
test_CD4$bin_overallScore <- ifelse(test_CD4$High_pondered > 0, 'High', 'Low')
test_CD8$bin_overallScore <- ifelse(test_CD8$High_pondered > 0, 'High', 'Low')
test_CD4$bin_overallScore <- factor(test_CD4$bin_overallScore, levels=c('Low', 'High'))
test_CD8$bin_overallScore <- factor(test_CD8$bin_overallScore, levels=c('Low', 'High'))

epivalcd4 <- round(Epi::twoby2(table(test_CD4$new_OS, test_CD4$bin_overallScore))$p.value[[2]], 3)
epivalcd8 <- round(Epi::twoby2(table(test_CD8$new_OS, test_CD8$bin_overallScore))$p.value[[2]], 3)

cowplot::plot_grid(
plot_HM(table(test_CD4$new_OS, test_CD4$bin_overallScore), 'CD4 signature'),
plot_HM(table(test_CD8$new_OS, test_CD8$bin_overallScore), 'CD8 signature'),
plot_HM(round(prop.table(as.table(table(test_CD4$new_OS, test_CD4$bin_overallScore)), 1)*100, 3), paste0('CD4 signature\nRelative Risk Pval:', epivalcd4)), 
plot_HM(round(prop.table(as.table(table(test_CD8$new_OS, test_CD8$bin_overallScore)), 1)*100, 3), paste0('CD8 signature\nRelative Risk Pval:', epivalcd8)), 
# plot_HM(table(test_CD4$new_OS, test_CD4$bin_overallScore), paste0('CD4 signature\nRelative Risk Pval:', epivalcd4)),
# plot_HM(table(test_CD8$new_OS, test_CD8$bin_overallScore), paste0('CD8 signature\nRelative Risk Pval:', epivalcd8)), 
ncol =2
)

dev.off()

# + viridis::scale_fill_viridis(discrete = TRUE, alpha=0.6)

library(ggplot2)
library(ggprism)

pdf('./Plots/Signature_HighScore_boxplot.pdf')

test_CD4_tmp <- test_CD4
test_CD8_tmp <- test_CD8
test_CD4_tmp$CD <-  'CD4'
test_CD8_tmp$CD <-  'CD8'
tmp <- rbind(test_CD4_tmp, test_CD8_tmp)

ggplot2::ggplot(tmp, ggplot2::aes(x = new_OS, y = High_pondered,)) + ggplot2::geom_boxplot(color='black') + ggplot2::geom_point(size = 0.6, alpha = 0.4) + ggplot2::theme_classic()  + ggplot2::labs(title='Wilcoxon signed-rank test over NR and CR', x='Overall Survival', y='CAR density score') + ggplot2::facet_wrap(~CD) + guides() +   theme_prism() + 
ggsignif::geom_signif(comparisons = list(c("CR", "NR")), map_signif_level = TRUE)+ ggplot2::theme(legend.position='none')

dev.off()

+ ggplot2::scale_fill_manual(values=c(c('#FCB357', '#d3d3d3')))

# BULK CIMA



normalize_BULK <- function(path){
    print(path)
    file_name <- sub('.RData','_Normalized.rds', path)
    if (file.exists(file_name)){
        normalized_data <- readRDS(file_name)
        return(normalized_data)
    }
    cd_type <- stringr::str_extract(path, '(CD[\\d]{1})(?=\\.RData)')
    load(file = path)
    countData<-as.matrix(MyData_filt[, grep("_High|_Low", colnames(MyData_filt))])
    rownames(countData)<-row.names(MyData_filt)
    countData[!is.finite(countData)] <- 0
    class(countData)<-"integer"
    colnames(countData)

    colData <- as.data.frame(colnames(countData))
    colData <- cbind(colData, annotation_col1[grep("_High|_Low", rownames(annotation_col1)),])
    colnames(colData) <- c("SampleName", "CellType", "Day")

    colData$Donor<- as.factor(sapply(colnames(countData), function(x){
    ifelse(grepl("D10", x),"D10",
            ifelse(grepl("D14", x), "D14", 
                    ifelse(grepl("D18",x), "D18", 
                        ifelse(grepl("D7", x), "D7", 
                                ifelse(grepl("D9", x), "D9", "D8")))))}))

    colData$CellType <- as.factor(colData$CellType) 
    colData$Interactions <- as.factor(paste(colData$Day, colData$CellType, sep="."))
    colData$Batch.Interactions <- as.factor(paste(colData$Donor, colData$Day, sep="."))
    print(paste("Both elements have compatible dimensions?...",ncol(countData)==nrow(colData)))

    # Relevel Sample_Group levels to get the control level as first level (prefered)
    colData$CellType<- relevel(colData$CellType, "Low")
    colData$Interactions<- relevel(colData$Interactions, "1st.Low")

    dds <- DESeq2::DESeqDataSetFromMatrix(countData = countData,
                                        colData = colData,
                                        design = ~Donor + CellType)

    # NOW::::Run the differential expression pipeline.
    dds <- DESeq2::DESeq(dds)
    # Build the results table. We can set the alpha signification level and 
    # a lfcThreshold.
    res <- DESeq2::results(dds, alpha = 0.05, parallel = TRUE)
    res <- DESeq2::results(dds, contrast=c("CellType","High","Low"), alpha = 0.05, lfcThreshold = 1)
    normalized_data <- DESeq2::counts(dds, normalized=T)
    saveRDS(normalized_data, file_name)
    return(normalized_data)
}

CAR_NAME <- 'CAR_pCCL_BCMA'

CD4_bulk_CIMA  <- normalize_BULK("/home/sevastopol/data/mcallejac/RNA_HighLow_ALL/results/HL_CD4.RData")
metadata_CD4_bulk_CIMA <- data.frame(Sample = colnames(CD4_bulk_CIMA), type = 'CD4')
results_CD4_bulk_CIMA <-  apply_signature(CD4_bulk_CIMA, CD4_signature, CD4_signature_genes, metadata_CD4_bulk_CIMA)

CD8_bulk_CIMA <- normalize_BULK("/home/sevastopol/data/mcallejac/RNA_HighLow_ALL/results/HL_CD8.RData")
metadata_CD8_bulk_CIMA <- data.frame(Sample = colnames(CD8_bulk_CIMA), type = 'CD8')
results_CD8_bulk_CIMA<-  apply_signature(CD8_bulk_CIMA, CD4_signature, CD4_signature_genes, metadata_CD8_bulk_CIMA)


results_CD4_bulk_CIMA$BinScore <- ifelse(results_CD4_bulk_CIMA$High_pondered >=0, 'High','Low')
results_CD8_bulk_CIMA$BinScore <- ifelse(results_CD8_bulk_CIMA$High_pondered >=0, 'High','Low')
results_CD4_bulk_CIMA$NamedScore <- ifelse(stringr::str_detect(as.character(results_CD4_bulk_CIMA$Sample), 'High'), 'High','Low')
results_CD8_bulk_CIMA$NamedScore <- ifelse(stringr::str_detect(as.character(results_CD8_bulk_CIMA$Sample), 'High'), 'High','Low')
results_bulk_CIMA <- rbind(results_CD8_bulk_CIMA, results_CD4_bulk_CIMA)

car_expr <- as.data.frame(rbind(t(CD8_bulk_CIMA[CAR_NAME, , drop=FALSE]), t(CD4_bulk_CIMA[CAR_NAME, , drop=FALSE])))
car_expr$Sample <- rownames(car_expr)

results_bulk_CIMA <- merge(results_bulk_CIMA, car_expr, by='Sample')

results_bulk_CIMA$NamedScore <- factor(results_bulk_CIMA$NamedScore , levels=c('Low', 'High'))

pdf('./Plots/BULK_CIMA.pdf')
ggplot2::ggplot(results_bulk_CIMA, ggplot2::aes(x=NamedScore, y=High_pondered, fill=NamedScore)) + ggplot2::geom_boxplot() + ggplot2::facet_wrap(~type) +
ggplot2::scale_fill_manual(values = c('#DBBE78', '#7F7F7F')) + ggplot2::theme_classic() + ggplot2::ggtitle('CAR surface score')

ggplot2::ggplot(results_bulk_CIMA, ggplot2::aes(x=NamedScore, y=CAR_pCCL_BCMA, fill=NamedScore)) + ggplot2::geom_boxplot() + ggplot2::facet_wrap(~type) +
ggplot2::scale_fill_manual(values = c('#DBBE78', '#7F7F7F')) + ggplot2::theme_classic()+ ggplot2::ggtitle('CAR expression')

ggplot2::ggplot(results_bulk_CIMA, ggplot2::aes(x=High_pondered, y=CAR_pCCL_BCMA, color=NamedScore, shape=type)) + ggplot2::geom_point(size =4) + 
ggplot2::scale_color_manual(values = c('#DBBE78', '#7F7F7F')) + ggplot2::theme_classic() + ggplot2::geom_vline(xintercept=0, colour='red')
dev.off()



results_bulk_CIMA_tmp <- results_bulk_CIMA[,c('Sample', 'NamedScore', 'type', 'CAR_pCCL_BCMA', 'High_pondered')]
results_bulk_CIMA_tmp$Sample <- stringr::str_extract(results_bulk_CIMA_tmp$Sample, '^D[\\da]{1,3}')

FACS_High <- data.frame(Sample = c('D7','D8a','D9','D10','D14','D18'),
            NamedScore = 'High',
            CD4 = c(26993,18169,22593,30468,38541,34778),
            CD8 = c(22675,21135,21643,22857,36487,31001))

FACS_Low <- data.frame(Sample = c('D7','D8a','D9','D10','D14','D18'),
            NamedScore = 'Low',
            CD4 = c(2272, 2286,2338,2073,3169,3043),
            CD8 = c(2298, 2299,2313,2059,3110,2826))


FACS_High_melted <-  setNames(reshape::melt(FACS_High, measure.vars = c('CD4','CD8')), c('Sample', 'NamedScore', 'type', 'FACS_value'))
FACS_High_melted <- merge(results_bulk_CIMA_tmp, FACS_High_melted, by=c('Sample', 'NamedScore', 'type'))


FACS_Low_melted <-  setNames(reshape::melt(FACS_Low, measure.vars = c('CD4','CD8')), c('Sample', 'NamedScore', 'type', 'FACS_value'))
FACS_Low_melted <- merge(results_bulk_CIMA_tmp, FACS_Low_melted, by=c('Sample', 'NamedScore', 'type'))


FACs <- rbind(FACS_High_melted, FACS_Low_melted)

pdf('./Plots/FACS_comparison.pdf')
cowplot::plot_grid(
ggplot2::ggplot(FACs, ggplot2::aes(x=CAR_pCCL_BCMA, y=FACS_value, color=NamedScore)) + ggplot2::geom_point() +
ggplot2::scale_color_manual(values = c('#DBBE78', '#7F7F7F')) + ggplot2::theme_bw() + ggplot2::theme(legend.position='bottom') + 
ggplot2::facet_wrap(~type, nrow=2)
,
ggplot2::ggplot(FACs, ggplot2::aes(x=High_pondered, y=FACS_value, color=NamedScore)) + ggplot2::geom_point() +
ggplot2::scale_color_manual(values = c('#DBBE78', '#7F7F7F')) + ggplot2::theme_bw() + ggplot2::theme(legend.position='bottom') + 
ggplot2::facet_wrap(~type, nrow=2)
,ncol=2)
dev.off()