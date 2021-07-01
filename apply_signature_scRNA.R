library(ggplot2)
##########################
##########################
##########################

data_path = '/home/sevastopol/data/gserranos/CART_HL/Data/antiCD19/GSE151511_RAW/'
folders <- list.files(data_path, pattern = '^ac[0-9]{2}')


CD4_signature_genes <-read.delim(file="./Data/signature/BatchK_CD4_SIGGenes_BASAL_BOTH_LowvsHigh.tsv",sep="\t",header=T)
CD8_signature_genes <-read.delim(file="./Data/signature/BatchK_CD8_SIGGenes_BASAL_BOTH_LowvsHigh.tsv",sep="\t",header=T)

CD4_signature<-read.delim(file="./Data/signature/BatchK_CD4_output_BASAL_BOTH_LowvsHigh.tsv",sep="\t",header=T)
CD8_signature<-read.delim(file="./Data/signature/BatchK_CD8_output_BASAL_BOTH_LowvsHigh.tsv",sep="\t",header=T)


apply_signature <- function(dataset, signature, genes){
    res.df <- as.data.frame(signature)
    res2.cd4 <- res.df[res.df$GeneID %in% genes$sigGenes_symbol, c('GeneID', 'log2FoldChange', 'padj')]
    # we keep the genes present on the signature list
    # bulkNormalized$gene_id <- rownames(bulkNormalized)
    dataset_bck <- dataset
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

    pdf('./Plots/Hist_signatures.pdf')
    ggplot2::ggplot(DT.cd4.high , ggplot2::aes(x=scale(High_pondered))) + ggplot2::geom_histogram (bins=50)
    ggplot2::ggplot(DT.cd4.low , ggplot2::aes(x=scale(Low_pondered))) + ggplot2::geom_histogram (bins=50)
    dev.off()

    all_signatures <- merge(DT.cd4.high, DT.cd4.low, by='Patient_ID')
    # stoped at line 202....
    set.seed(123)
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
    metadata <- data.frame(cell_id = colnames(dataset))
    metadata$CAR_HIGH_SCORE_LEVEL_HL <- car_exp.hl$CAR_level
    metadata$CAR_HIGH_SCORE_COD <- car_exp.hl$CAR_High_level_COD
    metadata$CAR_HIGH_SCORE_pondered_scaled <- car_exp.hl$V1
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
    metadata$CAR_LOW_SCORE_LEVEL_HL <- car_exp.hl$CAR_level
    metadata$CAR_LOW_SCORE_COD <- car_exp.hl$CAR_High_level_COD
    metadata$CAR_LOW_SCORE_pondered_scaled <- car_exp.hl$V1
    metadata$Overall_score <- (metadata$CAR_HIGH_SCORE_LEVEL_HL>0)*(metadata$CAR_HIGH_SCORE_LEVEL_HL +3) + metadata$CAR_LOW_SCORE_LEVEL_HL
    metadata$High_pondered_bin <- ifelse(all_signatures$High_pondered >0, 'High', 'Low')
    metadata$High_pondered <- all_signatures$High_pondered

    all_signatures_ <- all_signatures[, c('High_pondered', 'Patient_ID'), drop=FALSE]
    CAR_gene_expr <- as.data.frame(t(dataset_bck[grepl('CD19SCFV', rownames(dataset_bck)), ]))
    CAR_gene_expr$Patient_ID <- rownames(CAR_gene_expr)
    all_signatures_ <- merge(all_signatures_ , CAR_gene_expr, by='Patient_ID')
    return(metadata)
}



results_CD4 <- list()
results_CD8 <- list()
tables_CD4 <- matrix(seq(1:6), 6,1)
tables_CD8 <- matrix(seq(1:6), 6,1)
tables_percent_CD4 <- matrix(seq(1:6), 6,1)
tables_percent_CD8 <- matrix(seq(1:6), 6,1)
metadata <- data.frame(Patient.ID = sort(folders))
metadata$OS <- 'NR'
metadata$OS[metadata$Patient.ID %in% c('ac14','ac07','ac08','ac05','ac16','ac10','ac01','ac09','ac12')] <- 'CR'
car_high_score_cd4 <- list()
car_high_score_cd8 <- list()

pb <- progress::progress_bar$new(total=length(folders))
for (folder in folders){
    data_cd4 <- readRDS(paste0(data_path, folder, '/', folder, '_normalized_only_car_cd4.rds'))
    data_cd8 <- readRDS(paste0(data_path, folder, '/', folder, '_normalized_car_cd8.rds'))

    results_CD4[[folder]] <- apply_signature(data_cd4, CD4_signature, CD4_signature_genes)
    tables_CD4 <- cbind(tables_CD4, as.matrix(table(results_CD4[[folder]]$Overall_score)))
    tables_percent_CD4 <- cbind(tables_percent_CD4, as.matrix(table(results_CD4[[folder]]$Overall_score))/ sum(as.matrix(table(results_CD4[[folder]]$Overall_score)))* 100)
    ifelse(results_CD8[[folder]]$Overall_score >3, 'High', 'Low')
    results_CD8[[folder]] <- apply_signature(data_cd8, CD8_signature, CD8_signature_genes)
    tables_CD8 <- cbind(tables_CD8, as.matrix(table(results_CD8[[folder]]$Overall_score)))
    tables_percent_CD8 <- cbind(tables_percent_CD8, as.matrix(table(results_CD8[[folder]]$Overall_score))/ sum(as.matrix(table(results_CD8[[folder]]$Overall_score)))* 100)
    ifelse(results_CD8[[folder]]$Overall_score >3, 'High', 'Low')
    pb$tick()
    car_high_score_cd4[[folder]] <- results_CD4[[folder]]$High_pondered_bin
    car_high_score_cd8[[folder]] <- results_CD8[[folder]]$High_pondered_bin
}




car_high_score_df_cd4 <- as.data.frame(t(as.data.frame(lapply(car_high_score_cd4, function(x) as.numeric(prop.table(table(x))*100)))))
colnames(car_high_score_df_cd4) <- c('High', 'Low')
car_high_score_df_cd4$Patient.ID <- rownames(car_high_score_df_cd4)
car_high_score_df_cd4 <- merge(car_high_score_df_cd4, metadata, by = 'Patient.ID')

car_high_score_df_cd8 <- as.data.frame(t(as.data.frame(lapply(car_high_score_cd8, function(x) as.numeric(prop.table(table(x))*100)))))
colnames(car_high_score_df_cd8) <- c('High', 'Low')
car_high_score_df_cd8$Patient.ID <- rownames(car_high_score_df_cd8)
car_high_score_df_cd8 <- merge(car_high_score_df_cd8, metadata, by = 'Patient.ID')
plotter_cd4 <- reshape2::melt(car_high_score_df_cd4)
plotter_cd4$variable <- factor(plotter_cd4$variable, levels=c('Low', 'High'))
plotter_cd8 <- reshape2::melt(car_high_score_df_cd8)
plotter_cd8$variable <- factor(plotter_cd8$variable, levels=c('Low', 'High'))
pdf('./Plots/Signature_HIGH_TEST.pdf')
ggplot(plotter_cd4, aes(x=variable, y=value, fill=variable)) + geom_boxplot() + scale_fill_manual(values = c('#DBBE78', '#7F7F7F'))+ scale_alpha_manual(values=c(0.8)) + theme_classic() + facet_wrap(~OS)+ ggtitle('CD4')
ggplot(plotter_cd8, aes(x=variable, y=value, fill=variable)) + geom_boxplot() + scale_fill_manual(values = c('#DBBE78', '#7F7F7F'))+ scale_alpha_manual(values=c(0.8)) + theme_classic() + facet_wrap(~OS)+ ggtitle('CD8')
dev.off()

tables_percent_CD4 <- as.data.frame(t(tables_percent_CD4[,-1]))
tables_percent_CD4$Patient.ID <- folders
tables_percent_CD4 <- merge(tables_percent_CD4, metadata, by = 'Patient.ID')


tables_percent_bin_CD4 <- tables_percent_CD4[, c('Patient.ID', 'OS')]
tables_percent_bin_CD4$Low <- rowSums(tables_percent_CD4[, c('1', '2', '3')])
tables_percent_bin_CD4$High <- rowSums(tables_percent_CD4[, c('4', '5', '6')])

tables_percent_CD8 <- as.data.frame(t(tables_percent_CD8[,-1]))
tables_percent_CD8$Patient.ID <- folders
tables_percent_CD8 <- merge(tables_percent_CD8, metadata, by = 'Patient.ID')


tables_percent_bin_CD8 <- tables_percent_CD8[, c('Patient.ID', 'OS')]
tables_percent_bin_CD8$Low <- rowSums(tables_percent_CD8[, c('1', '2', '3')])
tables_percent_bin_CD8$High <- rowSums(tables_percent_CD8[, c('4', '5', '6')])

pdf('./Plots/Signature_SC_OnlyCar_CD4.pdf')
ggplot(reshape2::melt(tables_percent_CD4), aes(x=variable, y=value, group=Patient.ID, color = OS)) + geom_line() + theme_bw() + facet_wrap(~Patient.ID) + scale_color_manual(values = c('#3C77AF', '#8E221A')) + ggtitle('Percentage of cells per CAR density level')

boxplot_plotter <- reshape2::melt(tables_percent_bin_CD4)
boxplot_plotter$variable <- factor(boxplot_plotter$variable, levels=c('Low', 'High'))
cowplot::plot_grid(
    ggplot(reshape2::melt(tables_percent_CD4), aes(x=variable, y=value, group=Patient.ID, color = OS)) + geom_line() + theme_classic() + facet_wrap(~OS)+ scale_color_manual(values = c('#3C77AF', '#8E221A')) + theme(legend.position='bottom'), 
    ggplot(boxplot_plotter, aes(x=variable, y=value, fill=variable)) + geom_boxplot() + scale_fill_manual(values = c('#DBBE78', '#7F7F7F'))+ scale_alpha_manual(values=c(0.8)) + theme_classic() + facet_wrap(~OS)+ theme(legend.position='bottom') , 
ncol=2)

dev.off()


pdf('./Plots/Signature_SC_OnlyCar_CD8.pdf')
ggplot(reshape2::melt(tables_percent_CD8), aes(x=variable, y=value, group=Patient.ID, color = OS)) + geom_line() + theme_bw() + facet_wrap(~Patient.ID) + scale_color_manual(values = c('#3C77AF', '#8E221A')) + ggtitle('Percentage of cells per CAR density level')

boxplot_plotter <- reshape2::melt(tables_percent_bin_CD8)
boxplot_plotter$variable <- factor(boxplot_plotter$variable, levels=c('Low', 'High'))
cowplot::plot_grid(ggplot(reshape2::melt(tables_percent_CD8), aes(x=variable, y=value, group=Patient.ID, color = OS)) + geom_line() + theme_classic() + facet_wrap(~OS)+ scale_color_manual(values = c('#3C77AF', '#8E221A')) + theme(legend.position='bottom'), 
ggplot(boxplot_plotter, aes(x=variable, y=value, fill=variable)) + geom_boxplot() + scale_fill_manual(values = c('#DBBE78', '#7F7F7F'))+ scale_alpha_manual(values=c(0.8)) + theme_classic() + facet_wrap(~OS)+ theme(legend.position='bottom') , 
ncol=2)

dev.off()


# fisher.test(boxplot_plotter[boxplot_plotter$OS== 'CR' & boxplot_plotter$variable == 'Low','value' ], boxplot_plotter[boxplot_plotter$OS== 'CR' & boxplot_plotter$variable == 'High','value' ])