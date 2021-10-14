
library(ggplot2)
library(Seurat)
##########################
##########################
##########################

data_path = '/home/sevastopol/data/gserranos/CART_HL/Data/antiCD19/GSE151511_RAW/'
folders <- list.files(data_path, pattern = '^ac[0-9]{2}')


CD4_signature_genes <-read.delim(file="./Data/signature/BatchK_CD4_SIGGenes_BASAL_BOTH_LowvsHigh.tsv",sep="\t",header=T)
CD8_signature_genes <-read.delim(file="./Data/signature/BatchK_CD8_SIGGenes_BASAL_BOTH_LowvsHigh.tsv",sep="\t",header=T)

CD4_signature <- read.delim(file="./Data/signature/BatchK_CD4_output_BASAL_BOTH_LowvsHigh.tsv",sep="\t",header=T)
CD8_signature <- read.delim(file="./Data/signature/BatchK_CD8_output_BASAL_BOTH_LowvsHigh.tsv",sep="\t",header=T)


apply_signature <- function(dataset, signature, genes, sample_name, cargene = 'CD19SCFV'){
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


metadata <- data.frame(Sample = sort(folders))
metadata$OS <- 'NR'
metadata$OS[metadata$Sample %in% c('ac14','ac07','ac08','ac05','ac16','ac10','ac01','ac09','ac12')] <- 'CR'
pb <- progress::progress_bar$new(total=length(folders))
for (folder in folders){
    data_cd4 <- readRDS(paste0(data_path, folder, '/', folder, '_normalized_only_car_cd4.rds'))
    data_cd8 <- readRDS(paste0(data_path, folder, '/', folder, '_normalized_car_cd8.rds'))

    if (folder == folders[1]){
        results_CD4 <- apply_signature(data_cd4, CD4_signature, CD4_signature_genes, folder)
        results_CD8 <- apply_signature(data_cd8, CD8_signature, CD8_signature_genes, folder)

    }else{
        results_CD4 <- rbind(results_CD4, apply_signature(data_cd4, CD4_signature, CD4_signature_genes, folder))
        results_CD8 <- rbind(results_CD8, apply_signature(data_cd8, CD8_signature, CD8_signature_genes, folder))
    }
    pb$tick()

}

results_CD4 <- merge(results_CD4, metadata, by='Sample') 
results_CD8 <- merge(results_CD8, metadata, by='Sample') 

pdf('./Plots/Signature_HIGH_TEST.pdf')
plotter_cd4 <- results_CD4[, c('High_pondered_bin', 'Sample')]
plotter_cd4 <- t(table(plotter_cd4$High_pondered_bin, plotter_cd4$Sample))
plotter_cd4 <- as.data.frame(t(apply(plotter_cd4 , 1, FUN=function(x) x/sum(x)*100)))
plotter_cd4$Sample <- rownames(plotter_cd4)
plotter_cd4 <- merge(plotter_cd4, metadata, by='Sample')
plotter_cd4 <- reshape2::melt(plotter_cd4)
plotter_cd4$variable <- factor(plotter_cd4$variable, levels=c('Low', 'High'))

plotter_cd8 <- results_CD8[, c('High_pondered_bin', 'Sample')]
plotter_cd8 <- t(table(plotter_cd8$High_pondered_bin, plotter_cd8$Sample))
plotter_cd8 <- as.data.frame(t(apply(plotter_cd8 , 1, FUN=function(x) x/sum(x)*100)))
plotter_cd8$Sample <- rownames(plotter_cd8)
plotter_cd8 <- merge(plotter_cd8, metadata, by='Sample')
plotter_cd8 <- reshape2::melt(plotter_cd8)
plotter_cd8$variable <- factor(plotter_cd8$variable, levels=c('Low', 'High'))
cowplot::plot_grid(
    ggplot(plotter_cd4, aes(x=variable, y=value, fill=variable)) + geom_boxplot() + scale_fill_manual(values = c('#DBBE78', '#7F7F7F'))+ scale_alpha_manual(values=c(0.8)) + theme_classic() + facet_wrap(~OS) + ggtitle('CD4') + theme(legend.position='none'), 
    ggplot(plotter_cd8, aes(x=variable, y=value, fill=variable)) + geom_boxplot() + scale_fill_manual(values = c('#DBBE78', '#7F7F7F'))+ scale_alpha_manual(values=c(0.8)) + theme_classic() + facet_wrap(~OS) + ggtitle('CD8') + theme(legend.position='none'), 
ncol=2)
dev.off()


pdf('./Plots/Signature_HIGH_TEST_cd8.pdf')

plotter_cd8 <- results_CD8[, c('High_pondered_bin', 'Sample')]
plotter_cd8 <- t(table(plotter_cd8$High_pondered_bin, plotter_cd8$Sample))
plotter_cd8 <- as.data.frame(t(apply(plotter_cd8 , 1, FUN=function(x) x/sum(x)*100)))
plotter_cd8$Sample <- rownames(plotter_cd8)
plotter_cd8 <- merge(plotter_cd8, metadata, by='Sample')
plotter_cd8 <- reshape2::melt(plotter_cd8)
plotter_cd8$variable <- factor(plotter_cd8$variable, levels=c('Low', 'High'))
    ggplot(plotter_cd8, aes(x=variable, y=value, fill=variable)) + geom_boxplot() + scale_fill_manual(values = c('#DBBE78', '#7F7F7F'))+ scale_alpha_manual(values=c(0.8)) + theme_classic() + facet_wrap(~OS) + ggtitle('DENG single cell data') + theme(legend.position='none') + ggsignif::geom_signif(comparisons = list(c("High", "Low")), map_signif_level = TRUE, vjust=0.5) + ggprism::theme_prism() + labs(y='Cell percentage')
dev.off()


pdf('./Plots/Signature_SC_OnlyCar_CD4.pdf')
    plotter_cd4 <- table(results_CD4$Sample, results_CD4$Overall_score)
    plotter_cd4 <- as.data.frame(t(apply(plotter_cd4 , 1, FUN=function(x) x/sum(x)*100)))
    plotter_cd4$Sample <- rownames(plotter_cd4)
    plotter_cd4 <- merge(plotter_cd4, metadata, by = 'Sample')
    ggplot(reshape2::melt(plotter_cd4), aes(x=variable, y=value, group=Sample, color = OS)) + geom_line() + theme_bw() + facet_wrap(~Sample) + scale_color_manual(values = c('#3C77AF', '#8E221A')) + ggtitle('Percentage of cells per CAR density level')
    plotter <- plotter_cd4
    plotter$Low  <- rowSums(plotter_cd4[, c('1', '2', '3')])
    plotter$High <- rowSums(plotter_cd4[, c('4', '5', '6')])
    plotter <- reshape2::melt(plotter[, c('Sample', 'OS', 'Low', 'High')])
    plotter$variable <-  factor(plotter$variable, levels=c('Low', 'High'))
    cowplot::plot_grid(
        ggplot(reshape2::melt(plotter_cd4), aes(x=variable, y=value, group=Sample, color = OS)) + geom_line() + theme_bw() + facet_wrap(~OS) + scale_color_manual(values = c('#3C77AF', '#8E221A')) + theme(legend.position='none'),
        ggplot(plotter, aes(x=variable, y=value, fill=variable)) + geom_boxplot() + scale_fill_manual(values = c('#DBBE78', '#7F7F7F'))+ scale_alpha_manual(values=c(0.8)) + theme_classic() + facet_wrap(~OS)+ theme(legend.position='none') , 
    ncol=2)

dev.off()


pdf('./Plots/Signature_SC_OnlyCar_CD8.pdf')
    plotter_cd8 <- table(results_CD8$Sample, results_CD8$Overall_score)
    plotter_cd8 <- as.data.frame(t(apply(plotter_cd8 , 1, FUN=function(x) x/sum(x)*100)))
    plotter_cd8$Sample <- rownames(plotter_cd8)
    plotter_cd8 <- merge(plotter_cd8, metadata, by = 'Sample')
    ggplot(reshape2::melt(plotter_cd8), aes(x=variable, y=value, group=Sample, color = OS)) + geom_line() + theme_bw() + facet_wrap(~Sample) + scale_color_manual(values = c('#3C77AF', '#8E221A')) + ggtitle('Percentage of cells per CAR density level')
    plotter <- plotter_cd8
    plotter$Low  <- rowSums(plotter_cd8[, c('1', '2', '3')])
    plotter$High <- rowSums(plotter_cd8[, c('4', '5', '6')])
    plotter <- reshape2::melt(plotter[, c('Sample', 'OS', 'Low', 'High')])
    plotter$variable <-  factor(plotter$variable, levels=c('Low', 'High'))
    cowplot::plot_grid(
        ggplot(reshape2::melt(plotter_cd8), aes(x=variable, y=value, group=Sample, color = OS)) + geom_line() + theme_bw() + facet_wrap(~OS) + scale_color_manual(values = c('#3C77AF', '#8E221A')) + theme(legend.position='none'),
        ggplot(plotter, aes(x=variable, y=value, fill=variable)) + geom_boxplot() + scale_fill_manual(values = c('#DBBE78', '#7F7F7F'))+ scale_alpha_manual(values=c(0.8)) + theme_classic() + facet_wrap(~OS)+ theme(legend.position='none') , 
    ncol=2)

dev.off()


results_CD4$CD <- 'CD4'
results_CD8$CD <- 'CD8'
CAR_NAME <- 'FMC63-CD19SCFV'
car_exp_level <- setNames(rbind(results_CD4[, c('Sample', 'cell_id', 'High_pondered', 'Overall_score', CAR_NAME, 'OS', 'CD')], 
                       results_CD8[, c('Sample', 'cell_id', 'High_pondered', 'Overall_score',CAR_NAME, 'OS', 'CD')] ), 
                       c('Sample', 'cell_id', 'High_pondered', 'Overall_score','CAR_Exp', 'OS', 'CD'))

car_exp_level$Sample_OS <- paste0(car_exp_level$Sample, '_', car_exp_level$OS)

pdf('./Plots/Correlation_CarVsSig.pdf')
annotation_percentage <- as.data.frame.matrix(table(car_exp_level$Sample_OS, car_exp_level$High_pondered >=0))
annotation_percentage <- as.data.frame.matrix(t(apply(annotation_percentage, 1, FUN=function(x) x/sum(x)*100)))
annotation_percentage$Sample_OS <- rownames(annotation_percentage)
annotation_percentage$Hig_percent <- signif(annotation_percentage[, 'TRUE'], digits=3)
annotation_percentage$Sample <- 'CD4'
annotation_percentage$CD <- 'CD4'

ggplot(car_exp_level, aes(x = High_pondered, y = CAR_Exp, group=Sample, color=CD)) + geom_point(size= 0.8, alpha =0.6) + scale_color_manual(values = c('#3C77AF', '#8E221A')) + theme_bw() + geom_vline(xintercept = 0, color = "red") + ggtitle('CD4+CD8')+ facet_wrap(~Sample_OS, nrow=6) + geom_text(x = 100, y = 3.5, aes(label = Hig_percent), data = annotation_percentage)#+ geom_rug(aes(color=CD), outside = TRUE)
ggplot(car_exp_level, aes(x = Overall_score, y = CAR_Exp, group=Sample, color=CD)) + geom_point(size= 0.8, alpha =0.6) + scale_color_manual(values = c('#3C77AF', '#8E221A')) + theme_bw()  + ggtitle('CD4+CD8')+ facet_wrap(~Sample_OS, nrow=6)
ggplot(car_exp_level, aes(x = Overall_score, y = CAR_Exp, group=Sample, color=CD)) + geom_point(size= 0.8, alpha =0.6) + scale_color_manual(values = c('#3C77AF', '#8E221A')) + theme_bw() + ggtitle('CD4+CD8')+ facet_wrap(~Sample_OS, nrow=6)
dev.off()


pdf('./Plots/DENG_results_signature.pdf')
annotation_percentage <- as.data.frame.matrix(table(car_exp_level$Sample_OS, car_exp_level$High_pondered >=0))
annotation_percentage <- as.data.frame.matrix(t(apply(annotation_percentage, 1, FUN=function(x) x/sum(x)*100)))
annotation_percentage$Sample_OS <- rownames(annotation_percentage)
annotation_percentage$Hig_percent <- signif(annotation_percentage[, 'TRUE'], digits=3)
annotation_percentage$Sample <- 'CD4'
annotation_percentage$CD <- 'CD4'

ggplot(car_exp_level, aes(x = High_pondered, y = CAR_Exp, group=Sample, color=CD)) + geom_point(size= 0.8, alpha =0.6) + scale_color_manual(values = c('#3C77AF', '#8E221A')) + theme_bw() + geom_vline(xintercept = 0, color = "red") + ggtitle('CD4+CD8')+ facet_wrap(~Sample_OS, nrow=6) + geom_text(x = 100, y = 3.5, aes(label = Hig_percent), data = annotation_percentage)#+ geom_rug(aes(color=CD), outside = TRUE)
ggplot(car_exp_level, aes(x = OS, y= High_pondered, fill=OS)) + geom_boxplot() + scale_fill_manual(values = c('#DBBE78', '#7F7F7F'))+ scale_alpha_manual(values=c(0.8)) + theme_classic() + facet_wrap(~CD)+ theme(legend.position='none')
dev.off()


########################
#### Single cell by CIMA
########################

exprMatrix_cd4 <-read.csv(file="./Data/signature/ForGuille_exprMatrix_cd4.csv", row.names = "X")
exprMatrix_cd8 <-read.csv(file="./Data/signature/ForGuille_exprMatrix_cd8.csv", row.names = "X")
cd4_metadata <- read.csv(file="./Data/signature/ForGuille_metadata_cd4.csv", row.names = "X")
cd8_metadata <- read.csv(file="./Data/signature/ForGuille_metadata_cd8.csv", row.names = "X")

cd4_metadata$cell_id <- sub('-', '.', cd4_metadata$don_cells)
cd8_metadata$cell_id <- sub('-', '.', cd8_metadata$don_cells)


results_CD4_SC_FP <- apply_signature(exprMatrix_cd4, CD4_signature, CD4_signature_genes, 'CD4', 'CAR-pCCL-BCMA')
results_CD8_SC_FP <- apply_signature(exprMatrix_cd8, CD8_signature, CD8_signature_genes, 'CD8', 'CAR-pCCL-BCMA')

results_CD4_SC_FP <- merge(cd4_metadata, results_CD4_SC_FP, by='cell_id')
results_CD8_SC_FP <- merge(cd8_metadata, results_CD8_SC_FP, by='cell_id')
results_CD4_SC_FP$CD <- 'CD4'
results_CD8_SC_FP$CD <- 'CD8'

car_exp_level_SC_FP <- setNames(rbind(results_CD4_SC_FP[, c('Sample', 'cell_id', 'High_pondered', 'donor', 'Overall_score', 'CAR-pCCL-BCMA', 'CD')], 
                       results_CD8_SC_FP[, c('Sample', 'cell_id', 'High_pondered', 'donor','Overall_score','CAR-pCCL-BCMA', 'CD')] ), 
                       c('Sample', 'cell_id', 'High_pondered', 'donor','Overall_score','CAR_Exp', 'CD'))


pdf('./Plots/Correlation_CarVsSig_SCFP.pdf')
# ggplot(car_exp_level, aes(x = High_pondered, y = CAR_Exp, group=Sample, color=OS)) + geom_point(size= 0.8, alpha =0.6) + scale_color_manual(values = c('#3C77AF', '#8E221A')) + theme_bw() + geom_vline(xintercept = 0, color = "red") + ggtitle('CD4+CD8')+ facet_wrap(~Sample, nrow=6) 
annotation_percentage <- as.data.frame.matrix(table(car_exp_level_SC_FP$donor, car_exp_level_SC_FP$High_pondered >=0))
annotation_percentage <- as.data.frame.matrix(t(apply(annotation_percentage, 1, FUN=function(x) x/sum(x)*100)))
annotation_percentage$donor <- rownames(annotation_percentage)
annotation_percentage$Hig_percent <- signif(annotation_percentage[, 'TRUE'], digits=3)
annotation_percentage$Sample <- 'CD4'
annotation_percentage$CD <- 'CD4'
ggplot(car_exp_level_SC_FP, aes(x = High_pondered, y = CAR_Exp, group=Sample, color=CD)) + geom_point(size= 0.8, alpha =0.6) + scale_color_manual(values = c('#3C77AF', '#8E221A')) + theme_bw() + geom_vline(xintercept = 0, color = "red") + ggtitle('CD4+CD8')+ facet_wrap(~donor, nrow=6) + geom_text(x = 150, y = 5.5, aes(label = Hig_percent), data = annotation_percentage) #+ geom_rug(aes(color=CD), outside = TRUE)


ggplot(car_exp_level_SC_FP, aes(x = Overall_score, y = CAR_Exp, group=Sample, color=CD)) + geom_point(size= 0.8, alpha =0.6) + scale_color_manual(values = c('#3C77AF', '#8E221A')) + theme_bw() + ggtitle('CD4+CD8')+ facet_wrap(~donor, nrow=6)
dev.off()




######

car_exp_level_SC_FP$cell_id <- sub('\\.', '-', car_exp_level_SC_FP$cell_id)
rds_test <- readRDS('/home/sevastopol/data/mcallejac/JuanRo_SimiC/data/CART_HIGHLOW/Scores_Improved_Apr/HighLowCod_ctrl_integrated_seurat_cd4cd8_clusters.rds')
coords <- as.data.frame(rds_test@reductions$umap@cell.embeddings)
coords$cell_id <- rownames(coords)

clusters <- setNames(as.data.frame(rds_test$ClusterNames_0.8_by_JR), 'Cluster')
clusters$cell_id <- rownames(clusters)
coords <- merge(coords, car_exp_level_SC_FP[, c('cell_id', 'High_pondered', 'CAR_Exp', 'CD')], by='cell_id')
coords$BinScore <- ifelse(coords$High_pondered > 0, 'High', 'Low')

coords <- merge(coords, clusters, by='cell_id')
coords$Cluster <- stringr::str_remove(coords$Cluster, '^C?[\\d]{1,2}\\.')

pdf('./Plots/Umap_Score_CarExp.pdf', width=8)
cowplot::plot_grid(
ggplot(coords, aes(x=UMAP_1, y=UMAP_2, color= High_pondered, shape = CD)) + geom_point(alpha=0.6) + viridis::scale_color_viridis()  + theme_bw() + labs(subtitle = 'High-low score') + theme(legend.position='bottom') +guides(shape=FALSE) ,
ggplot(coords, aes(x=UMAP_1, y=UMAP_2, color= CAR_Exp, shape = CD)) + geom_point(alpha=0.6) + viridis::scale_color_viridis()  + theme_bw() + labs(subtitle = 'CAR expression')+ theme(legend.position='bottom', legend.key.width= unit(4, 'mm')) ,
ncol=2)

ggplot(coords, aes(x=UMAP_1, y=UMAP_2, color= BinScore, shape = CD)) + geom_point(alpha=0.6) + scale_color_manual(values=c('#73C272', '#3F5588'))  + theme_bw()+ labs(subtitle = 'High-low distribution')
ggplot(coords, aes(x=UMAP_1, y=UMAP_2, color= Cluster)) + geom_point(alpha=0.8) + hues::scale_color_iwanthue()  + theme_bw() + facet_wrap(~BinScore) + theme(legend.position='bottom')+ labs(subtitle = 'High-low distribution per cluster')

dev.off()

tmp <- setNames(as.data.frame(rds_test$donor), 'Donor')
tmp$cell_id <- rownames(tmp)
coords <- merge(coords, tmp, by='cell_id')

stats <- table(coords$Cluster, coords$BinScore)
stats_percent <- t(apply(stats, 1, FUN=function(x) round(x/sum(x) *100, digits =3)))
colnames(stats_percent) <- paste0(colnames(stats_percent), '_percent')
stats <- cbind(stats, stats_percent)
pdf('./Plots/HighLow_per_cluster.pdf')
print(gridExtra::grid.table(stats))
dev.off()


normData <- as.data.frame(t(rds_test@assays$SCT@scale.data))
get_markers <- function(data, genelist){
    tmp <- data[, genelist]
    tmp$cell_id <- rownames(tmp)
    return(tmp)
}

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

get_umap <- function(data, gene){
    gene_expr <- normData[, gene, drop=FALSE]
    gene_expr$cell_id <- rownames(gene_expr)
    coords_markers <- merge(coords, gene_expr, by='cell_id')
    plot <- ggplot(coords_markers, aes(x=UMAP_1, y=UMAP_2, color= get(gene))) + 
    geom_point(alpha=0.3, size = 0.8) + 
    scale_color_gradient(low="grey90", high ="blue", name = 'Expression') + theme_void() + 
    theme(legend.position='right', plot.title = element_text(hjust = 0.5), , legend.key.height = unit(5, 'mm'), legend.key.width = unit(2, 'mm')) + 
    ggtitle(gene)
    return(plot)
}


genelist <- c('CAR-pCCL-BCMA', 'CD4', 'CD8A','MYC', 'CD28', 'IL7R')
genelist %in% colnames(normData)

pdf('./Plots/Umap_Markers.pdf')
# legend <- cowplot::get_legend(get_umap(coords_markers, 'CD4')+ theme(legend.position='right'))
cowplot::plot_grid(
    get_umap(coords, 'CAR-pCCL-BCMA'),
    get_umap(coords, 'CD4'),
    get_umap(coords, 'CD8A'),
    get_umap(coords, 'MYC'),
    get_umap(coords, 'CD28'),
    get_umap(coords, 'IL7R'),
    ggplot(coords, aes(x=UMAP_1, y=UMAP_2, color= Cluster)) + geom_point(alpha=0.6, size = 0.5) + 
    hues::scale_color_iwanthue() + theme_void() + theme(legend.position='none', plot.title = element_text(hjust = 0.5)) + ggtitle('Populations'),
    # legend,
    ncol=2
)

dev.off()
