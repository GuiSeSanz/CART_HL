library(ggplot2)
library(Seurat)

get_annotation <- function(data){
    annotation <- data.frame(FACS_Level= stringr::str_extract(colnames(data), '(?<=_)[a-zA-Z]+(?=$|_)') )
    rownames(annotation) <- colnames(data)
    return(annotation)
}

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
    # set.seed(123)
    # pdf('./Plots/Hist_signatures.pdf')
    # cowplot::plot_grid(
    #     ggplot2::ggplot(all_signatures , ggplot2::aes(x=High_pondered)) + ggplot2::geom_histogram (bins=50, fill = '#7F7F7F', alpha = 0.6)+ geom_vline(xintercept = 0, 
    #             color = "red") + theme_classic(), 
    #     ggplot2::ggplot(all_signatures , ggplot2::aes(x=Low_pondered)) + ggplot2::geom_histogram (bins=50, fill = '#DBBE78', alpha = 0.6)+ geom_vline(xintercept = 0, 
    #             color = "red", size=1)+ theme_classic(), 
    #     ggplot2::ggplot(all_signatures , ggplot2::aes(x=scale(High_pondered))) + ggplot2::geom_histogram (bins=50, fill = '#7F7F7F', alpha = 0.6)+ geom_vline(xintercept = 0, 
    #             color = "red", size=1)+ theme_classic(), 
    #     ggplot2::ggplot(all_signatures , ggplot2::aes(x=scale(Low_pondered))) + ggplot2::geom_histogram (bins=50, fill = '#DBBE78', alpha = 0.6)+ geom_vline(xintercept = 0, 
    #             color = "red", size=1)+ theme_classic(), 
    # ncol=2)
    # dev.off()

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

get_markers <- function(data, genelist){
    tmp <- data[, genelist]
    tmp$cell_id <- rownames(tmp)
    return(tmp)
}

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

get_umap <- function(data, gene){
    gene_expr <- data[, gene, drop=FALSE]
    gene_expr$cell_id <- sub('-','.', rownames(gene_expr))
    coords_markers <- merge(coords, gene_expr, by='cell_id')
    plot <- ggplot(coords_markers, aes(x=UMAP_1, y=UMAP_2, color= get(gene))) + 
    geom_point(alpha=0.3, size = 0.8) + 
    # viridis::scale_color_viridis()+
    scale_color_gradient(low="grey90", high ="blue", name = 'Expression') + 
    guides(color = guide_colourbar(barwidth = 0.5, barheight = 2, label = TRUE)) + 
    theme_void() + ggtitle(gene)+
    theme(legend.position='none', plot.title = element_text(hjust = 0.5, size = 8))#, legend.key.height = unit(2, 'mm'), legend.key.width = unit(1, 'mm')) + 
    return(plot)
}

get_violin <- function(data, title){
    tmp <- ggplot(data, aes(x=Cluster, y=Value)) + geom_violin(aes(fill=Cluster)) + hues::scale_fill_iwanthue() + theme_classic() + ggtitle(title) + theme(legend.position='none') + labs(y = 'Value')
    return(tmp)
}

get_normalized_bulk <- function(path){
    cd_type <- stringr::str_extract(path, '(CD[\\d]{1})(?=\\.RData)')
    load(file = path)
    countData<-as.matrix(MyData_filt[, grep("_High|_Low", colnames(MyData_filt))])
    rownames(countData)<-row.names(MyData_filt)
    countData[!is.finite(countData)] <- 0
    class(countData)<-"integer"
    # colnames(countData)
    # Coldata table
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
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~Donor + CellType)
    dds <- DESeq2::estimateSizeFactors(dds)
    normalized_counts <- DESeq2::counts(dds, normalized=TRUE)
    normalized_counts <- as.data.frame(normalized_counts)
    return(normalized_counts)
}

get_heatmap <- function(data, ann_legend= TRUE){
    ann_color = list(Density = c(High='#7F7F7F', Low='#DBBE78'))
    annotation <- setNames(as.data.frame(ifelse(grepl('High',colnames(data)), 'High', 'Low')), 'Density')
    rownames(annotation) <- colnames(data)
    new_col_order <- c(base::grep('Low', colnames(data), value=TRUE), base::grep('High', colnames(data), value=TRUE))
    htmap <- pheatmap::pheatmap(data[, new_col_order], scale='row', fontsize=5, angle_col =45, show_rownames=FALSE, show_colnames = FALSE, treeheight_row = 5, treeheight_col = 5, cluster_col = FALSE, legend=TRUE, silent = TRUE, color = viridis::viridis(50), annotation_col=annotation, annotation_color = ann_color, annotation_legend = ann_legend, fontsize_row=2) #cellheight = 0.5, cellwidth = 8
    return(htmap)
}

get_Score_bulk <- function(norm_data, type){
    metadata_bulk_CIMA <- data.frame(Sample = colnames(norm_data), type = type)
    results_bulk_CIMA <-  apply_signature(norm_data, CD4_signature, CD4_signature_genes, metadata_bulk_CIMA)
    results_bulk_CIMA$Sample <- results_bulk_CIMA$Sample$Sample
    results_bulk_CIMA$type <- type
    results_bulk_CIMA$BinScore <- ifelse(results_bulk_CIMA$High_pondered >=0, 'High','Low')
    results_bulk_CIMA$NamedScore <- ifelse(stringr::str_detect(as.character(results_bulk_CIMA$Sample), 'High'), 'High','Low')
    return(results_bulk_CIMA)
}


get_correlation <- function(data, name, y_just = 35000){
    tmp <- ggplot(data, aes(x=get(name), y=FACS_value)) + geom_point(aes(color=NamedScore)) + scale_color_manual(values = c('#30A3CC', '#FCB357')) + theme_bw() + theme(legend.position='bottom') + geom_smooth(method=lm, se=FALSE, color="black", size = 0.5) + ggpubr::stat_regline_equation(label.y = y_just, aes(label=..rr.label..))
    return(tmp)
}

get_LOO_table <- function(data){
    rownames(data) <- ifelse(grepl('^Control', data$Sample), 'Control', paste0('Leaved_out_',data$Sample))
    data <- data[, !colnames(data) %in% c('Sample', 'type')]
    annot_colors <- list(FACS_Level = c(High='#30A3CC', Low='#FCB357'))
    annot <- get_annotation(data)
    tmp <- pheatmap::pheatmap(data, cluster_cols=TRUE, cluster_rows=FALSE, 
    fontsize=5, angle_col =45, show_rownames=TRUE, show_colnames = FALSE, treeheight_row = 5, treeheight_col = 5,
    scale='none', color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name ="RdYlBu")))(100), 
    annotation_col = annot, annotation_colors = annot_colors, cutree_cols = 2, display_numbers=TRUE, main='Raw score', 
    legend=FALSE, silent=TRUE)$gtable
    return(tmp)
}


check_genes <- function(gene_list, df){
    not_found_genes <- gene_list[which(!gene_list %in% colnames(normData))]
    if(length(not_found_genes)!=0){
        message(paste0('Left behind ', length(not_found_genes), ': ', paste0(not_found_genes, collapse='; ')))
        gene_list <- gene_list[which(sign_genes %in% colnames(normData))]
    }else{
        gene_list <- gene_list
    }

    return(gene_list)
}

get_umap_signature <- function(df, title){
    plot <- ggplot(df, aes(x=UMAP_1, y=UMAP_2, color= Signature)) + 
    geom_point(alpha=0.9, size = 0.8) + 
    viridis::scale_color_viridis()+
    # scale_color_gradient(low="grey90", high ="blue", name = 'Expression') + 
    guides(color = guide_colourbar(barwidth = 0.5, barheight = 2, label = TRUE)) + 
    theme_void() + ggtitle(title)+
    theme(legend.position='none', plot.title = element_text(hjust = 0.5, size = 8))#, legend.key.height = unit(2, 'mm'), legend.key.width = unit(1, 'mm')) + 
    return(plot)
}

get_expression_signature <- function(gene_list_name, df, coords){
    gene_list <- as.character(read.table(paste0(signatures_path, '/', gene_list_name))$V1)
    gene_list <- check_genes(gene_list, df)
    tmp <- df[, colnames(df) %in% gene_list]
    tmp$Signature <- rowMeans(tmp)
    tmp$cell_id <- sub('-', '\\.',rownames(tmp))
    tmp <- merge(tmp, coords, by='cell_id')
    title <- gsub('_', ' ', stringr::str_remove(gene_list_name, '\\.txt$'))
    plot <- get_umap_signature(tmp, title)
    return(plot)
}


########################
#### Single cell by CIMA
########################

# Data to calculate the signature

CD4_signature_genes <-read.delim(file="./Data/signature/BatchK_CD4_SIGGenes_BASAL_BOTH_LowvsHigh.tsv",sep="\t",header=T)
CD8_signature_genes <-read.delim(file="./Data/signature/BatchK_CD8_SIGGenes_BASAL_BOTH_LowvsHigh.tsv",sep="\t",header=T)
CD4_signature <- read.delim(file="./Data/signature/BatchK_CD4_output_BASAL_BOTH_LowvsHigh.tsv",sep="\t",header=T)
CD8_signature <- read.delim(file="./Data/signature/BatchK_CD8_output_BASAL_BOTH_LowvsHigh.tsv",sep="\t",header=T)



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

car_exp_level_SC_FP <- setNames(
                        rbind(results_CD4_SC_FP[, c('Sample', 'cell_id', 'High_pondered', 'donor', 'Overall_score', 'CAR-pCCL-BCMA', 'CD')], 
                              results_CD8_SC_FP[, c('Sample', 'cell_id', 'High_pondered', 'donor', 'Overall_score', 'CAR-pCCL-BCMA', 'CD')] ), 
                       c('Sample', 'cell_id', 'High_pondered', 'donor','Overall_score','CAR_Exp', 'CD'))


# Read the RDS to get the UMAP coords
rds_test <- readRDS('/home/sevastopol/data/mcallejac/JuanRo_SimiC/data/CART_HIGHLOW/Scores_Improved_Apr/HighLowCod_ctrl_integrated_seurat_cd4cd8_clusters.rds')
coords <- as.data.frame(rds_test@reductions$umap@cell.embeddings)
coords$cell_id <- sub('-', '.',rownames(coords))
clusters <- setNames(as.data.frame(rds_test$ClusterNames_0.8_by_JR), 'Cluster')
clusters$cell_id <- sub('-', '.',rownames(clusters))

coords <- merge(coords, car_exp_level_SC_FP[, c('cell_id', 'High_pondered', 'CAR_Exp', 'CD')], by='cell_id')
coords$BinScore <- ifelse(coords$High_pondered > 0, 'High', 'Low')
coords <- merge(coords, clusters, by='cell_id')
coords$Cluster <- stringr::str_remove(coords$Cluster, '^C?[\\d]{1,2}\\.')

normData <- as.data.frame(t(rds_test@assays$RNA@scale.data))


color_list_populations <- c('#A55E34', '#C6B2D3', '#D0342B', '#8E221A', '#2E6B34', '#BBDE93', '#AECDE1', '#3C77AF', '#ED9E9B', '#DA913D', '#821851', '#643F95', '#DBBE78', '#7F7F7F', '#000000')


Idents(rds_test) <- rds_test$seurat_clusters
all_genes <- rownames(rds_test)
rds_test <-  ScaleData(rds_test, features = all_genes)
rds_test <-  FindVariableFeatures(rds_test, selection.method = "vst", nfeatures = 2000)

Idents(rds_test) <- rds_test$seurat_clusters
rds_test$seurat_clusters <- factor(rds_test$seurat_clusters, levels=sort(unique(rds_test$seurat_clusters)))

markers <- FindAllMarkers(rds_test, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

marker_list <- list()
for (clust in unique(markers$cluster)){
    cluster_name <- grep(paste0('(^C?)',clust, '\\.'), unique(rds_test$ClusterNames_0.8_by_JR), value=TRUE)
    marker_list[[cluster_name]] <- markers[markers$cluster == clust,]
}
# WriteXLS::WriteXLS(marker_list, ExcelFileName= './Plots/Markers_cluster_scRNA_HL.xlsx', SheetNames=names(marker_list),  row.names=TRUE, BoldHeaderRow=TRUE)
library(dplyr)
all_markers_top20 <- markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)


# pdf('./Plots/Test.pdf')
#     DoHeatmap(rds_test, group.by='seurat_clusters', features=all_markers_top20$gene, angle = 0, group.colors = hues::iwanthue(length(levels(all_markers_top20$cluster))), size = 3, combine = TRUE) + NoLegend() + theme(axis.text.y = element_text(size = 3))
# dev.off()

htmap <- DoHeatmap(rds_test, group.by='seurat_clusters', features=all_markers_top20$gene, angle = 0, group.colors = hues::iwanthue(length(levels(all_markers_top20$cluster))), size = 3, combine = TRUE) + NoLegend() + theme(axis.text.y = element_text(size = 3))


c('CD4', 'CD8A', 'SELL', 'CCR7', 'TCF7', 'GZMA', 'GNLY', 'NKG7', 'HLA-DR','IL2RA','TIGIT','LAG3') %in% colnames(normData)
pdf('./Plots/Figure3.pdf')
    # legend <- cowplot::get_legend(get_umap(coords_markers, 'CD4')+ theme(legend.position='right'))
    cowplot::plot_grid(
        cowplot::plot_grid(
        ggplot(coords, aes(x=UMAP_1, y=UMAP_2, color= Cluster)) + geom_point(alpha=0.6, size = 0.5) + scale_color_manual(values=color_list_populations) + theme_void() + theme(legend.position='right', plot.title = element_text(hjust = 0.5)) + ggtitle('Populations') + guides(color = guide_legend(override.aes = list(size=3, alpha = 1))),
        cowplot::plot_grid(
            get_umap(normData, 'CD4'),
            get_umap(normData, 'CD8A'),
            get_umap(normData, 'SELL'),
            get_umap(normData, 'CCR7'),
            get_umap(normData, 'TCF7'),
            get_umap(normData, 'GZMA'),
            get_umap(normData, 'GNLY'),
            get_umap(normData, 'NKG7'),
            get_umap(normData, 'HLA-DRA'),
            get_umap(normData, 'IL2RA'),
            get_umap(normData, 'TIGIT'),
            get_umap(normData, 'LAG3'),
            # legend,
            ncol=4),
        nrow=2),

        cowplot::plot_grid(htmap, ncol=1),
    ncol=2, rel_widths = c(2, 1)
    )
dev.off()


pdf('./Plots/FigureS5.pdf')
    clusters_tmp  <- setNames(as.data.frame(rds_test$seurat_clusters) , 'Cluster')
    clusters_tmp$cell_id <- rownames(clusters_tmp)
    g2mScore  <- setNames(as.data.frame(rds_test$G2M.Score) , 'Value')
    g2mScore <- merge(g2mScore, clusters_tmp, by.x=0, by.y='cell_id')
    sScore    <- setNames(as.data.frame(rds_test$S.Score) , 'Value')
    sScore <- merge(sScore, clusters_tmp, by.x=0, by.y='cell_id')
    mitoRatio <- setNames(as.data.frame(rds_test$mitoRatio) , 'Value')
    mitoRatio <- merge(mitoRatio, clusters_tmp, by.x=0, by.y='cell_id')
    g2mScore  <- get_violin(g2mScore, 'G2M phase')
    sScore    <- get_violin(sScore, 'S phase')
    mitoRatio <- get_violin(mitoRatio, 'Ratio of mitochondrial genes')
    contribution <- merge(clusters_tmp, setNames(as.data.frame(rds_test$donor) , 'Donor'), by.x='cell_id', by.y=0)
    contribution2 <- setNames(as.data.frame(prop.table(table(contribution[contribution$Donor == 'd10', 'Cluster']))*100), c('Cluster', 'd10'))
    contribution2 <- merge(contribution2, setNames(as.data.frame(prop.table(table(contribution[contribution$Donor == 'd14', 'Cluster']))*100), c('Cluster', 'd14')), by='Cluster')
    contribution2 <-merge(contribution2,  setNames(as.data.frame(prop.table(table(contribution[contribution$Donor == 'd18', 'Cluster']))*100), c('Cluster', 'd18')), by='Cluster')
    contribution2 <- setNames(reshape2::melt(contribution2), c('Cluster', 'Donor', 'value'))
    contribution <- ggplot(contribution, aes(x= Cluster, fill=Donor))+  geom_bar(position="fill", colour="black") + theme_classic() + labs(title = "Contribution by donor", x= 'Cluster', y = '% of cells') + scale_fill_brewer(palette="RdBu")  #scale_fill_manual(values=c('#585123', '#f2a65a', '#772f1a'))
    contribution2 <- ggplot(contribution2, aes(x= Cluster, y=value, fill=Donor))+  geom_bar(position="fill", stat='identity', colour="black") + theme_classic() + labs(title = "Contribution by donor", x= 'Cluster', y = '% of cells') + scale_fill_brewer(palette="RdBu")  #scale_fill_manual(values=c('#585123', '#f2a65a', '#772f1a'))

    tmp <- setNames(as.data.frame(rds_test$donor), 'Donor')
    tmp$cell_id <- sub('-', '.', rownames(tmp))
    tmp <- merge(coords, tmp, by='cell_id')
    cowplot::plot_grid(   
        cowplot::plot_grid(
            ggplot(coords, aes(x=UMAP_1, y=UMAP_2, color= Cluster)) + geom_point(alpha=0.9, size = 0.5) + scale_color_manual(values=color_list_populations) + theme_void() + theme(legend.position='right', plot.title = element_text(hjust = 0.5)) + ggtitle('Populations') + guides(color = guide_legend(override.aes = list(size=3, alpha =1))), 
            # cowplot::plot_grid(
            #     ggplot(tmp, aes(x=UMAP_1, y=UMAP_2, color= Cluster)) + geom_point(alpha=0.9, size = 0.5) + scale_color_manual(values=color_list_populations) + theme_void() + theme(legend.position='none', plot.title = element_text(hjust = 0.5)) + ggtitle('Populations by donor') + guides(color = guide_legend(override.aes = list(size=5)))  + facet_wrap(~Donor)
            # ), 
            # contribution,
            contribution2,
            nrow=2, rel_heights=c(2,1)
        ),
        cowplot::plot_grid(
            g2mScore,
            sScore,
            mitoRatio,
            nrow=3)
    , ncol=2, rel_widths=c(3,2))
        
dev.off()


bin_high_low <- car_exp_level_SC_FP[, c('cell_id', 'High_pondered')]
bin_high_low$cell_id <- sub('\\.', '-', bin_high_low$cell_id)
bin_high_low$BinScore <- ifelse(bin_high_low$High_pondered >=0, 'High', 'Low')
BinScore <- bin_high_low$BinScore
names(BinScore) <- bin_high_low$cell_id
rds_test <- AddMetaData(rds_test, BinScore, col.name = 'BinScore')
Idents(rds_test) <- rds_test$BinScore

High_Low_DE_markers <- FindMarkers(rds_test, ident.1 = "High", ident.2 = "Low")
High_Low_DE_markers$gene <- rownames(High_Low_DE_markers)
# WriteXLS::WriteXLS(High_Low_DE_markers, ExcelFileName= './Plots/Markers_scRNA_HL.xlsx', SheetNames='High_Low',  row.names=TRUE, BoldHeaderRow=TRUE)
htmap <- DoHeatmap(rds_test, group.by='BinScore', features=High_Low_DE_markers$gene, angle = 0, group.colors = c('#30A3CC', '#FCB357'), size = 3, combine = TRUE) + NoLegend() + theme(axis.text.y = element_text(size = 3))

signatures_path <- '/home/sevastopol/data/gserranos/CART_HL/Data/signature/OtherSignatures'
signatures <- list.files(signatures_path)


pdf('./Plots/Figure4.pdf', width=7.2, height=10)
HighAndLow = FALSE
if(HighAndLow){
    clusters_tmp  <- setNames(as.data.frame(rds_test$seurat_clusters) , 'Cluster')
    clusters_tmp$cell_id <- sub('-','.',rownames(clusters_tmp))
    composition <- merge(coords[, c('cell_id', 'BinScore')], clusters_tmp, by='cell_id')
    composition <- ggplot(composition, aes(x= Cluster, fill=BinScore))+  geom_bar(position="fill") + theme_classic() + labs(title = "High CART distribution", x= 'Cluster', y = '% of cells') + scale_fill_manual(values=c('#30A3CC', '#FCB357')) 
}else{
    clusters_tmp  <- setNames(as.data.frame(rds_test$seurat_clusters) , 'Cluster')
    clusters_tmp$cell_id <- sub('-','.',rownames(clusters_tmp))
    composition <- merge(coords[, c('cell_id', 'BinScore')], clusters_tmp, by='cell_id')
    composition$Cluster <- as.character(composition$Cluster)
    composition <- table(composition$Cluster, composition$BinScore)
    composition <- as.data.frame.matrix(composition)
    composition$High_prop <- apply(composition, 1, FUN=function(x) (x[1]/sum(x))*100 )
    composition$Cluster <- factor(rownames(composition) , levels=c(sort(as.numeric(rownames(composition)))))
    composition <- ggplot(composition, aes(x=Cluster, y=High_prop))+  geom_bar(stat='identity', fill = "#30A3CC") + theme_classic() + labs(title = "High CART distribution", x= 'Cluster', y = '% of cells') + theme(panel.grid.major = element_line(colour="#f0f0f0"))
}
cowplot::plot_grid(
    cowplot::plot_grid(
        ggplot(coords, aes(x=UMAP_1, y=UMAP_2, color= BinScore, shape = CD)) + geom_point(alpha=0.6) + scale_color_manual(values=c('#30A3CC', '#FCB357'))  + theme_void()+ labs(subtitle = 'High-low distribution') + guides(color = guide_legend(override.aes = list(size=3, alpha = 1))),
        composition,
        cowplot::plot_grid(
            get_expression_signature("Genes_Activation.txt", normData, coords ),
            get_expression_signature("Genes_HLA.txt", normData, coords ),
            get_expression_signature("Genes_Polyfun.txt", normData, coords ),
            get_expression_signature("Genes_Prolif.txt", normData, coords ),
            get_expression_signature("Genes_Teff.txt", normData, coords ),
            get_expression_signature("Genes_Th17.txt", normData, coords ),
            get_expression_signature("Genes_Th1.txt", normData, coords ),
            get_expression_signature("Genes_Th2.txt", normData, coords ),
            get_expression_signature("Genes_Tonic.txt", normData, coords ),
            ncol=3
        ),
    ncol=1, rel_heights=c(2,1,2)), 
htmap,
ncol=2, rel_widths=c(2,1))
dev.off()





normalized_counts_CD4 <- get_normalized_bulk("/home/sevastopol/data/mcallejac/RNA_HighLow_ALL/results/HL_CD4.RData")
normalized_counts_CD8 <- get_normalized_bulk("/home/sevastopol/data/mcallejac/RNA_HighLow_ALL/results/HL_CD8.RData")



results_CD8_bulk_CIMA <- get_Score_bulk(normalized_counts_CD4, 'CD4')
results_CD4_bulk_CIMA <- get_Score_bulk(normalized_counts_CD8, 'CD8')
results_bulk_CIMA <- rbind(results_CD8_bulk_CIMA, results_CD4_bulk_CIMA)

CAR_NAME <- 'CAR_pCCL_BCMA'
car_expr <- as.data.frame(rbind(t(normalized_counts_CD8[CAR_NAME, , drop=FALSE]), t(normalized_counts_CD4[CAR_NAME, , drop=FALSE])))
car_expr$Sample <- rownames(car_expr)

results_bulk_CIMA <- merge(results_bulk_CIMA, car_expr, by='Sample')

results_bulk_CIMA$NamedScore <- factor(results_bulk_CIMA$NamedScore , levels=c('Low', 'High'))

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

# FACs_Vs_Expr <- ggplot(FACs, aes(x=CAR_pCCL_BCMA, y=FACS_value, color=NamedScore)) + geom_point() + scale_color_manual(values = c('#DBBE78', '#7F7F7F')) + theme_bw() + theme(legend.position='bottom') +  facet_wrap(~type, nrow=2)
# FACs_Vs_Score <- ggplot(FACs, aes(x=High_pondered, y=FACS_value, color=NamedScore)) + geom_point() + scale_color_manual(values = c('#DBBE78', '#7F7F7F')) + theme_bw() + theme(legend.position='bottom') + facet_wrap(~type, nrow=2)




FACs_Vs_Expr_CD4  <- get_correlation(FACs[FACs$type == 'CD4',], 'CAR_pCCL_BCMA') + labs(x = 'CAR expression')
FACs_Vs_Expr_CD8  <- get_correlation(FACs[FACs$type == 'CD8',], 'CAR_pCCL_BCMA') + labs(x = 'CAR expression')
FACs_Vs_Score_CD4 <- get_correlation(FACs[FACs$type == 'CD4',], 'High_pondered') + labs(x = 'Signature Score')
FACs_Vs_Score_CD8 <- get_correlation(FACs[FACs$type == 'CD8',], 'High_pondered') + labs(x = 'Signature Score')
CD4_htm <- get_heatmap(normalized_counts_CD4[rownames(normalized_counts_CD4) %in% CD4_signature_genes$sigGenes_symbol,])
CD8_htm <- get_heatmap(normalized_counts_CD8[rownames(normalized_counts_CD8) %in% CD8_signature_genes$sigGenes_symbol,])

Results_LOO_CD8 <- readRDS('./Data/signature/Results_LOO_CD8.rds')
Results_LOO_CD4 <- readRDS('./Data/signature/Results_LOO_CD4.rds')



LOO_CD4 <- get_LOO_table(Results_LOO_CD4)
LOO_CD8 <- get_LOO_table(Results_LOO_CD8)




pdf('./Plots/FigureS6.pdf')
    legend <- cowplot::get_legend(FACs_Vs_Score_CD4+ theme(legend.box.margin = margin(0, 0, 0, 0, 'mm')))
    cowplot::plot_grid(
        cowplot::plot_grid(
            FACs_Vs_Expr_CD4 + theme(legend.position='none'),
            FACs_Vs_Score_CD4 + theme(legend.position='none'),
            legend,
            LOO_CD4,
        nrow=4, rel_heights =c(2,2,0.2,3)),
        CD4_htm$gtable,
    ncol=2, rel_widths=c(2,1))

    legend <- cowplot::get_legend(FACs_Vs_Score_CD8+ theme(legend.box.margin = margin(0, 0, 0, 0, 'mm')))
    cowplot::plot_grid(
        cowplot::plot_grid(
            FACs_Vs_Expr_CD8 + theme(legend.position='none'),
            FACs_Vs_Score_CD8 + theme(legend.position='none'),
            legend,
            LOO_CD8,
        nrow=4, rel_heights =c(2,2,0.2,3)),
        CD8_htm$gtable,
    ncol=2, rel_widths=c(2,1))

dev.off()
