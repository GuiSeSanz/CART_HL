
library(rlang)
library(ggplot2)
library(Seurat)
library(Gviz)
library(rtracklayer)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
require(tidyverse)
require(biomaRt)
require(org.Hs.eg.db)

get_colorStats <- function(data, num_breaks){
    mat <- as.matrix(data)
    mat_breaks <- seq(min(mat), max(mat), length.out = 10)
    dat <- data.frame(values=unlist(data))
    dat_colors <- data.frame(
        xmin = mat_breaks[1:(length(mat_breaks)-1)],
        xmax = mat_breaks[2:length(mat_breaks)],
        ymin = 0,
        ymax = max(density(as.matrix(mat), bw = "SJ")$y),
        fill = rev(viridis::inferno(length(mat_breaks) - 1)),
        stringsAsFactors = FALSE
    )
    dat2 <- as.data.frame(table(cut(
    mat, mat_breaks
    )))
    dat2$fill <- viridis::inferno(nrow(dat2))

    quantile_breaks <- function(xs, n = 10) {
        breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
        breaks[!duplicated(breaks)]
    }
    mat_breaks_q <- quantile_breaks(mat, n = num_breaks)
    dat2_q <- as.data.frame(table(cut(
    mat, mat_breaks_q
    )))
    dat2_q$fill <- viridis::inferno(nrow(dat2_q))
    dat_colors_q <- data.frame(
        xmin = mat_breaks_q[1:(length(mat_breaks_q)-1)],
        xmax = mat_breaks_q[2:length(mat_breaks_q)],
        ymin = 0,
        ymax = max(density(mat, bw = "SJ")$y),
        fill = rev(viridis::inferno(length(mat_breaks_q) - 1)),
        stringsAsFactors = FALSE
    )

    pdf('./Plots/ColorStats.pdf')
    print(cowplot::plot_grid(
        ggplot() + geom_rect(
            data = dat_colors,
            mapping = aes(
            xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill
            ))
         + geom_density(
            data = dat,
            mapping = aes(values),
            bw = "SJ", color = "cyan"
        ) +
        theme_classic() +
        scale_fill_manual(values = dat_colors$fill) +
        theme(legend.position = "none") +
        labs(title = "Uniform breaks")
        ,
        ggplot() +
        geom_bar(
            data = dat2,
            mapping = aes(x = Var1, weight = Freq, fill = Var1),
            color = "black", size = 0.1
        ) +
        coord_flip() + theme_classic() +
        scale_fill_manual(values = dat2$fill) +
        theme(legend.position = "none") +
        labs(y = "data points", x = "breaks",
            title = "Number of data points per color"),
    ncol=1))

    print(cowplot::plot_grid(

        ggplot() +
        geom_rect(
            data = dat_colors_q,
            mapping = aes(
            xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill
            )) +
        geom_density(
            data = dat,
            mapping = aes(values),
            bw = "SJ", color = "cyan"
        ) + theme_classic() +
        scale_fill_manual(values = dat_colors_q$fill) +
        theme(legend.position = "none") +
        labs(title = "Quantile breaks")
        ,
        ggplot() +
        geom_bar(
            data = dat2_q,
            mapping = aes(x = Var1, weight = Freq, fill = Var1),
            color = "black", size = 0.1
        ) +
        coord_flip() + theme_classic() +
        scale_fill_manual(values = dat2_q$fill) +
        theme(legend.position = "none") +
        labs(y = "data points", x = "breaks",
            title = "Number of data points per color"),
        nrow=2))
        dev.off()
    return(mat_breaks_q)
}

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

get_umap <- function(data, gene, upplim=NULL ){
    gene_expr <- data[, gene, drop=FALSE]
    if(!is.null(upplim)){
        gene_expr[gene_expr>upplim] <- upplim
    }
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

get_violin <- function(data, title, cluster_color = color_list_populations){
    tmp <- ggplot(data, aes(x=Cluster, y=Value)) + geom_violin(aes(fill=Cluster)) + 
    theme_classic() + ggtitle(title) + theme(legend.position='none', axis.text=element_text(size=6)) + labs(y = 'Value') +
    scale_fill_manual(values=cluster_color)
    # hues::scale_fill_iwanthue() 
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

get_correlation <- function(data, name, y_just = log10(1e5)){ # y_just = 3500
    tmp <- ggplot(data, aes(x=get(name), y=FACS_value)) + geom_point(aes(color=NamedScore)) + scale_y_continuous(trans = 'log10')+ 
    scale_color_manual(values = c('#30A3CC', '#FCB357')) + theme_bw() + theme(legend.position='bottom') + 
    geom_smooth(method=lm, se=FALSE, color="black", size = 0.5) + 
    ggpubr::stat_regline_equation(label.y = y_just, aes(label=..rr.label..))
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
        gene_list <- gene_list[which(gene_list %in% colnames(normData))]
    }else{
        gene_list <- gene_list
    }

    return(gene_list)
}

get_umap_signature <- function(df, title){
    plot <- ggplot(df, aes(x=UMAP_1, y=UMAP_2, color= Signature)) + 
    geom_point(alpha=0.9, size = 0.8) + 
    viridis::scale_color_viridis()+
    # scale_color_gradientn(low="grey90", high ="blue", name = 'Expression') + 
    # guides(color = guide_colourbar(barwidth = 0.5, barheight = 2, label = TRUE)) + 
    ggprism::theme_prism() + ggtitle(title)+  labs(x = "UMAP 1", y = 'UMAP 2') +
    theme(legend.position='none', plot.title = element_text(hjust = 0.5, size = 12),  axis.text=element_blank(), axis.ticks=element_blank())#, legend.key.height = unit(2, 'mm'), legend.key.width = unit(1, 'mm')) + 
    return(plot)
}

get_expression_signature <- function(gene_list_name, df, coords, upplim=NULL, split_HL=FALSE){
    gene_list <- as.character(read.table(paste0(signatures_path, '/', gene_list_name))$V1)
    gene_list <- check_genes(gene_list, df)
    tmp <- df[, colnames(df) %in% gene_list]
    tmp$Signature <- rowMeans(tmp)
    tmp$cell_id <- sub('-', '\\.',rownames(tmp))
    tmp <- merge(tmp, coords, by='cell_id')
    title <- gsub('_', ' ', stringr::str_remove(gene_list_name, '\\.txt$'))
    if(!is.null(upplim)){
        upp_outliers <- which(tmp$Signature>upplim) 
        tmp[upp_outliers , 'Signature'] <- upplim
    }
    plot <- get_umap_signature(tmp, title)
    if(split_HL){
        plot <- plot + facet_wrap(~BinScore)
    }
    return(plot)
}

get_heatmap_atac <- function(plotter, title, col_order_as_is=TRUE){
    ann_color = list(Density = c(High = '#30A3CC', Low = '#FCB357'))
    annotation <- setNames(as.data.frame(ifelse(grepl('High',colnames(plotter)), 'High', 'Low')), 'Density')
    rownames(annotation) <- colnames(plotter)
    annotation_row <- setNames(as.data.frame(stringr::str_extract(rownames(plotter), '[\\w]+(?=_)')), 'Gene')
    rownames(annotation_row) <- rownames(plotter)
    color_genes <- hues::iwanthue(length(unique(annotation_row$Gene)))
    names(color_genes) <- unique(annotation_row$Gene)
    ann_color[['Gene']] <- color_genes
    if(!col_order_as_is){
        order_low <- plotter[, grep('Low', colnames(plotter))]
        order_low <- hclust(dist(t(order_low)))$order
        order_high <- plotter[, grep('High', colnames(plotter))]
        order_high <- hclust(dist(t(order_high)))$order
        order_high <- rev(order_high)
        names <- c(grep('Low', colnames(plotter), value=TRUE)[order_low], grep('High', colnames(plotter), value=TRUE)[order_high])
        plotter <- plotter[, names]
    }
    phm <- pheatmap::pheatmap(plotter, scale='row', fontsize=5, angle_col =45, show_rownames=FALSE, show_colnames = TRUE, treeheight_row = 5, treeheight_col = 5, cluster_row=FALSE, cluster_col = FALSE, legend=TRUE, silent = TRUE, color = viridis::viridis(50),annotation_col=annotation, annotation_row= annotation_row, annotation_color = ann_color, annotation_legend = TRUE, fontsize_row=2, main=title)
    return(phm)
}

get_heatmap_bulk <- function(plotter, title){
    ann_color = list(Density = c(High = '#30A3CC', Low = '#FCB357'))
    annotation <- setNames(as.data.frame(ifelse(grepl('High',colnames(plotter)), 'High', 'Low')), 'Density')
    rownames(annotation) <- colnames(plotter)
    phm <- pheatmap::pheatmap(plotter, scale='row', fontsize=5, angle_col =45, show_rownames=FALSE, show_colnames = TRUE, treeheight_row = 5, treeheight_col = 5, cluster_row=TRUE, cluster_col = TRUE, legend=TRUE, silent = TRUE, color = viridis::viridis(50),annotation_col=annotation, annotation_color = ann_color, annotation_legend = TRUE, fontsize_row=2, main=title)
    return(phm)
}

get_pca <- function(dataPCA, variation, color2plot){
    p = switch(color2plot, 
    'HighLow' = ggplot(data=dataPCA, aes(x = PC1, y = PC2, color=HighLow)) + geom_point(size=7) + xlab(paste0("PC1: ", variation[1], "% variance")) + ylab(paste0("PC2: ", variation[2], "% variance")) + theme_bw() +  scale_color_manual(values=c('High'='#30A3CC', 'Low'='#FCB357')) + theme(legend.position='none', axis.text=element_blank(), axis.ticks=element_blank(), panel.grid = element_blank()) ,
    'Donor'   = ggplot(data=dataPCA, aes(x = PC1, y = PC2, color=Donor)) + geom_point(size=7) + xlab(paste0("PC1: ", variation[1], "% variance")) + ylab(paste0("PC2: ", variation[2], "% variance")) + theme_bw()  + hues::scale_color_iwanthue() + theme(legend.position='none', axis.text=element_blank(), axis.ticks=element_blank(), panel.grid = element_blank()))
    return(p)
}

get_pie_Annot <- function(data){
    pie_colors <- c('#AECDE1', '#2D82AF', '#98D176', '#6E9E4B', '#F16667', '#F06C45', '#F7982D', '#D9A294', '#7D54A5', '#F0EB99', '#B15928')
    data_df <- as.data.frame(data@annoStat)
    data_df$Feature_Freq <- paste0(data_df$Feature, '  (' , round(data_df$Frequency/sum(data_df$Frequency)*100, 2), '%)') 
    data_df$ypos <- cumsum(data_df$Frequency)- 0.5*data_df$Frequency
    names(pie_colors) <- data_df$Feature_Freq
    gg <- ggplot2::ggplot(data_df, aes(x="", y=Frequency, fill=Feature_Freq)) +
    geom_bar(stat="identity", width=1, color="black") +
    coord_polar("y", start=0) +
    theme_void() + 
    theme(legend.position="right", text=element_text(family='serif'), legend.key.size = unit(0.15, 'cm')) +
    scale_fill_manual(values = c(pie_colors), name="Feature")
    return(gg)
}

get_peaks <- function(CHR, START, END, HIGH, LOW, filter_genes=NULL,HLstart=NULL, HLend=NULL, MIN =0, MAX=3){
    if(START > END){
        print('inversing start-end')
        tmp_start <- START
        tmp_end <- END
        START <- tmp_end
        END <- tmp_start
    }
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene 
    gtrack <- GenomeAxisTrack()
    CHR <- as.character(CHR) # chromosome number
    itrack <- IdeogramTrack(genome="hg38", 
                        chromosome=paste0("chr",CHR),
                        from =START, to=END)
    itrack@chromosome <- CHR

    # remove chr from chromosome naming
    # levels(itrack@bandTable$chrom) <- sub("^chr", "", levels(itrack@bandTable$chrom), ignore.case=T)
    HIGH_NAME <- stringr::str_extract(HIGH, '(?<=\\/)[Da0-9]+')
    LOW_NAME  <- stringr::str_extract(LOW, '(?<=\\/)[Da0-9]+')
    bw_high <- import.bw(HIGH, as="GRanges")
    bw_low  <- import.bw(LOW, as="GRanges")
    # START <- START - 20000
    # END <- END + 20000
    bw_high_track <- DataTrack(range=bw_high,
                        name = paste0('High\n', HIGH_NAME),
                        chromosome=CHR,
                        from =START,
                        to=END,  col.histogram=c('#36A3CC'), 
                        ylim = c(MIN, MAX))
    bw_low_track <- DataTrack(range=bw_low,
                            name = paste0('Low\n', LOW_NAME),
                            chromosome=CHR,
                            from =START,
                            to=END,  col.histogram=c('#DBBE78'),
                            ylim = c(MIN, MAX))

    bm <- useMart( biomart = "ENSEMBL_MART_ENSEMBL", 
                dataset = "hsapiens_gene_ensembl")
    if(!is.null(filter_genes)){
        biomTrack <- BiomartGeneRegionTrack(genome = "hg38", chromosome = CHR, 
                                        start = START, end = END,
                                        filters=list(external_gene_name=filter_genes),
                                        name = "ENSEMBL", biomart = bm,col.line = NULL, col= NULL, 
                                        fontface.group = 4)

    }else{
        biomTrack <- BiomartGeneRegionTrack(genome = "hg38", chromosome = CHR, 
                                            start = START, end = END,
                                            # filters=list(external_gene_name=c("TRPM3", 'KLF9')),
                                            name = "ENSEMBL", biomart = bm,col.line = NULL, col= NULL, 
                                            fontface.group = 4)
    }

    if(!is.null(HLstart) & !is.null(HLend)){
        ht <- HighlightTrack(trackList = c(bw_high_track,bw_low_track, biomTrack), 
        start = c(70200000), 
        width =100000, 
        chromosome = CHR,
        fill="#A8A8A8", 
        col='grey')
    }else{
        filename <- paste0('./Plots/',filter_genes[1] ,CHR, '_',START,'_',END,'_tmp.pdf')
        pdf(filename, height=4, width=7)
        # c( itrack, bw_high_track, bw_low_track, biomTrack, gtrack)
        Gviz::plotTracks(c(bw_high_track, bw_low_track, biomTrack),
        sizes = c(5,5,2),
        transcriptAnnotation="symbol",
        from =START,
        to=END,
        chromosome = CHR,
        fontcolor = "black",
        col.axis="black",
        fontsize=15,
        showTitle=TRUE,
        collapseTranscripts = 'longest', window="auto", 
        type="histogram", cex.title=0.7, fontsize=10,littleTicks = TRUE, add=TRUE)
        dev.off()
    }
    return(filename)

}

get_venn <- function(data){
    venn<- ggvenn::ggvenn(data,
                fill_color = hues::iwanthue(length(data)),
                stroke_size = 0.4,
                show_percentage = F,
                fill_alpha = 0.4,
                stroke_color = 'white',
                stroke_alpha = 1,
                stroke_linetype = 'solid',
                text_color = 'black',
                set_name_size = 4, 
                text_size = 3.5,
                label_sep = ','
                )+ theme(plot.title = element_text(hjust = 0.5))
    return(venn)
}

def_getHeatmap_signature <- function(dataset, title){
    GeneSets_HL <- readRDS('./Data/GeneSets_HL.rds')
    norm_matrix <- readRDS('./Data/Norm_counts_zscaled_rlog_HL.rds')
    sigGenes_symbol <-read.delim(file="./Data/signature/BatchK_CD8_SIGGenes_BASAL_BOTH_LowvsHigh.tsv",sep="\t",header=T)$sigGenes_symbol
    topVarGenes <- GeneSets_HL[[dataset]][complete.cases(GeneSets_HL[[dataset]]),]
    #print(topVarGenes)
    mat  <- norm_matrix[rownames(norm_matrix) %in% topVarGenes,]
    anno <- data.frame(Sample = colnames(norm_matrix))
    anno$CellType <- ifelse(stringr::str_detect(anno$Sample, 'High'), 'High', 'Low')
    anno$Donor <- stringr::str_extract(anno$Sample, '^D[0-9]+')
    rownames(anno) <- anno$Sample
    anno$Sample <- NULL
    anno.row <- data.frame(Sign = ifelse(rownames(mat) %in% sigGenes_symbol, "Sig", "NS" ))
    anno.row$Sign <- as.factor(anno.row$Sign)
    rownames(anno.row) <- rownames(mat)    
    cols <- colorRampPalette(RColorBrewer::brewer.pal(4, "Accent"))
    row_mycolors <- cols(length(unique(anno.row$sig)))
    names(row_mycolors) <- unique(anno.row$sig)
    colors_ann <- list(
        'CellType' = c(High = '#30A3CC', Low = '#FCB357'),
        'Sign' = c(Sig='#EE8FF7', NS='#F0978D'),
        'Donor' = c(setNames(viridis::viridis(length(unique(anno$Donor))), unique(anno$Donor)))
    )
  
    pm <- pheatmap::pheatmap(mat, scale = "row",
           treeheight_row=15, treeheight_col=20,
           #cellheight = 2, cellwidth = 1,
           color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
           cutree_rows = 2, 
           # cutree_cols = 3,
           fontsize_row = 6,
           legend=FALSE,
           annotation_col=anno,
           annotation_colors=colors_ann,
           annotation_legend=FALSE,
           annotation_row = anno.row, 
           show_colnames = FALSE,
           silent=TRUE,
           border_color='NA',
           main=title)
    return(pm)
}


get_plot_signature <- function(signature, data, cluster, signif=TRUE){
    genes_sign <- unique(as.character(read.table(paste0(signatures_path, '/', signature))$V1))
    tmp <- data[, colnames(data) %in% genes_sign]
    tmp$Signature_Score <- rowSums(tmp)
    tmp$cell_id <- sub('-','.',rownames(tmp))
    tmp <- merge(tmp, coords, by='cell_id')
    tmp <- tmp[tmp$Cluster == cluster,]
    p <- ggplot(tmp, aes(y=Signature_Score, x=BinScore, fill=BinScore)) + geom_boxplot(alpha = 0.8) + scale_fill_manual(values=c('High'='#30A3CC', 'Low'='#d3d3d3'))+   ggprism::theme_prism()+ labs(y= 'Signature Score') + theme(legend.position='none', axis.title.x = element_blank(), axis.title.y = element_text(size=8),plot.title = element_text(hjust = 0), plot.subtitle=element_text(hjust = 0)) 
    if (signif){
        p <- p + ggsignif::geom_signif(comparisons = list(c("High", "Low")), map_signif_level = TRUE, vjust=0.5)
    }
    return(p)
}



RdWhBl <- colorRampPalette(colors = c("blue", "white", "red")) (100)
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
rds_unfiltered <- readRDS('/home/sevastopol/data/mcallejac/JuanRo_SimiC/data/CART_UFILTER/integrated_seurat_clusters_afterpipeline.rds')
rds_test <- readRDS('/home/sevastopol/data/mcallejac/JuanRo_SimiC/data/CART_HIGHLOW/Scores_Improved_Apr/HighLowCod_ctrl_integrated_seurat_cd4cd8_clusters.rds')
coords <- as.data.frame(rds_test@reductions$umap@cell.embeddings)
coords$cell_id <- sub('-', '.',rownames(coords))
clusters <- setNames(as.data.frame(rds_test$ClusterNames_0.8_by_JR), 'Cluster')
clusters$cell_id <- sub('-', '.',rownames(clusters))
levels(clusters$Cluster)[levels(clusters$Cluster)=="21.CD4 Cytotoxic"] <- "C21.CD4 Cytotoxic"
coords <- merge(coords, car_exp_level_SC_FP[, c('cell_id', 'High_pondered', 'CAR_Exp', 'CD')], by='cell_id')
coords$BinScore <- ifelse(coords$High_pondered > 0, 'High', 'Low')
coords <- merge(coords, clusters, by='cell_id')
# coords$Cluster <- stringr::str_remove(coords$Cluster, '^C?[\\d]{1,2}\\.')

normData <- as.data.frame(t(rds_test@assays$RNA@scale.data))

color_list_populations <- c('#A55E34', '#C6B2D3', '#D0342B', '#8E221A', '#2E6B34', '#BBDE93', '#AECDE1', '#3C77AF', '#ED9E9B', '#DA913D', '#821851', '#643F95', '#DBBE78', '#7F7F7F', '#000000')
names(color_list_populations) <- levels(rds_test$ClusterNames_0.8_by_JR)

# color_all_clusters <- c('#A55E34', '#C6B2D3', '#D0342B', '#8E221A', '#2E6B34', '#BBDE93', '#AECDE1', '#3C77AF', '#ED9E9B', '#DA913D', '#821851', '#643F95', '#DBBE78', '#7F7F7F', '#000000')
# names(color_all_clusters) <- 

Idents(rds_test) <- rds_test$seurat_clusters
all_genes <- rownames(rds_test)
rds_test <-  ScaleData(rds_test, features = all_genes)
rds_test <-  FindVariableFeatures(rds_test, selection.method = "vst", nfeatures = 2000)

Idents(rds_test) <- rds_test$seurat_clusters
rds_test$seurat_clusters <- factor(rds_test$seurat_clusters, levels=sort(unique(rds_test$seurat_clusters)))

if (file.exists('./Data/Markers_cluster_scRNA_HL.rds')){
    markers <- readRDS('./Data/Markers_cluster_scRNA_HL.rds')
}else{
    markers <- FindAllMarkers(rds_test, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    saveRDS(markers, './Data/Markers_cluster_scRNA_HL.rds')
}

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



############# PCA for SC and BULK #############
#### CD8
# == atac ==
load('/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/RSession/cd8_2105_Csaw_fin.RData')
counts <- edgeR::getCounts(y)
cols_2_keep <- grep('High|Low', grep('d0', colnames(y), value=TRUE), value=TRUE)
counts <- counts[, cols_2_keep]
group <-factor(stringr::str_extract(colnames(counts), '^D[0-9a]+'))
dge <- edgeR::DGEList(counts=counts, group=group)
dge <- edgeR::calcNormFactors(dge, method = "TMM")
logCPM <- edgeR::cpm(dge, log=TRUE, prior.count=2)

combat_edata1 = sva::ComBat(dat=logCPM, batch=group, mod=NULL, par.prior=TRUE, prior.plots=FALSE)
pca <- prcomp(t(combat_edata1))
print(summary(pca))
pcaData_atac = as.data.frame(pca$x)
pcaData_atac$sample=rownames(pcaData_atac)
pcaData_atac$HighLow <- stringr::str_extract(pcaData_atac$sample, '(?<=_)[A-Za-z]+(?=_)')
pcaData_atac$Donor <- stringr::str_extract(pcaData_atac$sample, '^[A-Z0-9]+')
percentVar_atac = round(100 * (pca$sdev^2 / sum( pca$sdev^2 ) ))

# == Bulk ==

counts_norm <- read.table('/home/sevastopol/data/gserranos/CART_HL/Data/signature/BatchK_CD8_output_BASAL_BOTH_LowvsHigh.tsv', sep='\t', header=TRUE)
rownames(counts_norm) <- counts_norm$GeneID
counts_norm <- counts_norm[,grepl('^D', colnames(counts_norm))]
group <-factor(stringr::str_extract(colnames(counts_norm), '^D[0-9a]+'))

combat_edata2 = sva::ComBat(dat=as.matrix(counts_norm), batch=group, mod=NULL, par.prior=TRUE, prior.plots=FALSE)
pca <- prcomp(t(combat_edata2))
print(summary(pca))
pcaData_Bulk = as.data.frame(pca$x)
pcaData_Bulk$sample=rownames(pcaData_Bulk)
pcaData_Bulk$HighLow <- stringr::str_extract(pcaData_Bulk$sample, '(?<=_)[A-Za-z]+(?=_)')
pcaData_Bulk$Donor <- stringr::str_extract(pcaData_Bulk$sample, '^[A-Z0-9]+')
percentVar_Bulk = round(100 * (pca$sdev^2 / sum( pca$sdev^2 ) ))

pca_components <- pca$rotation[,1:2]
PC1 <- names(head(pca_components[order(pca_components[,'PC1'], decreasing = TRUE),1], 100)) 
PC2 <- names(head(pca_components[order(pca_components[,'PC2'], decreasing = TRUE),1], 100))
signatures_path <- '/home/sevastopol/data/gserranos/CART_HL/Data/signature/OtherSignatures'

pdf('./Plots/Venn_signaturesVsPCA.pdf')
cowplot::plot_grid(
    get_venn(list(PC1 = unique(PC1), PC2 = unique(PC2), signature = unique(as.character(read.table(paste0(signatures_path, '/', 'Genes_Activation.txt'))$V1)))),
    get_venn(list(PC1 = unique(PC1), PC2 = unique(PC2), signature = unique(as.character(read.table(paste0(signatures_path, '/', 'Genes_Tonic.txt'))$V1)))),
    get_venn(list(PC1 = unique(PC1), PC2 = unique(PC2), signature = unique(as.character(read.table(paste0(signatures_path, '/', 'Genes_Prolif.txt'))$V1))))
, labels= c('Activation', 'Tonic', 'Prolif'), nrow=2)
dev.off()


HLA_DRA <- get_peaks(6, 32432985, 32449203, '/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/BigWig/D14_High_CD8_d0.sort.rmdup.rmblackls.rmchr.norm.bw',  '/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/BigWig/D8a_Low_CD8_S8.sort.rmdup.rmblackls.rmchr.norm.bw', filter_genes = c('HLA-DRA'))
CTLA4   <- get_peaks(2, 203865466, 203876964, '/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/BigWig/D14_High_CD8_d0.sort.rmdup.rmblackls.rmchr.norm.bw',  '/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/BigWig/D9_Low_CD8_S14.sort.rmdup.rmblackls.rmchr.norm.bw', filter_genes = c('CTLA4'), MAX=4)
MCM3   <- get_peaks(6, 52254317, 52319190, '/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/BigWig/D14_High_CD8_d0.sort.rmdup.rmblackls.rmchr.norm.bw',  '/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/BigWig/D9_Low_CD8_S14.sort.rmdup.rmblackls.rmchr.norm.bw', filter_genes = c('MCM3'), MAX=5)

peaks_CD8 <- read.table('/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/results/cd8_Lowd0_vs_Highd0_csaw_denovo_trended-windows_allcounts_Annotated.txt', sep="\t", header=T, fill=TRUE)
peaks_CD8 <- peaks_CD8[grepl('^chr', peaks_CD8$seqnames),]
peaks_CD8 <- makeGRangesFromDataFrame(peaks_CD8)
annotation_Peaks <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
peaks_CD8 <- ChIPseeker::annotatePeak(peaks_CD8, tssRegion=c(-3000, 3000),TxDb=annotation_Peaks, annoDb="org.Hs.eg.db")

DE_Genes_CD8 <- read.table('/home/sevastopol/data/gserranos/CART_HL/Data/signature/BatchK_CD8_output_BASAL_BOTH_LowvsHigh.tsv', sep='\t', header=TRUE, stringsAsFactors=FALSE)[,1:6]
DE_peaks_CD8 <- read.delim("/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/results/cd8_Lowd0_vs_Highd0_csaw_denovo_trended_csaw-windows_significant_Annotated.txt", sep="\t", header=T)
ven_list_CD8 = list('Bulk RNA' = unique(DE_Genes_CD8[DE_Genes_CD8$padj < 0.05 ,'GeneID']), 'ATAC RNA' = unique(as.character(DE_peaks_CD8$SYMBOL)))

get_geom_point <- function(gene, signif=TRUE){
    tmp <- counts_norm[rownames(counts_norm) == gene, ,drop=FALSE]
    tmp <- reshape2::melt(tmp)
    tmp$HighLow <- ifelse(stringr::str_detect( as.character(tmp$variable),'High'), 'High', 'Low' )
    p <- ggplot(tmp, aes(x=HighLow, y=value, fill=HighLow)) + geom_boxplot(alfa=0.8) + ggprism::theme_prism() + scale_fill_manual(values=c('High'='#FCB357', 'Low'='#d3d3d3'))+   ggprism::theme_prism() + labs(y='Expression', title=gene) + theme(legend.position='none', axis.title.x = element_blank(), axis.title.y = element_text(size=10), axis.text.y = element_text(size=6), plot.title = element_text(hjust = 0, size = 10), plot.subtitle=element_text(hjust = 0))
    if (signif){
        p <- p + ggsignif::geom_signif(comparisons = list(c("High", "Low")), map_signif_level = TRUE, vjust=0.5)
    }
    return(p)
}

# old figure
    # pdf('./Plots/Figure2.pdf', width = 15)
    # lgnd <- cowplot::get_legend(get_pca(pcaData_Bulk, percentVar_Bulk, 'HighLow')+ theme(legend.position='bottom'))
    # cowplot::plot_grid(
    #     cowplot::plot_grid(
    #         get_pca(pcaData_Bulk, percentVar_Bulk, 'HighLow') + ggtitle('RNA-seq'), 
    #         lgnd, 
    #         get_pca(pcaData_atac, percentVar_atac, 'HighLow')+ ggtitle('ATAC-seq'), 
    #     nrow=3, rel_heights=c(1,0.2,1)),

    #     cowplot::plot_grid(
    #         cowplot::plot_grid(def_getHeatmap_signature('Genes_Activation', 'Activation')$gtable, 
    #         def_getHeatmap_signature('Genes_Tonic', 'Tonic signal')$gtable, 
    #         def_getHeatmap_signature('Genes_Prolif', 'Proliferation')$gtable, nrow=1), 
    #         cowplot::plot_grid( get_pie_Annot(peaks_CD8 ), get_venn(ven_list_CD8), rel_widths=c(4,2),
    #         nrow=1), 
    #     nrow=2, rel_heights=c(2, 1), align = "h"),

    #     cowplot::plot_grid(
    #         cowplot::plot_grid(
    #             get_geom_point('HLA-DRA'),
    #             get_geom_point('CD74'),
    #             get_geom_point('TNFRSF4'),
    #             get_geom_point('TNFRSF9'),
    #         ncol=2),
    #         cowplot::plot_grid(
    #             cowplot::ggdraw() + cowplot::draw_image(magick::image_read_pdf(MCM3, density = 600)), cowplot::ggdraw() + cowplot::draw_image(magick::image_read_pdf(CTLA4, density = 600)), cowplot::ggdraw() + cowplot::draw_image(magick::image_read_pdf(HLA_DRA, density = 600))
    #         , ncol=2), 
    #     nrow=2, rel_heights=c(3,3)),
    # ncol=3, rel_widths=c(1,4,3))
    # dev.off()
#

pdf('./Plots/Figure2.pdf', width=10, height=15)
lgnd <- cowplot::get_legend(get_pca(pcaData_Bulk, percentVar_Bulk, 'HighLow')+ theme(legend.position='right'))
cowplot::plot_grid(

    cowplot::plot_grid(
        get_pca(pcaData_Bulk, percentVar_Bulk, 'HighLow') + ggtitle('RNA-seq'), 
        get_pca(pcaData_atac, percentVar_atac, 'HighLow')+ ggtitle('ATAC-seq'), 
        lgnd, 
    ncol=3, rel_widths=c(1,1,0.5)),

    cowplot::plot_grid(
        def_getHeatmap_signature('Genes_Activation', 'Activation')$gtable, 
        def_getHeatmap_signature('Genes_Tonic', 'Tonic signal')$gtable, 
        def_getHeatmap_signature('Genes_Prolif', 'Proliferation')$gtable, 
    nrow=1), 
    cowplot::plot_grid(
        cowplot::plot_grid(
                get_geom_point('HLA-DRA'),
                get_geom_point('CD74'),
                get_geom_point('TNFRSF4'),
                get_geom_point('TNFRSF9'),
        ncol=2),
        cowplot::plot_grid( get_pie_Annot(peaks_CD8 ), get_venn(ven_list_CD8), nrow=2, rel_widths=c(1, 1.5), rel_heights=c(1,1.5))
    ,ncol=2),
    cowplot::plot_grid(
            cowplot::ggdraw() + cowplot::draw_image(magick::image_read_pdf(MCM3, density = 600)), 
            cowplot::ggdraw() + cowplot::draw_image(magick::image_read_pdf(CTLA4, density = 600)), 
            cowplot::ggdraw() + cowplot::draw_image(magick::image_read_pdf(HLA_DRA, density = 600))
        , ncol=3),
nrow=4, rel_heights=c(1.5, 3, 2, 2))
dev.off()



 

#### CD4 
# == atac ==
load('/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/RSession/cd4_2105_Csaw_fin.RData')
counts <- edgeR::getCounts(y)
cols_2_keep <- grep('High|Low', grep('d0', colnames(y), value=TRUE), value=TRUE)
counts <- counts[, cols_2_keep]
group <-factor(stringr::str_extract(colnames(counts), '^D[0-9a]+'))
dge <- edgeR::DGEList(counts=counts, group=group)
dge <- edgeR::calcNormFactors(dge, method = "TMM")
logCPM <- edgeR::cpm(dge, log=TRUE, prior.count=2)

combat_edata1 = sva::ComBat(dat=logCPM, batch=group, mod=NULL, par.prior=TRUE, prior.plots=FALSE)
pca <- prcomp(t(combat_edata1))
print(summary(pca))
pcaData_atac = as.data.frame(pca$x)
pcaData_atac$sample=rownames(pcaData_atac)
pcaData_atac$HighLow <- stringr::str_extract(pcaData_atac$sample, '(?<=_)[A-Za-z]+(?=_)')
pcaData_atac$Donor <- stringr::str_extract(pcaData_atac$sample, '^[A-Z0-9]+')
percentVar_atac = round(100 * (pca$sdev^2 / sum( pca$sdev^2 ) ))

# == Bulk ==

counts_norm <- read.table('/home/sevastopol/data/gserranos/CART_HL/Data/signature/BatchK_CD4_output_BASAL_BOTH_LowvsHigh.tsv', sep='\t', header=TRUE)
rownames(counts_norm) <- counts_norm$GeneID
counts_norm <- counts_norm[,grepl('^D', colnames(counts_norm))]
group <-factor(stringr::str_extract(colnames(counts_norm), '^D[0-9a]+'))

combat_edata2 = sva::ComBat(dat=as.matrix(counts_norm), batch=group, mod=NULL, par.prior=TRUE, prior.plots=FALSE)
pca <- prcomp(t(combat_edata2))
print(summary(pca))
pcaData_Bulk = as.data.frame(pca$x)
pcaData_Bulk$sample=rownames(pcaData_Bulk)
pcaData_Bulk$HighLow <- stringr::str_extract(pcaData_Bulk$sample, '(?<=_)[A-Za-z]+')
pcaData_Bulk$Donor <- stringr::str_extract(pcaData_Bulk$sample, '^[A-Z0-9]+')
percentVar_Bulk = round(100 * (pca$sdev^2 / sum( pca$sdev^2 ) ))

pca_components <- pca$rotation[,1:2]
PC1 <- names(head(pca_components[order(pca_components[,'PC1'], decreasing = TRUE),1], 100)) 
PC2 <- names(head(pca_components[order(pca_components[,'PC2'], decreasing = TRUE),1], 100))
signatures_path <- '/home/sevastopol/data/gserranos/CART_HL/Data/signature/OtherSignatures'

pdf('./Plots/Venn_signaturesVsPCA_CD4.pdf')
cowplot::plot_grid(
    get_venn(list(PC1 = unique(PC1), PC2 = unique(PC2), signature = unique(as.character(read.table(paste0(signatures_path, '/', 'Genes_Activation.txt'))$V1)))),
    get_venn(list(PC1 = unique(PC1), PC2 = unique(PC2), signature = unique(as.character(read.table(paste0(signatures_path, '/', 'Genes_Tonic.txt'))$V1)))),
    get_venn(list(PC1 = unique(PC1), PC2 = unique(PC2), signature = unique(as.character(read.table(paste0(signatures_path, '/', 'Genes_Prolif.txt'))$V1))))
, labels= c('Activation', 'Tonic', 'Prolif'), nrow=2)
dev.off()

# Heatmaps

HLA_DRA <- get_peaks(6, 32432985, 32449203, '/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/BigWig/D14_High_CD4_d0.sort.rmdup.rmblackls.rmchr.norm.bw',  '/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/BigWig/D8a_Low_CD4_S7.sort.rmdup.rmblackls.rmchr.norm.bw', filter_genes = c('HLA-DRA'))
CTLA4   <- get_peaks(2, 203865466, 203876964, '/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/BigWig/D14_High_CD4_d0.sort.rmdup.rmblackls.rmchr.norm.bw',  '/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/BigWig/D9_Low_CD4_S13.sort.rmdup.rmblackls.rmchr.norm.bw', filter_genes = c('CTLA4'), MAX=4)
MCM3   <- get_peaks(6, 52254317, 52319190, '/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/BigWig/D14_High_CD4_d0.sort.rmdup.rmblackls.rmchr.norm.bw',  '/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/BigWig/D9_Low_CD4_S13.sort.rmdup.rmblackls.rmchr.norm.bw', filter_genes = c('MCM3'), MAX=5)

peaks_CD4 <- read.table('/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/results/Lowd0_vs_Highd0_csaw_denovo_trended-windows_allcounts_Annotated.txt', sep="\t", header=T, fill=TRUE)
peaks_CD4 <- peaks_CD4[grepl('^chr', peaks_CD4$seqnames),]
peaks_CD4 <- makeGRangesFromDataFrame(peaks_CD4)
annotation_Peaks <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
peaks_CD4 <- ChIPseeker::annotatePeak(peaks_CD4, tssRegion=c(-3000, 3000),TxDb=annotation_Peaks, annoDb="org.Hs.eg.db")

DE_Genes_CD4 <- read.table('/home/sevastopol/data/gserranos/CART_HL/Data/signature/BatchK_CD4_output_BASAL_BOTH_LowvsHigh.tsv', sep='\t', header=TRUE, stringsAsFactors=FALSE)[,1:6]
DE_peaks_CD4 <- read.delim("/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/results/Lowd0_vs_Highd0_csaw_denovo_trended_csaw-windows_significant_Annotated.txt", sep="\t", header=T)
ven_list_CD4 = list('Bulk RNA' = unique(DE_Genes_CD4[DE_Genes_CD4$padj < 0.05 ,'GeneID']), 'ATAC RNA' = unique(as.character(DE_peaks_CD4$SYMBOL)))



pdf('./Plots/Figure2_CD4.pdf', width = 15)
lgnd <- cowplot::get_legend(get_pca(pcaData_Bulk, percentVar_Bulk, 'HighLow')+ theme(legend.position='bottom'))
cowplot::plot_grid(
    cowplot::plot_grid(
        get_pca(pcaData_Bulk, percentVar_Bulk, 'HighLow') + ggtitle('RNA-seq'), lgnd, get_pca(pcaData_atac, percentVar_atac, 'HighLow')+ ggtitle('ATAC-seq'), 
    nrow=3, rel_heights=c(1,0.2,1)),

    cowplot::plot_grid(
        cowplot::plot_grid(cowplot::ggdraw() + cowplot::draw_label("Activation"), cowplot::ggdraw() + cowplot::draw_label("Tonic signal"), cowplot::ggdraw() + cowplot::draw_label("Proliferation"),nrow=1),
        cowplot::plot_grid(def_getHeatmap_signature('Genes_Activation.txt')$gtable, def_getHeatmap_signature('Genes_Tonic.txt')$gtable, def_getHeatmap_signature('Genes_Prolif.txt')$gtable, nrow=1), 
        cowplot::plot_grid( get_pie_Annot(peaks_CD4 ), get_venn(ven_list_CD4), rel_widths=c(4,2),
        nrow=1), 
    nrow=3, rel_heights=c(0.1, 1, 1), align = "h"),

    cowplot::plot_grid(
        NULL, 
        cowplot::plot_grid(cowplot::ggdraw() + cowplot::draw_image(magick::image_read_pdf(MCM3, density = 600)), cowplot::ggdraw() + cowplot::draw_image(magick::image_read_pdf(CTLA4, density = 600)), cowplot::ggdraw() + cowplot::draw_image(magick::image_read_pdf(HLA_DRA, density = 600)), ncol=1), 
    nrow=2, rel_heights=c(2,3)),
ncol=3, rel_widths=c(1,1.5,1))
dev.off()








htmap <- DoHeatmap(rds_test, group.by='seurat_clusters', features=all_markers_top20$gene, angle = 0, group.colors = hues::iwanthue(length(levels(all_markers_top20$cluster))), size = 3, combine = TRUE) + NoLegend() + theme(axis.text.y = element_text(size = 3))
pdf('./Plots/Figure3.pdf')
    # legend <- cowplot::get_legend(get_umap(coords_markers, 'CD4')+ theme(legend.position='right'))
    cowplot::plot_grid(
        cowplot::plot_grid(
        ggplot(coords, aes(x=UMAP_1, y=UMAP_2, color= Cluster)) + geom_point(alpha=0.6, size = 0.5) + scale_color_manual(values=color_list_populations) + theme_void() + 
        theme(legend.position='right', plot.title = element_text(hjust = 0.5)) + ggtitle('Populations') + 
        guides(color = guide_legend(override.aes = list(size=3, alpha = 1, font=4))),
        cowplot::plot_grid(
            get_umap(normData, 'CD4'),
            get_umap(normData, 'TCF7'),
            get_umap(normData, 'HLA-DRA'),
            get_umap(normData, 'GZMA'),
            get_umap(normData, 'LAG3', 5),
            get_umap(normData, 'CD8A'),
            get_umap(normData, 'CCR7'),
            get_umap(normData, 'GATA3'),
            get_umap(normData, 'PRF1'),
            get_umap(normData, 'TIGIT', 5),
            # legend,
            nrow=2),
        nrow=2),

        cowplot::plot_grid(htmap, ncol=1),
    ncol=2, rel_widths = c(2, 1)
    )
dev.off()

""

cells_clust20 <- sub('\\.', '-',coords_unfiltered[coords_unfiltered$Cluster == '20','cell_id'])
test<-as.data.frame(rds_unfiltered@assays$SCT@counts)
test <- t(test[c('CD8A', 'CD4', 'CD3D'), cells_clust20])
test <- reshape2::melt(test)
pdf('./Plots/Cluster_20.pdf', width=15)
ggplot(test, aes(x=Var2, y=value)) + geom_violin(aes(fill=Var2),adjust = 1.8) + 
    theme_classic() + ggtitle('Cluster 20') + theme(legend.position='none', axis.text=element_text(size=6)) + labs(y = 'Value') + hues::scale_fill_iwanthue()
Idents(rds_unfiltered) <- rds_unfiltered$integrated_snn_res.0.8
VlnPlot(rds_unfiltered, features = c("CD8A", 'CD4', 'CD3D'),
    pt.size = 0, ncol = 3, cols = color_list_all_stripped[0:length(unique(coords_unfiltered$Cluster))])
dev.off()

pdf('./Plots/FigureS5.pdf',height=12, width=10)
    coords_unfiltered <- as.data.frame(rds_unfiltered@reductions$umap@cell.embeddings)
    coords_unfiltered$cell_id <- sub('-', '.',rownames(coords_unfiltered))
    coords_unfiltered$Cluster <- rds_unfiltered$integrated_snn_res.0.8
    clusters_tmp  <- setNames(as.data.frame(rds_unfiltered$seurat_clusters) , 'Cluster')
    clusters_tmp$cell_id <- rownames(clusters_tmp)
    # Set the colors, maintining the colors for the conserved clusters
    color_list_all <- c('#A55E34', '#C6B2D3', '#D0342B', '#8E221A', '#2E6B34', '#BBDE93', '#AECDE1', '#3C77AF', '#ED9E9B', '#DA913D', '#821851', '#643F95', '#DBBE78', '#7F7F7F', '#000000')
    names(color_list_all) <- stringr::str_extract(levels(rds_test$ClusterNames_0.8_by_JR), 'C?[0-9]{1,2}(?=\\.)')
    names(color_list_all)[15] <- 'C21'
    new_colors_names <- setdiff(paste0('C', seq(0,length(unique(coords_unfiltered$Cluster)))),names(color_list_all))
    new_colors <- hues::iwanthue(length(new_colors_names))
    names(new_colors) <- new_colors_names
    color_list_all_stripped <- c(color_list_all, new_colors)
    color_list_all_stripped <- color_list_all_stripped[order(stringr::str_extract(names(color_list_all_stripped),'C?[0-9]{1,2}(?=\\.)'))]
    names(color_list_all_stripped) <- stringr::str_extract(names(color_list_all_stripped), '[0-9]{1,2}')
    color_list_all_stripped <-  color_list_all_stripped[names(color_list_all_stripped)[order(as.numeric(names(color_list_all_stripped)))]]
    # end of colors
    g2mScore  <- setNames(as.data.frame(rds_unfiltered$G2M.Score) , 'Value')
    g2mScore <- merge(g2mScore, clusters_tmp, by.x=0, by.y='cell_id')
    sScore    <- setNames(as.data.frame(rds_unfiltered$S.Score) , 'Value')
    sScore <- merge(sScore, clusters_tmp, by.x=0, by.y='cell_id')
    mitoRatio <- setNames(as.data.frame(rds_unfiltered$mitoRatio) , 'Value')
    mitoRatio <- merge(mitoRatio, clusters_tmp, by.x=0, by.y='cell_id')
    g2mScore  <- get_violin(g2mScore, 'G2M phase', cluster_color = color_list_all_stripped)
    sScore    <- get_violin(sScore, 'S phase', cluster_color = color_list_all_stripped)
    mitoRatio <- get_violin(mitoRatio, 'Ratio of mitochondrial genes', cluster_color = color_list_all_stripped)
    contribution <- merge(clusters_tmp, setNames(as.data.frame(rds_unfiltered$donor) , 'Donor'), by.x='cell_id', by.y=0)
    contribution2 <- setNames(as.data.frame(prop.table(table(contribution[contribution$Donor == 'd10', 'Cluster']))*100), c('Cluster', 'd10'))
    contribution2 <- merge(contribution2, setNames(as.data.frame(prop.table(table(contribution[contribution$Donor == 'd14', 'Cluster']))*100), c('Cluster', 'd14')), by='Cluster')
    contribution2 <-merge(contribution2,  setNames(as.data.frame(prop.table(table(contribution[contribution$Donor == 'd18', 'Cluster']))*100), c('Cluster', 'd18')), by='Cluster')
    contribution2 <- setNames(reshape2::melt(contribution2), c('Cluster', 'Donor', 'value'))
    contribution <- ggplot(contribution, aes(x= Cluster, fill=Donor))+  geom_bar(position="fill", colour="black") + theme_classic() + labs(title = "Contribution by donor", x= 'Cluster', y = '% of cells') + scale_fill_brewer(palette="RdBu")  #scale_fill_manual(values=c('#585123', '#f2a65a', '#772f1a'))
    contribution2 <- ggplot(contribution2, aes(x= Cluster, y=value, fill=Donor))+  geom_bar(position="fill", stat='identity', colour="black") + theme_classic() + labs(title = "Contribution by donor", x= 'Cluster', y = '% of cells') + scale_fill_brewer(palette="RdBu")+ theme(axis.text=element_text(size=6))  #scale_fill_manual(values=c('#585123', '#f2a65a', '#772f1a'))
    cluster_20_CD8A <- setNames(as.data.frame(rds_unfiltered@assays$RNA@data[c('CD8A'),]), 'Value')
    cluster_20_CD8A <- merge(cluster_20_CD8A, clusters_tmp, by.x=0, by.y='cell_id')
    cluster_20_CD8A <- get_violin(cluster_20_CD8A, 'CD8A', cluster_color = color_list_all_stripped) + ylab('Expression Level')
    cluster_20_CD4 <- setNames(as.data.frame(rds_unfiltered@assays$RNA@data[c('CD4'),]), 'Value') 
    cluster_20_CD4 <- merge(cluster_20_CD4, clusters_tmp, by.x=0, by.y='cell_id')
    cluster_20_CD4 <- get_violin(cluster_20_CD4, 'CD4', cluster_color = color_list_all_stripped) + theme(axis.title.y=element_blank())
    cluster_20_CD3D <- setNames(as.data.frame(rds_unfiltered@assays$RNA@data[c('CD3D'),]), 'Value')
    cluster_20_CD3D <- merge(cluster_20_CD3D, clusters_tmp, by.x=0, by.y='cell_id')
    cluster_20_CD3D <- get_violin(cluster_20_CD3D, 'CD3D', cluster_color = color_list_all_stripped)+ theme(axis.title.y=element_blank())
    tmp <- setNames(as.data.frame(rds_unfiltered$donor), 'Donor')
    tmp$cell_id <- sub('-', '.', rownames(tmp))
    tmp <- merge(coords, tmp, by='cell_id')
    cowplot::plot_grid(  
        cowplot::plot_grid(   
            cowplot::plot_grid(
                ggplot(coords_unfiltered, aes(x=UMAP_1, y=UMAP_2, color= Cluster)) + geom_point(alpha=0.9, size = 0.5) + scale_color_manual(values=color_list_all_stripped[0:length(unique(coords_unfiltered$Cluster))]) + theme_void() + theme(legend.position='right', plot.title = element_text(hjust = 0.5)) + ggtitle('Populations') + guides(color = guide_legend(override.aes = list(size=3, alpha =1))), 
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
        , ncol=2, rel_widths=c(4,3)),
        cowplot::plot_grid( cluster_20_CD8A,
                            cluster_20_CD4,
                            cluster_20_CD3D,        
                        ncol=3),
    nrow=2, rel_heights=c(3,1))
        
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


rds_test_5Clusters <- rds_test
Idents(rds_test_5Clusters) <- rds_test_5Clusters$ClusterNames_0.8_by_JR
cluster_names <- c('C3.CD8 Memory', 'C8.CD8 Cytotoxic', 'C9.CD8 Cytotoxic (late)', 'C6.CD4 Activated', 'C17.CD4 Activated')
rds_test_5Clusters <- subset(rds_test_5Clusters, idents=cluster_names)
rds_test_5Clusters$ClusterNames_0.8_by_JR <- factor(rds_test_5Clusters$ClusterNames_0.8_by_JR, level=c(as.character(unique(rds_test_5Clusters$ClusterNames_0.8_by_JR))))
Idents(rds_test_5Clusters) <- rds_test_5Clusters$BinScore
High_Low_DE_markers <- FindMarkers(rds_test_5Clusters, ident.1 = "High", ident.2 = "Low")
High_Low_DE_markers$gene <- rownames(High_Low_DE_markers)
htmap <- DoHeatmap(rds_test_5Clusters, group.by='ClusterNames_0.8_by_JR', features=High_Low_DE_markers$gene, angle = 0, group.colors = c('#30A3CC', '#FCB357'), size = 3, combine = TRUE) + theme(axis.text.y = element_text(size = 3))


test <- as.data.frame(rds_test_5Clusters@assays$RNA@scale.data)
test <- test[rownames(test) %in% High_Low_DE_markers$gene, ]
annotation_HL <- as.data.frame(rds_test_5Clusters$BinScore)
annotation_Cl <- as.data.frame(rds_test_5Clusters$ClusterNames_0.8_by_JR)
annotation_col <- setNames(merge(annotation_HL, annotation_Cl, by=0), c('cell_id', 'HighLow', 'Cluster'))
col_order <- annotation_col[with(annotation_col, order(Cluster, HighLow)),]
rownames(annotation_col) <- annotation_col$cell_id
annotation_col$cell_id <- NULL
colors_ann <- list(
    'HighLow' = c(High = '#30A3CC', Low = '#FCB357'),
    'Cluster' = c(color_list_populations)
)

col_gaps <- c()
for (clust in levels(col_order$Cluster)){
    col_gaps <- c(col_gaps, which(col_order$Cluster == clust)[1]-1)
}



test[test >2.5] <- 2.5
test[test < -2.5] <- -2.5

colors_ann[['HighLow']]['Low'] <- '#d3d3d3'
htmap <- pheatmap::pheatmap(test[,c(col_order$cell_id) ],scale='none', fontsize=5, angle_col =45, show_rownames=TRUE, show_colnames = FALSE, treeheight_row = 0, treeheight_col = 5, cluster_row =TRUE, cluster_col = FALSE, gaps_col = col_gaps, legend=TRUE, annotation_legend=FALSE, silent = TRUE, color = PurpleAndYellow()[5:45]
, fontsize_row=2, annotation_col = annotation_col, annotation_colors=colors_ann)
# , cellwidth=0.03)

signatures_path <- '/home/sevastopol/data/gserranos/CART_HL/Data/signature/OtherSignatures'
signatures <- list.files(signatures_path)

pdf('./Plots/Figure4.pdf', width=7.5, height=10)
    HighAndLow = FALSE
    if(HighAndLow){
        clusters_tmp  <- setNames(as.data.frame(rds_test$seurat_clusters) , 'Cluster')
        clusters_tmp$cell_id <- sub('-','.',rownames(clusters_tmp))
        composition <- merge(coords[, c('cell_id', 'BinScore')], clusters_tmp, by='cell_id')
        composition <- ggplot(composition, aes(x= Cluster, fill=BinScore))+  geom_bar(position="fill") + ggprism::theme_prism() + labs(title = "High CART distribution", x= 'Cluster', y = '% of cells') + scale_fill_manual(values=c('#30A3CC', '#FCB357')) 
    }else{
        clusters_tmp  <- setNames(as.data.frame(rds_test$seurat_clusters) , 'Cluster')
        clusters_tmp$cell_id <- sub('-','.',rownames(clusters_tmp))
        composition <- merge(coords[, c('cell_id', 'BinScore')], clusters_tmp, by='cell_id')
        composition$Cluster <- as.character(composition$Cluster)
        composition <- table(composition$Cluster, composition$BinScore)
        composition <- as.data.frame.matrix(composition)
        composition$High_prop <- apply(composition, 1, FUN=function(x) (x[1]/sum(x))*100 )
        composition$Cluster <- factor(rownames(composition) , levels=c(sort(as.numeric(rownames(composition)))))
        composition <- ggplot(composition, aes(x=Cluster, y=High_prop))+  geom_bar(stat='identity', fill = "#30A3CC") + ggprism::theme_prism() + labs(title = "High CART distribution", x= 'Cluster', y = '% of cells') + theme(panel.grid.major = element_line(colour="#f0f0f0"))
    }
    cowplot::plot_grid(
        cowplot::plot_grid(
            ggplot() + 
            geom_point(coords[coords$BinScore == 'Low',], mapping = aes(x=UMAP_1, y=UMAP_2, color= BinScore, shape = CD),alpha=0.6) +
            geom_point(coords[coords$BinScore == 'High',], mapping = aes(x=UMAP_1, y=UMAP_2, color= BinScore, shape = CD),alpha=0.8) + 
            scale_color_manual(values=c('High'='#30A3CC', 'Low'='#d3d3d3'))  + ggprism::theme_prism()+ labs(subtitle = 'High-low distribution') + 
            guides(fill = guide_legend(override.aes = list(size=3, alpha = 1))) + labs(x = "UMAP 1", y = 'UMAP 2') + 
            theme(legend.position='bottom', plot.title = element_text(hjust = 0.5, size = 10), axis.text=element_blank(), axis.ticks=element_blank(), legend.spacing.x = unit(0.1, 'cm')) ,
            composition, 
        ncol=2),
        cowplot::plot_grid(
            cowplot::plot_grid(
                cowplot::plot_grid(
                    get_expression_signature("Genes_Activation.txt", normData, coords , 1.2, split_HL=TRUE) + theme(strip.text = element_text(size=8)),
                    get_expression_signature("Genes_Tonic.txt", normData, coords , 1.2, split_HL=TRUE)+ theme(strip.text = element_text(size=8)),
                ncol=1),
                cowplot::get_legend(get_expression_signature("Genes_Activation.txt", normData, coords , 1.2)+theme(legend.position='bottom')+ guides(color = guide_colourbar(barheight = 0.5))),
            ncol=1, rel_heights=c(2, 0.1)),
            cowplot::plot_grid(
                get_plot_signature("Genes_Activation.txt", normData, "C6.CD4 Activated") + labs(title = "C6.CD4 Activated", subtitle = "Activation genes") + theme(plot.title = element_text(size = 8),plot.subtitle = element_text(size = 6)),
                get_plot_signature("Genes_Activation.txt", normData, "C17.CD4 Activated", FALSE)+ labs(title = "C17.CD4 Activated", subtitle = "Activation genes")+ theme(plot.title = element_text(size = 8),plot.subtitle = element_text(size = 6)),
                get_plot_signature("Genes_Activation.txt", normData, "C9.CD8 Cytotoxic (late)")+ labs(title = "C9.CD8 Cytotoxic (late)", subtitle = "Activation genes")+ theme(plot.title = element_text(size = 8),plot.subtitle = element_text(size = 6)),
                get_plot_signature("Genes_Tonic.txt", normData, "C9.CD8 Cytotoxic (late)")+ labs(title = "C9.CD8 Cytotoxic (late)", subtitle = "Tonic signal")+ theme(plot.title = element_text(size = 8),plot.subtitle = element_text(size = 6)),
            ncol=2),
        ncol=2)    
    , nrow=2, rel_heights=c(2,3))
dev.off()


pdf('./Plots/Figure4_old.pdf', width=9, height=12)
HighAndLow = FALSE
if(HighAndLow){
    clusters_tmp  <- setNames(as.data.frame(rds_test$seurat_clusters) , 'Cluster')
    clusters_tmp$cell_id <- sub('-','.',rownames(clusters_tmp))
    composition <- merge(coords[, c('cell_id', 'BinScore')], clusters_tmp, by='cell_id')
    composition <- ggplot(composition, aes(x= Cluster, fill=BinScore))+  geom_bar(position="fill") + ggprism::theme_prism() + labs(title = "High CART distribution", x= 'Cluster', y = '% of cells') + scale_fill_manual(values=c('#30A3CC', '#FCB357')) 
}else{
    clusters_tmp  <- setNames(as.data.frame(rds_test$seurat_clusters) , 'Cluster')
    clusters_tmp$cell_id <- sub('-','.',rownames(clusters_tmp))
    composition <- merge(coords[, c('cell_id', 'BinScore')], clusters_tmp, by='cell_id')
    composition$Cluster <- as.character(composition$Cluster)
    composition <- table(composition$Cluster, composition$BinScore)
    composition <- as.data.frame.matrix(composition)
    composition$High_prop <- apply(composition, 1, FUN=function(x) (x[1]/sum(x))*100 )
    composition$Cluster <- factor(rownames(composition) , levels=c(sort(as.numeric(rownames(composition)))))
    composition <- ggplot(composition, aes(x=Cluster, y=High_prop))+  geom_bar(stat='identity', fill = "#30A3CC") + ggprism::theme_prism() + labs(title = "High CART distribution", x= 'Cluster', y = '% of cells') + theme(panel.grid.major = element_line(colour="#f0f0f0"))
}
cowplot::plot_grid(
    cowplot::plot_grid(
        cowplot::plot_grid(
            ggplot() + 
            geom_point(coords[coords$BinScore == 'Low',], mapping = aes(x=UMAP_1, y=UMAP_2, color= BinScore, shape = CD),alpha=0.6) +
            geom_point(coords[coords$BinScore == 'High',], mapping = aes(x=UMAP_1, y=UMAP_2, color= BinScore, shape = CD),alpha=0.8) + 
            scale_color_manual(values=c('High'='#30A3CC', 'Low'='#d3d3d3'))  + ggprism::theme_prism()+ labs(subtitle = 'High-low distribution') + 
            guides(fill = guide_legend(override.aes = list(size=3, alpha = 1))) + labs(x = "UMAP 1", y = 'UMAP 2') + 
            theme(legend.position='bottom', plot.title = element_text(hjust = 0.5, size = 10), axis.text=element_blank(), axis.ticks=element_blank()),
            composition, 
        nrow=2, rel_heights=c(1,0.6)),
        htmap$gtable,
        ncol=2, rel_widths=c(3,2)),
    cowplot::plot_grid(
        cowplot::plot_grid(
            cowplot::plot_grid(
                get_expression_signature("Genes_Activation.txt", normData, coords , 1.2, split_HL=TRUE) + theme(strip.text = element_text(size=8)),
                get_expression_signature("Genes_Tonic.txt", normData, coords , 1.2, split_HL=TRUE)+ theme(strip.text = element_text(size=8)),
            ncol=1),
            cowplot::get_legend(get_expression_signature("Genes_Activation.txt", normData, coords , 1.2)+theme(legend.position='bottom')+ guides(color = guide_colourbar(barheight = 0.5))),
        ncol=1, rel_heights=c(2, 0.2)),
        cowplot::plot_grid(
            get_plot_signature("Genes_Activation.txt", normData, "C6.CD4 Activated") + labs(title = "C6.CD4 Activated", subtitle = "Activation genes") + theme(plot.title = element_text(size = 8),plot.subtitle = element_text(size = 6)),
            get_plot_signature("Genes_Activation.txt", normData, "C17.CD4 Activated")+ labs(title = "C17.CD4 Activated", subtitle = "Activation genes")+ theme(plot.title = element_text(size = 8),plot.subtitle = element_text(size = 6)),
            get_plot_signature("Genes_Activation.txt", normData, "C9.CD8 Cytotoxic (late)")+ labs(title = "C9.CD8 Cytotoxic (late)", subtitle = "Activation genes")+ theme(plot.title = element_text(size = 8),plot.subtitle = element_text(size = 6)),
            get_plot_signature("Genes_Tonic.txt", normData, "C9.CD8 Cytotoxic (late)")+ labs(title = "C9.CD8 Cytotoxic (late)", subtitle = "Tonic signal")+ theme(plot.title = element_text(size = 8),plot.subtitle = element_text(size = 6)),
        ncol=2),
    ncol=2),
nrow=2)
dev.off()

# pdf('./Plots/Test.pdf')
# ggplot() + 
# geom_point(coords[coords$BinScore == 'Low',], mapping = aes(x=UMAP_1, y=UMAP_2, color= BinScore, shape = CD),alpha=0.6) +
# geom_point(coords[coords$BinScore == 'High',], mapping = aes(x=UMAP_1, y=UMAP_2, color= BinScore, shape = CD),alpha=0.8) + 
# scale_color_manual(values=c('High'='#30A3CC', 'Low'='#d3d3d3'))  + theme_void()+ labs(subtitle = 'High-low distribution') + 
# guides(fill = guide_legend(override.aes = list(size=3, alpha = 1))) + 
# theme(legend.position='bottom', plot.title = element_text(hjust = 0.5, size = 10))
# dev.off()

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


# Figure SIMIC
df_auc <- readRDS('/home/sevastopol/data/gserranos/CART_HL/SimiC/Data/SimiC_aucs.rds')
MinMax_clust <- readRDS('/home/sevastopol/data/gserranos/CART_HL/SimiC/Data/MinMax_clust.rds')
clusters_2_keep <- c( "C3.CD8 Memory","C4.CD4 Memory", "C6.CD4 Activated","C8.CD8 Cytotoxic",'C9.CD8 Cytotoxic (late)', "C12.CD4 Th2 helper", "C17.CD4 Activated")
Paula_list <- c('NR4A1', 'BATF', 'ARID5A', 'RFX5', 'SATB1', 'ATF5', 'ZBTB7B', 'EED' )

get_densities <- function(driver, cluster, title=NA){
    plotter2 <- df_auc[df_auc$cluster_id == cluster,]
    p <- ggplot(plotter2[plotter2$driver ==driver,], aes(x=value, fill=.id)) + 
    geom_density(alpha = 0.6, adjust = 1/8) + ggprism::theme_prism() + xlim(0,1)+
    scale_fill_manual( values=c('#30A3CC', '#bfbfbf')) +
    theme(legend.position = 'none',plot.title = element_text(hjust = 0, size = 12), axis.title.x = element_blank(), axis.title.y = element_blank())
    if(!is.na(title)){
        p <- p + ggtitle(title)
    }
    return(p)
}

get_curves <- function(driver_list){
    tmp <- MinMax_clust[rownames(MinMax_clust) %in% driver_list, c('C3.CD8 Memory', 'C8.CD8 Cytotoxic', 'C9.CD8 Cytotoxic (late)')]
    tmp$driver <- rownames(tmp)
    tmp <- reshape2::melt(tmp)
    p <- ggplot(tmp, aes(x=variable, y=value, group=driver, color=driver)) + geom_line() + geom_point() + ggprism::theme_prism() + ggthemes::scale_colour_economist() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    return(p)
}

get_ribbons <- function(driver_list){
    tmp <- df_auc[df_auc$cluster_id %in%  c('C3.CD8 Memory', 'C8.CD8 Cytotoxic', 'C9.CD8 Cytotoxic (late)') & df_auc$driver %in% driver_list,]
    p <- ggplot(tmp, aes(x=value, y=cluster_id, group=driver, color=driver)) + geom_line() + geom_ribbon() + ggprism::theme_prism() + ggthemes::scale_colour_economist() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    return(p)
}

plot_name <- function(name, angle=90){
    p <- cowplot::ggdraw() + cowplot::draw_label(name, angle=angle)
    return(p)
}

pdf('./Plots/SimiC_image.pdf', width=12, height=15)
cowplot::plot_grid(
    cowplot::plot_grid(plot_name('Transcription\nFactor',0), plot_name('C3.CD8 Memory', 0),plot_name('C8.CD8 Cytotoxic',0), plot_name('C9.CD8 Cytotoxic (late)',0),nrow=1,rel_widths=c(1,5,5,5)),
    cowplot::plot_grid(plot_name('NR4A1'), get_densities('NR4A1' , 'C3.CD8 Memory'), get_densities('NR4A1' , 'C8.CD8 Cytotoxic'), get_densities('NR4A1' , 'C9.CD8 Cytotoxic (late)'), nrow=1,rel_widths=c(1,5,5,5)),
    cowplot::plot_grid(plot_name('BATF'), get_densities('BATF'  , 'C3.CD8 Memory'), get_densities('BATF'  , 'C8.CD8 Cytotoxic'), get_densities('BATF'  , 'C9.CD8 Cytotoxic (late)'), nrow=1,rel_widths=c(1,5,5,5)),
    cowplot::plot_grid(plot_name('ARID5A'), get_densities('ARID5A', 'C3.CD8 Memory'), get_densities('ARID5A', 'C8.CD8 Cytotoxic'), get_densities('ARID5A', 'C9.CD8 Cytotoxic (late)'), nrow=1,rel_widths=c(1,5,5,5)),
    cowplot::plot_grid(plot_name('RFX5'), get_densities('RFX5'  , 'C3.CD8 Memory'), get_densities('RFX5'  , 'C8.CD8 Cytotoxic'), get_densities('RFX5'  , 'C9.CD8 Cytotoxic (late)'), nrow=1,rel_widths=c(1,5,5,5)),
    cowplot::plot_grid(plot_name('SATB1'), get_densities('SATB1' , 'C3.CD8 Memory'), get_densities('SATB1' , 'C8.CD8 Cytotoxic'), get_densities('SATB1' , 'C9.CD8 Cytotoxic (late)'), nrow=1,rel_widths=c(1,5,5,5)),
    cowplot::plot_grid(plot_name('ATF5'), get_densities('ATF5'  , 'C3.CD8 Memory'), get_densities('ATF5'  , 'C8.CD8 Cytotoxic'), get_densities('ATF5'  , 'C9.CD8 Cytotoxic (late)'), nrow=1,rel_widths=c(1,5,5,5)),
    cowplot::plot_grid(plot_name('ZBTB7B'), get_densities('ZBTB7B', 'C3.CD8 Memory'), get_densities('ZBTB7B', 'C8.CD8 Cytotoxic'), get_densities('ZBTB7B', 'C9.CD8 Cytotoxic (late)'), nrow=1,rel_widths=c(1,5,5,5)),
    cowplot::plot_grid(plot_name('EED'), get_densities('EED'   , 'C3.CD8 Memory'), get_densities('EED'   , 'C8.CD8 Cytotoxic'), get_densities('EED'   , 'C9.CD8 Cytotoxic (late)'), nrow=1,rel_widths=c(1,5,5,5)),
nrow=9, rel_heights=c(1,5,5,5,5,5,5,5,5))
dev.off()

pdf('./Plots/SimiC_image_2.pdf')
get_curves(Paula_list)
# get_ribbons(Paula_list)
dev.off()

plot_dims <- get_plot_dims(p)

pdf('./Plots/SimiC_image_HM.pdf' )
pheatmap::pheatmap(MinMax_clust,color=viridis::viridis(50, direction = 1, option = "C"), fontsize=5, angle_col =45, cellwidth=40, treeheight_col=20, treeheight_row=10)
pheatmap::pheatmap(MinMax_clust[, c('C3.CD8 Memory', 'C8.CD8 Cytotoxic', 'C9.CD8 Cytotoxic (late)')],color=viridis::viridis(50, direction = 1, option = "C"), fontsize=5, angle_col =45, cellwidth=40, treeheight_col=20, treeheight_row=10)
dev.off()


pm <- pheatmap::pheatmap(MinMax_clust[, c('C3.CD8 Memory', 'C8.CD8 Cytotoxic', 'C9.CD8 Cytotoxic (late)')],color=viridis::viridis(50, direction = 1, option = "C"), fontsize=5, angle_col =45, cellwidth=40, treeheight_col=20, treeheight_row=10, silent=TRUE, legend=TRUE)$gtable

pdf('./Plots/Figure5.pdf', 12, 15)
cowplot::plot_grid(
    cowplot::plot_grid(
        pm, get_curves(c('RFX5', 'SATB1','NR4A1')),
    ncol=2)
,
    cowplot::plot_grid(
        cowplot::plot_grid(
            plot_name('NR4A1'), get_densities('NR4A1' , 'C3.CD8 Memory', 'C3.CD8 Memory'), get_densities('NR4A1' , 'C8.CD8 Cytotoxic', 'C8.CD8 Cytotoxic'), get_densities('NR4A1' , 'C9.CD8 Cytotoxic (late)', 'C9.CD8 Cytotoxic (late)') ,NULL, 
        nrow=1,rel_widths=c(1,5,5,5,5)),
        cowplot::plot_grid(
            plot_name('RFX5'), get_densities('RFX5'  , 'C3.CD8 Memory'), get_densities('RFX5'  , 'C8.CD8 Cytotoxic'), get_densities('RFX5'  , 'C9.CD8 Cytotoxic (late)'),NULL, 
        nrow=1,rel_widths=c(1,5,5,5,5)),
        cowplot::plot_grid(
            plot_name('SATB1'), get_densities('SATB1' , 'C3.CD8 Memory'), get_densities('SATB1' , 'C8.CD8 Cytotoxic'), get_densities('SATB1' , 'C9.CD8 Cytotoxic (late)'),NULL,  
        nrow=1,rel_widths=c(1,5,5,5,5))
    ,nrow=3)
,ncol=1,nrow=2, rel_heights = 2,3)
dev.off()





pdf('./Plots/FigureS6_log.pdf')
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



# ATAC-seq
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(EnsDb.Hsapiens.v86)
edb <- EnsDb.Hsapiens.v86
annotation_Peaks <- TxDb.Hsapiens.UCSC.hg38.knownGene

# CD4
final.merged.peaks <- read.table("/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/results/Lowd0_vs_Highd0_csaw_denovo_trended_csaw-windows_all.txt", sep="\t", header=T)
FDR.thresh <- 0.05
final.merged.peaks.sig <- final.merged.peaks[final.merged.peaks$FDR < FDR.thresh, ]

final.merged.peaks$sig <- "n.s."
final.merged.peaks$sig[final.merged.peaks$FDR < FDR.thresh] <- "significant"

pdf('./Plots/CD4_DE_peaks_Low_Vs_High.pdf')
ggplot(final.merged.peaks, aes(x = logCPM, y = logFC, col = factor(sig, levels=c("n.s.", "significant")))) + 
  geom_point() + scale_color_manual(values = c("black", "red")) + 
#   geom_smooth() + # smoothed loess fit; can add span=0.5 to reduce computation load/time
  geom_hline(yintercept = 0) + labs(col = NULL)
ggplot(final.merged.peaks, aes(x = logCPM, y = logFC, col = factor(sig, levels=c("n.s.", "significant")))) + 
  geom_point() + scale_color_manual(values = c("black", "red")) + 
  geom_smooth() + # smoothed loess fit; can add span=0.5 to reduce computation load/time
  geom_hline(yintercept = 0) + labs(col = NULL)
dev.off()


# we need to use the gapped peaks
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(EnsDb.Hsapiens.v86)
edb <- EnsDb.Hsapiens.v86
annotation_Peaks <- TxDb.Hsapiens.UCSC.hg38.knownGene



path_High_CD8 <- '/home/sevastopol/data/gserranos/CART_HL/Data/High_low/ATACseq/HIGH/CD8/High_CD8.narrowPeak'
path_High_CD4 <- '/home/sevastopol/data/gserranos/CART_HL/Data/High_low/ATACseq/HIGH/CD4/High_CD4.narrowPeak'
path_Low_CD8  <- '/home/sevastopol/data/gserranos/CART_HL/Data/High_low/ATACseq/LOW/Low_CD8.narrowPeak'
path_Low_CD4  <- '/home/sevastopol/data/gserranos/CART_HL/Data/High_low/ATACseq/LOW/Low_CD4.narrowPeak'

peak_High_CD8 <- readPeakFile(path_High_CD8, header=FALSE)
peak_High_CD4 <- readPeakFile(path_High_CD4, header=FALSE)
peak_Low_CD8  <- readPeakFile(path_Low_CD8, header=FALSE)
peak_Low_CD4  <- readPeakFile(path_Low_CD4, header=FALSE)

peak_High_CD8 <- keepStandardChromosomes(peak_High_CD8, pruning.mode='coarse')
peak_High_CD4 <- keepStandardChromosomes(peak_High_CD4, pruning.mode='coarse')
peak_Low_CD8  <- keepStandardChromosomes(peak_Low_CD8, pruning.mode='coarse')
peak_Low_CD4  <- keepStandardChromosomes(peak_Low_CD4, pruning.mode='coarse')


pdf('./Plots/peakLocationGenome.pdf')
covplot(peak_High_CD8, weightCol='V10', title = 'High_CD8')
covplot(peak_High_CD4, weightCol='V10', title = 'High_CD4')
covplot(peak_Low_CD8,  weightCol='V10', title = 'Low_CD8')
covplot(peak_Low_CD4,  weightCol='V10', title = 'Low_CD4')
dev.off()

promoter <- getPromoters(TxDb=annotation_Peaks, upstream=3000, downstream=3000)


tagMatrix_High_CD8 <- getTagMatrix(peak_High_CD8, windows=promoter)
tagMatrix_High_CD4 <- getTagMatrix(peak_High_CD4, windows=promoter)
tagMatrix_Low_CD8  <- getTagMatrix(peak_Low_CD8, windows=promoter)
tagMatrix_Low_CD4  <- getTagMatrix(peak_Low_CD4, windows=promoter)


pdf('./Plots/TSSbindings.pdf', width= 4)
tagHeatmap(tagMatrix_High_CD8, xlim=c(-3000, 3000), color="red", title='High_CD8') 
tagHeatmap(tagMatrix_High_CD4, xlim=c(-3000, 3000), color="red", title='High_CD4') 
tagHeatmap(tagMatrix_Low_CD8 , xlim=c(-3000, 3000), color="red", title='Low_CD8') 
tagHeatmap(tagMatrix_Low_CD4 , xlim=c(-3000, 3000), color="red", title='Low_CD4') 
dev.off()


pdf('./Plots/AverageProfile.pdf')
High_CD8 <- plotAvgProf(tagMatrix_High_CD8, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
High_CD4 <- plotAvgProf(tagMatrix_High_CD4, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
Low_CD8 <- plotAvgProf(tagMatrix_Low_CD8, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
Low_CD4 <- plotAvgProf(tagMatrix_Low_CD4, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")

High_CD8 <- High_CD8$data
High_CD8$HighLow <- 'High'
Low_CD8 <- Low_CD8$data
Low_CD8$HighLow <- 'Low'
plotter_CD8 <- rbind(High_CD8, Low_CD8)
ggplot(plotter_CD8, aes(x=pos, y=value, color=HighLow)) + geom_line() + 
xlim(-3000, 3000) + ggtitle('CD8') +xlab("Genomic Region (5'->3')")+ ylab("Read Count Frequency") + theme_classic() + scale_color_manual(values=c('High'='#30A3CC', 'Low'='#FCB357'))

High_CD4 <- High_CD4$data
High_CD4$HighLow <- 'High'
Low_CD4 <- Low_CD4$data
Low_CD4$HighLow <- 'Low'
plotter_CD4<- rbind(High_CD4, Low_CD4)
ggplot(plotter_CD4, aes(x=pos, y=value, color=HighLow)) + geom_line() + 
xlim(-3000, 3000)+ ggtitle('CD4') + xlab("Genomic Region (5'->3')")+ ylab("Read Count Frequency") + theme_classic() + scale_color_manual(values=c('High'='#30A3CC', 'Low'='#FCB357'))

dev.off()




peakAnno_High_CD8 <- annotatePeak(peak_High_CD8, tssRegion=c(-3000, 3000),TxDb=annotation_Peaks, annoDb="org.Hs.eg.db")
peakAnno_High_CD4 <- annotatePeak(peak_High_CD4, tssRegion=c(-3000, 3000),TxDb=annotation_Peaks, annoDb="org.Hs.eg.db")
peakAnno_Low_CD8  <- annotatePeak(peak_Low_CD8, tssRegion=c(-3000, 3000),TxDb=annotation_Peaks, annoDb="org.Hs.eg.db")
peakAnno_Low_CD4  <- annotatePeak(peak_Low_CD4, tssRegion=c(-3000, 3000),TxDb=annotation_Peaks, annoDb="org.Hs.eg.db")


pdf('./Plots/GenAnnotationGR.pdf')
plotAnnoPie(peakAnno_High_CD8)
plotAnnoPie(peakAnno_High_CD4)
plotAnnoPie(peakAnno_Low_CD8)
plotAnnoPie(peakAnno_Low_CD4)
dev.off()



peakAnno_High_CD8_dfAnn <- as.data.frame(peakAnno_High_CD8)
peakAnno_High_CD4_dfAnn <- as.data.frame(peakAnno_High_CD4)
peakAnno_Low_CD8_dfAnn <- as.data.frame(peakAnno_Low_CD8)
peakAnno_Low_CD4_dfAnn <- as.data.frame(peakAnno_Low_CD4)


for(signature in signatures){
    message(signature)
    gene_list <- as.character(read.table(paste0(signatures_path, '/', signature))$V1)
    signature <- stringr::str_extract(signature, '[\\S]+(?=\\.)')
    pdf(paste0('./Plots/ATAC_',signature,'.pdf'))
   
    peakAnno_High_CD8_tmp <- makeGRangesFromDataFrame(peakAnno_High_CD8_dfAnn[peakAnno_High_CD8_dfAnn$SYMBOL %in% gene_list,], keep.extra.columns=TRUE)
    peakAnno_High_CD8_tmp <- getTagMatrix(peakAnno_High_CD8_tmp, windows=promoter)
    print(tagHeatmap(peakAnno_High_CD8_tmp, xlim=c(-3000, 3000), color="red", title=paste0('High_CD8_', signature)))
    
    peakAnno_High_CD4_tmp <- makeGRangesFromDataFrame(peakAnno_High_CD4_dfAnn[peakAnno_High_CD4_dfAnn$SYMBOL %in% gene_list,], keep.extra.columns=TRUE)
    peakAnno_High_CD4_tmp <- getTagMatrix(peakAnno_High_CD4_tmp, windows=promoter)
    print(tagHeatmap(peakAnno_High_CD4_tmp, xlim=c(-3000, 3000), color="red", title=paste0('High_CD4_', signature)))
    
    peakAnno_Low_CD8_tmp <- makeGRangesFromDataFrame(peakAnno_Low_CD8_dfAnn[peakAnno_Low_CD8_dfAnn$SYMBOL %in% gene_list,], keep.extra.columns=TRUE)  
    peakAnno_Low_CD8_tmp <- getTagMatrix(peakAnno_Low_CD8_tmp, windows=promoter)
    print(tagHeatmap(peakAnno_Low_CD8_tmp, xlim=c(-3000, 3000), color="red", title=paste0('Low_CD4_', signature)))
    
    peakAnno_Low_CD4_tmp <- makeGRangesFromDataFrame(peakAnno_Low_CD4_dfAnn[peakAnno_Low_CD4_dfAnn$SYMBOL %in% gene_list,])
    peakAnno_Low_CD4_tmp <- getTagMatrix(peakAnno_Low_CD4_tmp, windows=promoter)
    print(tagHeatmap(peakAnno_Low_CD4_tmp, xlim=c(-3000, 3000), color="red", title=paste0('Low_CD4_', signature)))

    tmp_High_CD8 <- plotAvgProf(peakAnno_High_CD8_tmp, xlim=c(-3000, 3000))
    tmp_Low_CD8  <- plotAvgProf(peakAnno_Low_CD8_tmp, xlim=c(-3000, 3000))

    tmp_High_CD8 <- tmp_High_CD8$data
    tmp_High_CD8$HighLow <- 'High'
    tmp_Low_CD8 <- tmp_Low_CD8$data
    tmp_Low_CD8$HighLow <- 'Low'
    plotter_CD8 <- rbind(tmp_High_CD8, tmp_Low_CD8)
    print(ggplot(plotter_CD8, aes(x=pos, y=value, color=HighLow)) + geom_line() + 
    xlim(-3000, 3000) + ggtitle(paste0('CD8  ',signature)) +xlab("Genomic Region (5'->3')")+ ylab("Read Count Frequency") + theme_classic() + scale_color_manual(values=c('High'='#30A3CC', 'Low'='#FCB357')))

    tmp_High_CD4 <- plotAvgProf(peakAnno_High_CD4_tmp, xlim=c(-3000, 3000))
    tmp_Low_CD4  <- plotAvgProf(peakAnno_Low_CD4_tmp, xlim=c(-3000, 3000))

    tmp_High_CD4 <- tmp_High_CD4$data
    tmp_High_CD4$HighLow <- 'High'
    tmp_Low_CD4 <- tmp_Low_CD4$data
    tmp_Low_CD4$HighLow <- 'Low'
    plotter_CD4 <- rbind(tmp_High_CD4, tmp_Low_CD4)
    print(ggplot(plotter_CD4, aes(x=pos, y=value, color=HighLow)) + geom_line() + 
    xlim(-3000, 3000)+ ggtitle(paste0('CD4  ', signature)) + xlab("Genomic Region (5'->3')")+ ylab("Read Count Frequency") + theme_classic() + scale_color_manual(values=c('High'='#30A3CC', 'Low'='#FCB357')))

    dev.off()
}



# Genome coverage with differential peaks
consensusPeaks_CD8 <- read.table('/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/results/cd8_Lowd0_vs_Highd0_csaw_denovo_trended_csaw-windows_significant.txt', sep="\t", header=T)
consensusPeaks_CD4 <- read.table('/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/results/Lowd0_vs_Highd0_csaw_denovo_trended_csaw-windows_significant.txt', sep="\t", header=T)

consensusPeaks_CD8 <- GenomicRanges::makeGRangesFromDataFrame(consensusPeaks_CD8, keep.extra.columns=TRUE)
consensusPeaks_CD4 <- GenomicRanges::makeGRangesFromDataFrame(consensusPeaks_CD4, keep.extra.columns=TRUE)

write.table(data.frame(consensusPeaks_CD4), file="/home/sevastopol/data/gserranos/CART_HL/Data/High_low/ATACseq/Diff_peaks_CD4.bed", quote=F, sep="\t", row.names=F, col.names=F)
write.table(data.frame(consensusPeaks_CD8), file="/home/sevastopol/data/gserranos/CART_HL/Data/High_low/ATACseq/Diff_peaks_CD8.bed", quote=F, sep="\t", row.names=F, col.names=F)


pdf('./Plots/peakLocationGenomeDifferential.pdf')
ChIPseeker::covplot(consensusPeaks_CD8, weightCol='logCPM', title = 'CD8')
ChIPseeker::covplot(consensusPeaks_CD4, weightCol='logCPM', title = 'CD4')
dev.off()



# TCR_Clonality
colors_ann <- list( 'HighLow' = c(High = '#30A3CC', Low = '#FCB357'))

bin_high_low$Patient <- stringr::str_extract(bin_high_low$cell_id, '^[0-9da]+')
bin_high_low$cell <- stringr::str_extract(bin_high_low$cell_id, '(?<=_)[A-Z0-9\\-]+')
folders <- list.files('/home/sevastopol/data/gserranos/CART_HL/Data/TCR_clonality/', pattern='day0')
plotter <- data.frame(Var1 = NULL, Freq= NULL, Sample = NULL, BinScore= NULL)
plotter_all <- data.frame(Var1 = NULL, Freq= NULL, Sample = NULL, BinScore= NULL)

for(folder in folders){
    print(folder)
    tcr <- read.table(paste0('/home/sevastopol/data/gserranos/CART_HL/Data/TCR_clonality/',folder,'/outs/filtered_contig_annotations.csv'), sep=',', header=TRUE)

    tcr_filtered <- unique(tcr[, c( 'barcode', 'raw_clonotype_id')])
    tcr_filtered <- tcr_filtered[tcr_filtered$raw_clonotype_id != '',]
    n_cells <- nrow(tcr_filtered)
    #check Top frequence Clonotype CDR3 Sequences
    print(max(table(tcr_filtered$raw_clonotype_id)))
    props <- as.data.frame(table(tcr_filtered$raw_clonotype_id))
    props <- props[props$Freq !=0, ]
    props$Freq <- (props$Freq/n_cells)*100
    props$Sample <- folder
    plotter_all <- rbind(plotter_all, props)

    bin_high_low_tmp <- bin_high_low[bin_high_low$cell %in% tcr_filtered$barcode,]
    tcr_filtered <- merge(tcr_filtered, bin_high_low_tmp[, c('cell', 'BinScore')], by.x='barcode', by.y='cell')
    tcr_filtered$raw_clonotype_id <- as.character(tcr_filtered$raw_clonotype_id)

    props <- as.data.frame(table(tcr_filtered$raw_clonotype_id))
    props$Freq <- (props$Freq/n_cells)*100
    props$Sample <- folder
    props <- merge(props, tcr_filtered[, c('raw_clonotype_id', 'BinScore')], by.x='Var1', by.y='raw_clonotype_id')
    plotter <- rbind(plotter, props)
}


pdf('./Plots/TCR_Clonality.pdf')
ggplot(plotter_all, aes(y=Freq, x=Sample)) + geom_point(alpha=0.7) + theme_classic()  + geom_jitter(width = 0.25, height = 0.001) + scale_color_manual(values=c('#30A3CC', '#FCB357'))
ggplot(plotter, aes(y=Freq, x=Sample, color=BinScore)) + geom_point(alpha=0.7) + theme_classic()  + geom_jitter(width = 0.25, height = 0.001) + scale_color_manual(values=c('#30A3CC', '#FCB357'))
ggplot(plotter, aes(y=Freq, x=Sample, fill=BinScore)) + geom_violin(alpha=0.7) + theme_classic()  + scale_fill_manual(values=c('#30A3CC', '#FCB357'))
dev.off()



all_folder_plotter <- data.frame(Clonotype=NULL, Freq=NULL, Sample=NULL)
for(folder in folders){
    print(folder)
    tcr <- read.table(paste0('/home/sevastopol/data/gserranos/CART_HL/Data/TCR_clonality/',folder,'/outs/filtered_contig_annotations.csv'), sep=',', header=TRUE)

    tcr_filtered <- unique(tcr[, c( 'barcode', 'raw_clonotype_id')])
    tcr_filtered <- tcr_filtered[grepl('^clonotype',tcr_filtered$raw_clonotype_id),]
    n_cells <- nrow(tcr_filtered)
    #check Top frequence Clonotype CDR3 Sequences
    print(max(table(tcr_filtered$raw_clonotype_id)))
    props <- as.data.frame(table(tcr_filtered$raw_clonotype_id))
    props <- props[order(-props$Freq),]
    props <- props[props$Freq !=0,]
    props_table <- setNames(as.data.frame(table(props$Freq)) , c('Clonotype_uniqueness', 'Frequency'))
    props_table$Percentage <- (props_table$Frequency/sum(props_table$Frequency)) *100
    props_table_High <- tcr_filtered[tcr_filtered$barcode %in%  bin_high_low[bin_high_low$BinScore == 'High', 'cell'],]
    props_table_High <- as.data.frame(table(props_table_High$raw_clonotype_id))
    props_table_High <- setNames(as.data.frame(table(props_table_High$Freq)) , c('Clonotype_uniqueness', 'Frequency'))
    props_table_High$Percentage <- (props_table_High$Frequency/sum(props_table_High$Frequency)) *100
    data_2_xlsx <- list('All'=props_table, 'High'=props_table_High)
    WriteXLS::WriteXLS(data_2_xlsx, ExcelFileName=paste0('./Plots/Clonality_',folder,'.xlsx'), SheetNames = names(data_2_xlsx),  col.names=TRUE, row.names=FALSE, BoldHeaderRow=TRUE)

    props_2_plot <- props[1:10,]
    props_2_plot$Var1 <- paste0('Clonotype', seq(1,nrow(props_2_plot)))
    other_clonotype <- tcr_filtered$raw_clonotype_id[!tcr_filtered$raw_clonotype_id %in% props_2_plot$Var1]

    props_2_plot <- rbind(props_2_plot, data.frame(Var1='Other', Freq=length(other_clonotype)))
    colnames(props_2_plot) <- c('Clonotype', 'Freq')
    props_2_plot$Sample <- folder
    all_folder_plotter <- rbind(all_folder_plotter, props_2_plot)


}
colors <- hues::iwanthue(10)
names(colors) <- grep('^Clonotype', props_2_plot$Clonotype, value=TRUE)

pdf('./Plots/TCR_all_Barplot.pdf')
colors <- c(colors, 'Other'='#d3d3d3')
ggplot(all_folder_plotter, aes(fill=Clonotype, y = Freq, x=Sample)) + geom_bar(position='fill', stat='identity') + scale_fill_manual(values=colors) + theme_classic()
dev.off()


#### Venn diagram

DE_peaks_CD4 <- read.delim("/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/results/Lowd0_vs_Highd0_csaw_denovo_trended_csaw-windows_significant_Annotated.txt", sep="\t", header=T)
DE_Genes_CD4 <- read.delim('/home/sevastopol/data/gserranos/CART_HL/Data/signature/BatchK_CD4_output_BASAL_BOTH_LowvsHigh.tsv',sep="\t",header=T)
# DE_Genes$GeneID[which(DE_Genes$GeneID %in% DE_peaks$SYMBOL)]

DE_peaks_CD8 <- read.delim("/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/results/cd8_Lowd0_vs_Highd0_csaw_denovo_trended_csaw-windows_significant_Annotated.txt", sep="\t", header=T)
DE_Genes_CD8 <- read.delim('/home/sevastopol/data/gserranos/CART_HL/Data/signature/BatchK_CD8_output_BASAL_BOTH_LowvsHigh.tsv',sep="\t",header=T)
# DE_Genes_CD8 <- read.delim('/home/sevastopol/data/gserranos/CART_HL/Data/genesrnaseq/BatchK_CD8_output_BASAL_BOTH_LowvsHigh.tsv',sep="\t",header=T)
# VennDiagram::venn.diagram(x = list(DE_Genes$GeneID, DE_peaks$SYMBOL), category.names = c('Bulk RNA', 'ATAC RNA'),
# filename = './Plots/Bulk_ATAC_venn.png', output=TRUE)

ven_list_CD4 = list('Bulk RNA' = unique(DE_Genes_CD4[DE_Genes_CD4$padj < 0.05 ,'GeneID']), 'ATAC RNA' = unique(DE_peaks_CD4$SYMBOL))
ven_list_CD8 = list('Bulk RNA' = unique(DE_Genes_CD8[DE_Genes_CD8$padj < 0.05 ,'GeneID']), 'ATAC RNA' = unique(DE_peaks_CD8$SYMBOL))
ven_list_CD4_flt = list('Bulk RNA' = unique(DE_Genes_CD4[DE_Genes_CD4$padj < 0.05 ,'GeneID']), 'ATAC RNA' = unique(DE_peaks_CD4[grepl('Promoter',DE_peaks_CD4$annotation), 'SYMBOL']))
ven_list_CD8_flt = list('Bulk RNA' = unique(DE_Genes_CD8[DE_Genes_CD8$padj < 0.05 ,'GeneID']), 'ATAC RNA' = unique(DE_peaks_CD8[grepl('Promoter',DE_peaks_CD8$annotation), 'SYMBOL']))


pdf('./Plots/Bulk_ATAC_venn.pdf')
cowplot::plot_grid(
ggvenn::ggvenn(ven_list_CD4,
    fill_color = hues::iwanthue(2),
    stroke_size = 0.4,
    set_name_size = 5,
    show_percentage = F,
    fill_alpha = 0.4,
    stroke_color = 'white',
    stroke_alpha = 1,
    stroke_linetype = 'solid',
    text_color = 'black',
    text_size = 6,
    label_sep = ','
) + ggtitle('CD4')+ theme(plot.title = element_text(hjust = 0.5)),
ggvenn::ggvenn(ven_list_CD8,
    fill_color = hues::iwanthue(2),
    stroke_size = 0.4,
    set_name_size = 5,
    show_percentage = F,
    fill_alpha = 0.4,
    stroke_color = 'white',
    stroke_alpha = 1,
    stroke_linetype = 'solid',
    text_color = 'black',
    text_size = 6,
    label_sep = ','
) + ggtitle('CD8')+ theme(plot.title = element_text(hjust = 0.5)),
ncol=2)
cowplot::plot_grid(

ggvenn::ggvenn(ven_list_CD4_flt,
    fill_color = hues::iwanthue(2),
    stroke_size = 0.4,
    set_name_size = 5,
    show_percentage = F,
    fill_alpha = 0.4,
    stroke_color = 'white',
    stroke_alpha = 1,
    stroke_linetype = 'solid',
    text_color = 'black',
    text_size = 6,
    label_sep = ','
) + ggtitle('CD4 Only Promoter')+ theme(plot.title = element_text(hjust = 0.5)),
ggvenn::ggvenn(ven_list_CD8_flt,
    fill_color = hues::iwanthue(2),
    stroke_size = 0.4,
    set_name_size = 5,
    show_percentage = F,
    fill_alpha = 0.4,
    stroke_color = 'white',
    stroke_alpha = 1,
    stroke_linetype = 'solid',
    text_color = 'black',
    text_size = 6,
    label_sep = ','
) + ggtitle('CD8 Only Promoter')+ theme(plot.title = element_text(hjust = 0.5)),
ncol=2)
dev.off()

write.table(intersect( unique(DE_Genes_CD4[DE_Genes_CD4$padj < 0.05 ,'GeneID']), unique(DE_peaks_CD4$SYMBOL) ), row.names=FALSE, col.names=FALSE, quote=FALSE, './Plots/CD4_DEG_Bulk_and_Atac_genes.txt')
write.table(intersect( unique(DE_Genes_CD8[DE_Genes_CD8$padj < 0.05 ,'GeneID']), unique(DE_peaks_CD8$SYMBOL) ), row.names=FALSE, col.names=FALSE, quote=FALSE, './Plots/CD8_DEG_Bulk_and_Atac_genes.txt')

# Heatmaps for SIMIC Bulk and ATAC
library("ChIPseeker")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("org.Hs.eg.db")
library(gridExtra)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene


TF_list <- c('NR4A1', 'BATF', 'ARID5A', 'RFX5', 'SATB1', 'ATF5', 'NR4A1', 'ZBTB7B', 'EED' )
simic_out <- read.table('/home/sevastopol/data/gserranos/CART_HL/SimiC/Data/data_2_nets_V2.tsv', sep='\t', header=TRUE)

for(CD in c('CD4', 'CD8')){
    message(CD)
    Bulk_signature <- read.delim(file=paste0('./Data/signature/BatchK_',CD,'_output_BASAL_BOTH_LowvsHigh.tsv'),sep="\t",header=T)
    results_atac <- read.table(paste0('/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/results/AVER_guille_',CD,'.txt'), sep='\t', header=TRUE)

    results_atac_ranges <- makeGRangesFromDataFrame(results_atac[,1:3])
    peakAnno <- annotatePeak(results_atac_ranges, tssRegion=c(-3000, 3000),
                            TxDb=txdb, annoDb="org.Hs.eg.db")
    patients_2_keep <- grep('Low|High', grep('d0',grep('^D[0-9]',names(results_atac), value=TRUE), value=TRUE), value=TRUE)
    results_atac_ann <- cbind(results_atac_ranges, as.data.frame(peakAnno),results_atac )

    pdf(paste0('./Plots/ATAC_',CD,'_SIMIC_phm.pdf'))
    for(TF in TF_list){
        message(TF)
        targets_tmp <- as.character(simic_out[simic_out$driver == TF, 'target'])
        results_atac_ann_tmp <- results_atac_ann[results_atac_ann$SYMBOL %in% targets_tmp, ]
        results_atac_ann_tmp$unique_ID <- paste0(results_atac_ann_tmp$SYMBOL, '_', results_atac_ann_tmp$start, '-',results_atac_ann_tmp$end)
        data_2_plot <- results_atac_ann_tmp[,patients_2_keep]
        rownames(data_2_plot) <- results_atac_ann_tmp$unique_ID 
        print(grid.arrange(arrangeGrob(grobs=list(get_heatmap_atac(data_2_plot, TF, FALSE)[[4]]))))

    }
    dev.off()

    pdf(paste0('./Plots/Bulk_',CD,'_SIMIC_phm.pdf'))
    for(TF in TF_list){
        targets_tmp <- as.character(simic_out[simic_out$driver == TF, 'target'])
        Bulk_signature_tmp <- Bulk_signature[Bulk_signature$GeneID %in% targets_tmp, ]
        data_2_plot <- Bulk_signature_tmp[, grepl('^D', colnames(Bulk_signature_tmp))]
        rownames(data_2_plot) <- Bulk_signature_tmp$GeneID
        print(grid.arrange(arrangeGrob(grobs=list(get_heatmap_bulk(data_2_plot, TF)[[4]]))))


    }
    dev.off()
}







# ATAC_peaks 



txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene 
gtrack <- GenomeAxisTrack()



all_bgwg <- list.files('/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/BigWig/')

high_bw <- '/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/BigWig/D14_High_CD8_d0.sort.rmdup.rmblackls.rmchr.norm.bw'
low_bw  <- '/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/BigWig/D14_Low_CD8_d0.sort.rmdup.rmblackls.rmchr.norm.bw'

HLA_DRA <- get_peaks(6, 32435305, 32449627, high_bw,  low_bw, filter_genes = c('HLA-DRA'))
CIITA   <- get_peaks(16, 10865209, 10942562, high_bw,  low_bw, filter_genes = c('CIITA'))
TNFRSF9 <- get_peaks(1, 7901385, 7955324, high_bw,  low_bw, filter_genes = c('TNFRSF9'))
MAP3K8  <- get_peaks(10, 30388299, 30507554, high_bw,  low_bw, filter_genes = c('MAP3K8'))
BATF3   <- get_peaks(1, 212662788, 212723568, high_bw,  low_bw, filter_genes = c('BATF3'))
TNFSF4  <- get_peaks(1, 173100132, 173366893, high_bw,  low_bw, filter_genes = c('TNFSF4'))
CTLA4   <- get_peaks(2, 203862006, 203878398, high_bw,  low_bw, filter_genes = c('CTLA4'))
HELLS   <- get_peaks(10, 94493840, 94634098, high_bw,  low_bw, filter_genes = c('HELLS'))

HLA_DRA <- get_peaks(6, 32432985, 32449203, '/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/BigWig/D14_High_CD8_d0.sort.rmdup.rmblackls.rmchr.norm.bw',  '/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/BigWig/D8a_Low_CD8_S8.sort.rmdup.rmblackls.rmchr.norm.bw', filter_genes = c('HLA-DRA'))
CTLA4   <- get_peaks(2, 203865466, 203876964, '/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/BigWig/D14_High_CD8_d0.sort.rmdup.rmblackls.rmchr.norm.bw',  '/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/BigWig/D9_Low_CD8_S14.sort.rmdup.rmblackls.rmchr.norm.bw', filter_genes = c('CTLA4'), MAX=4)
MCM3   <- get_peaks(6, 52254317, 52319190, '/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/BigWig/D14_High_CD8_d0.sort.rmdup.rmblackls.rmchr.norm.bw',  '/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/BigWig/D9_Low_CD8_S14.sort.rmdup.rmblackls.rmchr.norm.bw', filter_genes = c('MCM3'), MAX=5)



# CHR <- 16
# START<- 10865209
# END <- 10942562
# HIGH <- '/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/BigWig/D14_High_CD8_d0.sort.rmdup.rmblackls.rmchr.norm.bw'
# LOW <- '/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/BigWig/D14_Low_CD8_d0.sort.rmdup.rmblackls.rmchr.norm.bw'

# file.remove(CIITA)
