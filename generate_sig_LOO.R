

get_annotation <- function(data){
    annotation <- data.frame(FACS_Level= stringr::str_extract(colnames(data), '(?<=_)[a-zA-Z]+(?=$|_)') )
    rownames(annotation) <- colnames(data)
    return(annotation)
}


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
    DT.high <- cbind(DT.high, ann_markers[, c('GeneID', 'log2FoldChange', 'padj')])
    DT.cd4.high <- setNames(as.data.frame(colSums(DT.high[, !colnames(DT.high) %in% c('GeneID', 'log2FoldChange', 'padj')], )), c('High_pondered'))
    DT.cd4.high$Patient_ID <- rownames(DT.cd4.high)
   
    annotation$High_pondered <- DT.cd4.high$High_pondered
    return(annotation)
}


leave_one_out_sample <- function(counts, colData, normalized_counts){
    for (sample in colnames(counts)){
        print(paste0('Leaving out:  ', sample) )
        countData_tmp <- counts[, colnames(counts) != sample]
        colData_tmp <- colData[colData$SampleName != sample,]

        dds <- DESeq2::DESeqDataSetFromMatrix(countData = countData_tmp,
                                            colData = colData_tmp,
                                            design = ~Donor + CellType)


        # NOW::::Run the differential expression pipeline.
        dds <- DESeq2::DESeq(dds)
        # Build the results table. We can set the alpha signification level and 
        # a lfcThreshold.
        res <- DESeq2::results(dds, alpha = 0.05, parallel = TRUE)
        res <- DESeq2::results(dds, contrast=c("CellType","High","Low"), alpha = 0.05, lfcThreshold = 1)
        table(res$padj < 0.05)

        caseB <- "Low"
        caseA <- "High"
        res <- DESeq2::results(dds, contrast=c("CellType","High","Low"), alpha = 0.05, lfcThreshold = 1)

        res.shr <- DESeq2::lfcShrink(dds=dds, contrast=c("CellType","High","Low"), type="ashr")

        resOrdered <- res.shr[order(res.shr$pvalue),]
        resOrdered <- merge(as.data.frame(resOrdered), as.data.frame(DESeq2::counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
        names(resOrdered)[1] <- "GeneID"
        row.names(resOrdered) <- resOrdered$GeneID

        # write.table(data.frame(resOrdered),
        #             file = file.path(FIGURES_DIR,paste0("BatchK_CD4_output_BASAL_BOTH_", caseB, "vs", caseA, ".tsv")),
        #             sep='\t',
        #             row.names=F,
        #             quote=F)

        sigGenes_symbol <- rownames(resOrdered[abs(resOrdered$log2FoldChange)>1  
                                            & resOrdered$padj < 0.05 & !(is.na(resOrdered$padj)), ])
        top50Genes <- sigGenes_symbol[1:50]
        sigGenes_symbol <- sigGenes_symbol[!(is.na(sigGenes_symbol))]
        top50Genes <- top50Genes[!(is.na(top50Genes))]
        # write.table(data.frame(sigGenes_symbol),
        #             file = file.path(FIGURES_DIR,paste0("BatchK_CD4_SIGGenes_BASAL_BOTH_", caseB, "vs", caseA, ".tsv")),
        #             sep='\t',
        #             row.names=F,
        #             quote=F)

        res.df <- data.frame(resOrdered)
        res2.cd4 <- res.df[res.df$GeneID %in% sigGenes_symbol, c('GeneID', 'log2FoldChange', 'padj')]
        dataset <- normalized_counts[rownames(normalized_counts) %in% sigGenes_symbol, ]

        logplusone <- function(x) {log(x + 0.5)}
        l <- as.data.frame(apply(dataset, 2, logplusone))
        # head(l[,1:5])
        zscore <- as.data.frame(scale(l))

        zscore$GeneID <- rownames(zscore)
        ann_markers <- merge(zscore, res2.cd4, by='GeneID')

        DT.high <- ann_markers[, !colnames(ann_markers) %in% c('GeneID', 'log2FoldChange', 'padj'), drop=FALSE] * ann_markers[, 'log2FoldChange']
        DT.high <- cbind(DT.high, ann_markers[, c('GeneID', 'log2FoldChange', 'padj')])
        DT.cd4.high <- sum(DT.high[, !colnames(DT.high) %in% c('GeneID', 'log2FoldChange', 'padj')])
        DT.cd4.high_2 <- data.frame(t(unlist(colSums(DT.high[, !colnames(DT.high) %in% c('GeneID', 'log2FoldChange', 'padj')]))))
        
        if (!exists('results') ){
            results <- cbind(data.frame(Sample=sample, type=cd_type), DT.cd4.high_2)
        }else{
            results <- rbind(results, cbind(data.frame(Sample=sample, type=cd_type), DT.cd4.high_2))
        }
        # results <- rbind(results, data.frame(Sample=sample, Score= DT.cd4.high, type=cd_type ))
    }
    return(results)
}


leave_one_out <- function(counts, colData, normalized_counts){
    samples <-unique(stringr::str_extract(colnames(counts), '^D[\\d]{1,2}'))
    for (sample in samples){
        print(paste0('Leaving out:  ', sample) )
        countData_tmp <- countData[, !grepl(paste0('^',sample) ,colnames(counts) )]
        colData_tmp <- colData[!grepl(paste0('^',sample) ,colData$SampleName),]

        dds <- DESeq2::DESeqDataSetFromMatrix(countData = countData_tmp,
                                            colData = colData_tmp,
                                            design = ~Donor + CellType)


        # NOW::::Run the differential expression pipeline.
        dds <- DESeq2::DESeq(dds)
        # Build the results table. We can set the alpha signification level and 
        # a lfcThreshold.
        res <- DESeq2::results(dds, alpha = 0.05, parallel = TRUE)
        res <- DESeq2::results(dds, contrast=c("CellType","High","Low"), alpha = 0.05, lfcThreshold = 1)
        table(res$padj < 0.05)

        caseB <- "Low"
        caseA <- "High"
        res <- DESeq2::results(dds, contrast=c("CellType","High","Low"), alpha = 0.05, lfcThreshold = 1)

        res.shr <- DESeq2::lfcShrink(dds=dds, contrast=c("CellType","High","Low"), type="ashr")

        resOrdered <- res.shr[order(res.shr$pvalue),]
        resOrdered <- merge(as.data.frame(resOrdered), as.data.frame(DESeq2::counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
        names(resOrdered)[1] <- "GeneID"
        row.names(resOrdered) <- resOrdered$GeneID
        # write.table(data.frame(resOrdered),
        #             file = file.path(FIGURES_DIR,paste0("BatchK_CD4_output_BASAL_BOTH_", caseB, "vs", caseA, ".tsv")),
        #             sep='\t',
        #             row.names=F,
        #             quote=F)

        sigGenes_symbol <- rownames(resOrdered[abs(resOrdered$log2FoldChange)>1  
                                            & resOrdered$padj < 0.05 & !(is.na(resOrdered$padj)), ])
        top50Genes <- sigGenes_symbol[1:50]
        sigGenes_symbol <- sigGenes_symbol[!(is.na(sigGenes_symbol))]
        top50Genes <- top50Genes[!(is.na(top50Genes))]
        # write.table(data.frame(sigGenes_symbol),
        #             file = file.path(FIGURES_DIR,paste0("BatchK_CD4_SIGGenes_BASAL_BOTH_", caseB, "vs", caseA, ".tsv")),
        #             sep='\t',
        #             row.names=F,
        #             quote=F)

        res.df <- data.frame(resOrdered)
        res2.cd4 <- res.df[res.df$GeneID %in% sigGenes_symbol, c('GeneID', 'log2FoldChange', 'padj')]
        dataset <- normalized_counts[rownames(normalized_counts) %in% sigGenes_symbol, ]

        logplusone <- function(x) {log(x + 0.5)}
        l <- as.data.frame(apply(dataset, 2, logplusone))
        # head(l[,1:5])
        zscore <- as.data.frame(scale(l))

        zscore$GeneID <- rownames(zscore)
        ann_markers <- merge(zscore, res2.cd4, by='GeneID')

        DT.high <- ann_markers[, !colnames(ann_markers) %in% c('GeneID', 'log2FoldChange', 'padj'), drop=FALSE] * ann_markers[, 'log2FoldChange']
        DT.high <- cbind(DT.high, ann_markers[, c('GeneID', 'log2FoldChange', 'padj')])
        DT.cd4.high <- sum(DT.high[, !colnames(DT.high) %in% c('GeneID', 'log2FoldChange', 'padj')])
        DT.cd4.high_2 <- data.frame(t(unlist(colSums(DT.high[, !colnames(DT.high) %in% c('GeneID', 'log2FoldChange', 'padj')]))))
        
        if (!exists('results') ){
            results <- cbind(data.frame(Sample=sample, type=cd_type), DT.cd4.high_2)
        }else{
            results <- rbind(results, cbind(data.frame(Sample=sample, type=cd_type), DT.cd4.high_2))
        }
        # results <- rbind(results, data.frame(Sample=sample, Score= DT.cd4.high, type=cd_type ))
    }
    return(results)
}



CAR_NAME <- 'CAR_pCCL_BCMA'
paths = c("/home/sevastopol/data/mcallejac/RNA_HighLow_ALL/results/HL_CD4.RData", 
"/home/sevastopol/data/mcallejac/RNA_HighLow_ALL/results/HL_CD8.RData")
for (path in paths){
    print(path)
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
    tmp <- leave_one_out(countData, colData, normalized_counts)
    # Add the control (score from unmodified signature)
    switch(cd_type, 
        CD8 = {
            signature_genes <-read.delim(file="./Data/signature/BatchK_CD8_SIGGenes_BASAL_BOTH_LowvsHigh.tsv",sep="\t",header=T)
            signature<-read.delim(file="./Data/signature/BatchK_CD8_output_BASAL_BOTH_LowvsHigh.tsv",sep="\t",header=T)

        }, 
        CD4 ={
            signature_genes <-read.delim(file="./Data/signature/BatchK_CD4_SIGGenes_BASAL_BOTH_LowvsHigh.tsv",sep="\t",header=T)
            signature<-read.delim(file="./Data/signature/BatchK_CD4_output_BASAL_BOTH_LowvsHigh.tsv",sep="\t",header=T)

        }
    )
    annotation <- data.frame(Sample = colnames(normalized_counts))
    control <- apply_signature(normalized_counts, signature, signature_genes, annotation)
    control <- setNames(control, c('V1', 'V2'))
    control <- rbind(control, data.frame(V1=c('Sample', 'type'), V2=c('Control', cd_type)))
    rownames(control) <- control$V1
    control <- t(control[, 'V2', drop=FALSE])
    tmp <- rbind(tmp, control)
    tmp[, grepl('^D', colnames(tmp))] <- sapply(tmp[, grepl('^D', colnames(tmp))], as.numeric)
    assign(paste0('results_', cd_type),  tmp)  
  
}



rownames(results_CD8) <- ifelse(grepl('^Control', results_CD8$Sample), 'Control', paste0('Leaved_out_',results_CD8$Sample))
results_CD8 <- results_CD8[, !colnames(results_CD8) %in% c('Sample', 'type')]
rownames(results_CD4) <-ifelse(grepl('^Control', results_CD4$Sample), 'Control', paste0('Leaved_out_',results_CD4$Sample))
results_CD4 <- results_CD4[, !colnames(results_CD4) %in% c('Sample', 'type')]

results_CD4_bin <- ifelse(results_CD4 >0, 1,0)

annot_colors <- list(FACS_Level = c( High = '#C63F32' , Low = '#5074AF'))
annot <- get_annotation(results_CD4)

pdf('./Plots/Signature_LOO_CD4.pdf', width=10)
pheatmap::pheatmap(results_CD4, cluster_cols=TRUE, cluster_rows=FALSE, scale='row', color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name =
       "RdYlBu")))(100), annotation_col = annot, annotation_colors = annot_colors, cutree_cols = 2, display_numbers=TRUE, main='Scaled score', cellwidth = 30)
pheatmap::pheatmap(results_CD4, cluster_cols=TRUE, cluster_rows=FALSE, scale='none', color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name =
       "RdYlBu")))(100), annotation_col = annot, annotation_colors = annot_colors, cutree_cols = 2, display_numbers=TRUE, main='Raw score', cellwidth = 35)
dev.off()


annot <- get_annotation(results_CD8)

pdf('./Plots/Signature_LOO_CD8.pdf', width=10)
pheatmap::pheatmap(results_CD8, cluster_cols=TRUE, cluster_rows=FALSE, scale='row', color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name =
       "RdYlBu")))(100), annotation_col = annot, annotation_colors = annot_colors, cutree_cols = 2, display_numbers=TRUE, main='Scaled score', cellwidth = 30)
pheatmap::pheatmap(results_CD8, cluster_cols=TRUE, cluster_rows=FALSE, scale='none', color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name =
       "RdYlBu")))(100), annotation_col = annot, annotation_colors = annot_colors, cutree_cols = 2, display_numbers=TRUE, main='Raw score', cellwidth = 35)
dev.off()


results$BinScore <- ifelse(results$Score >0, 'High', 'Low')

results$ScoreSignif <- signif(results$Score, digits = 3)
write.table(results, '/home/sevastopol/data/gserranos/CART_HL/Data/LOO_results.tsv',sep='\t', quote=FALSE, row.names=FALSE)

pdf('./Plots/Tets.pdf')
ggplot2::ggplot(results, ggplot2::aes(x=BinScore, y=Score, fill=BinScore)) + ggplot2::geom_boxplot() + ggplot2::facet_wrap(~type) +
ggplot2::scale_fill_manual(values = c('#DBBE78', '#7F7F7F')) + ggplot2::theme_classic()
dev.off()
