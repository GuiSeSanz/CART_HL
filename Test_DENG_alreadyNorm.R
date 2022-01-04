
library(Seurat)
library(AnnotationHub)
library(ggplot2)
library(ggraph)
library(ensembldb)
library(harmony)

get_expression_signature <- function(gene_list_name, df, coords = coords, upplim=NULL, split_HL=FALSE){
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

    gene_list <- as.character(read.table(paste0(signatures_path, '/', gene_list_name))$V1)
    gene_list <- check_genes(gene_list, df)
    tmp <- df[, colnames(df) %in% gene_list]
    tmp$Signature <- rowMeans(tmp)
    # tmp$cell_id <- sub('-', '\\.',rownames(tmp))
    tmp <- merge(tmp, coords, by=0)
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

get_umap_signature <- function(df, title){
    plot <- ggplot(df, aes(x=UMAP_1, y=UMAP_2, color= Signature)) + 
    geom_point(alpha=0.9, size = 0.8) + 
    # scale_color_gradientn(low="grey90", high ="blue", name = 'Expression') + 
    # guides(color = guide_colourbar(barwidth = 0.5, barheight = 2, label = TRUE)) + 
    ggprism::theme_prism() + ggtitle(title)+  labs(x = "UMAP 1", y = 'UMAP 2') +
    theme(legend.position='none', plot.title = element_text(hjust = 0.5, size = 12),  axis.text=element_blank(), axis.ticks=element_blank())#, legend.key.height = unit(2, 'mm'), legend.key.width = unit(1, 'mm')) + 
    if(is.numeric(df$Signature)){
          plot <- plot + viridis::scale_color_viridis()
    }else{
      plot <- plot + scale_color_manual(values=c('High'='#30A3CC', 'Low'='#d3d3d3'))
    }
    return(plot)
}

get_umap_HL <- function(df, title){
    plot <- ggplot() + 
    geom_point(df[df$Signature == 'Low',], mapping = aes(x=UMAP_1, y=UMAP_2, color= Signature),alpha=0.6) +
    geom_point(df[df$Signature == 'High',], mapping = aes(x=UMAP_1, y=UMAP_2, color= Signature),alpha=0.8) + 
    geom_point(alpha=0.9, size = 0.8) + 
    ggprism::theme_prism() + ggtitle(title) +  labs(x = "UMAP 1", y = 'UMAP 2') +
    theme(legend.position='none', plot.title = element_text(hjust = 0.5, size = 12),  axis.text=element_blank(), axis.ticks=element_blank()) + scale_color_manual(values=c('High'='#30A3CC', 'Low'='#d3d3d3'))
    return(plot)
}

get_umap <- function(data, gene, upplim=NULL ){
    gene_expr <- data[, gene, drop=FALSE]
    if(!is.null(upplim)){
        gene_expr[gene_expr>upplim] <- upplim
    }
    gene_expr$cell_id <- sub('-','.', rownames(gene_expr))
    coords_markers <- merge(coords, gene_expr, by=0)
    coords_markers$orderRank <- rank(coords_markers[,gene], ties.method="first")
    plot <- ggplot(coords_markers, aes(x=UMAP_1, y=UMAP_2, color= get(gene)), order=orderRank) + 
    geom_point(alpha=0.9, size = 0.6) + 
    # viridis::scale_color_viridis()+
    scale_color_gradient(low="grey90", high ="blue", name = 'Expression') + 
    guides(color = guide_colourbar(barwidth = 0.5, barheight = 2, label = TRUE)) + 
    theme_void() + ggtitle(gene)+
    theme(legend.position='none', plot.title = element_text(hjust = 0.5, size = 8))
	#, legend.key.height = unit(2, 'mm'), legend.key.width = unit(1, 'mm')) + 
    return(plot)
}

plot_stats <- function(obj, pdf_name){
  pdf(paste0(PLOT_PATH, pdf_name,'.pdf'))
    # shows the distribution of the transcripts per cells
    print(ggplot(obj@meta.data,aes(color=Phenotype, x=nUMI, fill= Phenotype)) + 
        geom_density(alpha = 0.2) + 
        scale_x_log10() + 
        theme_classic() +
        ylab("Cell density") +
        geom_vline(xintercept = 500) + ggtitle('UMISpercell'))
    print(ggplot(obj@meta.data, aes(color=Phenotype, x=nGene, fill= Phenotype)) + 
        geom_density(alpha = 0.2) + 
        theme_classic() +
        scale_x_log10() + 
        geom_vline(xintercept = 300) + ggtitle('GENESpercell'))
    # Visualize the distribution of genes detected per cell via boxplot
    print(ggplot(obj@meta.data, aes(x=Phenotype, y=log10(nGene), fill=Phenotype)) + 
        geom_boxplot() + 
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
        theme(plot.title = element_text(hjust=0.5, face="bold")) +
        ggtitle("NCells vs NGenes"))
    # correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
    p <- ggplot(obj@meta.data, aes(x=nUMI, y=nGene, color=mitoRatio)) + 
        geom_point() + 
        scale_colour_gradient(low = "gray90", high = "black", limits=c(0,1)) +
        stat_smooth(method=lm) +
        scale_x_log10() + 
        scale_y_log10() + 
        theme_classic() +
        geom_vline(xintercept = 500) +
        geom_hline(yintercept = 250) + ggtitle(paste0('UMIvsGENESpercell  Ncell:: ', ncol(obj)))
    print(p)
    print(p + facet_wrap(~Phenotype))
    print(VlnPlot(obj, features = c("nGene", "nUMI", "mitoRatio", "RPSRatio"), ncol = 4, pt.size=0))
  dev.off()
}

get_CellCycleGenes_per_PCs <- function(PCAs){
  plotter <- data.frame(PC=NULL, CellCycleGenes=NULL)
  for(pc in colnames(PCAs)[1:20]){
    b <- names(sort(abs(PCAs[,pc]), decreasing=TRUE)[1:50])
    plotter <- rbind(plotter,data.frame(PC=pc, CellCycleGenes=sum(b %in% c(g2m_genes, s_genes))))
  }
  p <- ggplot(plotter, aes(x=PC, y=CellCycleGenes, fill=PC)) + geom_bar(alpha=0.8, stat="identity") + ggprism::theme_prism()  + theme(legend.position='none', axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + scale_fill_manual(values=c(ggthemes::tableau_color_pal('Classic 20')(20) )) + ylim(0,25)
  return(p)
}

how_many_PCs <- function(obj, pdf_name){
  # determine the correct PC number
  # Determine percent of variation associated with each PC
  pct <- obj[["pca"]]@stdev / sum(obj[["pca"]]@stdev) * 100
  # Calculate cumulative percents for each PC
  cumu <- cumsum(pct)
  # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
  co1 <- which(cumu > 90 & pct < 5)[1]
  # Determine the difference between variation of PC and subsequent PC
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  # last point where change of % of variation is more than 0.1%.
  # Minimum of the two calculation
  pcs <- min(co1, co2)
  plot_df <- data.frame(pct = pct, 
           cumu = cumu, 
           rank = 1:length(pct))
  # Elbow plot to visualize 
  pdf(paste0(PLOT_PATH ,pdf_name,'.pdf'))
  print(ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw())
  dev.off()
  return(pcs)
}

get_cell_distribution <- function(dataset){
	n_clusters <- nrow(dataset)
	if( n_clusters < 20){
		colors <- ggthemes::tableau_color_pal('Classic 20')(n_clusters)
	}else{
		colors <- colorRampPalette(ggthemes::tableau_color_pal('Classic 20')(20))(n_clusters)
	}
	p <- ggplot(dataset, aes(x=Cluster, y=Ncells, fill=Cluster)) + geom_bar(alpha=0.8, stat="identity") + 
	ggprism::theme_prism()  + theme(legend.position='none') + scale_fill_manual(values=c(colors))
	return(p)
}

get_cells_per_cluster <- function(obj){
  n_cells <- FetchData(obj, vars = c(paste0('SCT_snn_res.', seq(0.2,1.2, by=0.2)))) 
  n_cells <- apply(n_cells, 2, table)
  results <- data.frame(Res=NULL, Cluster=NULL, Ncells=NULL)
  for (name in names(n_cells)){
    tmp <- setNames(as.data.frame(t(n_cells[[name]])), c('Res', 'Cluster', 'Ncells'))
    tmp$Res <- name
    results <- rbind(results, tmp)
  }
  return(results)
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


# cell_cycle_genes <- read.csv('/home/sevastopol/data/gserranos/CART_HL/Data/cell_cycle_genes_hsa.csv')
# # get the names of the genes corresponding to the different cell proliferation states
# ah <- AnnotationHub() # Connect to AnnotationHub
# ahDb <- query(ah, pattern = c("Homo sapiens", "EnsDb"), ignore.case = TRUE) # Access the Ensembl database for organism
# id <- tail(rownames(mcols(ahDb)),n=1) # Acquire the latest annotation files
# edb <- ah[[id]] # Download the appropriate Ensembldb database
# annotations <- genes(edb, return.type = "data.frame") # Extract gene-level information from database
# # Select annotations of interest
# annotations <- annotations %>% dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)
# # Get gene names for Ensembl IDs for each gene
# cell_cycle_markers <- dplyr::left_join(cell_cycle_genes, annotations, by = c("geneID" = "gene_id"))
# s_genes <- cell_cycle_markers[cell_cycle_markers$phase == 'S', 'gene_name']# Acquire the S phase genes
# g2m_genes <- cell_cycle_markers[cell_cycle_markers$phase == 'G2/M', 'gene_name'] # Acquire the G2M phase genes   
# # dbDisconnect()



s_genes   <- c("UBR7","RFC2","RAD51","MCM2","TIPIN","MCM6","UNG","POLD3","WDR76",
"CLSPN","CDC45","CDC6","MSH2","MCM5","POLA1","MCM4","RAD51AP1","GMNN","RPA2",
"CASP8AP2","HELLS","E2F8","GINS2","PCNA","NASP","BRIP1","DSCC1","DTL","CDCA7",
"CENPU","ATAD2","CHAF1B","USP1","SLBP","RRM1","FEN1","RRM2","EXO1","CCNE2",
"TYMS","BLM","PRIM1","UHRF1")
g2m_genes <- c("NCAPD2","ANLN","TACC3","HMMR","GTSE1","NDC80","AURKA","TPX2",
"BIRC5","G2E3","CBX5","RANGAP1","CTCF","CDCA3","TTK","SMC4","ECT2","CENPA",
"CDC20","NEK2","CENPF","TMPO","HJURP","CKS2","DLGAP5","PIMREG","TOP2A","PSRC1",
"CDCA8","CKAP2","NUSAP1","KIF23","KIF11","KIF20B","CENPE","GAS2L3","KIF2C",
"NUF2","ANP32E","LBR","MKI67","CCNB2","CDC25C","HMGB2","CKAP2L","BUB1","CDK1",
"CKS1B","UBE2C","CKAP5","AURKB","CDCA2","TUBB4B","JPT1")


###################################################
# USING THE DATA ALREDY NORMALIZED
###################################################


CD4_signature_genes <-read.delim(file="./Data/signature/BatchK_CD4_SIGGenes_BASAL_BOTH_LowvsHigh.tsv",sep="\t",header=T)
CD8_signature_genes <-read.delim(file="./Data/signature/BatchK_CD8_SIGGenes_BASAL_BOTH_LowvsHigh.tsv",sep="\t",header=T)
CD4_signature <- read.delim(file="./Data/signature/BatchK_CD4_output_BASAL_BOTH_LowvsHigh.tsv",sep="\t",header=T)
CD8_signature <- read.delim(file="./Data/signature/BatchK_CD8_output_BASAL_BOTH_LowvsHigh.tsv",sep="\t",header=T)

data_path = '/home/sevastopol/data/gserranos/CART_HL/Data/antiCD19/GSE151511_RAW/'
SAMPLE_NAMES <- list.files(data_path, pattern = '^ac[0-9]{2}')


# SAMPLE <- 'ac01'

all_seurat <- list()
for (SAMPLE in SAMPLE_NAMES){
	message(SAMPLE)
	PLOT_PATH <- paste0(getwd(), '/Plots/DENG_Tests/',SAMPLE,'/')
	dir.create(PLOT_PATH, showWarnings = FALSE)

	data_cd4 <- readRDS(paste0(data_path, SAMPLE, '/',SAMPLE,'_normalized_only_car_cd4.rds'))
	data_cd8 <- readRDS(paste0(data_path, SAMPLE, '/',SAMPLE,'_normalized_car_cd8.rds'))
	seurat_CD4 <- CreateSeuratObject(counts = data_cd4, min.features = 100, project = 'CD4') 
	seurat_CD8 <- CreateSeuratObject(counts = data_cd8, min.features = 100, project = 'CD8') 

	merged_seurat <- merge(x = seurat_CD4, y = seurat_CD8 , add.cell.id = c('CD4', 'CD8'),  project='CART_HL')

	merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)
	merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^MT-")
	merged_seurat$RPSRatio  <- PercentageFeatureSet(object = merged_seurat, pattern = "^RP[SL]")
	merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100 
	merged_seurat$RPSRatio  <- merged_seurat@meta.data$RPSRatio / 100 
	merged_seurat$Phenotype <- merged_seurat$orig.ident
	metadata <- setNames(merged_seurat@meta.data, c('orig.ident', 'nUMI', 'nGene' ,'log10GenesPerUMI', 'mitoRatio', 'RPSRatio', 'Phenotype'))
	merged_seurat@meta.data <- metadata
	plot_stats(merged_seurat, paste0(SAMPLE, 'Prefilter'))
	# Filter out low quality reads using selected thresholds - these will change with experiment
	filtered_seurat <- subset(x = merged_seurat, 
	            subset= 
	            (nUMI >= 500) & 
	            (nGene >= 250) & 
	            (nGene <= 5000) &
	            (log10GenesPerUMI > 0.80) & 
	            (mitoRatio < 0.20)) 
	# Output a logical vector for every gene on whether the more than zero counts per cell
	# Extract counts
	counts <- GetAssayData(object = filtered_seurat, slot = "counts")
	# Output a logical vector for every gene on whether the more than zero counts per cell
	nonzero <- counts > 0
	# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
	keep_genes <- Matrix::rowSums(nonzero) >= 10
	# Only keeping those genes expressed in more than 10 cells
	filtered_counts <- counts[keep_genes, ]
	# Reassign to filtered Seurat object
	filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)
	plot_stats(filtered_seurat, paste0(SAMPLE, 'Postfilter'))

	filtered_seurat <- FindVariableFeatures(filtered_seurat, min.features = 100)
	filtered_seurat <- ScaleData(filtered_seurat)
	message('RUN PCA')
	filtered_seurat <- RunPCA(object = filtered_seurat)
	# pdf(paste0(PLOT_PATH ,SAMPLE, '_dims2umap.pdf'), width=12, height=10)
	#   ElbowPlot(filtered_seurat)
	#   DimHeatmap(filtered_seurat, dims = 1:10, cells = 500, balanced = TRUE)
	#   DimHeatmap(filtered_seurat, dims = 11:20, cells = 500, balanced = TRUE)
	#   VizDimLoadings(filtered_seurat, dims = 1:10, reduction = "pca")
	#   VizDimLoadings(filtered_seurat, dims = 11:20, reduction = "pca")
	# dev.off()
	message('get PCA dims')
	PCA_dims <- how_many_PCs(filtered_seurat, paste0(SAMPLE, 'NumPCs'))
	message('get Neighbors')
	filtered_seurat <- FindNeighbors(object = filtered_seurat, dims = 1:PCA_dims)
	# Determine the clusters for various resolutions                                
	filtered_seurat <- FindClusters(object = filtered_seurat, resolution = c(0.2,0.4, 0.6, 0.8, 1.0, 1.2))
	  message('RUN UMAP')
	filtered_seurat <- RunUMAP(filtered_seurat, reduction = "pca", dims = 1:PCA_dims)

	# pdf(paste0(PLOT_PATH ,SAMPLE, '_Umap_Stats.pdf'), height = 8)
	#   clustree::clustree(filtered_seurat@meta.data, prefix = "SCT_snn_res.")
	#   DimPlot(filtered_seurat,  reduction = "umap", group.by = "orig.ident",  cols = c(ggthemes::tableau_color_pal('Classic 10 Medium')(3))) 
	#   features = c('nGene', 'nUMI', 'mitoRatio', 'RPSRatio', 'S.Score', 'G2M.Score')
	#   FeaturePlot(filtered_seurat, features = features, order = TRUE) & NoLegend()
	#   features = c('CD8A', 'CD3E', 'CD4', 'CD3D', 'IL7R', 'CCR7', 'SELL', 'CD69', 'FMC63-CD19SCFV')
	#   FeaturePlot(filtered_seurat, features = features, order = TRUE, min.cutoff = 'q10')& NoLegend()
	# dev.off()

	signature_cd4 <- apply_signature(data_cd4, CD4_signature, CD4_signature_genes, 'CD4', 'FMC63-CD19SCFV')
	signature_cd8 <- apply_signature(data_cd8, CD8_signature, CD8_signature_genes, 'CD8', 'FMC63-CD19SCFV')
	signature_cd8 <- signature_cd8[,c('cell_id','High_pondered')]
	rownames(signature_cd8) <- paste0('CD8_',signature_cd8$cell_id)
	signature_cd4 <- signature_cd4[,c('cell_id','High_pondered')]
	rownames(signature_cd4) <- paste0('CD4_',signature_cd4$cell_id)

	signature <- rbind(signature_cd8, signature_cd4)
	signature <- signature[colnames(filtered_seurat),]

	filtered_seurat$signature <- signature[, 'High_pondered', drop=F]
	filtered_seurat$signature_bin <- ifelse(filtered_seurat$signature > 0, 'High', 'Low')

	# pdf(paste0(PLOT_PATH ,SAMPLE, '_HL.pdf'), height = 8)
	#   DimPlot(filtered_seurat,  reduction = "umap", group.by = "signature_bin",  cols = c(ggthemes::tableau_color_pal('Classic 10 Medium')(3))) 
	#   FeaturePlot(filtered_seurat, features = 'signature', order = TRUE)& NoLegend()
	# dev.off()
	filtered_seurat$sample <- SAMPLE
	all_seurat[[SAMPLE]] <- filtered_seurat
}



for (idx in seq(1, length(SAMPLE_NAMES), by = 2)){
	print(idx)
	name1 <- SAMPLE_NAMES[idx]
	name2 <- SAMPLE_NAMES[idx+1]
	if (idx == 1) {
		all_merged <- merge(x = all_seurat[[name1]], y = all_seurat[[name2]] , add.cell.id = c(name1, name2),  project='CART_HL')
	} else {
		all_merged <- merge(x = all_merged, y = all_seurat[[name1]] , add.cell.id2 =  name1,  project='CART_HL')
		all_merged <- merge(x = all_merged, y = all_seurat[[name2]] , add.cell.id2 =  name2,  project='CART_HL')
	}

}


all_merged <- FindVariableFeatures(all_merged, nfeatures = 2000)
all_merged <- ScaleData(all_merged)
all_merged <- RunPCA(object = all_merged)
PCA_dims <- how_many_PCs(all_merged, paste0(SAMPLE, 'NumPCs'))
all_merged <- FindNeighbors(object = all_merged, dims = 1:PCA_dims)
# Determine the clusters for various resolutions
all_merged <- FindClusters(object = all_merged, resolution = c(0.2,0.4, 0.6, 0.8, 1.0, 1.2))
all_merged <- RunUMAP(all_merged, reduction = "pca", dims = 1:PCA_dims)
pdf('/home/sevastopol/data/gserranos/CART_HL/Plots/DENG_Tests/AllUmap_Stats.pdf', height = 8)
	clustree::clustree(all_merged@meta.data, prefix = "RNA_snn_res.")
	DimPlot(all_merged,  reduction = "umap", group.by = "sample") 
	features = c('nGene', 'nUMI', 'mitoRatio', 'RPSRatio')
	FeaturePlot(all_merged, features = features, order = TRUE) & NoLegend()
	features = c('CD8A', 'CD3E', 'CD4', 'CD3D', 'IL7R', 'CCR7', 'SELL', 'CD69', 'FMC63-CD19SCFV')
	FeaturePlot(all_merged, features = features, order = TRUE, min.cutoff = 'q10')& NoLegend()
dev.off()



# seurat integration
features <- SelectIntegrationFeatures(object.list = all_seurat)
all_seurat.anchors <- FindIntegrationAnchors(object.list = all_seurat, anchor.features = features)
all_seurat.combined <- IntegrateData(anchorset = all_seurat.anchors)

# harmony integration
library(harmony)
all_merged_int <- RunHarmony(all_merged, "sample")
harmony_embeddings <- Embeddings(all_merged_int, 'harmony')
all_merged_int <- RunUMAP(all_merged_int, reduction = "harmony",  dims = 1:PCA_dims)

# saveRDS(all_merged_int, paste0(PLOT_PATH, 'all_merged_int.rds'))
all_merged_int <- readRDS('/home/sevastopol/data/gserranos/CART_HL/Plots/DENG_Tests/ac24/all_merged_int.rds')

 pdf('/home/sevastopol/data/gserranos/CART_HL/Plots/DENG_Tests/Test.pdf', height = 8)
    DimPlot(object = all_merged_int, reduction = "umap", pt.size = .1, group.by = "sample")
    features = c('CD8A', 'CD3E', 'CD4', 'CD3D', 'IL7R', 'CCR7', 'SELL', 'CD69', 'FMC63-CD19SCFV')
    FeaturePlot(all_merged_int, features = features, order = TRUE, min.cutoff = 'q10')& NoLegend()
dev.off()

normData <- as.data.frame(t(all_merged_int@assays$RNA@data))
coords <- as.data.frame(all_merged_int@reductions$umap@cell.embeddings)
signatures_path <- '/home/sevastopol/data/gserranos/CART_HL/Data/signature/OtherSignatures'
tmp <- merge(coords, setNames(as.data.frame(all_merged_int$signature), 'Signature'), by=0)

tmp_bin <- tmp
tmp_bin$Signature <- ifelse(tmp_bin$Signature  > 0, 'High', 'Low')

coords_HL <- merge(coords, setNames(tmp_bin[, c('Row.names', 'Signature')], c('Row.names', 'BinScore')), by.x=0, by.y='Row.names')
rownames(coords_HL) <- coords_HL$Row.names

# pdf('/home/sevastopol/data/gserranos/CART_HL/Plots/DENG_Tests/DENG_signatures.pdf', height = 8)
# cowplot::plot_grid(
# get_expression_signature("Genes_Activation.txt", normData, coords , 1.2, split_HL=F) + theme(strip.text = element_text(size=8)),
# get_expression_signature("Genes_Tonic.txt", normData, coords , 1.2, split_HL=F)+ theme(strip.text = element_text(size=8)),
# get_umap_signature(tmp , 'HL signature'),
# get_umap_signature(tmp_bin , 'HL bin') + theme(legend.position = "bottom"),
# ncol=2)
# dev.off()



pdf('/home/sevastopol/data/gserranos/CART_HL/Plots/DENG_Tests/DENG_signatures_expanded.pdf', height = 12, width=10)
cowplot::plot_grid(
	cowplot::plot_grid(
		get_umap_signature(tmp , 'HL signature') + theme(axis.title.x = element_blank()),
		get_umap_signature(tmp_bin , 'HL bin') + theme(legend.position = "top", 
												 legend.margin=margin(0,0,0,0),
												 legend.box.margin=margin(-10,-10,-10,-10), 
												 axis.title = element_blank()),
	ncol=2),
	cowplot::plot_grid(
		get_expression_signature("Genes_Activation.txt", normData, coords , 1.2, split_HL=F) + 
								theme(strip.text = element_text(size=8),
									  axis.title.x = element_blank()),
		get_expression_signature("Genes_Tonic.txt", normData, coords , 1.2, split_HL=F) + 
								theme(strip.text = element_text(size=8), 
									  axis.title = element_blank()),
		ncol=2
	),
	cowplot::plot_grid(
		get_expression_signature("Genes_Activation.txt", normData, coords = coords_HL , 1.2, split_HL=T) + 
								theme(strip.text = element_text(size=8)),
		get_expression_signature("Genes_Tonic.txt", normData, coords = coords_HL , 1.2, split_HL=T) + 
								theme(strip.text = element_text(size=8),
									  axis.title.y = element_blank()),
		ncol=2
	),
	cowplot::plot_grid(
		get_umap(normData, 'CD4'),
		get_umap(normData, 'PDCD1'),
		get_umap(normData, 'HLA-DRA'),
		get_umap(normData, 'GZMA'),
		get_umap(normData, 'LAG3', 5),
		get_umap(normData, 'CD8A'),
		get_umap(normData, 'CCR7'),
		get_umap(normData, 'GATA3'),
		get_umap(normData, 'PRF1'),
		get_umap(normData, 'TIGIT', 5),
	nrow=2),
nrow=4)
dev.off()



pdf('/home/sevastopol/data/gserranos/CART_HL/Plots/DENG_Tests/FigureS9.pdf', height = 12, width=10)
legend <- cowplot::get_legend(get_expression_signature("Genes_Activation.txt", normData, coords = coords_HL , 1.2, split_HL=T) + 
															theme(legend.position = 'right'))
cowplot::plot_grid(
	cowplot::plot_grid(
		NULL,
		get_umap_HL(tmp_bin , 'HL bin') + theme(legend.position = "top", 
												 legend.margin=margin(0,0,0,0),
												 legend.box.margin=margin(-10,-10,-10,-10), 
												 axis.title = element_blank()),
		legend,
	ncol=3, rel_widths = c(0.2, 1, 0.2)),
	cowplot::plot_grid(
		get_expression_signature("Genes_Activation.txt", normData, coords = coords_HL , 1.2, split_HL=T) + 
								theme(strip.text = element_text(size=8)),
		get_expression_signature("Genes_Tonic.txt", normData, coords = coords_HL , 1.2, split_HL=T) + 
								theme(strip.text = element_text(size=8),
									  axis.title.y = element_blank()),
		ncol=2
	),
	cowplot::plot_grid(
		get_umap(normData, 'CD4'),
		get_umap(normData, 'PDCD1'),
		get_umap(normData, 'HLA-DRA'),
		get_umap(normData, 'GZMA'),
		get_umap(normData, 'LAG3', 5),
		get_umap(normData, 'CD8A'),
		get_umap(normData, 'CCR7'),
		get_umap(normData, 'GATA3'),
		get_umap(normData, 'PRF1'),
		get_umap(normData, 'TIGIT', 5),
	nrow=2),
nrow=3, rel_heights = c(1, 0.8,0.8))
dev.off()





# pdf('/home/sevastopol/data/gserranos/CART_HL/Plots/DENG_Tests/aaa.pdf')
# df <- as.data.frame(diamonds[order(diamonds$price, decreasing=TRUE), ])
# ggplot(data = df,aes(x=factor(cut),y=carat,colour=price)) + 
# geom_point(position=position_jitter(width=.4)) +
# scale_colour_gradientn(colours=c("grey20","orange","orange3"))
# dev.off()
