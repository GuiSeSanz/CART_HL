library(Seurat)
library(R.utils)
library(ggplot2)
library(gridExtra)






data_path = '/home/sevastopol/data/gserranos/CART_HL/Data/antiCD19/GSE151511_RAW/'
folders <- list.files(data_path, pattern = '^ac[0-9]{2}')


create_softlink <- function(folder){
    files_in_folder <- list.files(paste0(data_path, folder))
    if(!'barcodes.tsv' %in% files_in_folder){
        origin <- paste0(data_path, folder, '/', folder, '_barcodes.tsv')
        destination <- paste0(data_path, folder, '/', 'barcodes.tsv')
        command <- paste('ln -s', origin, destination, sep=' ')
        system(command)
    }
    if(!'matrix.mtx' %in% files_in_folder){
        origin <- paste0(data_path, folder, '/', folder, '_matrix.mtx')
        destination <- paste0(data_path, folder, '/', 'matrix.mtx')
        command <- paste('ln -s', origin, destination, sep=' ')
        system(command)
    }
    if(!'genes.tsv' %in% files_in_folder){
        origin <- paste0(data_path, folder, '/', folder, '_genes.tsv')
        destination <- paste0(data_path, folder, '/', 'genes.tsv')
        command <- paste('ln -s', origin, destination, sep=' ')
        system(command)
    }
    
}


cell_cycle_genes <- read.csv("./Data/cell_cycle_genes_hsa.csv") 
ah <- AnnotationHub::AnnotationHub() # Connect to AnnotationHub
ahDb <- AnnotationHub::query(ah, pattern = c("Homo sapiens", "EnsDb"), ignore.case = TRUE) # Access the Ensembl database for organism
id <- names(ahDb@.db_uid)[length(ahDb@.db_uid)]
edb <- ah[[id]]
annotations <- genes(edb, return.type = "data.frame")
annotations <- annotations[, c('gene_id', 'gene_name', 'seq_name', 'gene_biotype', 'description')]
cell_cycle_markers <- merge(cell_cycle_genes, annotations, by.x = "geneID" , by.y = "gene_id")
s_genes <- cell_cycle_markers[cell_cycle_markers$phase == 'S', 'gene_name']
g2m_genes <- cell_cycle_markers[cell_cycle_markers$phase == 'G2/M', 'gene_name']
dbDisconnect()


STATS <- data.frame(PATIENT=NULL,PercCAR=NULL, CD4=NULL, CD8=NULL, TOTAL=NULL )
for (folder in folders){
    print(folder)
    tmp_plot_path = paste(getwd(), 'Plots', 'antiCD19', folder, sep ='/')
    dir.create(tmp_plot_path, showWarnings=FALSE)

    create_softlink(folder)

    seurat_data <- Read10X(data.dir = paste0(data_path ,folder))
    seurat_obj <- CreateSeuratObject(counts = seurat_data,
                                   min.features = 100,
                                   project = folder)
    # metadata <- seurat_obj@meta.data
    # metadata$don_cells <- rownames(metadata) 
    # metadata <- setNames(metadata, c('seq_folder', 'nUMI', 'nGene'))
    # seurat_obj@meta.data <- metadata
    seurat_obj@meta.data$log10GenesPerUMI <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)
    seurat_obj$mitoRatio <- PercentageFeatureSet(object = seurat_obj, pattern = "^MT-")
    seurat_obj$RPSRatio <- PercentageFeatureSet(object =   seurat_obj, pattern = "^RP[SL]")
    seurat_obj$mitoRatio <- seurat_obj$mitoRatio / 100 
    seurat_obj$RPSRatio <- seurat_obj$RPSRatio / 100 

    # pdf(paste(tmp_plot_path, 'QC_1.pdf', sep = '/'))
    #     FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "mitoRatio")
    #     FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    #     FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "RPSRatio")
    # dev.off()

    count.min <- quantile(seurat_obj@meta.data$nCount_RNA,c(0.025))
    count.max <- quantile(seurat_obj@meta.data$nCount_RNA,c(0.975))
    feat.min <- quantile(seurat_obj@meta.data$nFeature_RNA,c(0.025))
    feat.max <- quantile(seurat_obj@meta.data$nFeature_RNA,c(0.975))

    # print("nUMIs = nCount_RNA")
    # print(summary(seurat_obj@meta.data$nCount_RNA))
    # print(c(count.min, count.max))
    # print("nGene = nFeature_RNA")
    # print(summary(seurat_obj@meta.data$nFeature_RNA))
    # print(c(feat.min, feat.max))
    
    DF_nCounts_nUmi <- data.frame(matrix(ncol = 6, nrow = 0))
    DF_nFeatureRNA_nGene <- data.frame(matrix(ncol = 6, nrow = 0))

    DF_nCounts_nUmi_maxmin <- data.frame(matrix(ncol = 2, nrow = 0))
    DF_nFeatureRNA_nGene_maxmin <- data.frame(matrix(ncol = 2, nrow = 0))

    DF_POST_nCounts_nUmi <- data.frame(matrix(ncol = 6, nrow = 0))
    DF_POST_nFeatureRNA_nGene <- data.frame(matrix(ncol = 6, nrow = 0))

    summaryofnCounts_nUmi <- summary(seurat_obj@meta.data$nCount_RNA)
    DF_nCounts_nUmi <- rbind(DF_nCounts_nUmi, summaryofnCounts_nUmi) 
    
    summaryofnFeatures_nGene <- summary(seurat_obj@meta.data$nFeature_RNA)
    DF_nFeatureRNA_nGene <- rbind(DF_nFeatureRNA_nGene, summaryofnFeatures_nGene) 
    
    mixmax <- c(count.min, count.max)
    DF_nCounts_nUmi_maxmin <- rbind(DF_nCounts_nUmi_maxmin, mixmax) 
    
    mixmax <- c(feat.min, feat.max)
    DF_nFeatureRNA_nGene_maxmin <- rbind(DF_nFeatureRNA_nGene_maxmin, mixmax)


    # pdf(paste(tmp_plot_path , "nGene_nUMI_mitoratio_RPSratio_PRENORM.pdf", sep='/'), width=11, height=8.5, onefile = T) 
    #     print(VlnPlot(seurat_obj, features = c("nGene", "nUMI", "mitoRatio", "RPSRatio"), ncol = 3,pt.size = -1))

    #     plot1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "mitoRatio")
    #     plot2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    #     plot3 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "RPSRatio")

    #     print(CombinePlots(plots = list(plot1, plot2, plot3)))
    #     print(plot(seurat_obj@meta.data$nCount_RNA,seurat_obj@meta.data$nFeature_RNA,pch=16,cex=0.7,bty="n"))
    #     print(abline(h=c(count.min,count.max),v=c(feat.min,feat.max),lty=2,lwd=1,col="red"))
    #     layout(matrix(1:3,ncol=3))
    #     print(plot(density(seurat_obj@meta.data$nCount_RNA),bty="n"))
    #     print(plot(density(seurat_obj@meta.data$nFeature_RNA),bty="n"))
    #     print(plot(density(seurat_obj@meta.data$mitoRatio),bty="n"))
    #     print(plot(density(seurat_obj@meta.data$RPSRatio),bty="n"))
    # dev.off()
    

  filtered_seurat <- subset(x = seurat_obj, 
                            subset = 
                              (nCount_RNA >= count.min) &
                              (nCount_RNA <= count.max) &
                              (nFeature_RNA >= feat.min) &
                              (nFeature_RNA <= feat.max) &
                              ## nFeature_RNA > 611 & 
                              ## nFeature_RNA < 6530 & 
                              ## nCount_RNA > 867 & 
                              ## nCount_RNA < 32441 &
                              # (nUMI >= 500) & 
                              # (nGene >= 300) & 
                              (log10GenesPerUMI > 0.80) & 
                              (mitoRatio < 0.1) & 
                              (RPSRatio < 0.4))
  counts <- GetAssayData(object = filtered_seurat, slot = "counts")
  nonzero <- counts > 0L 
  keep_genes <- rowSums(as.matrix(nonzero)) >= 10 
  filtered_counts <- counts[keep_genes, ] 
  clean_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data) 
  
    # pdf(paste(tmp_plot_path , "nGene_nUMI_mitoratio_RPSratio_PostNORM.pdf", sep='/'), width=11, height=8.5, onefile = T) 
    #     print(VlnPlot(clean_seurat, features = c("nGene", "nUMI", "mitoRatio", "RPSRatio"), ncol = 3,pt.size = -1))

    #     plot1 <- FeatureScatter(clean_seurat, feature1 = "nCount_RNA", feature2 = "mitoRatio")
    #     plot2 <- FeatureScatter(clean_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    #     plot3 <- FeatureScatter(clean_seurat, feature1 = "nCount_RNA", feature2 = "RPSRatio")

    #     print(CombinePlots(plots = list(plot1, plot2, plot3)))
    #     print(plot(clean_seurat@meta.data$nCount_RNA,clean_seurat@meta.data$nFeature_RNA,pch=16,cex=0.7,bty="n"))
    #     print(abline(h=c(count.min,count.max),v=c(feat.min,feat.max),lty=2,lwd=1,col="red"))
    #     layout(matrix(1:3,ncol=3))
    #     print(plot(density(clean_seurat@meta.data$nCount_RNA),bty="n"))
    #     print(plot(density(clean_seurat@meta.data$nFeature_RNA),bty="n"))
    #     print(plot(density(clean_seurat@meta.data$mitoRatio),bty="n"))
    #     print(plot(density(clean_seurat@meta.data$RPSRatio),bty="n"))
    # dev.off()
    


    clean_seurat <- NormalizeData(clean_seurat, verbose = TRUE)
    clean_seurat <- CellCycleScoring(clean_seurat, g2m.features=g2m_genes, s.features=s_genes)
    clean_seurat <- SCTransform(clean_seurat, vars.to.regress = c("mitoRatio", "nFeature_RNA"))
    # clean_seurat <- RunPCA(object = clean_seurat)
    #filter by cd4, cd8, cd8A and cargene


    # clean_seurat <- NormalizeData(object = clean_seurat)
    # clean_seurat <- FindVariableFeatures(object = clean_seurat)
    # clean_seurat <- ScaleData(object = clean_seurat)
    # clean_seurat <- RunPCA(object = clean_seurat)
    NormalizedData <- as.data.frame(GetAssayData(clean_seurat))
    CAR_NAME <- 'FMC63-CD19SCFV'
    NormalizedData_only_car <- NormalizedData[, NormalizedData[c(CAR_NAME), ]>0]
    NormalizedData_only_car_cd4 <- NormalizedData_only_car[, NormalizedData_only_car[c('CD4'), ]>0]
    NormalizedData_only_car_cd8 <- NormalizedData_only_car[, NormalizedData_only_car[c('CD8A'), ]>0 & NormalizedData_only_car[c('CD8B'), ]>0]
    # remove samples present in both groups
    samples_2_remove <- intersect(colnames(NormalizedData_only_car_cd4),colnames(NormalizedData_only_car_cd8))
    NormalizedData_only_car_cd4 <- NormalizedData_only_car_cd4[, !colnames(NormalizedData_only_car_cd4) %in% samples_2_remove]
    NormalizedData_only_car_cd8 <- NormalizedData_only_car_cd8[, !colnames(NormalizedData_only_car_cd8) %in% samples_2_remove]
    stats <- prop.table(table( NormalizedData[c(CAR_NAME), ]>0))*100
    STATS <- rbind(STATS, data.frame(PATIENT=folder,PercCAR=stats['TRUE'], CD4=ncol(NormalizedData_only_car_cd4), CD8=ncol(NormalizedData_only_car_cd8), TOTAL=ncol(NormalizedData_only_car)))

    path2save <- paste0(data_path, folder, '/', folder,'_normalized.rds')
    message(paste0('Writing data to ', path2save))
    saveRDS(NormalizedData, path2save)

    path2save <- paste0(data_path, folder, '/', folder,'_normalized_only_car.rds')
    message(paste0('Writing data to ', path2save))
    saveRDS(NormalizedData_only_car, path2save)

    path2save <- paste0(data_path, folder, '/', folder,'_normalized_only_car_cd4.rds')
    message(paste0('Writing data to ', path2save))
    saveRDS(NormalizedData_only_car_cd4, path2save)

    path2save <- paste0(data_path, folder, '/', folder,'_normalized_car_cd8.rds')
    message(paste0('Writing data to ', path2save))
    saveRDS(NormalizedData_only_car_cd8, path2save)

    # clean_seurat <- FindNeighbors(object = clean_seurat,  dims = 1:17)
    # clean_seurat <- FindClusters(object = clean_seurat, resolution = 0.6)
    # clean_seurat <- RunUMAP(clean_seurat, reduction = "pca", dims = 1:17)
    # pdf(paste(tmp_plot_path , "UMAP_clusters.pdf", sep='/'))
    #   DimPlot(clean_seurat, reduction = "umap")
    # dev.off()

    # Markers <- FindAllMarkers(clean_seurat, only.pos=TRUE)
    # Markers_sep <-  split(Markers, f=Markers$cluster)
    # path2save <- paste0(data_path, folder, '/', folder,'_Markers_per_cluster.rds')
    # saveRDS(Markers_sep, path2save)
    

    ## and RunPCA(object = clean_seurat) are already in explore1
}


print(STATS)
rownames(STATS) <- STATS$PATIENT