#!/usr/bin/env Rscript
#===============================================================================
#' Author: Maria E. Calleja
#' Date: 2019/09/30
#' https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html
#' As well as JP tutorial. 
#' This script is over JuanRo's samples, following Satija Lab's vignettes.
#' https://github.com/satijalab/seurat/wiki/Seurat
#===============================================================================
## PACKAGES.
library(Seurat)
library(Matrix)
library(dplyr)
library(cowplot)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(RCurl)
library(AnnotationHub)
library(ensembldb)
library(gridExtra)
library(grid)
library(lattice)
library(rlang)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(scCATCH)
library(SingleR)
library(viridis)
library(ggplot2)
library(RColorBrewer)
library(hues)
options(stringsAsFactors = FALSE)
#===============================================================================

#===============================================================================
## GLOBAL VARIABLES.
#Local
#PARENT_DIR<-"Cluster/Touch/MECC"
PARENT_DIR<-"/home/mecc/cluster_juanro"
#PARENT_DIR<-"/home/jrrodriguez"
#PARENT_DIR<-"/Users/bioinfo104/cluster_juanro"
PROJECT_DIR<-file.path(PARENT_DIR,"SC_HighLow/data")


BS_D10_RAW_DATA_DIR<-file.path(PARENT_DIR,"SC_HighLow/data/edited/GEX_D10_BCMA_BBzB_CAR_d0/outs/filtered_feature_bc_matrix/")
ST_D10_RAW_DATA_DIR<-file.path(PARENT_DIR,"SC_HighLow/data/edited/GEX_D10_BCMA_BBzB_CAR_d7/outs/filtered_feature_bc_matrix/")
BS_D14_RAW_DATA_DIR<-file.path(PARENT_DIR,"SC_HighLow/data/edited/sc5_D14_BCMA_d0/outs/filtered_feature_bc_matrix/")
ST_D14_RAW_DATA_DIR<-file.path(PARENT_DIR,"SC_HighLow/data/edited/sc5_D14_BCMA_d7/outs/filtered_feature_bc_matrix/")
BS_D18_RAW_DATA_DIR<-file.path(PARENT_DIR,"SC_HighLow/data/edited/sc5_D18_BCMA_d0/outs/filtered_feature_bc_matrix/")
ST_D18_RAW_DATA_DIR<-file.path(PARENT_DIR,"SC_HighLow/data/edited/sc5_D18_BCMA_d7/outs/filtered_feature_bc_matrix/")


FIGURES_DIR<-file.path(PARENT_DIR,"SC_HighLow/temp_feb")
RESULTS_DIR<-file.path(PARENT_DIR,"SC_HighLow/temp_feb")
#source("SC_HighLow/src/Explore1.R")
#load((file =file.path(PARENT_DIR,"SC_HighLow/RSession/7ene21.RData")))
#===============================================================================
## FUNCTIONS.
##
## FUNCTIONS.
ScQCPlots <- function(metadata, nUMI, nGene, log10GenesPerUMI, mitoRatio, sample, Plotname, directoriodeseado){
  nUMI <- enquo(nUMI)
  nGene <- enquo(nGene)
  log10GenesPerUMI <- enquo(log10GenesPerUMI)
  mitoRatio <- enquo(mitoRatio)
  sample <- enquo(sample)
  FIGURES_DIR <- directoriodeseado
  
  # Visualize the number of cell counts per cell
  gg.cel.counts<- metadata %>% 
    ggplot(aes(x=sample, fill=sample)) + 
    geom_bar() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5, face="bold")) +
    ggtitle("NCells")
  # Visualize the number UMIs/transcripts per cellmetadata %>% 
  gg.umis.per.cell<- metadata %>% 
    ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    ylab("Cell density") +
    geom_vline(xintercept = 500)
  # Genes detected per cell
  gg.genes.per.cell<- metadata %>% 
    ggplot(aes(color=sample, x=nGene, fill= sample)) + 
    geom_density(alpha = 0.2) + 
    theme_classic() +
    scale_x_log10() + 
    geom_vline(xintercept = 300)
  # Visualize the distribution of genes detected per cell via boxplot
  gg.genes.per.cell.bx<- metadata %>% 
    ggplot(aes(x=sample, y=log10(nGene), fill=sample)) + 
    geom_boxplot() + 
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5, face="bold")) +
    ggtitle("NCells vs NGenes")
  # UMIs vs. genes detected
  gg.umis.vs.genes<- metadata %>% 
    ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
    geom_point() + 
    stat_smooth(method=lm) +
    scale_x_log10() + 
    scale_y_log10() + 
    theme_classic() +
    geom_vline(xintercept = 500) +
    geom_hline(yintercept = 250) +
    facet_wrap(~sample)
  # Mitochondrial counts ratio
  gg.mt.ct.ratio<- metadata %>% 
    ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    geom_vline(xintercept = 0.2)
  #Novelty
  gg.novelty<- metadata %>%
    ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
    geom_density(alpha = 0.2) +
    theme_classic() +
    geom_vline(xintercept = 0.8)
  
  ggplot_list <- Filter(function(x) is(x, "ggplot"), mget(ls()))
  ggsave(file=file.path(FIGURES_DIR,paste0(Plotname,"multipage_plots_in_function.pdf")), 
         do.call(arrangeGrob, c(ggplot_list, ncol=2)), width=11, height=8.5)
  #dev.off()
  message("Sc QC Plot ha terminado")
  return(ggplot_list)
}
#===============================================================================
## MAIN.
################################################################################
## All steps together 
################################################################################
# Download cell cycle genes for organism at https://github.com/hbc/tinyatlas/tree/master/cell_cycle. Read it in with:
cc_file <- getURL("https://raw.githubusercontent.com/hbc/tnyatlas/master/cell_cycle/Homo_sapiens.csv") 
cc_file <- file.path(PARENT_DIR,"SC_HighLow/temp/cell_cycle_genes_hsa.csv") 
cell_cycle_genes <- read.csv(cc_file)
ah <- AnnotationHub() # Connect to AnnotationHub
ahDb <- query(ah, pattern = c("Homo sapiens", "EnsDb"), ignore.case = TRUE) # Access the Ensembl database for organism
id <- ahDb %>% mcols() %>% rownames() %>%tail(n = 1) # Acquire the latest annotation files
edb <- ah[[id]] # Download the appropriate Ensembldb database
annotations <- genes(edb, return.type = "data.frame") # Extract gene-level information from database
# Select annotations of interest
annotations <- annotations %>% dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)
# Get gene names for Ensembl IDs for each gene
cell_cycle_markers <- dplyr::left_join(cell_cycle_genes, annotations, by = c("geneID" = "gene_id"))
s_genes <- cell_cycle_markers %>% dplyr::filter(phase == "S") %>% pull("gene_name") # Acquire the S phase genes
g2m_genes <- cell_cycle_markers %>% dplyr::filter(phase == "G2/M") %>% pull("gene_name") # Acquire the G2M phase genes   
dbDisconnect()

# Create each individual Seurat object for every sample
GEX_Directories <- list(BS_D10_RAW_DATA_DIR=BS_D10_RAW_DATA_DIR, 
                        ST_D10_RAW_DATA_DIR=ST_D10_RAW_DATA_DIR,
                        BS_D14_RAW_DATA_DIR=BS_D14_RAW_DATA_DIR,
                        ST_D14_RAW_DATA_DIR=ST_D14_RAW_DATA_DIR,
                        BS_D18_RAW_DATA_DIR=BS_D18_RAW_DATA_DIR,
                        ST_D18_RAW_DATA_DIR=ST_D18_RAW_DATA_DIR)
for (ind in seq_along(GEX_Directories)){
  nm<- names(GEX_Directories)[ind]
  print(paste(strsplit(nm,"_")[[1]][1:2], collapse="_"))
  seurat_data <- Read10X(data.dir = file.path(GEX_Directories[ind]))
  seurat_obj <- CreateSeuratObject(counts = seurat_data,
                                   min.features = 100,
                                   project = paste(strsplit(nm,"_")[[1]][1:2], collapse="_"))
  assign(paste(strsplit(nm,"_")[[1]][1:2], collapse="_"), seurat_obj)}

rm(seurat_data,seurat_obj, GEX_Directories, nm, ind) 

merged_seurat_d10 <- merge(x = BS_D10, 
                       y = ST_D10, 
                       add.cell.id = c("ctrl", "stim")) 

merged_seurat_d14 <- merge(x = BS_D14, 
                       y = ST_D14, 
                       add.cell.id = c("ctrl", "stim")) 

merged_seurat_d18 <- merge(x = BS_D18, 
                           y = ST_D18, 
                           add.cell.id = c("ctrl", "stim")) 

# merged_seurats <- list(merged_seurat_d10=merged_seurat_d10,
#                        merged_seurat_d14=merged_seurat_d14,
#                        merged_seurat_d18=merged_seurat_d18)

rm(BS_D10, ST_D10, BS_D14, ST_D14, BS_D18, ST_D18)

clean_seurat <- merge(x = merged_seurat_d10, y = list(merged_seurat_d14, merged_seurat_d18), add.cell.ids = c("d10","d14","d18"))

metadata <- clean_seurat@meta.data 
metadata$don_cells <- rownames(metadata) 
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident, nUMI = nCount_RNA, nGene = nFeature_RNA) 
metadata$donor <- NA 
metadata$donor[which(str_detect(metadata$don_cells, "^d10_"))] <- "d10"
metadata$donor[which(str_detect(metadata$don_cells, "^d14_"))] <- "d14"
metadata$donor[which(str_detect(metadata$don_cells, "^d18_"))] <- "d18"

metadata$donor_by <- NA 
metadata$donor_by[which(str_detect(metadata$don_cells, "^d10_ctrl_"))] <- "d10_ctrl"
metadata$donor_by[which(str_detect(metadata$don_cells, "^d14_ctrl_"))] <- "d14_ctrl"
metadata$donor_by[which(str_detect(metadata$don_cells, "^d18_ctrl_"))] <- "d18_ctrl"
metadata$donor_by[which(str_detect(metadata$don_cells, "^d10_stim_"))] <- "d10_stim"
metadata$donor_by[which(str_detect(metadata$don_cells, "^d14_stim_"))] <- "d14_stim"
metadata$donor_by[which(str_detect(metadata$don_cells, "^d18_stim_"))] <- "d18_stim"

clean_seurat@meta.data <- metadata 

seurat_raw <- clean_seurat
length(rownames(seurat_raw@meta.data))

rm(merged_seurat_d10, merged_seurat_d14, merged_seurat_d18, metadata, ah, ahDb)
gc(gc())
# 
seurat_list <- SplitObject(seurat_raw, split.by = "donor_by")
rm(clean_seurat)
#seurat_list <- seurat_list[c("d10_ctrl", "d10_stim", "d14_ctrl", "d14_stim", "d18_ctrl", "d18_stim")]

#options(future.globals.maxSize = 4000 * 1024^2)
rm(seurat_raw)
for (i in 1:length(seurat_list)) {
  nm <- names(seurat_list[i])
  print(nm)
  print(dim(seurat_list[[i]]@assays$RNA))
  print(sum(as.data.frame(seurat_list[[i]]@assays$RNA@counts) == 0)/ prod(dim(as.data.frame(seurat_list[[i]]@assays$RNA@counts))))
  #https://stackoverflow.com/questions/43663314/how-could-i-calculate-the-sparsity-of-a-data-frame-in-r
}


CTR_STIM_Seurats <- vector("list", 6)
DF_nCounts_nUmi <- data.frame(matrix(ncol = 6, nrow = 0))
DF_nFeatureRNA_nGene <- data.frame(matrix(ncol = 6, nrow = 0))

DF_nCounts_nUmi_maxmin <- data.frame(matrix(ncol = 2, nrow = 0))
DF_nFeatureRNA_nGene_maxmin <- data.frame(matrix(ncol = 2, nrow = 0))

DF_POST_nCounts_nUmi <- data.frame(matrix(ncol = 6, nrow = 0))
DF_POST_nFeatureRNA_nGene <- data.frame(matrix(ncol = 6, nrow = 0))

for (i in 1:length(seurat_list)) {
  nm <- names(seurat_list[i])
  seurat_list[[i]]@meta.data$log10GenesPerUMI <- log10(seurat_list[[i]]@meta.data$nGene) / log10(seurat_list[[i]]@meta.data$nUMI)
  seurat_list[[i]]$mitoRatio <- PercentageFeatureSet(object = seurat_list[[i]], pattern = "^MT-")
  seurat_list[[i]]$RPSRatio <- PercentageFeatureSet(object = seurat_list[[i]], pattern = "^RP[SL]")
  seurat_list[[i]]$mitoRatio <- seurat_list[[i]]$mitoRatio / 100 
  seurat_list[[i]]$RPSRatio <- seurat_list[[i]]$RPSRatio / 100 
  
  # rb.genes <- rownames(alldata)[grep("^RP[SL]",rownames(alldata))]
  # percent.ribo <- colSums(C[rb.genes,])/Matrix::colSums(C)*100
  # alldata <- AddMetaData(alldata, percent.ribo, col.name = "percent.ribo")
  
  metadata <- seurat_list[[i]]@meta.data 
  metadata$cells <- rownames(metadata) 
  metadata$sample <- NA
  metadata$sample[which(str_detect(metadata$cells, "_ctrl_"))] <- "ctrl"
  metadata$sample[which(str_detect(metadata$cells, "_stim_"))] <- "stim"
  seurat_list[[i]]@meta.data <- metadata 
  
  #Viene de "Seurat_Anañysis Organoides"
  plot1 <- FeatureScatter(seurat_list[[i]], feature1 = "nUMI", feature2 = "mitoRatio")
  plot2 <- FeatureScatter(seurat_list[[i]], feature1 = "nUMI", feature2 = "nGene")
  plot3 <- FeatureScatter(seurat_list[[i]], feature1 = "nUMI", feature2 = "RPSRatio")
  
  count.min <- quantile(seurat_list[[i]]@meta.data$nUMI,c(0.025))
  count.max <- quantile(seurat_list[[i]]@meta.data$nUMI,c(0.975))
  feat.min <- quantile(seurat_list[[i]]@meta.data$nGene,c(0.025))
  feat.max <- quantile(seurat_list[[i]]@meta.data$nGene,c(0.975))
  
  print("nUMIs = nCount_RNA")
  print(summary(seurat_list[[i]]@meta.data$nUMI))
  print(c(count.min, count.max))
  print("nGene = nFeature_RNA")
  print(summary(seurat_list[[i]]@meta.data$nGene))
  print(c(feat.min, feat.max))
  
  summaryofnCounts_nUmi <- summary(seurat_list[[i]]@meta.data$nUMI)
  DF_nCounts_nUmi <- rbind(DF_nCounts_nUmi, summaryofnCounts_nUmi) 
  
  summaryofnFeatures_nGene <- summary(seurat_list[[i]]@meta.data$nGene)
  DF_nFeatureRNA_nGene <- rbind(DF_nFeatureRNA_nGene, summaryofnFeatures_nGene) 
  
  mixmax <- c(count.min, count.max)
  DF_nCounts_nUmi_maxmin <- rbind(DF_nCounts_nUmi_maxmin, mixmax) 
  
  mixmax <- c(feat.min, feat.max)
  DF_nFeatureRNA_nGene_maxmin <- rbind(DF_nFeatureRNA_nGene_maxmin, mixmax)
  
  #plot4 <- plot(seurat_list[[i]]@meta.data$nUMI,seurat_list[[i]]@meta.data$nGene,pch=16,cex=0.7,bty="n")
  #plot4 <- plot4 + abline(h=c(count.min,count.max),v=c(feat.min,feat.max),lty=2,lwd=1,col="red") #viene de los cuantiles... 
  
  pdf(file.path(FIGURES_DIR,paste0("nGene_nUMI_mitoratio_RPSratio",nm,"_22032021_PRENORM.pdf")), width=11, height=8.5, onefile = T) 
  print(VlnPlot(seurat_list[[i]], features = c("nGene", "nUMI", "mitoRatio", "RPSRatio"), ncol = 3,pt.size = -1))
  print(CombinePlots(plots = list(plot1, plot2, plot3)))
  print(plot(seurat_list[[i]]@meta.data$nUMI,seurat_list[[i]]@meta.data$nGene,pch=16,cex=0.7,bty="n"))
  print(abline(h=c(count.min,count.max),v=c(feat.min,feat.max),lty=2,lwd=1,col="red"))
  layout(matrix(1:3,ncol=3))
  print(plot(density(seurat_list[[i]]@meta.data$nUMI),bty="n"))
  print(plot(density(seurat_list[[i]]@meta.data$nGene),bty="n"))
  print(plot(density(seurat_list[[i]]@meta.data$mitoRatio),bty="n"))
  print(plot(density(seurat_list[[i]]@meta.data$RPSRatio),bty="n"))
  dev.off()
  
  givename <- paste0("PRE_FILTERING_",nm,"_22032021")
  ScQCPlots(metadata, nUMI, nGene, log10GenesPerUMI, mitoRatio, sample, givename, directoriodeseado = FIGURES_DIR)
  
  filtered_seurat <- subset(x = seurat_list[[i]], 
                            subset = 
                              (nUMI >= count.min) &
                              (nUMI <= count.max) &
                              (nGene >= feat.min) &
                              (nGene <= feat.max) &
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
  
  print(nm)
  print("Filtered_seurat")
  print(dim(filtered_seurat@assays$RNA))
  print("Filtered_counts")
  print(dim(filtered_counts))
  print("clean_seurat")
  print(dim(clean_seurat))
  
  summaryofnCounts_nUmi <- summary(clean_seurat@meta.data$nUMI)
  DF_POST_nCounts_nUmi <- rbind(DF_POST_nCounts_nUmi, summaryofnCounts_nUmi) 
  
  summaryofnFeatures_nGene <- summary(clean_seurat@meta.data$nGene)
  DF_POST_nFeatureRNA_nGene <- rbind(DF_POST_nFeatureRNA_nGene, summaryofnFeatures_nGene) 
  
  plot1 <- FeatureScatter(clean_seurat, feature1 = "nUMI", feature2 = "mitoRatio")
  plot2 <- FeatureScatter(clean_seurat, feature1 = "nUMI", feature2 = "nGene")
  plot3 <- FeatureScatter(clean_seurat, feature1 = "nUMI", feature2 = "RPSRatio")
  
  pdf(file.path(FIGURES_DIR,paste0("nGene_nUMI_mitoratio_RPSratio",nm,"_22032021_POSNORM.pdf")), width=11, height=8.5, onefile = T) 
  print(VlnPlot(clean_seurat, features = c("nGene", "nUMI", "mitoRatio", "RPSRatio"), ncol = 3,pt.size = -1))
  print(CombinePlots(plots = list(plot1, plot2, plot3)))
  layout(matrix(1:3,ncol=3))
  print(plot(density(clean_seurat@meta.data$nUMI),bty="n"))
  print(plot(density(clean_seurat@meta.data$nGene),bty="n"))
  print(plot(density(clean_seurat@meta.data$mitoRatio),bty="n"))
  print(plot(density(clean_seurat@meta.data$RPSRatio),bty="n"))
  dev.off()
  rm(counts, filtered_counts, filtered_seurat, nonzero, keep_genes, count.max, count.min, feat.max, feat.min)
  
  CTR_STIM_Seurats[i] <- clean_seurat
  names(CTR_STIM_Seurats)[i] <- nm #https://stackoverflow.com/questions/38643000/naming-list-elements-in-r
  
  # Re-assess QC metrics
  metadata_clean <- clean_seurat@meta.data #Filtered subset: I'm naming clean_merged_seurat referring to a "personalized for this exp"
  save(clean_seurat, file=file.path(RESULTS_DIR,paste0("clean_merged_seurat",nm,"_22032021_Fangorn.RData",collapse = "_")))
  givename <- paste0("POST-FILTERING",nm,"_22032021",collapse = "_")
  ScQCPlots(metadata_clean, nUMI, nGene, log10GenesPerUMI, mitoRatio, sample, givename, directoriodeseado = FIGURES_DIR)
  
  # set.seed(98)
  # seurat_raw <- clean_seurat # Create new object
  # control_cell_ids <- rownames(seurat_raw@meta.data[which(seurat_raw@meta.data$sample == "ctrl"), ])
  # stim_cell_ids <- rownames(seurat_raw@meta.data[which(seurat_raw@meta.data$sample == "stim"), ])
  # seurat_control <- subset(seurat_raw, cells = control_cell_ids) 
  # seurat_stim <- subset(seurat_raw, cells = stim_cell_ids)
  
  #nm<-strsplit(nm,"_")[[1]][3]
  print(nm)
  rm(metadata, metadata_clean, clean_seurat, givename, plot1, plot2, plot3)
  
  # CTR_STIM_Seurats <- list( a = seurat_control, b = seurat_stim)  
  # names(CTR_STIM_Seurats) <- c(paste0("seurat_control_",nm), paste0("seurat_stim_",nm))
  # assign(paste0("seurat_control_",as.character(nm)),CTR_STIM_Seurats[[paste0("seurat_control_",as.character(nm))]], envir = .GlobalEnv )
  # assign(paste0("seurat_stim_",as.character(nm)),CTR_STIM_Seurats[[paste0("seurat_stim_",as.character(nm))]], envir = .GlobalEnv )
  # #set.seed(98)
  # #lapply(seq_along(CTR_STIM_Seurats), explore1, seurat_list=CTR_STIM_Seurats, directoriodeseado=FIGURES_DIR)
  # 
  # rm(seurat_control, seurat_stim)
  
}
colnames(DF_nCounts_nUmi) <- c("Min.","1st Qu.","Median","Mean", "3rd Qu.", "Max.") 
rownames(DF_nCounts_nUmi) <- names(seurat_list)
colnames(DF_nFeatureRNA_nGene) <- c("Min.","1st Qu.","Median","Mean",  "3rd Qu.", "Max.") 
rownames(DF_nFeatureRNA_nGene) <- names(seurat_list)

colnames(DF_POST_nCounts_nUmi) <- c("Min.","1st Qu.","Median","Mean", "3rd Qu.", "Max.") 
rownames(DF_POST_nCounts_nUmi) <- names(seurat_list)
colnames(DF_POST_nFeatureRNA_nGene) <- c("Min.","1st Qu.","Median","Mean",  "3rd Qu.", "Max.") 
rownames(DF_POST_nFeatureRNA_nGene) <- names(seurat_list)

colnames(DF_nCounts_nUmi_maxmin) <- c("q.min", "q.max") 
rownames(DF_nCounts_nUmi_maxmin) <- names(seurat_list)
colnames(DF_nFeatureRNA_nGene_maxmin) <- c("q.min", "q.max") 
rownames(DF_nFeatureRNA_nGene_maxmin) <- names(seurat_list)

DF_nCounts_nUmi <- cbind(DF_nCounts_nUmi, DF_nCounts_nUmi_maxmin, DF_POST_nCounts_nUmi)
DF_nFeatureRNA_nGene <- cbind(DF_nFeatureRNA_nGene, DF_nFeatureRNA_nGene_maxmin, DF_POST_nFeatureRNA_nGene)

write.csv(DF_nCounts_nUmi, file.path(RESULTS_DIR,"DF_nCounts_nUmi.csv"), quote = FALSE, row.names = TRUE)
write.csv(DF_nFeatureRNA_nGene, file.path(RESULTS_DIR,"DF_nFeatures_nGene.csv"), quote = FALSE, row.names = TRUE)

for (i in 1:length(CTR_STIM_Seurats)) {
  nm <- names(CTR_STIM_Seurats[i])
  
  print("PRE")
  print(nm)
  print(dim(seurat_list[[i]]@assays$RNA))
  print(sum(as.data.frame(seurat_list[[i]]@assays$RNA@counts) == 0)/ prod(dim(as.data.frame(seurat_list[[i]]@assays$RNA@counts))))
  print("")
  
  print("POST")
  print(nm)
  print(dim(CTR_STIM_Seurats[[i]]@assays$RNA))
  print(sum(as.data.frame(CTR_STIM_Seurats[[i]]@assays$RNA@counts) == 0)/ prod(dim(as.data.frame(CTR_STIM_Seurats[[i]]@assays$RNA@counts))))
}

# do.call(cbind, lapply(seurat_list[[i]]@meta.data, summary))
# s <- with(seurat_list[[i]]@meta.data, tapply(nUMI, summary))
# head(seurat_list[[i]]@meta.data)
# DF_nCounts_nUmi <- do.call(rbind, summaryofnCounts_nUmi) 

# list.files(RESULTS_DIR)
# load(file=file.path(PARENT_DIR,"/SC_HighLow/temp_oct/clean_merged_seuratd10_ctrl_22032021_Fangorn.RData"))
# clean_seurat_d10 <- clean_seurat
# load(file=file.path(PARENT_DIR,"/SC_HighLow/temp_oct/clean_merged_seuratd10_stim_22032021_Fangorn.RData"))
# clean_seurat_d10 <- clean_seurat
# load(file=file.path(PARENT_DIR,"/SC_HighLow/temp_oct/clean_merged_seuratd14_ctrl_22032021_Fangorn.RData"))
# clean_seurat_d10 <- clean_seurat
# load(file=file.path(PARENT_DIR,"/SC_HighLow/temp_oct/clean_merged_seuratd14_stim_22032021_Fangorn.RData"))
# clean_seurat_d10 <- clean_seurat
# load(file=file.path(PARENT_DIR,"/SC_HighLow/temp_oct/clean_merged_seuratd18_ctrl_22032021_Fangorn.RData"))
# clean_seurat_d10 <- clean_seurat
# load(file=file.path(PARENT_DIR,"/SC_HighLow/temp_oct/clean_merged_seuratd18_stim_22032021_Fangorn.RData"))
# clean_seurat_d10 <- clean_seurat
# 
# rm(clean_seurat)
# 
# clean_seurat <- merge(x = clean_seurat_d10, y = list(clean_seurat_d14, clean_seurat_d18), add.cell.ids = c("d10","d14","d18"))
# 
# metadata <- clean_seurat@meta.data 
# metadata$don_cells <- rownames(metadata) 
# metadata <- metadata %>%
#   dplyr::rename(seq_folder = orig.ident, nUMI = nCount_RNA, nGene = nFeature_RNA) 
# metadata$donor <- NA 
# metadata$donor[which(str_detect(metadata$don_cells, "^d10_"))] <- "d10"
# metadata$donor[which(str_detect(metadata$don_cells, "^d14_"))] <- "d14"
# metadata$donor[which(str_detect(metadata$don_cells, "^d18_"))] <- "d18"
# 
# metadata$donor_by <- NA 
# metadata$donor_by[which(str_detect(metadata$don_cells, "^d10_ctrl_"))] <- "d10_ctrl"
# metadata$donor_by[which(str_detect(metadata$don_cells, "^d14_ctrl_"))] <- "d14_ctrl"
# metadata$donor_by[which(str_detect(metadata$don_cells, "^d18_ctrl_"))] <- "d18_ctrl"
# metadata$donor_by[which(str_detect(metadata$don_cells, "^d10_stim_"))] <- "d10_stim"
# metadata$donor_by[which(str_detect(metadata$don_cells, "^d14_stim_"))] <- "d14_stim"
# metadata$donor_by[which(str_detect(metadata$don_cells, "^d18_stim_"))] <- "d18_stim"
# 
# clean_seurat@meta.data <- metadata 
# 
# seurat_raw <- clean_seurat
# length(rownames(seurat_raw@meta.data))
# 
# rm(clean_seurat_d10, clean_seurat_d14, clean_seurat_d18, metadata, ah, ahDb)
# gc(gc())
# seurat_list <- SplitObject(seurat_raw, split.by = "donor_by")
# rm(clean_seurat)
# #seurat_list <- seurat_list[c("d10_ctrl", "d10_stim", "d14_ctrl", "d14_stim", "d18_ctrl", "d18_stim")]
# 
# #options(future.globals.maxSize = 4000 * 1024^2)
# rm(seurat_raw, metadata)

save.image (file =file.path("/home/mecc/cluster_juanro/SC_HighLow/temp_mar/23mar21.RData"))
rm(list= ls(pattern = "DF_*"))
rm(mixmax, summaryofnCounts_nUmi, summaryofnFeatures_nGene)
rm(seurat_list)
gc()
gc(gc())
saveRDS(CTR_STIM_Seurats, file.path(PARENT_DIR,"SC_HighLow/temp_mar/clean_seurats_24mar.rds"))

###Inicio_encluster_CAR_v4
source(file.path(PARENT_DIR,"SC_HighLow/temp_mar/Combined/explore1.R"))
source(file.path(PARENT_DIR,"SC_HighLow/temp_mar/Combined/explore2.R"))
FIGURES_DIR<-file.path(PARENT_DIR,"SC_HighLow/temp_mar/Combined/")
RESULTS_DIR<-file.path(PARENT_DIR,"SC_HighLow/temp_mar/Combined/")

clean_seurat_d10 <- merge(x = CTR_STIM_Seurats$d10_ctrl, y = CTR_STIM_Seurats$d10_stim)
head(CTR_STIM_Seurats$d10_ctrl@meta.data)
head(clean_seurat_d10@meta.data)
table(clean_seurat_d10@meta.data$sample)
clean_seurat_d14 <- merge(x = CTR_STIM_Seurats$d14_ctrl, y = CTR_STIM_Seurats$d14_stim)
clean_seurat_d18 <- merge(x = CTR_STIM_Seurats$d18_ctrl, y = CTR_STIM_Seurats$d18_stim)

clean_seurat <- merge(x = clean_seurat_d10, y = list(clean_seurat_d14, clean_seurat_d18))
length(rownames(clean_seurat@meta.data))

rm(clean_seurat_d10, clean_seurat_d14, clean_seurat_d18)
gc(gc())

saveRDS(clean_seurat, file.path(PARENT_DIR,"SC_HighLow/temp_mar/Combined/clean_seurat_merged_24mar.rds"))

CTR_STIM_Seurats <- list(clean_seurat= clean_seurat)
lapply(seq_along(CTR_STIM_Seurats), explore1, seurat_list=CTR_STIM_Seurats, directoriodeseado=FIGURES_DIR)
## NormalizeData(object = clean_seurat), FindVariableFeatures(object = clean_seurat), ScaleData(object = clean_seurat)
## and RunPCA(object = clean_seurat) are already in explore1
set.seed(123)
clean_seurat <- FindNeighbors(object = clean_seurat,  dims = 1:17)
set.seed(123)
clean_seurat <- FindClusters(object = clean_seurat, resolution = c(0.8, 1.0))
Resolution_Identities<- c( "RNA_snn_res.0.8","RNA_snn_res.1")
CTR_STIM_Seurats <- list(clean_seurat= clean_seurat)
CTR_STIM_dims <- c(clean_seurat= "1:18")
lapply(seq_along(CTR_STIM_Seurats), explore2, seurat_list=CTR_STIM_Seurats, seurat_dims_list=CTR_STIM_dims, Resolution_Identities=Resolution_Identities,
       directoriodeseado=FIGURES_DIR )
rm(Resolution_Identities)
Idents(object = clean_seurat) <- "RNA_snn_res.1"
set.seed(123)
clean_seurat <- RunTSNE(object = clean_seurat)
cl_plot <- DimPlot(object = clean_seurat, reduction = "tsne")
cl_plot <- cl_plot + ggtitle("t-SNE")
cl_plot_2 <- DimPlot(object = clean_seurat, reduction = "tsne", split.by = "donor_by", group.by = "donor_by")
cl_plot_2 <- cl_plot_2 + ggtitle("t-SNE grouped by donor_by")
cl_plot_3 <- DimPlot(object = clean_seurat, reduction = "tsne", split.by = "donor")
cl_plot_3 <- cl_plot_3 + ggtitle("t-SNE split by donor")
set.seed(123)
clean_seurat <- RunUMAP(clean_seurat, reduction = "pca", dims = 1:17)
cl_plot_4 <- DimPlot(clean_seurat,reduction = "umap")
                     # label = TRUE,label.size = 6,
                     #split.by = "donor_by")
cl_plot_4 <- cl_plot_4 + ggtitle("UMAP")

cl_plot_5 <- DimPlot(clean_seurat,reduction = "umap",
                     # label = TRUE,label.size = 6,
                     split.by = "donor_by", group.by = "donor_by")
cl_plot_5 <- cl_plot_5 + ggtitle("UMAP")
cl_plot_6 <- DimPlot(clean_seurat,reduction = "umap",
                     # label = TRUE,label.size = 6,
                     split.by = "donor")
cl_plot_6 <- cl_plot_6 + ggtitle("UMAP")

pdf(file.path(FIGURES_DIR,"CleanSeurat_Tsne_UMAP_RNA_snn_res.1.pdf"),
    width=11, height=8.5, onefile = T)
print(cl_plot)
print(cl_plot_2)
print(cl_plot_3)
print(cl_plot_4)
print(cl_plot_5)
print(cl_plot_6)
dev.off()

rm(list = ls(pattern = "cl_*"))
rm(CTR_STIM_dims)
CTR_STIM_Seurats <- readRDS(file.path(PARENT_DIR,"SC_HighLow/temp_mar/clean_seurats_24mar.rds"))

FIGURES_DIR<-file.path(PARENT_DIR,"SC_HighLow/temp_mar/CAR_EXP/")
RESULTS_DIR<-file.path(PARENT_DIR,"SC_HighLow/temp_mar/CAR_EXP/")

for (i in 1:length(CTR_STIM_Seurats)) {
  CTR_STIM_Seurats[[i]] <- NormalizeData(CTR_STIM_Seurats[[i]], verbose = TRUE)
  CTR_STIM_Seurats[[i]] <- CellCycleScoring(CTR_STIM_Seurats[[i]], g2m.features=g2m_genes, s.features=s_genes)
  CTR_STIM_Seurats[[i]] <- SCTransform(CTR_STIM_Seurats[[i]], vars.to.regress = c("mitoRatio", "nGene"))
}
# CAR_v5, CAR_v6
DefaultAssay(CTR_STIM_Seurats[[i]])
saveRDS(CTR_STIM_Seurats, file.path(PARENT_DIR,"SC_HighLow/temp_mar/SCTransformed_seurats_24mar.rds"))

for (i in 1:length(CTR_STIM_Seurats)) {
  nm<-names(CTR_STIM_Seurats)[[i]]
  prueba_1 <- FetchData(CTR_STIM_Seurats[[i]], vars = c("CAR-pCCL-BCMA"), slot = "data")
  DefaultAssay(CTR_STIM_Seurats[[i]]) <- "RNA"
  prueba_2 <- FetchData(CTR_STIM_Seurats[[i]], vars = c("CAR-pCCL-BCMA"), slot = "data")
  
  pdf(file.path(FIGURES_DIR,paste0("Box_plot_car_exp_per_norm_",i,".pdf")),
      width=11, height=8.5, onefile = T)
  print(hist(prueba_1[,"CAR-pCCL-BCMA"], main = "SCT Transformation"))
  print(hist(prueba_2[,"CAR-pCCL-BCMA"], main = "Log Normalized"))
  print(boxplot(CTR_STIM_Seurats[[i]]@assays$SCT@data["CAR-pCCL-BCMA",], CTR_STIM_Seurats[[i]][["RNA"]]@data["CAR-pCCL-BCMA",]))
  dev.off()
  
  CTR_STIM_Seurats[[i]][["CAR_NORMALIZED_EXPRES"]] <- CTR_STIM_Seurats[[i]][["RNA"]]@data["CAR-pCCL-BCMA",]
  car_exp <- CTR_STIM_Seurats[[i]]@meta.data$CAR_NORMALIZED_EXPRES
  set.seed(123)
  car_p33 <- quantile(car_exp[car_exp>0], probs = c(0.33))
  set.seed(123)
  car_p66 <- quantile(car_exp[car_exp>0], probs = c(0.66))
  car_exp.hl <- as.data.frame(car_exp)
  car_exp.hl <- mutate(car_exp.hl, CAR_level = ifelse((car_exp>0 & car_exp<=car_p33), 1,
                                                      ifelse((car_exp>car_p33 & car_exp<=car_p66) ,2,
                                                             ifelse( (car_exp > car_p66), 3, 0 ))))
  car_exp.hl$CAR_level_COD <- plyr::mapvalues(x = car_exp.hl$CAR_level,
                                              from = c(0 , 1 , 2 , 3 ),
                                              to = c("No_CAR","Low","Med","High"))
  CTR_STIM_Seurats[[i]][["CAR_NORM_EXPRESS_HL"]]<- car_exp.hl$CAR_level
  CTR_STIM_Seurats[[i]][["CAR_NORM_level_COD"]]<- car_exp.hl$CAR_level_COD
  prueba_1 <- FetchData(CTR_STIM_Seurats[[i]],
                        vars = c("CAR_NORMALIZED_EXPRES", "CAR_NORM_level_COD")) %>%group_by(CAR_NORM_level_COD)
  bxplotcar_levels<- ggplot(prueba_1, aes(x=CAR_NORM_level_COD, y=CAR_NORMALIZED_EXPRES)) +
    geom_boxplot(outlier.shape=NA) +
    geom_jitter(position=position_jitter(width=.01, height=0))
  pdf(file.path(FIGURES_DIR,paste0("Box_plot_car_normalized_levels_",i,".pdf")),
      width=11, height=8.5, onefile = T)
  print(bxplotcar_levels)
  dev.off()
  
  CTR_STIM_Seurats[[i]][["CAR_SCT_EXPRES"]] <- CTR_STIM_Seurats[[i]][["SCT"]]@data["CAR-pCCL-BCMA",]
  car_exp <- CTR_STIM_Seurats[[i]]@meta.data$CAR_SCT_EXPRES
  set.seed(123)
  car_p33 <- quantile(car_exp[car_exp>0], probs = c(0.33))
  set.seed(123)
  car_p66 <- quantile(car_exp[car_exp>0], probs = c(0.66))
  car_exp.hl <- as.data.frame(car_exp)
  car_exp.hl <- mutate(car_exp.hl, CAR_level = ifelse((car_exp>0 & car_exp<=car_p33), 1,
                                                      ifelse((car_exp>car_p33 & car_exp<=car_p66) ,2,
                                                             ifelse( (car_exp > car_p66), 3, 0 ))))
  car_exp.hl$CAR_level_COD <- plyr::mapvalues(x = car_exp.hl$CAR_level,
                                              from = c(0 , 1 , 2 , 3 ),
                                              to = c("No_CAR","Low","Med","High"))
  CTR_STIM_Seurats[[i]][["CAR_SCT_EXPRESS_HL"]]<- car_exp.hl$CAR_level
  CTR_STIM_Seurats[[i]][["CAR_SCT_level_COD"]]<- car_exp.hl$CAR_level_COD
  prueba_1 <- FetchData(CTR_STIM_Seurats[[i]],
                        vars = c("CAR_SCT_EXPRES", "CAR_SCT_level_COD")) %>%group_by(CAR_SCT_level_COD)
  bxplotcar_levels<- ggplot(prueba_1, aes(x=CAR_SCT_level_COD, y=CAR_SCT_EXPRES)) +
    geom_boxplot(outlier.shape=NA) +
    geom_jitter(position=position_jitter(width=.01, height=0))
  pdf(file.path(FIGURES_DIR,paste0("Box_plot_car_SCT_levels_",i,".pdf")),
      width=11, height=8.5, onefile = T)
  print(bxplotcar_levels)
  dev.off()
  saveRDS(CTR_STIM_Seurats[[i]], file=file.path(RESULTS_DIR,paste0("clean_merged_seurat_SCT_namedCAR_",nm,"_230321.rds")))
  
  ridge_plot_sct <- RidgePlot(CTR_STIM_Seurats[[i]], features = c("CAR_SCT_EXPRES"), group.by = "CAR_SCT_level_COD")
  ridge_plot_norm <- RidgePlot(CTR_STIM_Seurats[[i]], features = c("CAR_NORMALIZED_EXPRES"), group.by = "CAR_NORM_level_COD")
  
  pdf(file.path(FIGURES_DIR,paste0("RidgePlots_car_levels_",i,".pdf")),
      width=11, height=8.5, onefile = T)
  print(ridge_plot_sct)
  print(ridge_plot_norm)
  dev.off()
  
  DefaultAssay(CTR_STIM_Seurats[[i]]) <- "SCT"
  #rm(car_exp, car_p33, car_p66, car_exp.hl, prueba_1, prueba_2, joined_tibble, bxplotcar_levels)
}
rm(car_exp, car_p33, car_p66, car_exp.hl, prueba_1, prueba_2, joined_tibble, bxplotcar_levels, ridge_plot_norm, ridge_plot_sct)

FIGURES_DIR<-file.path(PARENT_DIR,"SC_HighLow/temp_mar/Integrated/")
RESULTS_DIR<-file.path(PARENT_DIR,"SC_HighLow/temp_mar/Integrated/")

integ_features <- SelectIntegrationFeatures(object.list = CTR_STIM_Seurats, 
                                            nfeatures = 3000) 
CTR_STIM_Seurats <- PrepSCTIntegration(object.list = CTR_STIM_Seurats,
                                  anchor.features = integ_features)
# Find best integration anchors - can take a while to run
integ_anchors <- FindIntegrationAnchors(object.list = CTR_STIM_Seurats, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)
save.image (file =file.path(RESULTS_DIR,"23mar21_INTEGRATED_RData.RData"))

# Integrate across conditions
seurat_integrated <- IntegrateData(anchorset = integ_anchors,
                                   normalization.method = "SCT")

# Save integrated seurat object
saveRDS(seurat_integrated, file.path(PARENT_DIR,"SC_HighLow/temp_mar/Integrated/integrated_seurat_24mar.rds"))
rm(integ_anchors, CTR_STIM_Seurats)
gc(gc())

FIGURES_DIR<-file.path(PARENT_DIR,"SC_HighLow/temp_mar/Integrated/CAR_EXP")
RESULTS_DIR<-file.path(PARENT_DIR,"SC_HighLow/temp_mar/Integrated/CAR_EXP")

DefaultAssay(seurat_integrated)
#JUST ONCE
DefaultAssay(seurat_integrated) <- "integrated"
prueba_1 <- FetchData(seurat_integrated, vars = c("CAR-pCCL-BCMA"), slot = "data")
DefaultAssay(seurat_integrated) <- "RNA"
prueba_2 <- FetchData(seurat_integrated, vars = c("CAR-pCCL-BCMA"), slot = "data")
DefaultAssay(seurat_integrated) <- "SCT"
prueba_3 <- FetchData(seurat_integrated, vars = c("CAR-pCCL-BCMA"), slot = "data")

pdf(file.path(FIGURES_DIR,paste0("Box_plot_car_exp_per_norm_INTEGRATED.pdf")),
    width=11, height=8.5, onefile = T)
print(hist(prueba_1[,"CAR-pCCL-BCMA"], main = "Integrated Transformation"))
print(hist(prueba_2[,"CAR-pCCL-BCMA"], main = "Log Normalized"))
print(hist(prueba_3[,"CAR-pCCL-BCMA"], main = "SCT Transformation"))
print(boxplot(seurat_integrated[["integrated"]]@data["CAR-pCCL-BCMA",],
              seurat_integrated[["RNA"]]@data["CAR-pCCL-BCMA",],
              seurat_integrated@assays$SCT@data["CAR-pCCL-BCMA",],
              names= c("integrated", "log_norm","SCT")))
dev.off()
 
seurat_integrated[["CAR_NORMALIZED_EXPRES"]] <- seurat_integrated[["RNA"]]@data["CAR-pCCL-BCMA",]
car_exp <- seurat_integrated@meta.data$CAR_NORMALIZED_EXPRES
set.seed(123)
car_p33 <- quantile(car_exp[car_exp>0], probs = c(0.33))
set.seed(123)
car_p66 <- quantile(car_exp[car_exp>0], probs = c(0.66))
car_exp.hl <- as.data.frame(car_exp)
car_exp.hl <- mutate(car_exp.hl, CAR_level = ifelse((car_exp>0 & car_exp<=car_p33), 1,
                                                    ifelse((car_exp>car_p33 & car_exp<=car_p66) ,2,
                                                           ifelse( (car_exp > car_p66), 3, 0 ))))
car_exp.hl$CAR_level_COD <- plyr::mapvalues(x = car_exp.hl$CAR_level,
                                            from = c(0 , 1 , 2 , 3 ),
                                            to = c("No_CAR","Low","Med","High"))
seurat_integrated[["CAR_NORM_EXPRESS_HL"]]<- car_exp.hl$CAR_level
seurat_integrated[["CAR_NORM_level_COD"]]<- car_exp.hl$CAR_level_COD
prueba_1 <- FetchData(seurat_integrated,
                      vars = c("CAR_NORMALIZED_EXPRES", "CAR_NORM_level_COD")) %>%group_by(CAR_NORM_level_COD)
bxplotcar_levels<- ggplot(prueba_1, aes(x=CAR_NORM_level_COD, y=CAR_NORMALIZED_EXPRES)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(position=position_jitter(width=.01, height=0))
pdf(file.path(FIGURES_DIR,paste0("Box_plot_car_normalized_levels_INTEGRATED.pdf")),
    width=11, height=8.5, onefile = T)
print(bxplotcar_levels)
dev.off()

seurat_integrated[["CAR_SCT_EXPRES"]] <- seurat_integrated[["SCT"]]@data["CAR-pCCL-BCMA",]
car_exp <- seurat_integrated@meta.data$CAR_SCT_EXPRES
set.seed(123)
car_p33 <- quantile(car_exp[car_exp>0], probs = c(0.33))
set.seed(123)
car_p66 <- quantile(car_exp[car_exp>0], probs = c(0.66))
car_exp.hl <- as.data.frame(car_exp)
car_exp.hl <- mutate(car_exp.hl, CAR_level = ifelse((car_exp>0 & car_exp<=car_p33), 1,
                                                    ifelse((car_exp>car_p33 & car_exp<=car_p66) ,2,
                                                           ifelse( (car_exp > car_p66), 3, 0 ))))
car_exp.hl$CAR_level_COD <- plyr::mapvalues(x = car_exp.hl$CAR_level,
                                            from = c(0 , 1 , 2 , 3 ),
                                            to = c("No_CAR","Low","Med","High"))
seurat_integrated[["CAR_SCT_EXPRESS_HL"]]<- car_exp.hl$CAR_level
seurat_integrated[["CAR_SCT_level_COD"]]<- car_exp.hl$CAR_level_COD
prueba_1 <- FetchData(seurat_integrated,
                      vars = c("CAR_SCT_EXPRES", "CAR_SCT_level_COD")) %>%group_by(CAR_SCT_level_COD)
bxplotcar_levels<- ggplot(prueba_1, aes(x=CAR_SCT_level_COD, y=CAR_SCT_EXPRES)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(position=position_jitter(width=.01, height=0))
pdf(file.path(FIGURES_DIR,paste0("Box_plot_car_SCT_levels_INTEGRATED.pdf")),
    width=11, height=8.5, onefile = T)
print(bxplotcar_levels)
dev.off()
saveRDS(seurat_integrated, file.path(PARENT_DIR,"SC_HighLow/temp_mar/Integrated/CAR_EXP/integrated_seurat_24mar_CAR_COD_ALLINT.rds"))

ridge_plot_sct <- RidgePlot(seurat_integrated, features = c("CAR_SCT_EXPRES"), group.by = "CAR_SCT_level_COD")
ridge_plot_norm <- RidgePlot(seurat_integrated, features = c("CAR_NORMALIZED_EXPRES"), group.by = "CAR_NORM_level_COD")
violin_plot_sct <- VlnPlot(seurat_integrated, features = c("CAR_SCT_EXPRES"), group.by = "CAR_SCT_level_COD", pt.size = -1)
violin_plot_norm <- VlnPlot(seurat_integrated, features = c("CAR_NORMALIZED_EXPRES"), group.by = "CAR_NORM_level_COD", pt.size = -1)

pdf(file.path(FIGURES_DIR,paste0("Ridge_vln_Plots_car_levels_INTEGRATED.pdf")),
    width=11, height=8.5, onefile = T)
print(ridge_plot_sct)
print(ridge_plot_norm)
print(violin_plot_sct)
print(violin_plot_norm)
dev.off()
DefaultAssay(seurat_integrated) <-  "integrated"

n_cells <- FetchData(seurat_integrated, 
                     vars = c("ident", "orig.ident", "CAR_NORM_level_COD","CAR_SCT_level_COD")) %>%
  dplyr::count(ident, orig.ident, CAR_NORM_level_COD,CAR_SCT_level_COD) %>%
  tidyr::spread(ident, n)

# View table
View(n_cells)

write.table(n_cells, file=file.path(FIGURES_DIR,paste0("N_cells_CAR_COD_BYEXPRESSION_per_integratedsamples.txt")), col.names=TRUE, row.names=TRUE, sep="\t", quote=FALSE)
rm(car_exp, car_p33, car_p66, car_exp.hl, prueba_1, prueba_2, prueba_3, bxplotcar_levels, 
   ridge_plot_norm, ridge_plot_sct, violin_plot_norm, violin_plot_sct, n_cells)
gc(gc())

FIGURES_DIR<-file.path(PARENT_DIR,"SC_HighLow/temp_mar/Integrated/RNA_Plots/")
RESULTS_DIR<-file.path(PARENT_DIR,"SC_HighLow/temp_mar/Integrated/RNA_Plots/")

# Run PCA
set.seed(123)
seurat_integrated <- RunPCA(object = seurat_integrated)

# Plot PCA
set.seed(123)
PCAPlot(seurat_integrated, split.by = "sample")

# Run UMAP
set.seed(123)
seurat_integrated <- RunUMAP(seurat_integrated,
                             dims = 1:40,
                             reduction = "pca")

ElbowPlot(object = seurat_integrated, 
          ndims = 40)

# Determine the K-nearest neighbor graph
set.seed(123)
seurat_integrated <- FindNeighbors(object = seurat_integrated, 
                                   dims = 1:40)
# Determine the clusters for various resolutions                                
set.seed(123)
seurat_integrated <- FindClusters(object = seurat_integrated,
                                  resolution = c(0.2,0.4,0.6,0.8))

save.image (file =file.path(RESULTS_DIR,"23mar21_INTEGRATED_PCA_UMAP_ELBOW_Neighbors_RData.RData"))

# Explore resolutions
seurat_integrated@meta.data %>% 
  View()

gc(gc())
# Assign identity of clusters, Please choose accordingly and change the downstream in "FetchData..." 
Resolution_Identity <- "integrated_snn_res.0.6"
Resolution_Identity <- "integrated_snn_res.0.8"
# Resolution_Identity <- "integrated_snn_res.0.2"

#Resolution_Identity <- "integrated_snn_res.0.4"

Idents(object = seurat_integrated) <- Resolution_Identity
table(Idents(seurat_integrated))
# Run UMAP/TSNE
set.seed(123)
seurat_integrated <- RunUMAP(seurat_integrated, 
                             dims = 1:40,
                             reduction = "pca")
# Calculation of t-SNE
set.seed(123)
seurat_integrated <- RunTSNE(object = seurat_integrated,
                             dims = 1:40,
                             reduction = "pca")

# Plot UMAP                             
DimPlot(seurat_integrated)                             

DimPlot(seurat_integrated,
        split.by = "sample")

DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6)

DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6, 
        split.by = "sample")


# Plotting t-SNE
plot1 <- DimPlot(object = seurat_integrated,label = TRUE,reduction = "tsne", label.size=6)
plot1 <- plot1 + ggtitle("t-SNE")
plot2 <- DimPlot(object = seurat_integrated,label = TRUE,reduction = "tsne", label.size=6, split.by="sample")
plot2 <- plot2 + ggtitle("t-SNE")
plot3 <- plot2 + NoLegend()
plot4 <- DimPlot(seurat_integrated,reduction = "tsne",label = TRUE,label.size = 2, split.by = "donor")
plot4 <- plot4 + ggtitle("t-SNE")

# Plot the UMAP
plot5 <- DimPlot(seurat_integrated,reduction = "umap",label = TRUE,label.size = 6)
plot5 <- plot5 + ggtitle("UMAP")
plot6 <- DimPlot(seurat_integrated,reduction = "umap",label = TRUE,label.size = 6, split.by="sample")
plot6 <- plot6 + ggtitle("UMAP")
plot6 <- plot6 + NoLegend()
plot7 <- DimPlot(seurat_integrated,reduction = "umap",label = TRUE,label.size = 2, split.by = "donor")
plot7 <- plot7 + ggtitle("UMAP")


pdf(file.path(FIGURES_DIR,paste0("Tsne_UMAP_INTEGRATE",Resolution_Identity,".pdf")), 
    width=11, height=8.5, onefile = T)
print(plot1)
print(plot2)
print(plot3)
print(plot4)
print(plot5)
print(plot6)
print(plot7)
dev.off()

rm(plot1,plot2, plot3, plot4, plot5, plot6, plot7)
# Extract identity and sample information from seurat object to determine the number of cells per cluster per sample
n_cells <- FetchData(seurat_integrated, 
                     vars = c("ident", "orig.ident")) %>%
  dplyr::count(ident, orig.ident) %>%
  tidyr::spread(ident, n)

# View table
View(n_cells)
write.table(n_cells, file=file.path(FIGURES_DIR,paste0("N_cells_per_cluster_samples_",Resolution_Identity,".txt")), col.names=TRUE, row.names=TRUE, sep="\t", quote=FALSE)

##Trying https://bioinformatics.stackexchange.com/questions/11197/stacked-barplot-for-single-cell-analysis
# DF <- seurat_integrated@meta.data
# DF_Count <- DF %>%group_by(donor_by,Phase,integrated_snn_res.0.4) %>% dplyr::count() %>% ungroup() %>% group_by(donor_by,integrated_snn_res.0.4) %>%mutate(Freq = n/sum(n)*100)
# Plotcellcycle <- ggplot(DF_Count, aes(x = integrated_snn_res.0.4, y = Freq, fill = Phase))+
#   geom_col()+
#   #geom_text(aes(label = paste(round(Freq, 2),"%")),position = position_stack(vjust = 0.5))+
#   facet_wrap(~donor_by)
# pdf(file.path(FIGURES_DIR,paste0("CellCyclePhase_bydonor_bycluster_INTEGRATED_",Resolution_Identity,".pdf")),
#     width=11, height=8.5, onefile = T)
# plot(Plotcellcycle)
# dev.off()
# 
# DF_Count <- DF %>%group_by(sample,Phase,integrated_snn_res.0.4) %>% dplyr::count() %>% ungroup() %>% group_by(sample,integrated_snn_res.0.4) %>%mutate(Freq = n/sum(n)*100)
# Plotcellcycle <- ggplot(DF_Count, aes(x = integrated_snn_res.0.4, y = Freq, fill = Phase))+
#   geom_col()+
#   #geom_text(aes(label = paste(round(Freq, 2),"%")),position = position_stack(vjust = 0.5))+
#   facet_wrap(~sample)
# pdf(file.path(FIGURES_DIR,paste0("CellCyclePhase_bydonor_by_ctr_stim_INTEGRATED_",Resolution_Identity,".pdf")),
#     width=11, height=8.5, onefile = T)
# plot(Plotcellcycle)
# dev.off()



DF <- seurat_integrated@meta.data
DF_Count <- DF %>%group_by(donor_by,Phase,integrated_snn_res.0.6) %>% dplyr::count() %>% ungroup() %>% group_by(donor_by,integrated_snn_res.0.6) %>%mutate(Freq = n/sum(n)*100)
Plotcellcycle <- ggplot(DF_Count, aes(x = integrated_snn_res.0.6, y = Freq, fill = Phase))+
  geom_col()+
  #geom_text(aes(label = paste(round(Freq, 2),"%")),position = position_stack(vjust = 0.5))+
  facet_wrap(~donor_by)
pdf(file.path(FIGURES_DIR,paste0("CellCyclePhase_bydonor_bycluster_INTEGRATED_",Resolution_Identity,".pdf")),
    width=11, height=8.5, onefile = T)
plot(Plotcellcycle)
dev.off()

DF_Count <- DF %>%group_by(sample,Phase,integrated_snn_res.0.6) %>% dplyr::count() %>% ungroup() %>% group_by(sample,integrated_snn_res.0.6) %>%mutate(Freq = n/sum(n)*100)
Plotcellcycle <- ggplot(DF_Count, aes(x = integrated_snn_res.0.6, y = Freq, fill = Phase))+
  geom_col()+
  #geom_text(aes(label = paste(round(Freq, 2),"%")),position = position_stack(vjust = 0.5))+
  facet_wrap(~sample)
pdf(file.path(FIGURES_DIR,paste0("CellCyclePhase_bydonor_by_ctr_stim_INTEGRATED_",Resolution_Identity,".pdf")),
    width=11, height=8.5, onefile = T)
plot(Plotcellcycle)
dev.off()




# Explore whether clusters segregate by cell cycle phase
plot1 <- DimPlot(seurat_integrated, label = TRUE, split.by = "Phase")  + NoLegend()
metrics <-  c("nUMI", "nGene", "S.Score", "G2M.Score", "mitoRatio", "RPSRatio")
plot2 <- FeaturePlot(seurat_integrated, reduction = "umap", features = metrics, pt.size = 0.4, 
                     sort.cell = TRUE, min.cutoff = 'q10', label = TRUE)
plot3 <- FeaturePlot(seurat_integrated, reduction = "umap", features = c("nUMI", "nGene", "mitoRatio"), pt.size = 0.4, 
                     sort.cell = TRUE, min.cutoff = 'q10', label = TRUE, split.by = "sample")
plot4 <- FeaturePlot(seurat_integrated, reduction = "umap", features = c("S.Score", "G2M.Score"), pt.size = 0.4, 
                     sort.cell = TRUE, min.cutoff = 'q10', label = TRUE, split.by = "sample")
plot5 <- FeaturePlot(seurat_integrated, reduction = "umap", features = c("RPSRatio"), pt.size = 0.4, 
                     sort.cell = TRUE, min.cutoff = 'q10', label = TRUE, split.by = "sample")
pdf(file.path(FIGURES_DIR,paste0("UMAP_CellCyclePhase_othervariations_INTEGRATED_",Resolution_Identity,".pdf")),
    width=11, height=8.5, onefile = T)
print(plot1)
print(plot2)
print(plot3)
print(plot4)
print(plot5)
dev.off()
rm(plot1, plot2, plot3, plot4, plot5, Plotcellcycle,n_cells)


pdf(file.path(FIGURES_DIR,paste0("ViolinPlot_KnownMarkers_INTEGRATED_",Resolution_Identity,".pdf")),
    width=11, height=8.5, onefile = T)
print(VlnPlot(seurat_integrated, features = c("CD3D", "mitoRatio"), pt.size = -1))
print(VlnPlot(seurat_integrated, features = c("S.Score", "G2M.Score"), pt.size = -1))
print(VlnPlot(seurat_integrated, features = c("RPSRatio"), pt.size = -1))
dev.off()

# Select the RNA counts slot to be the default assay
DefaultAssay(seurat_integrated)
DefaultAssay(seurat_integrated) <- "RNA"

# Normalize RNA data for visualization purposes DONE JUST ONCE
# seurat_integrated <- NormalizeData(seurat_integrated, verbose = FALSE)
# all_genes <- rownames(x = seurat_integrated)
# seurat_integrated <- ScaleData(seurat_integrated, features = all_genes)

EMjuanro <- c("CD4","CD8A","IL7R","TCF7","CCR7","HLA-DRA","CD40LG","GZMA","PRF1","GNLY","TIGIT","LAG3","OAS1","ISG15","MX1",
              "ALDOA","PGK1","GPI","FOXP3","CTLA4","IL2RA")

plotheatm <- DotPlot(seurat_integrated,features = EMjuanro) + RotatedAxis()
pdf(file.path(FIGURES_DIR,paste0("DotPlot_SelectedGENES_perCluster_INTEGRATED_",Resolution_Identity,".pdf")),
    width=11, height=8.5, onefile = T)
print(plotheatm)
dev.off()
rm(plotheatm, EMjuanro, DF, DF_Count)

pdf(file.path(FIGURES_DIR,paste0("UMAPS_For_KnownMarkers_INTEGRATED_",Resolution_Identity,".pdf")),
    width=11, height=8.5, onefile = T)
print(FeaturePlot(seurat_integrated, reduction = "umap", features = c("CD79A", "CD19","CD14", "CAR-pCCL-BCMA","CD3D", "PTPRC"), label = T))
print(FeaturePlot(seurat_integrated, reduction = "umap", features = c("CD4", "CD8A", "EOMES", "TBX21", "FOXP3"), label = T))
print(FeaturePlot(seurat_integrated, reduction = "umap", features = c("SELL","CCR7","CXCR3", "CD27", "CD28", "CD127"), label = T))

print(FeaturePlot(seurat_integrated, reduction = "umap", features = c("IL7R","KLF2","TCF7"), label = T))
print(FeaturePlot(seurat_integrated, reduction = "umap", features = c("GZMA", "GZMB", "GZMK", "GZMH", "GNLY", "NKG7"), label = T))
print(FeaturePlot(seurat_integrated, reduction = "umap", features = c("HLA-DRA", "IL2RA", "IL3RA", "CD69", "CD40LG"), label = T))

print(FeaturePlot(seurat_integrated, reduction = "umap", features = c("E2F4", "E2F7", "CDK1"), label = T))
print(FeaturePlot(seurat_integrated, reduction = "umap", features = c("TIGIT", "LAG3","PDCD1","CTLA4"), label = T))
print(FeaturePlot(seurat_integrated, reduction = "umap", features = c("TOX", "TOX2", "IRF8", "IRF4", "NR4A2", "ENTPD1"), label = T))
dev.off()

pdf(file.path(FIGURES_DIR,paste0("ViolinPlot_KnownMarkers_INTEGRATED_",Resolution_Identity,".pdf")),
    width=11, height=8.5, onefile = T)
print(VlnPlot(seurat_integrated, features = c("CD79A", "CD19","CD14", "CAR-pCCL-BCMA","CD3D", "PTPRC"), pt.size = -1))
print(VlnPlot(seurat_integrated, features = c("CD4", "CD8A", "EOMES", "TBX21", "FOXP3"), pt.size = -1))
print(VlnPlot(seurat_integrated, features = c("SELL","CCR7","CXCR3", "CD27", "CD28", "CD127"), pt.size = -1))

print(VlnPlot(seurat_integrated, features = c("IL7R","KLF2","TCF7"), pt.size = -1))
print(VlnPlot(seurat_integrated, features = c("GZMA", "GZMB", "GZMK", "GZMH", "GNLY", "NKG7"), pt.size = -1))
print(VlnPlot(seurat_integrated, features = c("HLA-DRA", "IL2RA", "IL3RA", "CD69", "CD40LG"), pt.size = -1))

print(VlnPlot(seurat_integrated, features = c("E2F4", "E2F7", "CDK1"), pt.size = -1))
print(VlnPlot(seurat_integrated, features = c("TIGIT", "LAG3","PDCD1","CTLA4"), pt.size = -1))
print(VlnPlot(seurat_integrated, features = c("TOX", "TOX2", "IRF8", "IRF4", "NR4A2", "ENTPD1"), pt.size = -1))
dev.off()

pdf(file.path(FIGURES_DIR,paste0("UMAPS_For_Selected_Markers_INTEGRATED_",Resolution_Identity,".pdf")),
    width=11, height=8.5, onefile = T)
print(FeaturePlot(seurat_integrated, reduction = "umap", features = c("CD4", "CD8A", "GZMB", "CAR-pCCL-BCMA"), 
                  ncol=2, cols=c("blue", "cyan", "red"), 
                  min.cutoff = "q33", max.cutoff = "q66"))
dev.off()

pdf(file.path(FIGURES_DIR,paste0("UMAPS_For_KnownMarkers_INTEGRATED_SPLIT",Resolution_Identity,".pdf")),
    width=11, height=8.5, onefile = T)
print(FeaturePlot(seurat_integrated, reduction = "umap", features = c("CD79A", "CD19","CD14", "CAR-pCCL-BCMA","CD3D", "PTPRC"), label = T, split.by = "sample"))
print(FeaturePlot(seurat_integrated, reduction = "umap", features = c("CD4", "CD8A", "EOMES", "TBX21", "FOXP3"), label = T, split.by = "sample"))
print(FeaturePlot(seurat_integrated, reduction = "umap", features = c("SELL","CCR7","CXCR3", "CD27", "CD28", "CD127"), label = T, split.by = "sample"))

print(FeaturePlot(seurat_integrated, reduction = "umap", features = c("IL7R","KLF2","TCF7"), label = T, split.by = "sample"))
print(FeaturePlot(seurat_integrated, reduction = "umap", features = c("GZMA", "GZMB", "GZMK", "GZMH", "GNLY", "NKG7"), label = T, split.by = "sample"))
print(FeaturePlot(seurat_integrated, reduction = "umap", features = c("HLA-DRA", "IL2RA", "IL3RA", "CD69", "CD40LG"), label = T, split.by = "sample"))

print(FeaturePlot(seurat_integrated, reduction = "umap", features = c("E2F4", "E2F7", "CDK1"), label = T, split.by = "sample"))
print(FeaturePlot(seurat_integrated, reduction = "umap", features = c("TIGIT", "LAG3","PDCD1","CTLA4"), label = T, split.by = "sample"))
print(FeaturePlot(seurat_integrated, reduction = "umap", features = c("TOX", "TOX2", "IRF8", "IRF4", "NR4A2", "ENTPD1"), label = T, split.by = "sample"))

dev.off()

pdf(file.path(FIGURES_DIR,paste0("ViolinPlot_KnownMarkers_INTEGRATED_SPLIT",Resolution_Identity,".pdf")),
    width=11, height=8.5, onefile = T)
print(VlnPlot(seurat_integrated, features = c("CD79A", "CD19","CD14", "CAR-pCCL-BCMA","CD3D", "PTPRC"), pt.size = -1, split.by = "sample"))
print(VlnPlot(seurat_integrated, features = c("CD4", "CD8A", "EOMES", "TBX21", "FOXP3"), pt.size = -1, split.by = "sample"))
print(VlnPlot(seurat_integrated, features = c("SELL","CCR7","CXCR3", "CD27", "CD28", "CD127"), pt.size = -1, split.by = "sample"))

print(VlnPlot(seurat_integrated, features = c("IL7R","KLF2","TCF7"), pt.size = -1, split.by = "sample"))
print(VlnPlot(seurat_integrated, features = c("GZMA", "GZMB", "GZMK", "GZMH", "GNLY", "NKG7"), pt.size = -1, split.by = "sample"))
print(VlnPlot(seurat_integrated, features = c("HLA-DRA", "IL2RA", "IL3RA", "CD69", "CD40LG"), pt.size = -1, split.by = "sample"))

print(VlnPlot(seurat_integrated, features = c("E2F4", "E2F7", "CDK1"), pt.size = -1, split.by = "sample"))
print(VlnPlot(seurat_integrated, features = c("TIGIT", "LAG3","PDCD1","CTLA4"), pt.size = -1, split.by = "sample"))
print(VlnPlot(seurat_integrated, features = c("TOX", "TOX2", "IRF8", "IRF4", "NR4A2", "ENTPD1"), pt.size = -1, split.by = "sample"))
dev.off()

pdf(file.path(FIGURES_DIR,paste0("UMAPS_For_Selected_Markers_INTEGRATED_SPLIT",Resolution_Identity,".pdf")),
    width=11, height=8.5, onefile = T)
print(FeaturePlot(seurat_integrated, reduction = "umap", features = c("CD4", "CD8A", "GZMB", "CAR-pCCL-BCMA"), ncol=2, cols=c("blue", "cyan", "red"), 
                  min.cutoff = "q33", max.cutoff = "q66",  split.by = "sample"))
dev.off()

source(file.path(PARENT_DIR,"SC_HighLow/temp_mar/Integrated/StackedVlnPlot.R"))
features<- c("CD4","CD8A","CAR-pCCL-BCMA")
plot111 <- StackedVlnPlot(obj = seurat_integrated, features = features)
pdf(file.path(FIGURES_DIR,paste0("StackedViolin_For_KnownMarkers_",Resolution_Identity,".pdf")),
    width=11, height=8.5, onefile = T)
print(plot111)
dev.off()
rm(plot111)
########## NOSE SI AL FINAL AQUI HACER LO DE CAR_EXP POR seurat_splitted.... 

##########

# Find markers for every cluster compared to all remaining cells, 
# Since we have samples representing different conditions in our dataset, our best option is to find conserved markers.
DefaultAssay(seurat_integrated) <- "RNA"

#annotations <- read.csv("cluster_juanro/SC_HighLow/temp/annotation.csv")

Idents(object = seurat_integrated) <- Resolution_Identity
# Create function to get conserved markers for any given cluster
get_conserved <- function(cluster){
  FindConservedMarkers(seurat_integrated,
                       ident.1 = cluster,
                       grouping.var = "sample",
                       only.pos = TRUE) %>%
    rownames_to_column(var = "gene") %>%
    left_join(y = unique(annotations[, c("gene_name", "description")]),
              by = c("gene" = "gene_name")) %>%
    cbind(cluster_id = cluster, .)
}
# Iterate function across desired clusters
#conserved_markers <- map_dfr(c(0:14), get_conserved) #For resolution 0.4
#conserved_markers <- map_dfr(c(0:18), get_conserved) #For resolution 0.6
conserved_markers <- map_dfr(c(0:22), get_conserved) #For resolution 0.8

# Extract top 10 markers per cluster
top10 <- conserved_markers %>% 
  mutate(avg_fc = (ctrl_avg_logFC + stim_avg_logFC) /2) %>% 
  group_by(cluster_id) %>% 
  top_n(n = 10, 
        wt = avg_fc)

write.table(top10, file=file.path(FIGURES_DIR,paste0("Top_10_ConserverdMarkers_perCluster_persample",Resolution_Identity,".txt")), col.names=TRUE, row.names=TRUE, sep="\t", quote=FALSE)
write.table(conserved_markers, file=file.path(FIGURES_DIR,paste0("ConserverdMarkers_perCluster_persample",Resolution_Identity,".txt")), col.names=TRUE, row.names=TRUE, sep="\t", quote=FALSE)

#Resolution_Identity <- "integrated_snn_res.0.4"
Idents(object = seurat_integrated) <- Resolution_Identity
table(Idents(seurat_integrated))

conserved_markers

plot4 <- DimPlot(seurat_integrated,reduction = "umap",label = TRUE,label.size = 6)
plot4 <- plot4 + ggtitle("UMAP")

clu <- 1
FeaturePlot(object = seurat_integrated, features = top10[top10$cluster_id == clu, "gene"] %>% pull(gene))

head(conserved_markers %>% pull(gene))
head(conserved_markers)
#Six is again a mitocondrial cluster
#CTR_STIM_clusts <- c(0:5,7:14) #For resolution 0.4
#CTR_STIM_clusts <- c(0:7, 9:18)#For resolution 0.6
CTR_STIM_clusts <- c(0:6, 8:17, 19:22)#For resolution 0.8
clusts <-CTR_STIM_clusts
print(clusts)
  
for (clu in clusts){
  pdf(file.path(FIGURES_DIR,paste0("FP_Top12Markers_Cluster_",clu,"_",Resolution_Identity,".pdf")), ###Need to attach the chose one
      width=11, height=8.5, onefile = T)
  print(plot4)
  print(FeaturePlot(object = seurat_integrated, features = top10[top10$cluster_id == clu, "gene"] %>% pull(gene)))
  dev.off()
  pdf(file.path(FIGURES_DIR,paste0("VP_Top12Markers_Cluster_",clu,"_",Resolution_Identity,".pdf")), ###Need to attach the chose one
      width=11, height=8.5, onefile = T)
  print(plot4)
  print(VlnPlot(object = seurat_integrated, features = top10[top10$cluster_id == clu, "gene"] %>% pull(gene), pt.size = -1))
  dev.off()
  #head(ann_markers %>% pull(gene))
  
  #assign(paste0("Markers_for_",clu,"_",phase), ann_markers %>% dplyr::filter(cluster==clu) %>%pull(gene),envir = .GlobalEnv)
  Gene.list <- conserved_markers %>% dplyr::filter(cluster_id==clu) %>%pull(gene)
  Translated.genes.list <- bitr(Gene.list,fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  
  EnrichGo.MF <- enrichGO(Translated.genes.list$ENTREZID,OrgDb= org.Hs.eg.db,ont= "MF",
                          pAdjustMethod = "BH",pvalueCutoff  = 0.05,
                          qvalueCutoff  = 0.05,readable= TRUE)
  EnrichGo.BP <- enrichGO(Translated.genes.list$ENTREZID,OrgDb= org.Hs.eg.db,ont= "BP",
                          pAdjustMethod = "BH",pvalueCutoff  = 0.05,
                          qvalueCutoff  = 0.05,readable= TRUE)
  EnrichGo.CC <- enrichGO(Translated.genes.list$ENTREZID,OrgDb= org.Hs.eg.db,ont= "CC",
                          pAdjustMethod = "BH",pvalueCutoff  = 0.05,
                          qvalueCutoff  = 0.05,readable= TRUE)
  
  EnrichKEGG <- enrichKEGG(Translated.genes.list$ENTREZID,organism = "hsa",pAdjustMethod = "BH",
                           pvalueCutoff  = 0.05)
  EnrichPathway <- enrichPathway(Translated.genes.list$ENTREZID)
  
  ####### ON THE NEXT ROUND I'LL HAVE TO ADD "named as" @ the end... 
  pdf(file.path(FIGURES_DIR,paste0("Cluster_",as.character(clu),"_",Resolution_Identity,".pdf")), width = 11, height = 8.5)
  print(dotplot(EnrichGo.MF, title = paste("GO.Molecular Function of Cluster: ",clu)))
  print(dotplot(EnrichGo.BP, title = paste("GO.Biological Process of Cluster: ",clu)))
  print(dotplot(EnrichGo.CC, title = paste("GO.Cellular Component of Cluster: ",clu)))
  print(dotplot(EnrichKEGG, title = paste("KEGG of Cluster: ",clu)))
  print(dotplot(EnrichPathway, title = paste("Pathway of Cluster: ",clu)))
  dev.off()
  
  write.table(EnrichGo.MF, file.path(FIGURES_DIR,paste0("GO_MF_Cluster",clu,"_",Resolution_Identity,".tsv")), col.names=TRUE, row.names=TRUE, sep="\t", quote=FALSE)
  write.table(EnrichGo.BP, file.path(FIGURES_DIR,paste0("GO_BP_Cluster",clu,"_",Resolution_Identity,".tsv")), col.names=TRUE, row.names=TRUE, sep="\t", quote=FALSE)
  write.table(EnrichGo.CC, file.path(FIGURES_DIR,paste0("GO_CC_Cluster",clu,"_",Resolution_Identity,".tsv")), col.names=TRUE, row.names=TRUE, sep="\t", quote=FALSE)
  write.table(EnrichKEGG, file.path(FIGURES_DIR,paste0("KEGG_Cluster",clu,"_",Resolution_Identity,".tsv")), col.names=TRUE, row.names=TRUE, sep="\t", quote=FALSE)
  write.table(EnrichPathway, file.path(FIGURES_DIR,paste0("Pathway_Cluster",clu,"_",Resolution_Identity,".tsv")), col.names=TRUE, row.names=TRUE, sep="\t", quote=FALSE)
}

rm(conserved_markers, EnrichGo.BP, EnrichGo.CC, EnrichGo.MF, EnrichKEGG, EnrichPathway, 
   plot4, top10, Translated.genes.list, Gene.list, clu, clusts, CTR_STIM_clusts)
gc()
gc(gc())
saveRDS(seurat_integrated, file.path(PARENT_DIR,"SC_HighLow/temp_mar/Integrated/integrated_seurat_clusters_afterpipeline_31mar.rds"))

################################################################################
## Saving and etcs 
################################################################################
save.image (file =file.path(RESULTS_DIR,"24MARCH21.RData"))
savehistory(file =file.path(RESULTS_DIR,"24MARCH21.Rhistory"))
sink(file =file.path(RESULTS_DIR,"24MARCH21.txt"))
toLatex(sessionInfo())
sink(NULL)

save.image (file =file.path(RESULTS_DIR,"31MARCH21.RData"))
savehistory(file =file.path(RESULTS_DIR,"31MARCH21.Rhistory"))
sink(file =file.path(RESULTS_DIR,"31MARCH21.txt"))
toLatex(sessionInfo())
sink(NULL)


rm(list=ls(pattern = "^BS_"))
rm(list=ls(pattern = "^ST_"))
rm(all_genes, integ_features, metrics, nm, i, id, features, edb)


