#!/usr/bin/env Rscript
#===============================================================================
#' Author: Maria E. Calleja copied from Gserrano SimiC_post.R (literally)
#' Date: 2020/Jan
#===============================================================================
# salloc -c 16 --mem=50G -p medium
# module load Python
# module load R/3.5.1-foss-2018b
# module load Seurat/3.0.2-foss-2018b-R-3.5.1
# module load Seurat-wrappers/20190818-foss-2018b-R-3.5.1
module load binutils/2.30-GCCcore-7.3.0
#===============================================================================

## PACKAGES.
#library(reticulate)
library(reticulate)
reticulate::use_python("/usr/bin/python3")
library(Seurat)
library(cowplot)
library(dplyr)
library(future)
library(viridis)
library(ggplot2)
library(plyr)
library(reshape2)
library(gridExtra)
library(ggridges)
library(stringr)
library(hues)
library(pheatmap)
library(Rtsne)

options(stringsAsFactors = FALSE)
#===============================================================================
## GLOBAL VARIABLES.
#Local or Cluster
#/home/sevastopol/data/mcallejac/JuanRo_SimiC/data/CART_UFILTER/
PARENT_DIR<-"/home/sevastopol/data/mcallejac"
PROJECT_DIR<-file.path(PARENT_DIR,"JuanRo_SimiC/data/CART_UFILTER/")
FIGURES_DIR<-file.path(PARENT_DIR,"JuanRo_SimiC/data/CART_UFILTER/Plots/")
#===============================================================================
## FUNCTIONS.
##
#===============================================================================
## MAIN.
################################################################################
setwd(FIGURES_DIR)
list.files()

get_plot_dims <- function(heat_map)
{
  plot_height <- sum(sapply(heat_map$gtable$heights, grid::convertHeight, "in"))
  plot_width  <- sum(sapply(heat_map$gtable$widths, grid::convertWidth, "in"))
  dev.off()
  return(list(height = plot_height, width = plot_width))
}


#data_root_dir <- "/home/sevastopol/data/gserranos/SimiC/Simic_Validation/Data/"
data_root_dir <- PROJECT_DIR
#plot_root_dir <- "/home/sevastopol/data/gserranos/SimiC/Simic_Validation/Plots/Test/"
plot_root_dir <- FIGURES_DIR

#file_idx <- 'ctrl_CAR_INT'
file_root <- 'CART_UFILTER_INT'

#file_root <- 'ClonalKinetics_filtered'
params <- '_L10.001_L20.01_'
file_idx <- paste0(file_root, params)
plasma <- viridis(50, direction = 1, option = "C")


# The phenotypes studied, be aware that must have the same order to the Simic input 
phenotypes <-  c('ctrl', 'stim')

# binaryze the phenotypes with base 0
assignment <- as.character(seq(0,length(phenotypes)-1))
names(assignment) <- phenotypes

# Load the weigths from SimiC
weights_file <- paste0(data_root_dir,file_idx,"Ws_BIC_filtered_BIC.pickle")
#weights_file <- paste0(data_root_dir,file_idx,"Ws_filtered_BIC.pickle")
SimiC_weights <- py_load_object(filename =weights_file)

# filter low R2 targets
unselected_targets <- list()
for (phenotype in phenotypes){
    unselected_targets[[phenotype]] <- SimiC_weights$query_targets[which(SimiC_weights$adjusted_r_squared[[assignment[[phenotype]]]] < 0.7)]

}

pdf(paste0(plot_root_dir ,file_idx,'Simic_Filtered_R2_targets_hist.pdf'))
for (phenotype in phenotypes){
    selectedTargets <- which(SimiC_weights$adjusted_r_squared[[assignment[[phenotype]]]] > 0.7)
    hist(SimiC_weights$adjusted_r_squared[[assignment[[phenotype]]]], col='grey', breaks=100, 
    xlab = 'Adjusted R2',
    main = paste0('Phenotype: ',phenotype ,'\n Targets selected: ', length(selectedTargets), ', mean R2: ', mean(SimiC_weights$adjusted_r_squared[[assignment[[phenotype]]]][selectedTargets])))
}
dev.off()


SimiC_weights_df <- data.frame(driver=NULL, target = NULL, value = NULL, .id = NULL, stringsAsFactors=F)
for(phenotype in phenotypes){
    tmp <-SimiC_weights$weight_dic[[assignment[[phenotype]]]][-nrow(SimiC_weights$weight_dic[[assignment[[phenotype]]]]),]
    rownames(tmp)<-SimiC_weights$TF_ids
    colnames(tmp)<-SimiC_weights$query_targets

    tmp <- setNames(melt(tmp), c('driver', 'target', 'value'))
    tmp$.id <- phenotype
    SimiC_weights_df <- rbind(SimiC_weights_df, tmp)
}


TF_2_remove <- c()
pdf(paste0(plot_root_dir ,file_idx,'Simic_TF_weigths.pdf'), onefile = TRUE, width=20)
for(drv in unique(SimiC_weights_df$driver)){
  message(drv)
  tmp_plotter <- SimiC_weights_df[SimiC_weights_df$driver == drv,]
  # remove the unselected targets for each phenotype
  tmp_plotter_filtered <- data.frame(driver=NULL,  target=NULL,  value=NULL,  .id=NULL , stringsAsFactors = F)
  for (phenotype in phenotypes){
      tmp_plotter_filtered <- rbind(tmp_plotter_filtered, tmp_plotter[tmp_plotter$.id == phenotype & !tmp_plotter$target %in% c(unselected_targets[[phenotype]]) ,])
  }
  #
  tmp_plotter_filtered <- tmp_plotter_filtered[order(abs(tmp_plotter_filtered$value), decreasing=T),]
  bests <- unique(tmp_plotter_filtered$target)[1:100]
  tmp_plotter_filtered <- tmp_plotter_filtered[tmp_plotter_filtered$target %in% bests,]
  tmp_plotter_filtered$target <- factor(tmp_plotter_filtered$target, levels = unique(tmp_plotter_filtered$target))

  if(sum(abs(tmp_plotter_filtered$value)) == 0){
    TF_2_remove <- c(TF_2_remove, drv)
    drv <- paste0(drv, '(deleted)')
  }
  p <- ggplot(tmp_plotter_filtered, aes(x=target, y=value, fill=.id)) + 
    geom_bar(stat='identity', position='dodge', color='black') + 
    scale_fill_iwanthue() + 
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=1, size=4)) + ggtitle(drv) 
  print(p)
}
dev.off()


pdf(paste0(plot_root_dir ,file_idx,'Simic_Targets_weigths.pdf'), onefile = TRUE, width=20)
plot_counter <- 1
plot_list <- list()
for(tgt in sort(as.character(unique(SimiC_weights_df$target)))){
  show(tgt)
  # get the target genes for each TF
  # if the target is discarded on all the phenotypes
  if(tgt %in% Reduce(intersect, unselected_targets)){
    next
  # if the target is discarded on at least one of the phenotypes
  }else if(tgt %in% unique(unlist(unselected_targets))){
    tmp_plotter <- data.frame(driver=NULL, target =NULL, value=NULL, .id=NULL)
    for (phenotype in phenotypes){
      if ( !tgt %in% unselected_targets[[phenotype]]){
        show(phenotype)
         tmp_plotter <- rbind(tmp_plotter, SimiC_weights_df[SimiC_weights_df$target == tgt & SimiC_weights_df$.id == phenotype,])
      }
    }  
  }else{
    tmp_plotter <- SimiC_weights_df[SimiC_weights_df$target == tgt,]
  }
  # tmp_plotter$value <- scale(tmp_plotter$value, center=FALSE)
  tmp_plotter <- tmp_plotter[order(abs(tmp_plotter$value), decreasing=T),]
  assign(paste0('p', plot_counter), ggplot(tmp_plotter, aes(x=reorder(driver, -abs(value)), y=value, fill=.id,  palette = "jco")) + 
    geom_bar(stat='identity', position='dodge', color='black') + 
    scale_fill_iwanthue() + 
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=1, size=4)) + ggtitle(tgt) )
  if(plot_counter == 4){
    grid.arrange(p1, p2, p3, p4, nrow=2, ncol=2)
    plot_list <- list()
    plot_counter <- 1
    print('plotted')
  }else{
    plot_counter <- plot_counter +1
  }
}
dev.off()


# Load the AUC data per phenotype
AUCs <- list()
for(i in seq_along(phenotypes)){
  index <- as.character(i-1)
  show(index)
  tmp <- read.table(paste0(data_root_dir,file_idx,"AUCs_filtered_BIC_",index,"_BIS.csv"), header=T, sep='\t')
  rownames(tmp) <- tmp$X
  tmp <- tmp[,!colnames(tmp) %in% c('X')]
  AUCs[[index]] <- tmp

}

#readRDS(file=file.path(PARENT_DIR,"JuanRo_SimiC/data/clean_merged_seuratintegrated_CAR_INTEGRATED_named_070121.rds"))
seurat_data <- readRDS(file=file.path(PARENT_DIR,"JuanRo_SimiC/data/CART_UFILTER/integrated_seurat_clusters_afterpipeline.rds"))

cell_populations <- c('ctrl', 'stim')
Idents(seurat_data) <- seurat_data$sample
seurat_data <- subset(seurat_data, idents=cell_populations)
#table(seurat_data@meta.data$sample)
#ctrl  stim 
#24453 19528 
#table(Idents(seurat_data))
#0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
#5033 4140 3550 3504 3315 2778 2765 2652 2537 2357 1976 1917 1715 1482 1235 1091 
#16   17   18 
#894  527  513 

labels_2_keep <- c(0,2:7,9:10,12,15:17)  #In this case, i was extracting 01 dividing cells, 8 mito genes, 
## 11 dividing cells, 13,14, 18 cd3 negative. This correspond to to integrated_snn_res.0.6 clusters... 
cells_2_keep_1 <- names(seurat_data$integrated_snn_res.0.6 [seurat_data$integrated_snn_res.0.6 %in% labels_2_keep])
labels_2_keep <- c("ctrl", "stim")
cells_2_keep_2 <- names(seurat_data$sample [seurat_data$sample %in% labels_2_keep])
cells_2_keep <- intersect(cells_2_keep_1, cells_2_keep_2)
length(cells_2_keep)

seurat_data <- subset(seurat_data, cells = cells_2_keep)
#table(seurat_data@meta.data$sample)
#ctrl  stim 
#17101 15056
#table(Idents(seurat_data))
#0    2    3    4    5    6    7    9   10   12   15   16   17 
#5033 3550 3504 3315 2778 2765 2652 2357 1976 1715 1091  894  527 

Idents(seurat_data) <- "integrated_snn_res.0.8"
table(Idents(seurat_data))
#No cluster 13, thanks... 
#0    1    2    3    4    5    6    7    8    9   10   11   12   14   15   16 
#4218 3709    3 3283 3166 2630 2612  135 2353 2140 1926    3 1714   19 1079 1063 
#17   18   19   20   21   22 
#933   71  605    1  457   37 

seurat_data <- subset(seurat_data, idents = c(2,7,11,14,18,20,22), invert = TRUE)
table(Idents(seurat_data))
table(Idents(seurat_data))
#0    1    3    4    5    6    8    9   10   12   15   16   17   19   21 
#4218 3709 3283 3166 2630 2612 2353 2140 1926 1714 1079 1063  933  605  457 

#mt_genes <- grep("^MT-", rownames(seurat_data@assays$RNA@counts), value = TRUE)
#rps_genes <- grep("^RP[SL]", rownames(seurat_data@assays$RNA@counts), value = TRUE)
#genes_2_discard <- union(mt_genes, rps_genes)
#hemato_data_raw <- hemato_data_raw[ !rownames(hemato_data_raw) %in% genes_2_discard, ]
#seurat_data <- subset(seurat_data, features = mt_genes, invert = TRUE) (esto claro que no funciona...)
#seurat_data <- GetAssayData(seurat_data)[!mt_genes,]  ###VOYY A OBVIARLO POR EL MOMENTO###

# remove undesired genes
#mt_genes <- names(which(str_detect(rownames(seurat_data_raw), "^MT-")))
#mt_genes <- grep("^MT-", rownames(seurat_data_raw), value = TRUE)

Resolution_Identity <- "integrated_snn_res.0.8"
Idents(object = seurat_data) <- Resolution_Identity
table(Idents(seurat_data))
CTR_STIM_clusts_cd4_cd8 <- c(0,1,3,4,5,6,8,9,10,12,15,16,17,19,21 )
named_clusters_withNames_cd4_cd8 <-   plyr::mapvalues(c(0,1,3,4,5,6,8,9,10,12,15,16,17,19,21 ),
                                              from = c("0","1","3",
                                                       "4","5","6",
                                                       "8", "9", 
                                                       "10","12","15",
                                                       "16","17","19", "21"),
                                              to = c("C0.CD4 Early Memory","C1.CD4 Th2 helper", "C3.CD8 Memory",
                                                     "C4.CD4 Memory","C5.CD4 Early Memory","C6.CD4 Activated",
                                                     "C8.CD8 Cytotoxic", "C9.CD8 Cytotoxic (late)",
                                                     "C10.CD4 Cytotoxic", "C12.CD4 Th2 helper", "C15.CD4 Glycolisis",
                                                     "C16.CD4 INF response","C17.CD4 Activated", "C19.CD4 Treg", "21.CD4 Cytotoxic"))

clusts <-CTR_STIM_clusts_cd4_cd8
print(clusts)
print(named_clusters_withNames_cd4_cd8)
named_clusters_Idents <- c(named_clusters_withNames_cd4_cd8)
names(named_clusters_Idents) <- levels(seurat_data)
seurat_data <- RenameIdents(object = seurat_data, named_clusters_Idents)
seurat_data[["ClusterNames_0.8_by_JR"]] <- Idents(object = seurat_data) 
Idents(seurat_data)
table(Idents(seurat_data))

# set the phenotype and cluster per TF
Idents(seurat_data) <- seurat_data$sample
AUCs_by_state<-list()
state_specific_AUCs<-NULL
for(phenotype in phenotypes){
  cell_names_phenotype <- colnames(subset(seurat_data, idents = phenotype))
  AUCs_cell_i<-AUCs[[assignment[phenotype]]]
  AUCs_cell_i$.id <- phenotype
  AUCs_by_state[[phenotype]]<-AUCs_cell_i[rownames(AUCs_cell_i) %in%  cell_names_phenotype,]
  state_specific_AUCs<-rbind(state_specific_AUCs, AUCs_cell_i[rownames(AUCs_cell_i) %in%  cell_names_phenotype,])
}

df <- do.call('rbind', AUCs_by_state)
df$cell_id<-rownames(state_specific_AUCs)

#load('/home/sevastopol/data/gserranos/SimiC/Simic_Validation/Data/clonalKinetics_subseted_noCARexpres.Robj') #subseted_object_forsimic

clusters_df <- setNames(as.data.frame(seurat_data$ClusterNames_0.8_by_JR),"cluster_id")
clusters_df$cell_id <- rownames(clusters_df)

df_w_cluster <- merge(df,clusters_df,by="cell_id")
df_auc <- melt(df_w_cluster ,  id.vars = c('.id','cell_id', 'cluster_id'), variable.name = 'driver', stringsAsFactors =F)

# check the population per phenotype and cluster
tmp <- with(df_w_cluster, table(cluster_id, .id))
tmp <- melt(tmp)

pdf(paste0(plot_root_dir ,file_idx,'Simic_Population_phenotype.pdf'))
  ggplot(tmp, aes(x=cluster_id, y=value, fill=.id)) + 
  geom_bar(stat='identity' , position='dodge',  color='black') + 
  geom_text(aes(label=value), position=position_dodge(width=0.9), vjust=-0.25) +
  scale_fill_iwanthue() + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=1, size=4))
dev.off()


##### keep the clusters with >10 representatnts in both populations

clusters_2_keep <- as.data.frame.matrix(with(df_auc, table(cluster_id, .id)))
clusters_2_keep$cluster_id <- rownames(clusters_2_keep)
clusters_2_keep <- clusters_2_keep[rowSums(clusters_2_keep > 10) == ncol(clusters_2_keep), 'cluster_id']

df_auc_common <- df_auc[df_auc$cluster_id %in% clusters_2_keep, ]

# calculate the regulation dissimilarity by the total variation score
MinMax_clust<-NULL
for(cluster in clusters_2_keep){
  tmp <- df_auc_common[df_auc_common$cluster_id==cluster,]
  MinMax_val <- NULL
  for (tf in unique(df_auc_common$driver)){
    n_breaks = 100
    Wauc_dist <- list()
    for (phenotype in phenotypes){
      Wauc_dist[[phenotype]] <- hist(tmp[tmp$.id ==phenotype & tmp$driver ==tf, 'value'], breaks = c(seq(0,1, 1/n_breaks)), plot=FALSE)$density
    }
    mat <-  do.call('rbind', Wauc_dist)
    mat <- mat[complete.cases(mat),]
    minmax_diff <- apply(na.omit(mat), 2, max) - apply(na.omit(mat), 2, min)
    variant <- sum(abs(minmax_diff)) / n_breaks
    variant <- variant/ sum(rowSums(mat)!=0)
    MinMax_val <- append(MinMax_val,variant)
  }
  MinMax_val <- setNames(as.data.frame(MinMax_val), c(cluster))
  rownames(MinMax_val) <-  unique(df_auc_common$driver)
  if (is.null(MinMax_clust)){
    MinMax_clust <- MinMax_val
  }else{
    MinMax_clust<-cbind(MinMax_clust,MinMax_val)
  }
}



# plot the densities of the AUC and the score per TF

for (clust in clusters_2_keep){
  show(clust)
  plotter2 <- df_auc[df_auc$cluster_id == clust,]
  pdf(paste0(plot_root_dir ,file_idx,'Simic_Auc_Cluster_',clust,'.pdf'), width = 15, onefile = TRUE)
  plot_counter <- 1
  for (tf in unique(plotter2$driver)){
    assign( paste0('p', plot_counter), 
    ggplot(plotter2[plotter2$driver ==tf,], aes(x=value, fill=.id)) + 
    geom_density(alpha = 0.6, adjust = 1/8) + theme_classic() + 
    scale_fill_iwanthue() +
    theme(legend.position = 'top')+ geom_rug() + 
    ggtitle(paste0(tf, '   ', MinMax_clust[rownames(MinMax_clust) == tf, clust])) )
    if(plot_counter == 2){
      grid.arrange(p1, p2, ncol=2)
      plot_counter <- 1
    }else{
      plot_counter <- plot_counter +1
    }
  }
  grid.arrange(p1, p2 , ncol=2)
  dev.off()
}


clust_order_asc<-names(sort(apply(MinMax_clust, 2, mean)))
MinMax_clust<-MinMax_clust[,clust_order_asc]
MinMax_df <- as.data.frame(MinMax_clust)
MinMax_df$driver <- rownames(MinMax_df)
MinMax_df <- melt(MinMax_df,variable.name = "cluster_id")
MinMax_clust <- MinMax_clust[!rownames(MinMax_clust) %in% TF_2_remove,]


# Plot the heatmap of the regulatory dissimilarity score
p <- pheatmap(MinMax_clust,color=plasma, fontsize=5, angle_col =45, cellwidth=40)
p <- pheatmap(MinMax_clust,color=plasma, angle_col =45, cellwidth=40, cellheight=10)
plot_dims <- get_plot_dims(p)
      
pdf(paste0(plot_root_dir ,file_idx,'Simic_HeatMaps_RegDissScore.pdf'), height = plot_dims$height, width = plot_dims$width )
print(p)
dev.off()

cd8_clusts <- c("C8.CD8 Cytotoxic","C9.CD8 Cytotoxic (late)", "C3.CD8 Memory")
p <- pheatmap(MinMax_clust[, cd8_clusts],color=plasma, angle_col =45, cellheight = 8)
plot_dims <- get_plot_dims(p)
pdf(paste0(plot_root_dir ,file_idx,'Simic_HeatMaps_RegDissScore_CD8.pdf.pdf'), height = plot_dims$height, width = plot_dims$width )
print(p)
dev.off()

# for (clust in clusters_2_keep){
#   show(clust)
#   plotter2 <- df_auc[df_auc$cluster_id == clust,]

#   pdf(paste0(plot_root_dir ,file_idx,'Simic_Auc_Cluster_',clust,'Combined.pdf'), width = 15, onefile = TRUE)
#   ggplot(plotter2, aes(x=value, fill=.id)) + 
#   geom_density(alpha = 0.6, adjust = 1/8) + theme_classic() + 
#   scale_fill_iwanthue() +
#   theme(legend.position = 'top')+ geom_rug() 
#   dev.off()
  
#   pdf(paste0(plot_root_dir ,file_idx,'Simic_Auc_Clsuter_',clust,'.pdf'), width = 15, onefile = TRUE)
#   plot_counter <- 1
#   for (tf in unique(plotter2$driver)){
#     assign( paste0('p', plot_counter), 
#     ggplot(plotter2[plotter2$driver ==tf,], aes(x=value, fill=.id)) + 
#     geom_density(alpha = 0.6, adjust = 1/8) + theme_classic() + 
#     scale_fill_iwanthue() +
#     theme(legend.position = 'top')+ geom_rug() + 
#     ggtitle(paste0(tf, '   ', MinMax_clust[rownames(MinMax_clust) == tf, clust])) )
#     if(plot_counter == 2){
#       grid.arrange(p1, p2, ncol=2)
#       plot_counter <- 1
#     }else{
#       plot_counter <- plot_counter +1
#     }
#   }
#   grid.arrange(p1, p2 , ncol=2)
#   dev.off()
# }

###I've added this to capture the AUC differences between the two conditions, disregarding the cluster. 
MinMax_clust<-NULL
tmp <- df_auc_common
MinMax_val <- NULL
for (tf in unique(df_auc_common$driver)){
  n_breaks = 100
  Wauc_dist <- list()
  for (phenotype in phenotypes){
    Wauc_dist[[phenotype]] <- hist(tmp[tmp$.id ==phenotype & tmp$driver ==tf, 'value'], breaks = c(seq(0,1, 1/n_breaks)), plot=FALSE)$density
  }
  mat <-  do.call('rbind', Wauc_dist)
  mat[is.na(mat)] <- 0
  minmax_diff <- apply(na.omit(mat), 2, max) - apply(na.omit(mat), 2, min)
  variant <- sum(abs(minmax_diff)) / n_breaks
  variant <- variant/ sum(rowSums(mat)!=0)
  MinMax_val <- append(MinMax_val,variant)
}
MinMax_clust<-cbind(MinMax_clust,MinMax_val)
rownames(MinMax_clust) <- unique(df_auc_common$driver)

plotter2 <- df_auc
pdf(paste0(plot_root_dir ,file_idx,'Simic_Auc_ALLCLUST.pdf'), width = 15, onefile = TRUE)
plot_counter <- 1
for (tf in unique(plotter2$driver)){
  assign( paste0('p', plot_counter), 
  ggplot(plotter2[plotter2$driver ==tf,], aes(x=value, fill=.id)) + 
  geom_density(alpha = 0.6, adjust = 1/8) + theme_classic() + 
  scale_fill_iwanthue() +
  theme(legend.position = 'top')+ geom_rug() + 
  ggtitle(paste0(tf, '   ', MinMax_clust[rownames(MinMax_clust) == tf])) )
  if(plot_counter == 2){
    grid.arrange(p1, p2, ncol=2)
    plot_counter <- 1
  }else{
    plot_counter <- plot_counter +1
  }
}
grid.arrange(p1, p2 , ncol=2)
dev.off()
###I've added the latter and it finishes here. 


################################################################################
## Saving and etcs 
################################################################################
save.image (file =file.path(PARENT_DIR,"JuanRo_SimiC/RSession/210413_CART_UFILTER.RData"))
savehistory(file =file.path(PARENT_DIR,"JuanRo_SimiC/RSession/210413_CART_UFILTER.Rhistory"))
sink(file =file.path(PARENT_DIR,"JuanRo_SimiC/RSession/210413_CART_UFILTER.txt"))
toLatex(sessionInfo())
sink(NULL)


################################################################################
## Re load for JuanRos petitions
################################################################################
load(file =file.path(PARENT_DIR,"JuanRo_SimiC/RSession/210413_CART_UFILTER.RData"))

#cd8_clusts <- c("C8.CD8 Cytotoxic","C9.CD8 Cytotoxic (late)", "C3.CD8 Memory")

clusters_2_keep <- as.data.frame.matrix(with(df_auc, table(cluster_id, .id)))
clusters_2_keep$cluster_id <- rownames(clusters_2_keep)
clusters_2_keep <- clusters_2_keep[rowSums(clusters_2_keep > 10) == ncol(clusters_2_keep), 'cluster_id']
clusters_2_keep <- clusters_2_keep[clusters_2_keep %in% cd8_clusts]

df_auc_common <- df_auc[df_auc$cluster_id %in% clusters_2_keep, ]

# calculate the regulation dissimilarity by the total variation score
MinMax_clust<-NULL
for(cluster in clusters_2_keep){
  tmp <- df_auc_common[df_auc_common$cluster_id==cluster,]
  MinMax_val <- NULL
  for (tf in unique(df_auc_common$driver)){
    n_breaks = 100
    Wauc_dist <- list()
    for (phenotype in phenotypes){
      Wauc_dist[[phenotype]] <- hist(tmp[tmp$.id ==phenotype & tmp$driver ==tf, 'value'], breaks = c(seq(0,1, 1/n_breaks)), plot=FALSE)$density
    }
    mat <-  do.call('rbind', Wauc_dist)
    mat <- mat[complete.cases(mat),]
    minmax_diff <- apply(na.omit(mat), 2, max) - apply(na.omit(mat), 2, min)
    variant <- sum(abs(minmax_diff)) / n_breaks
    variant <- variant/ sum(rowSums(mat)!=0)
    MinMax_val <- append(MinMax_val,variant)
  }
  MinMax_clust<-cbind(MinMax_clust,MinMax_val)
}
colnames(MinMax_clust)<-unique(df_auc_common$cluster_id)
rownames(MinMax_clust)<-unique(df_auc_common$driver)

clust_order_asc<-names(sort(apply(MinMax_clust, 2, mean)))
MinMax_clust<-MinMax_clust[,clust_order_asc]
MinMax_df <- as.data.frame(MinMax_clust)
MinMax_df$driver <- rownames(MinMax_df)
MinMax_df <- melt(MinMax_df,variable.name = "cluster_id")
MinMax_clust <- MinMax_clust[!rownames(MinMax_clust) %in% TF_2_remove,]


# Plot the heatmap of the regulatory dissimilarity score
p <- pheatmap(MinMax_clust,color=plasma, fontsize=5, angle_col =45, cellwidth=40)
p <- pheatmap(MinMax_clust,color=plasma, angle_col =45, cellwidth=40, cellheight=10)
plot_dims <- get_plot_dims(p)
      
pdf(paste0(plot_root_dir ,file_idx,'Simic_HeatMaps_RegDissScore_CD8.pdf'), height = plot_dims$height, width = plot_dims$width )
print(p)
dev.off()

save.image (file =file.path(PARENT_DIR,"JuanRo_SimiC/RSession/210420_CART_UFILTER.RData"))
#load("/home/sevastopol/data/mcallejac/JuanRo_SimiC/RSession/210420_CART_UFILTER.RData")
#MinMax_clust and df_auc_common

get_ridges <- function(tf_name, dataset){
    plotter <- dataset[, tf_name, drop=FALSE]
    plotter$cell_id <- rownames(plotter)
    plotter <- merge(plotter, ann[, c('cell_id', 'population', 'timePoint')], by = 'cell_id')
    plotter <- setNames(plotter, c('cell_id','value', 'population', 'Phenotype'))
    plotter$population <- factor(plotter$population, levels=c('CD8 Memory', 'CD8 Effector', 'CD8 Exhausted'))
    plotter$Phenotype <- factor(plotter$Phenotype, levels=c('IP', 'D12'))
    tmp <- ggplot(plotter,aes(x=value, y=population , fill = Phenotype)) +
    ggridges::geom_density_ridges(alpha=0.6, adjust=1/2) +
    scale_fill_iwanthue() +
    theme_classic()+ ggtitle(tf_name)+ theme(legend.position = "none", axis.text.y = element_text(angle = 45, vjust = 0.5, hjust=1), axis.title.x = element_blank(), axis.title.y = element_blank())
    return(tmp)
}





all_data <- as.data.frame(seurat_data@assays$RNA@scale.data)
tsne_out <- Rtsne(as.matrix(t(all_data)))
plotter <-  as.data.frame(tsne_out$Y)
rownames(plotter) <- colnames(all_data)
plotter$cell_id <- colnames(all_data)


genes_2_plot <- c('RUNX3', 'MYC')
pdf(paste0(plot_root_dir ,file_idx,'Simic_tSNE_selection.pdf'), height =8, width=18)
for (gene in genes_2_plot){
  plotter_tmp <- plotter[rownames(plotter) %in% df_auc[df_auc$.id %in% phenotypes, 'cell_id'],]
  plotter_tmp <- merge(plotter_tmp, df_auc[df_auc$driver == gene, c('cell_id', 'value', '.id', 'cluster_id')], by='cell_id')
  plotter_tmp$cluster_id <- factor(plotter_tmp$cluster_id, levels=clusters_2_keep)

  ggplot() +
  geom_point(plotter_tmp[,-6], mapping=aes(x=V1, y=V2, shape = .id), color = "grey", alpha = 0.4) + 
  geom_point(plotter_tmp, mapping=aes(x=V1, y=V2, color = value, shape = .id)) + 
  theme_classic() + geom_point(size = 1.5) + 
  scale_color_viridis(option='inferno') + 
  facet_wrap(~cluster_id, ncol=3, nrow=1) + labs(x= '', y = gene )
}

dev.off()
