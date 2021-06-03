#!/usr/bin/env Rscript
#===============================================================================
#' Author: Maria E. Calleja
#' Date: 2021
#' This script is over JuanRo's samples, High-Low signature
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
library(viridis)
library(ggplot2)
library(RColorBrewer)
library(hues)
library(VennDiagram)
options(stringsAsFactors = FALSE)
#===============================================================================
## GLOBAL VARIABLES.
PARENT_DIR<-"/home/mecc/cluster_juanro"
#PARENT_DIR<-"/home/jrrodriguez"
#PARENT_DIR<-"/Users/bioinfo104/cluster_juanro"
PROJECT_DIR<-file.path(PARENT_DIR,"SC_HighLow/data")
#===============================================================================
## FUNCTIONS.
#===============================================================================
## MAIN.
#===============================================================================
# Load signatures and tidy them up. 
CD4_signature_genes <-read.delim(file="./Data/signature/BatchK_CD4_SIGGenes_BASAL_BOTH_LowvsHigh.tsv",sep="\t",header=T)
CD8_signature_genes <-read.delim(file="./Data/signature/BatchK_CD8_SIGGenes_BASAL_BOTH_LowvsHigh.tsv",sep="\t",header=T)

CD4_signature <- as.list(CD4_signature_genes)
CD8_signature <- as.list(CD8_signature_genes)

CD4_signature<-read.delim(file="./Data/signature/BatchK_CD4_output_BASAL_BOTH_LowvsHigh.tsv",sep="\t",header=T)
CD8_signature<-read.delim(file="./Data/signature/BatchK_CD8_output_BASAL_BOTH_LowvsHigh.tsv",sep="\t",header=T)

# write.csv(exprMatrix_cd4, file.path("/home/mecc/Desktop/ForGuille_exprMatrix_cd4.csv"), quote = FALSE, row.names = TRUE)
# write.csv(exprMatrix_cd8, file.path("/home/mecc/Desktop/ForGuille_exprMatrix_cd8.csv"), quote = FALSE, row.names = TRUE)
# write.csv(seurat_integrated_cd4@meta.data, file.path("/home/mecc/Desktop/ForGuille_metadata_cd4.csv"), quote = FALSE, row.names = TRUE)
# write.csv(seurat_integrated_cd8@meta.data, file.path("/home/mecc/Desktop/ForGuille_metadata_cd8.csv"), quote = FALSE, row.names = TRUE)

exprMatrix_cd4<-read.csv(file="./Data/signature/ForGuille_exprMatrix_cd4.csv")
exprMatrix_cd8<-read.csv(file="./Data/signature/ForGuille_exprMatrix_cd8.csv")
cd4_metadata<- read.csv(file="./Data/signature/ForGuille_metadata_cd4.csv", row.names = "X")
cd8_metadata <- read.csv(file="./Data/signature/ForGuille_metadata_cd8.csv", row.names = "X")

res.df <- as.data.frame(CD4_signature)
res.df <- as_tibble(res.df)

res2.cd4 <- res.df %>% 
  dplyr::filter(GeneID %in% CD4_signature_genes$sigGenes_symbol) %>% 
  dplyr::select(GeneID, log2FoldChange, padj) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(GeneID) 
head(res2.cd4)

rm(res.df)
res.df <- as.data.frame(CD8_signature)
res.df <- as_tibble(res.df)

res2.cd8 <- res.df %>% 
  dplyr::filter(GeneID %in% CD8_signature_genes$sigGenes_symbol) %>% 
  dplyr::select(GeneID, log2FoldChange, padj) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(GeneID)
head(res2.cd8)

rm(res.df)

head(res2.cd4)
res2.cd4.highness <- res2.cd4 %>% dplyr::filter(log2FoldChange > 0) 
res2.cd4.lowness <- res2.cd4 %>% dplyr::filter(log2FoldChange < 0) 
res2.cd8.highness <- res2.cd8 %>% dplyr::filter(log2FoldChange > 0) 
res2.cd8.lowness <- res2.cd8 %>% dplyr::filter(log2FoldChange < 0) 

# exprMatrix_cd4 <-as.matrix(seurat_integrated_cd4@assays$RNA@data)
# exprMatrix_cd4 <-as.matrix(exprMatrix_cd4)
dim(exprMatrix_cd4)
head(exprMatrix_cd4[,1:5])

#exprMatrix_cd4 <- rownames_to_column(as.data.frame(exprMatrix_cd4))
colnames(exprMatrix_cd4)[1] <- "GeneID"
head(exprMatrix_cd4[,1:5])
exprMatrix_cd4 <- as_tibble(exprMatrix_cd4)
head(exprMatrix_cd4[, 1:5])
exprMatrix_cd4 <- exprMatrix_cd4 %>% 
  dplyr::filter(GeneID %in% CD4_signature_genes$sigGenes_symbol)

dim(exprMatrix_cd4)
head(exprMatrix_cd4[, 1:6])
head(exprMatrix_cd4[, ncol(exprMatrix_cd4)-4:ncol(exprMatrix_cd4)])

logplusone <- function(x) {log(x + 0.5)}
l <- apply(exprMatrix_cd4[2:ncol(exprMatrix_cd4)], 2, logplusone)
head(l[,1:5])
zcore <- scale(l)
dim(zcore)
head(zcore[,1:5])
zcore <- as.data.frame(zcore)
rownames(zcore) <- exprMatrix_cd4$GeneID
zcore$GeneID <- exprMatrix_cd4$GeneID
# colnames(zcore)[1] <- "GeneID"
zcore <-  as_tibble(zcore)
head(zcore[, 1:5])
ann_markers <- merge(zcore,res2.cd4[, c("GeneID","log2FoldChange", "padj")],by =  "GeneID")


head(ann_markers[, 1:5])
head(ann_markers[, ncol(ann_markers)-4:ncol(ann_markers)])
head(ann_markers[, 12860:12863])
head(ann_markers[, 12862:ncol(ann_markers)])
dim(ann_markers)

DT.high <- ann_markers %>% mutate(across(2:12862, ~ .x * log2FoldChange)) 
DT.low <- ann_markers %>% mutate(across(2:12862, ~ .x * -log2FoldChange)) 

head(DT.high[, 1:5])
head(DT.high[, (ncol(DT.high)-4):ncol(DT.high)])
head(DT.low[, 1:5])
head(DT.low[, (ncol(DT.low)-4):ncol(DT.low)])

DT.cd4.high <- DT.high %>% bind_rows(summarise(.,
                      across(where(is.numeric), sum),
                      across(where(is.character), ~"Total")))
tail(DT.cd4.high[,1:5])
DT.cd4.high %>% tidyr::pivot_longer(!GeneID, names_to= "cells", values_to = "result" )
# DT.cd4.high %>%  
#   dplyr::filter(GeneID  == "Total") %>% 
#   tidyr::pivot_longer(!GeneID, names_to= "cells", values_to = "result" ) %>% 
#   ggplot2::ggplot(ggplot2::aes(x = result)) + ggplot2::geom_histogram (bins=50)

DF.DT.cd4.high <- DT.cd4.high %>%  
  dplyr::filter(GeneID  == "Total") %>% 
  tidyr::pivot_longer(!GeneID, names_to= "cells", values_to = "High_pondered" ) 

head(DF.DT.cd4.high)
# cd4_metadata <- seurat_integrated_cd4@meta.data
cd4_metadata$cells <- gsub("-",".",cd4_metadata$cells)
head(cd4_metadata[,1:15])

cd4_metadata <- inner_join(x = DF.DT.cd4.high [, c("cells", "High_pondered")],
                           y = cd4_metadata[, 1:ncol(cd4_metadata)],
                           by = c("cells" = "cells")) %>% unique()
cd4_metadata <- as.data.frame(cd4_metadata)
head(cd4_metadata[,1:5])
rownames(cd4_metadata) <- cd4_metadata[,1]
head(cd4_metadata[,1:5])
cd4_metadata <- cd4_metadata[,-1]
head(cd4_metadata)
ncol(cd4_metadata)
### 

###

DT.cd4.low <- DT.low %>% bind_rows(summarise(.,
                                               across(where(is.numeric), sum),
                                               across(where(is.character), ~"Total")))
tail(DT.cd4.low[,1:5])
DT.cd4.low %>% tidyr::pivot_longer(!GeneID, names_to= "cells", values_to = "result" )
DT.cd4.low %>%  
  dplyr::filter(GeneID  == "Total") %>% 
  tidyr::pivot_longer(!GeneID, names_to= "cells", values_to = "result" ) %>% 
  ggplot(aes(x = result)) + geom_histogram (bins=50)

DF.DT.cd4.low <- DT.cd4.low %>%  
  dplyr::filter(GeneID  == "Total") %>% 
  tidyr::pivot_longer(!GeneID, names_to= "cells", values_to = "Low_pondered" ) 

dim(DF.DT.cd4.low)
head(cd4_metadata[,1:5])
head(cd4_metadata)
cd4_metadata <- rownames_to_column(as.data.frame(cd4_metadata))

cd4_metadata <- inner_join(x = DF.DT.cd4.low [, c("cells", "Low_pondered")],
                           y = cd4_metadata[, 1:ncol(cd4_metadata)],by = c("cells" = "rowname")) %>% unique()
cd4_metadata <- as.data.frame(cd4_metadata)
head(cd4_metadata[,1:5])
rownames(cd4_metadata) <- cd4_metadata[,1]
head(cd4_metadata[,1:5])
cd4_metadata <- cd4_metadata[,-1]
head(cd4_metadata)
ncol(cd4_metadata)

#####

car_exp <- scale(cd4_metadata$High_pondered)
set.seed(123)
car_p33 <- quantile(car_exp[car_exp>0], probs = c(0.25))
set.seed(123)
car_p66 <- quantile(car_exp[car_exp>0], probs = c(0.75))
quantile(car_exp[car_exp>0], probs = c(0.99))
car_exp.hl <- as.data.frame(car_exp)
car_exp.hl <- mutate(car_exp.hl, CAR_level = ifelse((car_exp>0 & car_exp<=car_p33), 1,
                                                    ifelse((car_exp>car_p33 & car_exp<=car_p66) ,2,
                                                           ifelse( (car_exp > car_p66), 3, 0 ))))
car_exp.hl$CAR_High_level_COD <- plyr::mapvalues(x = car_exp.hl$CAR_level,
                                            from = c(0 , 1 , 2 , 3 ),
                                            to = c("Negative","Low","Med_to_high","High"))
head(car_exp.hl)
table(car_exp.hl$CAR_High_level_COD)
table(car_exp.hl$CAR_level)
summary(car_exp.hl$V1)

cd4_metadata$CAR_HIGH_SCORE_LEVEL_HL <- car_exp.hl$CAR_level
cd4_metadata$CAR_HIGH_SCORE_COD <- car_exp.hl$CAR_High_level_COD
cd4_metadata$CAR_HIGH_SCORE_pondered_scaled <- car_exp.hl$V1
head(cd4_metadata)
#####

#####

car_exp <- scale(cd4_metadata$Low_pondered)
set.seed(123)
car_p33 <- quantile(car_exp[car_exp>0], probs = c(0.25))
set.seed(123)
car_p66 <- quantile(car_exp[car_exp>0], probs = c(0.75))
quantile(car_exp[car_exp>0], probs = c(0.99))
car_exp.hl <- as.data.frame(car_exp)
car_exp.hl <- mutate(car_exp.hl, CAR_level = ifelse((car_exp>0 & car_exp<=car_p33), 1,
                                                    ifelse((car_exp>car_p33 & car_exp<=car_p66) ,2,
                                                           ifelse( (car_exp > car_p66), 3, 0 ))))
car_exp.hl$CAR_High_level_COD <- plyr::mapvalues(x = car_exp.hl$CAR_level,
                                                 from = c(0 , 1 , 2 , 3 ),
                                                 to = c("Negative","Low","Med_to_high","High"))
head(car_exp.hl)
table(car_exp.hl$CAR_High_level_COD)
table(car_exp.hl$CAR_level)
summary(car_exp.hl$V1)

cd4_metadata$CAR_LOW_SCORE_LEVEL_HL <- car_exp.hl$CAR_level
cd4_metadata$CAR_LOW_SCORE_COD <- car_exp.hl$CAR_High_level_COD
cd4_metadata$CAR_LOW_SCORE_pondered_scaled <- car_exp.hl$V1
#####
head(cd4_metadata)



#####

colnames(exprMatrix_cd8)[1] <- "GeneID"
head(exprMatrix_cd8[,1:5])
exprMatrix_cd8 <- as_tibble(exprMatrix_cd8)
head(exprMatrix_cd8[, 1:5])
exprMatrix_cd8 <- exprMatrix_cd8 %>% 
  dplyr::filter(GeneID %in% CD8_signature_genes$sigGenes_symbol)

dim(exprMatrix_cd8)
logplusone <- function(x) {log(x + 0.5)}
l <- apply(exprMatrix_cd8[2:ncol(exprMatrix_cd8)], 2, logplusone)
head(l[,1:5])
zcore <- scale(l)
dim(zcore)
head(zcore[,1:5])

rownames(zcore) <- exprMatrix_cd8$GeneID
zcore <- rownames_to_column(as.data.frame(zcore))
colnames(zcore)[1] <- "GeneID"
zcore <-  as_tibble(zcore)
head(zcore[, 1:5])
ann_markers <- inner_join(x = zcore,
                          y = res2.cd8[, c("GeneID","log2FoldChange", "padj")],by = c("GeneID" = "GeneID")) %>% unique()

head(ann_markers[, 1:5])
head(ann_markers[, ncol(ann_markers)-4:ncol(ann_markers)])
head(ann_markers[, 4140:ncol(ann_markers)])
dim(ann_markers)

DT.high <- ann_markers %>% mutate(across(2:4140, ~ .x * log2FoldChange)) 
DT.low <- ann_markers %>% mutate(across(2:4140, ~ .x * -log2FoldChange)) 

head(DT.high[, 1:5])
head(DT.high[, (ncol(DT.high)-4):ncol(DT.high)])
head(DT.low[, 1:5])
head(DT.low[, (ncol(DT.low)-4):ncol(DT.low)])

DT.cd8.high <- DT.high %>% bind_rows(summarise(.,
                                               across(where(is.numeric), sum),
                                               across(where(is.character), ~"Total")))
tail(DT.cd8.high[,1:5])
DT.cd8.high %>% pivot_longer(!GeneID, names_to= "cells", values_to = "result" )
DT.cd8.high %>%  
  dplyr::filter(GeneID  == "Total") %>% 
  pivot_longer(!GeneID, names_to= "cells", values_to = "result" ) %>% 
  ggplot(aes(x = result)) + geom_histogram (bins=50)

DF.DT.cd8.high <- DT.cd8.high %>%  
  dplyr::filter(GeneID  == "Total") %>% 
  pivot_longer(!GeneID, names_to= "cells", values_to = "High pondered" ) 

dim(DF.DT.cd8.high)
cd8_metadata$cells <- gsub("-",".",cd8_metadata$cells)
head(cd8_metadata[,1:15])

cd8_metadata <- inner_join(x = DF.DT.cd8.high [, c("cells", "High pondered")],
                           y = cd8_metadata[, 1:ncol(cd8_metadata)],
                           by = c("cells" = "cells")) %>% unique()
cd8_metadata <- as.data.frame(cd8_metadata)
head(cd8_metadata[,1:5])
rownames(cd8_metadata) <- cd8_metadata[,1]
head(cd8_metadata[,1:5])
cd8_metadata <- cd8_metadata[,-1]
head(cd8_metadata)
ncol(cd8_metadata)
###

###
DT.cd8.low <- DT.low %>% bind_rows(summarise(.,
                                             across(where(is.numeric), sum),
                                             across(where(is.character), ~"Total")))
tail(DT.cd8.low[,1:5])
DT.cd8.low %>% pivot_longer(!GeneID, names_to= "cells", values_to = "result" )
DT.cd8.low %>%  
  dplyr::filter(GeneID  == "Total") %>% 
  pivot_longer(!GeneID, names_to= "cells", values_to = "result" ) %>% 
  ggplot(aes(x = result)) + geom_histogram (bins=50)

DF.DT.cd8.low <- DT.cd8.low %>%  
  dplyr::filter(GeneID  == "Total") %>% 
  pivot_longer(!GeneID, names_to= "cells", values_to = "Low pondered" ) 

dim(DF.DT.cd8.low)

head(cd8_metadata[,1:5])
head(cd8_metadata)
cd8_metadata <- rownames_to_column(as.data.frame(cd8_metadata))

cd8_metadata <- inner_join(x = DF.DT.cd8.low [, c("cells", "Low pondered")],
                           y = cd8_metadata[, 1:ncol(cd8_metadata)],by = c("cells" = "rowname")) %>% unique()
cd8_metadata <- as.data.frame(cd8_metadata)
head(cd8_metadata[,1:5])
rownames(cd8_metadata) <- cd8_metadata[,1]
head(cd8_metadata[,1:5])
cd8_metadata <- cd8_metadata[,-1]
head(cd8_metadata)
ncol(cd8_metadata)
############

#####
car_exp <- scale(cd8_metadata$`High pondered`)
set.seed(123)
car_p33 <- quantile(car_exp[car_exp>0], probs = c(0.25))
set.seed(123)
car_p66 <- quantile(car_exp[car_exp>0], probs = c(0.75))
quantile(car_exp[car_exp>0], probs = c(0.99))
car_exp.hl <- as.data.frame(car_exp)
car_exp.hl <- mutate(car_exp.hl, CAR_level = ifelse((car_exp>0 & car_exp<=car_p33), 1,
                                                    ifelse((car_exp>car_p33 & car_exp<=car_p66) ,2,
                                                           ifelse( (car_exp > car_p66), 3, 0 ))))
car_exp.hl$CAR_High_level_COD <- plyr::mapvalues(x = car_exp.hl$CAR_level,
                                                 from = c(0 , 1 , 2 , 3 ),
                                                 to = c("Negative","Low","Med_to_high","High"))
head(car_exp.hl)
table(car_exp.hl$CAR_High_level_COD)
table(car_exp.hl$CAR_level)
summary(car_exp.hl$V1)

cd8_metadata$CAR_HIGH_SCORE_LEVEL_HL <- car_exp.hl$CAR_level
cd8_metadata$CAR_HIGH_SCORE_COD<- car_exp.hl$CAR_High_level_COD 
cd8_metadata$CAR_HIGH_SCORE_pondered_scaled<- car_exp.hl$V1

##
car_exp <- scale(cd8_metadata$`Low pondered`)
set.seed(123)
car_p33 <- quantile(car_exp[car_exp>0], probs = c(0.25))
set.seed(123)
car_p66 <- quantile(car_exp[car_exp>0], probs = c(0.75))
quantile(car_exp[car_exp>0], probs = c(0.99))
car_exp.hl <- as.data.frame(car_exp)
car_exp.hl <- mutate(car_exp.hl, CAR_level = ifelse((car_exp>0 & car_exp<=car_p33), 1,
                                                    ifelse((car_exp>car_p33 & car_exp<=car_p66) ,2,
                                                           ifelse( (car_exp > car_p66), 3, 0 ))))
car_exp.hl$CAR_High_level_COD <- plyr::mapvalues(x = car_exp.hl$CAR_level,
                                                 from = c(0 , 1 , 2 , 3 ),
                                                 to = c("Negative","Low","Med_to_high","High"))
head(car_exp.hl)
table(car_exp.hl$CAR_High_level_COD)
table(car_exp.hl$CAR_level)
summary(car_exp.hl$V1)

cd8_metadata$CAR_LOW_SCORE_LEVEL_HL <- car_exp.hl$CAR_level
cd8_metadata$CAR_LOW_SCORE_COD<- car_exp.hl$CAR_High_level_COD 
cd8_metadata$CAR_LOW_SCORE_pondered_scaled<- car_exp.hl$V1

################################################################################
## Saving and etcs 
################################################################################
save.image ("ANNOTATION_FORGUILLE.RData")
savehistory("ANNOTATION_FORGUILLE.Rhistory")
sink("ANNOTATION_FORGUILLE.txt")
toLatex(sessionInfo())
sink(NULL)

