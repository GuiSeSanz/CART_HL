
# CD4_High
cat /home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/PeakCalling_MACS/D10_High_CD4_d0_peaks.broadPeak /home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/PeakCalling_MACS/D14_High_CD4_d0_peaks.broadPeak /home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/PeakCalling_MACS/D18_High_CD4_d0_peaks.broadPeak > ConsensusCD4_High_d0.broadPeak
bedtools sort  -i ConsensusCD4_High_d0.broadPeak > ConsensusCD4_High_d0_sorted.broadPeak
bedtools merge -i ConsensusCD4_High_d0_sorted.broadPeak > ConsensusCD4_High_d0_sorted_merged.broadPeak

# CD4_Low
cat /home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/PeakCalling_MACS/D10_Low_CD4_d0_peaks.broadPeak /home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/PeakCalling_MACS/D14_Low_CD4_d0_peaks.broadPeak /home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/PeakCalling_MACS/D18_Low_CD4_d0_peaks.broadPeak > ConsensusCD4_Low_d0.broadPeak
bedtools sort  -i ConsensusCD4_Low_d0.broadPeak > ConsensusCD4_Low_d0_sorted.broadPeak
bedtools merge -i ConsensusCD4_Low_d0_sorted.broadPeak > ConsensusCD4_Low_d0_sorted_merged.broadPeak


# CD8_High
cat /home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/PeakCalling_MACS/D10_High_CD8_d0_peaks.broadPeak /home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/PeakCalling_MACS/D14_High_CD8_d0_peaks.broadPeak /home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/PeakCalling_MACS/D18_High_CD8_d0_peaks.broadPeak > ConsensusCD8_High_d0.broadPeak
bedtools sort  -i ConsensusCD8_High_d0.broadPeak > ConsensusCD8_High_d0_sorted.broadPeak
bedtools merge -i ConsensusCD8_High_d0_sorted.broadPeak > ConsensusCD8_High_d0_sorted_merged.broadPeak

# CD8_Low
cat /home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/PeakCalling_MACS/D10_Low_CD8_d0_peaks.broadPeak /home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/PeakCalling_MACS/D14_Low_CD8_d0_peaks.broadPeak /home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/PeakCalling_MACS/D18_Low_CD8_d0_peaks.broadPeak > ConsensusCD8_Low_d0.broadPeak
bedtools sort  -i ConsensusCD8_Low_d0.broadPeak > ConsensusCD8_Low_d0_sorted.broadPeak
bedtools merge -i ConsensusCD8_Low_d0_sorted.broadPeak > ConsensusCD8_Low_d0_sorted_merged.broadPeak







library(ggplot2)
library(edgeR)
library(DESeq2)

load('/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/RSession/cd8_2105_Csaw_fin.RData')
dds <- edgeR::getCounts(y)
dds_norm <- edgeR::cpm(y,  normalized.lib.sizes = TRUE)
samples_2_keep <- grep('High|Low', grep('d0',colnames(dds), value = TRUE), value = TRUE)
dds <- dds[, samples_2_keep]
pca <- prcomp(t(dds))
print(summary(pca))

pcaData = as.data.frame(pca$x)
pcaData$sample=rownames(pcaData)
pcaData$HighLow <- stringr::str_extract(pcaData$sample, '(?<=_)[A-Za-z]+(?=_)')
pcaData$Donor <- stringr::str_extract(pcaData$sample, '^[A-Z0-9]+')
percentVar = round(100 * (pca$sdev^2 / sum( pca$sdev^2 ) ))

dds_obj <- DESeqDataSetFromMatrix( countData = dds , colData = as.data.frame(colnames(dds)), design = design)
rld <- rlog(dds_obj, blind = T)
head(assay(rld),3)
colData(rld)
#sample distances
sampleDists <- dist(t(assay(rld)))
sampleDists
sampleDistMatrix <- as.matrix( sampleDists )

pdf('./Plots/ATAC_pca_dist.pdf')
ggplot(data=pcaData, aes(x = PC1, y = PC2, color=sample)) + geom_point(size=3)  + xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance")) + theme_classic() + hues::scale_color_iwanthue()
ggplot(data=pcaData, aes(x = PC1, y = PC2, color=Donor)) + geom_point(size=3)   + xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance")) + theme_classic() + hues::scale_color_iwanthue()
ggplot(data=pcaData, aes(x = PC1, y = PC2, color=HighLow)) + geom_point(size=3) + xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance")) + theme_classic() +  scale_color_manual(values=c('High'='#30A3CC', 'Low'='#FCB357'))
# pheatmap::pheatmap(sampleDistMatrix,
#          clustering_distance_rows = sampleDists,
#          clustering_distance_cols = sampleDists,
#          col = viridis::viridis(12, direction=-1))
dev.off()


#with limma
counts <- edgeR::getCounts(y)
samples_2_keep <- grep('High|Low', grep('d0',colnames(counts), value = TRUE), value = TRUE)
counts <- counts[, samples_2_keep]

group <- as.data.frame(factor(stringr::str_extract(colnames(counts), '^D[0-9a]+')))
rownames(group) <- colnames(counts)
HighLow <- factor(stringr::str_extract(colnames(counts), '(?=_)[\\w]+(?<=_)'))
group$HighLow <- HighLow
colnames(group) <- c('Donor', 'HighLow')
coldata <- group 

dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ HighLow)
vsd <- vst(dds, blind=FALSE)
mat <- assay(vsd)
mm <- model.matrix(~HighLow, colData(vsd))
mat <- limma::removeBatchEffect(mat, batch=vsd$batch, design=mm)
assay(vsd) <- mat
pdf('./Plots/ATAC_pca_limma.pdf')
plotPCA(vsd, intgroup = 'HighLow')
dev.off()



# library(edgeR)
counts <- edgeR::getCounts(y)
cols_2_keep <- grep('High|Low', grep('d0', colnames(y), value=TRUE), value=TRUE)
counts <- counts[, cols_2_keep]
# group <- as.data.frame(factor(stringr::str_extract(colnames(counts), '^D[0-9a]+')))
group <-factor(stringr::str_extract(colnames(counts), '^D[0-9a]+'))
# rownames(group) <- colnames(counts)
# HighLow <- factor(stringr::str_extract(colnames(counts), '(?=_)[\\w]+(?<=_)'))
# group$HighLow <- HighLow
# colnames(group) <- c('group', 'HighLow')
dge <- DGEList(counts=counts, group=group)

dge <- calcNormFactors(dge, method = "TMM")

# model.matrix(~group$High_low)
logCPM <- cpm(dge, log=TRUE, prior.count=2)
logCPM_bc <- removeBatchEffect(logCPM,  group=group)
# plotMDS(logCPM_bc)
pca <- prcomp(t(logCPM_bc))
print(summary(pca))

pcaData = as.data.frame(pca$x)
pcaData$sample=rownames(pcaData)
pcaData$HighLow <- stringr::str_extract(pcaData$sample, '(?<=_)[A-Za-z]+(?=_)')
pcaData$Donor <- stringr::str_extract(pcaData$sample, '^[A-Z0-9]+')
percentVar = round(100 * (pca$sdev^2 / sum( pca$sdev^2 ) ))



# pdf('./Plots/ATAC_pca_dist_BACTHRemoved.pdf')
# ggplot(data=pcaData, aes(x = PC1, y = PC2, color=sample)) + geom_point(size=3) + xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance")) + theme_classic() + hues::scale_color_iwanthue()
# ggplot(data=pcaData, aes(x = PC1, y = PC2, color=Donor)) + geom_point(size=3) + xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance")) + theme_classic() + hues::scale_color_iwanthue()
# ggplot(data=pcaData, aes(x = PC1, y = PC2, color=HighLow)) + geom_point(size=3) + xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance")) + theme_classic() +  scale_color_manual(values=c('High'='#30A3CC', 'Low'='#FCB357'))
# dev.off()


combat_edata1 = sva::ComBat(dat=logCPM, batch=group, mod=NULL, par.prior=TRUE, prior.plots=FALSE)
pca <- prcomp(t(combat_edata1))
print(summary(pca))
pcaData = as.data.frame(pca$x)
pcaData$sample=rownames(pcaData)
pcaData$HighLow <- stringr::str_extract(pcaData$sample, '(?<=_)[A-Za-z]+(?=_)')
pcaData$Donor <- stringr::str_extract(pcaData$sample, '^[A-Z0-9]+')
percentVar = round(100 * (pca$sdev^2 / sum( pca$sdev^2 ) ))
pdf('./Plots/ATAC_pca_dist_BACTHRemovedCombat.pdf')
    ggplot(data=pcaData, aes(x = PC1, y = PC2, color=sample)) + geom_point(size=3) + xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance")) + theme_classic() + hues::scale_color_iwanthue()
    ggplot(data=pcaData, aes(x = PC1, y = PC2, color=Donor)) + geom_point(size=3) + xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance")) + theme_classic() + hues::scale_color_iwanthue()
    ggplot(data=pcaData, aes(x = PC1, y = PC2, color=HighLow)) + geom_point(size=3) + xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance")) + theme_classic() +  scale_color_manual(values=c('High'='#30A3CC', 'Low'='#FCB357'))
dev.off()


# RNASEQ_COMBAT 
counts_norm <- read.table('/home/sevastopol/data/gserranos/CART_HL/Data/signature/BatchK_CD8_output_BASAL_BOTH_LowvsHigh.tsv', sep='\t', header=TRUE)
counts_norm <- counts_norm[,grepl('^D', colnames(counts_norm))]
group <-factor(stringr::str_extract(colnames(counts_norm), '^D[0-9a]+'))

combat_edata2 = sva::ComBat(dat=as.matrix(counts_norm), batch=group, mod=NULL, par.prior=TRUE, prior.plots=FALSE)
pca <- prcomp(t(combat_edata2))
print(summary(pca))
pcaData = as.data.frame(pca$x)
pcaData$sample=rownames(pcaData)
pcaData$HighLow <- stringr::str_extract(pcaData$sample, '(?<=_)[A-Za-z]+(?=_)')
pcaData$Donor <- stringr::str_extract(pcaData$sample, '^[A-Z0-9]+')
percentVar = round(100 * (pca$sdev^2 / sum( pca$sdev^2 ) ))
pdf('./Plots/BulkRNA_pca_dist_BACTHRemovedCombat.pdf')
    ggplot(data=pcaData, aes(x = PC1, y = PC2, color=sample)) + geom_point(size=3) + xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance")) + theme_classic() + hues::scale_color_iwanthue()
    ggplot(data=pcaData, aes(x = PC1, y = PC2, color=Donor)) + geom_point(size=3) + xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance")) + theme_classic() + hues::scale_color_iwanthue()
    ggplot(data=pcaData, aes(x = PC1, y = PC2, color=HighLow)) + geom_point(size=3) + xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance")) + theme_classic() +  scale_color_manual(values=c('High'='#30A3CC', 'Low'='#FCB357'))
dev.off()





HighLow <- factor(stringr::str_extract(colnames(counts), '(?=_)[\\w]+(?<=_)'))
design_br <- cbind(group, HighLow)
rownames(design_br) <- colnames(counts)
br_data <- removeBatchEffect(logCPM, group, design=design_br)
pca <- prcomp(t(br_data))
print(summary(pca))
pcaData = as.data.frame(pca$x)
pcaData$sample=rownames(pcaData)
pcaData$HighLow <- stringr::str_extract(pcaData$sample, '(?<=_)[A-Za-z]+(?=_)')
pcaData$Donor <- stringr::str_extract(pcaData$sample, '^[A-Z0-9]+')
percentVar = round(100 * (pca$sdev^2 / sum( pca$sdev^2 ) ))
pdf('./Plots/ATAC_pca_dist_BACTHRemovedAnotherDesign.pdf')
    ggplot(data=pcaData, aes(x = PC1, y = PC2, color=sample)) + geom_point(size=3) + xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance")) + theme_classic() + hues::scale_color_iwanthue()
    ggplot(data=pcaData, aes(x = PC1, y = PC2, color=Donor)) + geom_point(size=3) + xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance")) + theme_classic() + hues::scale_color_iwanthue()
    ggplot(data=pcaData, aes(x = PC1, y = PC2, color=HighLow)) + geom_point(size=3) + xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance")) + theme_classic() +  scale_color_manual(values=c('High'='#30A3CC', 'Low'='#FCB357'))
dev.off()


library(DESeq2)
raw_counts <- edgeR::getCounts(y)

cell_cycle <- c("G2","G1","G2",rep("G1",3),rep("S",3),rep("G2",2),"G1","G2","G1",rep("S",3),"G2")
treatment <- c(rep("purified",3),rep("total",8),"purified","total",rep("purified",5))

coldata <- as.data.frame(cbind(colnames(counts),cell_cycle,treatment))
colnames(coldata)=c("sample","cell_cycle","treatment")

group <- as.data.frame(factor(stringr::str_extract(colnames(counts), '^D[0-9a]+')))
rownames(group) <- colnames(counts)
HighLow <- factor(stringr::str_extract(colnames(counts), '(?=_)[\\w]+(?<=_)'))
group$HighLow <- HighLow
colnames(group) <- c('Donor', 'HighLow')
coldata <- group 
dds <- DESeqDataSetFromMatrix( countData = counts , colData = coldata,
                               design = ~ Donor + HighLow)
#assign the controls 
dds$cell_cycle <- relevel (dds$cell_cycle , "G1")
dds$treatment <- relevel (dds$treatment , "total")








library(csaw)
####### CD4

MACS_PATH = '/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/PeakCalling_MACS/'
list_MACS_files <- grep('High|Low', list.files(MACS_PATH, pattern ='.broadPeak'), value=TRUE)
list_MACS_files <- grep('_d0', list_MACS_files, value=TRUE)
list_MACS_files <- grep('CD4', list_MACS_files, value=TRUE)
list_MACS_files <- paste0(MACS_PATH, list_MACS_files)

peaks <- lapply(list_MACS_files, function(x) {read.table(x, sep = "\t")[,1:3]})
peak_names <- stringr::str_extract(list_MACS_files, '(?=D)[\\w]+')
names(peaks) <- peak_names



new_colnames <- c("chrom", "start", "end")
peaks <- lapply(peaks, setNames, nm = new_colnames)
peaks <- lapply(peaks, GRanges)




High.d0   <- Reduce(intersect, peaks [grepl("High",   names(peaks)) & grepl("d0", names(peaks))])
Low.d0    <- Reduce(intersect, peaks [grepl("Low",    names(peaks)) & grepl("d0", names(peaks))])

lapply(peaks, length)
all.peaks <- union(High.d0, Low.d0)
length(all.peaks)

BAMS_MACS_PATH = '/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/BAMS_MACS/'

bamFiles_CD4_MACS <- grep('High|Low', list.files(path=file.path(BAMS_MACS_PATH),pattern = "*.nsorted.fixed.bam$"), value=TRUE)
bamFiles_CD4_MACS <- grep('_d0', bamFiles_CD4_MACS, value=TRUE)                    
bamFiles_CD4_MACS <- grep('CD4', bamFiles_CD4_MACS, value=TRUE)  

pe.bams <- paste0(BAMS_MACS_PATH, bamFiles_CD4_MACS)

# define read parameters
standard.chr <- paste0("chr", c(1:22)) # only use standard chromosomes
##param <- readParam(max.frag=1000, pe="both", discard=blacklist, restrict=standard.chr)
param <- readParam(max.frag=1000, pe="both", restrict=standard.chr)

##############################
# count reads in windows specified by MACS2                                      
peak.counts <- regionCounts(pe.bams, all.peaks, param=param)

library(edgeR)
peak.abundances <- aveLogCPM(asDGEList(peak.counts)) 
peak.counts.filt <- peak.counts[peak.abundances > -3, ] # only use peaks logCPM > -3

# count BAM reads in, e.g. 300 bp windows
counts <- windowCounts(pe.bams, width=300, param=param) # set width as desired from the fragment length distribution analyses

# filter uninteresting features (windows) by local enrichment
# local background estimator: 2kb neighborhood
neighbor <- suppressWarnings(resize(rowRanges(counts), width=2000, fix="center")) # change width parameter as desired
wider <- regionCounts(pe.bams, regions=neighbor, param=param) # count reads in neighborhoods

filter.stat <- filterWindows(counts, wider, type="local") 
counts.local.filt <- counts[filter.stat$filter > log2(3),] # threshold of 3-fold increase in enrichment over 2kb neighborhood abundance; change as desired

# csaw de novo peaks by local enrichment, csaw loess-normalization
counts.local.loess <- counts.local.filt
counts.local.loess <- normOffsets(counts.local.loess, se.out=TRUE)

working.windows <- counts.local.loess # csaw de novo peaks by local enrichment, for trended biases

# setup design matrix
# see edgeR manual for more information
y <- asDGEList(working.windows)

colnames(y$counts)  <- names(peaks)
rownames(y$samples) <- names(peaks)
y$samples$group <- stringr::str_extract





# Deeptools
# From the R sesion 
# write.table(data.frame(all.peaks), file="/home/sevastopol/data/gserranos/CART_HL/Data/High_low/ATACseq/Test_Guille/DeepTools/All_peaks.bed", quote=F, sep="\t", row.names=F, col.names=F)

#  ++++HIGH++++ 
# CD4
/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/BigWig/D10_High_CD4_d0.sort.rmdup.rmblackls.rmchr.norm.bw
/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/BigWig/D14_High_CD4_d0.sort.rmdup.rmblackls.rmchr.norm.bw
/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/BigWig/D18_High_CD4_d0.sort.rmdup.rmblackls.rmchr.norm.bw
# CD8
/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/BigWig/D10_High_CD8_d0.sort.rmdup.rmblackls.rmchr.norm.bw
/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/BigWig/D14_High_CD8_d0.sort.rmdup.rmblackls.rmchr.norm.bw
/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/BigWig/D18_High_CD8_d0.sort.rmdup.rmblackls.rmchr.norm.bw
#  ++++LOW++++ 
# CD4
/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/BigWig/D10_Low_CD4_d0.sort.rmdup.rmblackls.rmchr.norm.bw
/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/BigWig/D14_Low_CD4_d0.sort.rmdup.rmblackls.rmchr.norm.bw
/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/BigWig/D18_Low_CD4_d0.sort.rmdup.rmblackls.rmchr.norm.bw
# CD8
/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/BigWig/D10_Low_CD8_d0.sort.rmdup.rmblackls.rmchr.norm.bw
/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/BigWig/D14_Low_CD8_d0.sort.rmdup.rmblackls.rmchr.norm.bw
/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/BigWig/D18_Low_CD8_d0.sort.rmdup.rmblackls.rmchr.norm.bw



 computeMatrix scale-regions -S /home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/BigWig/D10_High_CD4_d0.sort.rmdup.rmblackls.rmchr.norm.bw /home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/BigWig/D14_High_CD4_d0.sort.rmdup.rmblackls.rmchr.norm.bw /home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/BigWig/D18_High_CD4_d0.sort.rmdup.rmblackls.rmchr.norm.bw -R All_peaks.bed  --beforeRegionStartLength 3000  --regionBodyLength 5000  --afterRegionStartLength 3000 --skipZeros -o CD4_HIGH.mat.gz

 computeMatrix scale-regions -S /home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/BigWig/D10_High_CD8_d0.sort.rmdup.rmblackls.rmchr.norm.bw /home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/BigWig/D10_High_CD8_d0.sort.rmdup.rmblackls.rmchr.norm.bw /home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/BigWig/D14_High_CD8_d0.sort.rmdup.rmblackls.rmchr.norm.bw /home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/BigWig/D18_High_CD8_d0.sort.rmdup.rmblackls.rmchr.norm.bw -R All_peaks.bed  --beforeRegionStartLength 3000  --regionBodyLength 5000  --afterRegionStartLength 3000 --skipZeros -o CD8_HIGH.mat.gz


MACS_PATH = '/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/PeakCalling_MACS/'
list_MACS_files <- grep('High|Low', list.files(MACS_PATH, pattern ='.broadPeak'), value=TRUE)
list_MACS_files <- grep('_d0', list_MACS_files, value=TRUE)
# list_MACS_files <- grep('CD4', list_MACS_files, value=TRUE)
list_MACS_files <- paste0(MACS_PATH, list_MACS_files)

peaks <- lapply(list_MACS_files, function(x) {read.table(x, sep = "\t")[,1:3]})
peak_names <- stringr::str_extract(list_MACS_files, '(?=D)[\\w]+')
names(peaks) <- peak_names



new_colnames <- c("chrom", "start", "end")
peaks <- lapply(peaks, setNames, nm = new_colnames)
peaks <- lapply(peaks, GRanges)

# Specific independent annotation
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(EnsDb.Hsapiens.v86)
edb <- EnsDb.Hsapiens.v86
annotation_Peaks <- TxDb.Hsapiens.UCSC.hg38.knownGene

High.CD8 <- Reduce(intersect, peaks [grepl("High",   names(peaks)) & grepl("CD8", names(peaks))])
High.CD4 <- Reduce(intersect, peaks [grepl("High",   names(peaks)) & grepl("CD4", names(peaks))])
Low.CD8  <- Reduce(intersect, peaks [grepl("Low",    names(peaks)) & grepl("CD8", names(peaks))])
Low.CD4  <- Reduce(intersect, peaks [grepl("Low",    names(peaks)) & grepl("CD4", names(peaks))])



peakAnno_High_CD8 <- annotatePeak(High.CD8, tssRegion=c(-3000, 3000),TxDb=annotation_Peaks, annoDb="org.Hs.eg.db")
peakAnno_High_CD4 <- annotatePeak(High.CD4, tssRegion=c(-3000, 3000),TxDb=annotation_Peaks, annoDb="org.Hs.eg.db")
peakAnno_Low_CD8  <- annotatePeak(Low.CD8, tssRegion=c(-3000, 3000),TxDb=annotation_Peaks, annoDb="org.Hs.eg.db")
peakAnno_Low_CD4  <- annotatePeak(Low.CD4, tssRegion=c(-3000, 3000),TxDb=annotation_Peaks, annoDb="org.Hs.eg.db")

promoter <- getPromoters(TxDb=annotation_Peaks, upstream=3000, downstream=3000)

High.CD8_Mat <- plotAvgProf(getTagMatrix(High.CD8, windows=promoter), xlim=c(-3000, 3000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
High.CD4_Mat <- plotAvgProf(getTagMatrix(High.CD4, windows=promoter), xlim=c(-3000, 3000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
Low.CD8_Mat  <- plotAvgProf(getTagMatrix(Low.CD8, windows=promoter), xlim=c(-3000, 3000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
Low.CD4_Mat  <- plotAvgProf(getTagMatrix(Low.CD4, windows=promoter), xlim=c(-3000, 3000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")

High.CD8_Mat <- High.CD8_Mat$data
High.CD8_Mat$HighLow <- 'High'
Low.CD8_Mat <- Low.CD8_Mat$data
Low.CD8_Mat$HighLow <- 'Low'
plotter_CD8 <- rbind(High.CD8_Mat, Low.CD8_Mat)
High.CD4_Mat <- High.CD4_Mat$data
High.CD4_Mat$HighLow <- 'High'
Low.CD4_Mat  <- Low.CD4_Mat$data
Low.CD4_Mat$HighLow <- 'Low'
plotter_CD4<- rbind(High.CD4_Mat, Low.CD4_Mat)

#Consensus
pdf('./Plots/AverageProfile.pdf')
ggplot(plotter_CD8, aes(x=pos, y=value, color=HighLow)) + geom_line() + 
xlim(-3000, 3000) + ggtitle('CD8') +xlab("Genomic Region (5'->3')")+ ylab("Read Count Frequency") + theme_classic() + scale_color_manual(values=c('High'='#30A3CC', 'Low'='#FCB357'))

ggplot(plotter_CD4, aes(x=pos, y=value, color=HighLow)) + geom_line() + 
xlim(-3000, 3000)+ ggtitle('CD4') + xlab("Genomic Region (5'->3')")+ ylab("Read Count Frequency") + theme_classic() + scale_color_manual(values=c('High'='#30A3CC', 'Low'='#FCB357'))
dev.off()





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
    theme(legend.position="right", text=element_text(family='serif'), legend.key.size = unit(0.2, 'cm')) +
    scale_fill_manual(values = c(pie_colors), name="Feature")
    return(gg)
}

pdf('./Plots/GenAnnotation.pdf', width =8.3, height = 5.8)
legend <- cowplot::get_legend(get_pie_Annot(peakAnno_High_CD8) + theme(legend.position='bottom', legend.key.size = unit(0.3, 'cm')))
cowplot::plot_grid(
    get_pie_Annot(peakAnno_High_CD8) ,
    get_pie_Annot(peakAnno_High_CD4) ,
    get_pie_Annot(peakAnno_Low_CD8 ) ,
    get_pie_Annot(peakAnno_Low_CD4 ) , ncol=2 ,
    labels = c('High CD8', 'High CD4', 'Low CD8', 'Low CD4')

)

dev.off()


peaks_CD4 <- read.table('/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/results/Lowd0_vs_Highd0_csaw_denovo_trended-windows_allcounts_Annotated.txt', sep="\t", header=T, fill=TRUE)
peaks_CD8 <- read.table('/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/results/cd8_Lowd0_vs_Highd0_csaw_denovo_trended-windows_allcounts_Annotated.txt', sep="\t", header=T, fill=TRUE)
peaks_CD4 <- makeGRangesFromDataFrame(peaks_CD4)
peaks_CD8 <- makeGRangesFromDataFrame(peaks_CD8)

annotation_Peaks <- TxDb.Hsapiens.UCSC.hg38.knownGene
peaks_CD4 <- annotatePeak(peaks_CD4, tssRegion=c(-3000, 3000),TxDb=annotation_Peaks, annoDb="org.Hs.eg.db")
peaks_CD8 <- annotatePeak(peaks_CD8, tssRegion=c(-3000, 3000),TxDb=annotation_Peaks, annoDb="org.Hs.eg.db")


pdf('./Plots/GenAnnotationDE_High_Low.pdf', width =8.3, height = 5.8)
cowplot::plot_grid(
    get_pie_Annot(peaks_CD4) ,
    get_pie_Annot(peaks_CD8 ), 
ncol =2, labels =c('CD4', 'CD8'))

dev.off()






source('trackplot.R')








library(Gviz)
library(rtracklayer)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
require(tidyverse)
require(biomaRt)
require(org.Hs.eg.db)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene 
gtrack <- GenomeAxisTrack()

chr_no <- "9" # chromosome number
chr_start <- 70119695 # start of region
chr_end <- 70660299 # end of region

# gtTrack <- GeneRegionTrack(txdb,
#                            chromosome=chr_no, # chromosome number
#                            start=chr_start, # start of region
#                            end=chr_end, # end of region
#                            transcriptAnnotation="symbol", # symbol is the gene symbol
#                            fontsize.group=20, # free to adjust font size
#                            stacking="dense"
#                            )

itrack <- IdeogramTrack(genome="hg38", 
                       chromosome=paste0("chr",chr_no), # specify chromosome in ucsc naming
                       from =chr_start,
                       to=chr_end)
 
itrack@chromosome <- chr_no

# remove chr from chromosome naming
levels(itrack@bandTable$chrom) <- sub("^chr", "", levels(itrack@bandTable$chrom), ignore.case=T)
 
bigwigs <- c('/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/BigWig/D10_High_CD8_d0.sort.rmdup.rmblackls.rmchr.norm.bw', '/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/BigWig/D10_Low_CD8_d0.sort.rmdup.rmblackls.rmchr.norm.bw')
bw_high <- import.bw('/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/BigWig/D10_High_CD8_d0.sort.rmdup.rmblackls.rmchr.norm.bw', as="GRanges")
bw_low <- import.bw('/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/BigWig/D10_Low_CD8_d0.sort.rmdup.rmblackls.rmchr.norm.bw', as="GRanges")

# change chr name to without chr (this depends on the data)
bw_high@seqnames@values <- as.factor(stringr::str_replace_all(bw_high@seqnames@values, "chr", ""))
bw_low@seqnames@values <- as.factor(stringr::str_replace_all(bw_low@seqnames@values, "chr", ""))

bw_high@seqinfo@seqnames <- stringr::str_replace_all(bw_high@seqinfo@seqnames, "chr", "")
bw_low@seqinfo@seqnames <- stringr::str_replace_all(bw_low@seqinfo@seqnames, "chr", "") 



convertensembl <- function(x){
  
  hs <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  convertedgene <- getBM(attributes = c("ensembl_transcript_id_version", "external_gene_name"),
                         filters = "ensembl_transcript_id_version",
                         values = x@range@elementMetadata@listData$symbol,
                         mart = hs)

    x@range@elementMetadata@listData$gene
  for(i in 1:nrow(convertedgene)){
    x@range@elementMetadata@listData$gene <- x@range@elementMetadata@listData$symbol %>% 
      str_replace_all(convertedgene[i,1], convertedgene[i,2]) 
    
  }
  return(x)
}

bw_high_track <- DataTrack(range=bw_high,
                     chromosome=chr_no,
                     name = 'High',
                     from =chr_start,
                     to=chr_end,
                     ylim=c(0,1.2),
                     col.histogram=c('#36A3CC'),
                     name = 'High'
)
# track 2 h3k
bw_low_track <- DataTrack(range=bw_low,
                         chromosome=chr_no,
                         from =chr_start,
                         name = 'Low',
                         to=chr_end,
                         ylim=c(0,1.2),
                         col.histogram=c('#DBBE78'), 
                         name = 'Low'
)

# convert ensembl id to gene name
gtTrack <- convertensembl(gtTrack)

bm <- useMart( biomart = "ENSEMBL_MART_ENSEMBL", 
              dataset = "hsapiens_gene_ensembl")

biomTrack <- BiomartGeneRegionTrack(genome = "hg38", chromosome = chr_no, 
                                    start = chr_start, end = chr_end,
                                    filters=list(external_gene_name=c("TRPM3", 'KLF9')),
                                    name = "ENSEMBL", biomart = bm,col.line = NULL, col= NULL, 
                                    fontface.group = 4)


ht <- HighlightTrack(trackList = c(bw_high_track,bw_low_track, biomTrack), 
start = c(70200000), 
width =100000, 
chromosome = chr_no,
fill="#A8A8A8", 
col='grey')

pdf('./Plots/Test2.pdf')
plot_list <- list(
plotTracks(c( itrack, bw_high_track, bw_low_track, biomTrack, gtrack),
           transcriptAnnotation="symbol",
           from =chr_start,
           to=chr_end,
           fontcolor = "black",
           col.axis="black",
           fontsize=15,
           showTitle=TRUE,
           collapseTranscripts = 'longest', window="auto", 
           type="histogram", cex.title=0.5, fontsize=10, cex.title=0.7,littleTicks = TRUE
),

plotTracks(c( itrack, ht, gtrack),
           transcriptAnnotation="symbol",
           from =chr_start,
           to=chr_end,
           fontcolor = "black",
           col.axis="black",
           fontsize=15,
           showTitle=TRUE,
           collapseTranscripts = 'longest', window="auto", 
           type="histogram", cex.title=0.5, fontsize=10, cex.title=0.7,littleTicks = TRUE
))




pdf('./Plots/Test2.pdf')
grid.newpage()
# 2x2 layout
pushViewport(viewport(layout=grid.layout(2, 2)))
# 1,1 first plot
pushViewport(viewport(layout.pos.col=1,layout.pos.row=1))
Gviz::plotTracks(c( itrack, bw_high_track, bw_low_track, biomTrack, gtrack),
           transcriptAnnotation="symbol",
           from =chr_start,
           to=chr_end,
           fontcolor = "black",
           col.axis="black",
           fontsize=15,
           showTitle=TRUE,
           collapseTranscripts = 'longest', window="auto", 
           type="histogram", cex.title=0.5, fontsize=10, cex.title=0.7,littleTicks = TRUE, add=TRUE)
popViewport()
# 2,2 second plot
pushViewport(viewport(layout.pos.col=2,layout.pos.row=2))
Gviz::plotTracks(c( itrack, ht, gtrack),
           transcriptAnnotation="symbol",
           from =chr_start,
           to=chr_end,
           fontcolor = "black",
           col.axis="black",
           fontsize=15,
           showTitle=TRUE,
           collapseTranscripts = 'longest', window="auto", 
           type="histogram", cex.title=0.5, fontsize=10, cex.title=0.7,littleTicks = TRUE, add=TRUE)
popViewport()
popViewport()

dev.off()

pdf('./Plots/Test3.pdf')
a <- cowplot::ggdraw() + cowplot::draw_image(magick::image_read_pdf("./Plots/Test2.pdf", density = 600))
cowplot::plot_grid(a,a, ncol=1)
dev.off()




bw_high <- '/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/BigWig/D14_High_CD8_d0.sort.rmdup.rmblackls.rmchr.norm.bw'
bw_low <- '/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/BigWig/D14_Low_CD8_d0.sort.rmdup.rmblackls.rmchr.norm.bw'



chr_no <- "9" # chromosome number
chr_start <- 70119695 # start of region
chr_end <- 70660299 # end of region


a <- get_peaks(14, 75466135, 75572235,   '/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/BigWig/D14_High_CD8_d0.sort.rmdup.rmblackls.rmchr.norm.bw', '/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/data/edited/BigWig/D14_Low_CD8_d0.sort.rmdup.rmblackls.rmchr.norm.bw')
a <- get_peaks(9, 70119695, 70660299,  bw_high, bw_low)

file.remove(a)

get_peaks <- function(CHR, START, END, HIGH, LOW, HLstart=NULL, HLend=NULL){
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
    bw_high <- import.bw(HIGH, as="GRanges")
    bw_low <- import.bw(LOW, as="GRanges")
    # change chr name to without chr (this depends on the data)
    # bw_high@seqnames@values <- as.factor(stringr::str_replace_all(bw_high@seqnames@values, "chr", ""))
    # bw_low@seqnames@values <- as.factor(stringr::str_replace_all(bw_low@seqnames@values, "chr", ""))
    # bw_high@seqinfo@seqnames <- stringr::str_replace_all(bw_high@seqinfo@seqnames, "chr", "")
    # bw_low@seqinfo@seqnames <- stringr::str_replace_all(bw_low@seqinfo@seqnames, "chr", "") 

    bw_high_track <- DataTrack(range=bw_high,
                        name = 'High',
                        chromosome=CHR,
                        from =START,
                        to=END,  col.histogram=c('#36A3CC'))
    bw_low_track <- DataTrack(range=bw_low,
                            name = 'Low',
                            chromosome=CHR,
                            from =START,
                            to=END,  col.histogram=c('#DBBE78'))

    bm <- useMart( biomart = "ENSEMBL_MART_ENSEMBL", 
                dataset = "hsapiens_gene_ensembl")

    biomTrack <- BiomartGeneRegionTrack(genome = "hg38", chromosome = CHR, 
                                        start = START, end = END,
                                        # filters=list(external_gene_name=c("TRPM3", 'KLF9')),
                                        name = "ENSEMBL", biomart = bm,col.line = NULL, col= NULL, 
                                        fontface.group = 4)

    if(!is.null(HLstart) & !is.null(HLend)){
        ht <- HighlightTrack(trackList = c(bw_high_track,bw_low_track, biomTrack), 
        start = c(70200000), 
        width =100000, 
        chromosome = CHR,
        fill="#A8A8A8", 
        col='grey')
    }else{
        filename <- paste0('./Plots/',CHR, '_',START,'_',END,'_tmp.pdf')
        pdf(filename)
        Gviz::plotTracks(c( itrack, bw_high_track, bw_low_track, biomTrack, gtrack),
        transcriptAnnotation="symbol",
        from =START,
        to=END,
        chromosome = CHR,
        fontcolor = "black",
        col.axis="black",
        fontsize=15,
        showTitle=TRUE,
        collapseTranscripts = 'longest', window="auto", 
        type="histogram", cex.title=0.5, fontsize=10, cex.title=0.7,littleTicks = TRUE, add=TRUE)
        dev.off()
    }
    return(filename)

}

