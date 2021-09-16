
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
ggplot(data=pcaData, aes(x = PC1, y = PC2, color=sample)) + geom_point(size=3) + xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance")) + theme_classic() + hues::scale_color_iwanthue()
ggplot(data=pcaData, aes(x = PC1, y = PC2, color=Donor)) + geom_point(size=3) + xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance")) + theme_classic() + hues::scale_color_iwanthue()
ggplot(data=pcaData, aes(x = PC1, y = PC2, color=HighLow)) + geom_point(size=3) + xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance")) + theme_classic() +  scale_color_manual(values=c('High'='#30A3CC', 'Low'='#FCB357'))
pheatmap::pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = viridis::viridis(12, direction=-1))
dev.off()


# library(edgeR)
load('/home/sevastopol/data/mcallejac/ATAC_HighLow_ALL/RSession/cd8_2105_Csaw_fin.RData')
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


pdf('./Plots/ATAC_pca_dist_BACTHRemoved2.pdf')
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
