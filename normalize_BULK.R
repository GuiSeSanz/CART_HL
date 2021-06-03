library(DESeq2)


# read the data 

counts <- readxl::read_excel('./Data/Counts_rnaseq_bulk_public.xlsx')
metadata <- read.table('./Data/METADATA_JRR.csv', sep=',' , header=TRUE)


counts_df <- as.data.frame(counts)
colnames(counts_df) <- counts_df[2,]
colnames(counts_df)[1] <- 'gene_name'
counts_df <- counts_df[3:nrow(counts_df),]

duplicated <- counts_df$gene_name[which(duplicated(counts_df$gene_name))]
counts_df <- counts_df[!counts_df$gene_name %in% duplicated, ]
rownames(counts_df) <- counts_df$gene_name
counts_df <- counts_df[,2:ncol(counts_df)]

rownames(metadata) <- metadata$Patient.ID
metadata$binary_response <- factor(ifelse( metadata$OS == 'CR' | metadata$OS == 'PRTD', 1, 0))

cts_matrix <- as.matrix(counts_df)
class(cts_matrix) <- 'numeric'

stopifnot(all(rownames(metadata) %in% colnames(cts_matrix)))
stopifnot(all(rownames(metadata) == colnames(cts_matrix)))

dds <- DESeqDataSetFromMatrix(countData = cts_matrix,
                              colData = metadata,
                              design = ~ binary_response)



featureData <- data.frame(gene=rownames(cts_matrix))
mcols(dds) <- DataFrame(mcols(dds), featureData)     

# Pre-filtering
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds$condition <- factor(dds$binary_response, levels = c("0","1"))
dds <- DESeq(dds)

res_2 <- results(dds, contrast=c("binary_response","0","1"))
resultsNames(res_2)

res <- results(dds)
res05 <- results(dds, alpha=0.05)
sum(res05$padj < 0.01, na.rm=TRUE)

pdf('./Plots/Volcano_0.5.pdf')
plotMA(res05)

plotMA( lfcShrink(dds, coef="binary_response_1_vs_0"))

dev.off()

res05_df <- as.data.frame(res05)
res05_df <- res05_df[!is.na(res05_df$padj), ]

# compare the DE genes with signature
signature_genes_CD4 <- read.table('./Data/signature/BatchK_CD4_SIGGenes_BASAL_BOTH_LowvsHigh.tsv', sep='\t', header=TRUE)
signature_genes_CD8 <- read.table('./Data/signature/BatchK_CD8_SIGGenes_BASAL_BOTH_LowvsHigh.tsv', sep='\t', header=TRUE)


futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
VennDiagram::venn.diagram(
    x=list(signature_genes_CD4$sigGenes_symbol, signature_genes_CD8$sigGenes_symbol, rownames(res05_df[res05_df$padj < 0.05,]) ),
    category.names = c('signature_genes_CD4', 'signature_genes_CD8', 'DE_genes_0.05'), filename='./Plots/VennDiagram_DeVsSignatures.png',
    output=TRUE, alpha = 0.50, fill = hues::iwanthue(3), margin = 0.05
)

dds <- estimateSizeFactors(dds)
normalized_data <- counts(dds, normalized=T)
saveRDS(normalized_data, './Data/Normalized_Counts_Bulk.rds')



