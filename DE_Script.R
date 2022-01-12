#Running the DEseq analysis
library("DESeq2")

#samples and condition 
samples <- c("SRR11068823", "SRR11068824","SRR11068825", "SRR11068826", "SRR11068827", "SRR11068828")
condition <- c(rep("untreated", times = 3), rep("treated", times = 3))

#make a experimental table
experiment.data <- data.frame(samples = samples, condition = condition)
experiment.data

#read in count table
counts <- as.matrix(read.csv("correct_counts.counts", sep = "\t", col.names = samples))
#counts <- as.matrix(read.csv("count_table.txt", sep = "\t", col.names = samples))
head(counts)
#make dseq matrix
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = experiment.data,
                              design = ~ condition)

dds$condition <- factor(dds$condition, levels = c("untreated","treated"))
#It will be convenient to make sure that Control is the first level in the treatment factor, so that the
#default log2 fold changes are calculated as treatment over control and not the other way around. The
#function relevel achieves this:
dds$condition <- relevel(dds$condition, ref = "untreated")
dds_s <- DESeq(dds)
res <- results(dds_s)
res


##### plots #######
#plot MA plot
plotMA(res, ylim=c(-2,2))

#shrinking for visualization and ranking
resLFC <- lfcShrink(dds_s, coef="condition_treated_vs_untreated", type="apeglm")
plotMA(resLFC, ylim=c(-2,2))

#plot count
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")

# plot PCA which needs some transformation of the data
vsd <- vst(dds_s, blind=FALSE)
rld <- rlog(dds_s, blind=FALSE)

plotPCA(vsd, intgroup=c("condition"))

#dispersion plot
plotDispEsts(dds_s)

#clustering
sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) 

library(pheatmap)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

#volcano plot (see: https://github.com/kevinblighe/EnhancedVolcano)
library(EnhancedVolcano)
EnhancedVolcano(res, lab = rownames(res), x = 'log2FoldChange', y = 'pvalue')

###### Filter results on FC en p-vals #######
res
res_filtered <- res[!is.na(res$padj) & res$padj <= 0.1,] #
res_filtered

res_up <- as.data.frame(res_filtered[res_filtered$log2FoldChange >= 2,]) #1
res_up

res_down <- as.data.frame(res_filtered[res_filtered$log2FoldChange <= -2,]) #1
res_down

###### Annotation #######
annotation <- read.csv("Annotation_Celegans.txt", sep = "\t", row.names = "Your.Input")
rownames(annotation)
annotation2 <- read.csv("mart_export.txt", sep = ",", row.names = "UniProtKB.TrEMBL.ID")

head(annotation2)
#upregulated
res_up_anno <- merge(res_up, annotation, by='row.names')
res_up_anno_sorted <- res_up_anno[order(res_up_anno$log2FoldChange, decreasing = TRUE),]
head(res_up_anno_sorted)

#downregulated
res_down_anno <- merge(res_down, annotation, by = 'row.names')
res_down_anno_sorted <- res_down_anno[order(res_down_anno$log2FoldChange, decreasing = FALSE),]
head(res_down_anno_sorted)

#write to new file
concat_res = rbind(res_up_anno_sorted, res_down_anno_sorted)

write.table(concat_res, "DEGs_Celegans.txt", append = FALSE, sep = "\t", dec = ",",
            row.names = FALSE, col.names = TRUE)
#overzicht
genesOfInterest <- read.csv(file.choose(), sep = "\t")
