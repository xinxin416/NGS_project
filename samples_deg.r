library(DESeq2)
library(xlsx)
library("BiocParallel")
library("RColorBrewer")
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
register(MulticoreParam(4))
options(stringsAsFactors = F)

# load and create data
cnts <- read.csv("gene_count_matrix.csv",header = T,row.names = 1)
colData <- read.xlsx("group_list.xlsx",sheetIndex = 1)
type <- paste(colData$condition,colData$batch,colData$treat_time,sep = "_")
colData <- data.frame(colData,type)

# create dds object and run deseq
colData$treat_time <- factor(colData$treat_time)
dds <- DESeqDataSetFromMatrix(cnts[,c(1:2,5:6)], colData[c(1:2,5:6),], ~treat_time)
dds2 <- DESeq(dds,parallel = T,BPPARAM=MulticoreParam(4))
dds2
res <- results(dds2,alpha = 0.05)
summary(res)


#treat 4h vs control 4h 
colData$condition <- factor(colData$condition)
dds <- DESeqDataSetFromMatrix(cnts[,c(5,6,10,11)], colData[c(5,6,10,11),], ~condition)
dds2 <- DESeq(dds,parallel = T,BPPARAM=MulticoreParam(4))
dds2
res <- results(dds2,alpha = 0.05)
summary(res)
sum(res$padj < 0.05, na.rm=TRUE)

# treat 24h vs control 24h 
colData$condition <- factor(colData$condition)
dds <- DESeqDataSetFromMatrix(cnts[,c(5:6,17:19)], colData[c(5,6,17:19),], ~condition)
dds2 <- DESeq(dds,parallel = T,BPPARAM=MulticoreParam(4))
dds2
res <- results(dds2,alpha = 0.05)
summary(res)
sum(res$padj < 0.05, na.rm=TRUE)

# get deseq2 diff gene list, transform IDs 
geneList <- rownames(subset(res, padj < 0.05 & pvalue<0.05 & (log2FoldChange > 2 | log2FoldChange < -2)))
geneList <- sort(geneList, decreasing = TRUE)
gene.df <- bitr(geneList, fromType = "ENSEMBL",
                toType = c("ENTREZID", "SYMBOL"),
                OrgDb = org.Hs.eg.db)

# GO over-representation test
ego_mf <- enrichGO(gene          = gene.df$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "MF",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego_mf)
barplot(ego_mf,showCategory = 20,title = "EnrichmentGO_treat_24h_MF")

ego_bp <- enrichGO(gene          = gene.df$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego_bp)
barplot(ego_bp,showCategory = 20,title = "EnrichmentGO_treat_24h_cc")
#dotplot(ego_bp,showCategory = 20,title = "EnrichmentGO_treat4h_BP")
