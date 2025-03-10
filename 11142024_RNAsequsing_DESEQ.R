################# Description #################
## RNA-Seq mRNA differential gene expression analysis using DESeq2 package
# The source of this script is from a tutorial
# The tutorial can be found at: https://darkomedin-datascience.medium.com/oncology-bioinformatics-project-tutorial-differential-gene-expression-and-biomarker-discovery-d3ac07db5652
# The design data can be found at: https://www.ebi.ac.uk/gxa/experiments/E-GEOD-50760/Experiment%20Design
# The Raw Counts can be found at: https://www.ebi.ac.uk/gxa/experiments/E-GEOD-50760/Downloads

###############################################

# Package Install: 
install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("ggplot2")
BiocManager::install("pheatmap")
BiocManager::install("edgeR")
BiocManager::install("ExpressionAtlas")
BiocManager::install("EnhancedVolcano", force = TRUE)


library(edgeR)
library(ExpressionAtlas)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(tidyverse)
library(tidyr)

# Data acquisition:

design <- read.delim("Data/E-GEOD-50760-experiment-design.tsv")
raw.counts <- read.delim("Data/E-GEOD-50760-raw-counts.tsv")

# Wrangling 
design<- design[1:37,]

dataset <- raw.counts[,c("Gene.Name",design$Run)]
genetable <- raw.counts[,design$Run]

Metadata = data.frame(id= design$Run, type = design[,2:ncol(design)])
colnames(Metadata)[2:ncol(Metadata)]<- sub("type.Sample.Characteristic.","",colnames(Metadata)[2:ncol(Metadata)])
metadata
#DESeq analysis

dds = DESeqDataSetFromMatrix(genetable,Metadata,~biopsy.site.)
dds<- DESeq(dds)


#Results: 
res = results(object = dds, contrast = c('biopsy.site.','primary tumor', 'normal'),
              pAdjustMethod = "holm",alpha = 0.000001)
row.names(res) <- dataset$Gene.Name
DESeq2::summary(res)


# create volcano plots using enhancedplots: 

volc <- EnhancedVolcano(res,
                lab = dataset$Gene.Name,
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 10e-7,
                FCcutoff = 2.5,
                xlim = c(-5.7,5.7),
                ylim = c(0,-log10(10.2E-12)),
                pointSize = 1.3,
                labSize = 2.6,
                title = "The results",
                subtitle = "Differnetial Expression Analysis",
                caption = "log2fc cutoff = 1.333; p-value cutoff = 10E-06",
                legendPosition = "right",
                legendLabSize = 14,
                col = c('lightblue','orange','blue','red2'),
                colAlpha = 0.6,
                drawConnectors = TRUE,
                hline = c(10e-8),
                widthConnectors = 0.5)

ggsave(filename = "Signif.png", plot = volc, units = "in",dpi = 600, width = 15,height = 10)

# save results
resord = as.data.frame(res)
finaltable1= resord[order(resord$padj),]

write.table(finaltable1, file = 'finaltable.csv',sep = ',', dec = ".", row.names = TRUE)
write.csv(finaltable1, file = 'finaltable1.csv')
