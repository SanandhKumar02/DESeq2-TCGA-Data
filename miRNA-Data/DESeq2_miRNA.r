if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")
BiocManager::install("EnhancedVolcano")
BiocManager::install("DESeq2")
library(TCGAbiolinks)
library(SummarizedExperiment)
library(EnhancedVolcano)
library(ggplot2)
library(DESeq2)

#Function for Identifiying Upregulated and Downregulated genes

get_upregulated <- function(df){
    key <- intersect(rownames(df)[which(df$log2FoldChange>=1)],
              rownames(df)[which(df$pvalue<=0.05)])
    
    results <- as.data.frame((df)[which(rownames(df) %in% key),])
    return(results)
  }

get_downregulated <- function(df){
  key <- intersect(rownames(df)[which(df$log2FoldChange<=-1)],
            rownames(df)[which(df$pvalue<=0.05)])
  
  results <- as.data.frame((df)[which(rownames(df) %in% key),])
  return(results)
}

#Downlading and Data wrangling

mirna_query <- GDCquery(project = "TCGA-PRAD",
                        data.category = "Transcriptome Profiling",
                        data.type = "miRNA Expression Quantification",
                        experimental.strategy = "miRNA-Seq")


GDCdownload(mirna_query, method = "api", files.per.chunk = 100)
miR_df <- GDCprepare(mirna_query)

#Remove the column not needed
rownames(miR_df) <- miR_df$miRNA_ID
miR_df <- miR_df[,-1]
number_cols <- ncol(miR_df)
subset <- seq(from = 1, to = number_cols, by = 3)
miR_df <- miR_df[, subset]

#Deleting read_count
colnames(miR_df) <- gsub(".*_","",colnames(miR_df))

#Matching Meta-Data
miR_meta <- mirna_query[[1]][[1]]
miR_meta <- miR_meta[,c("cases", "sample_type")]
miR_meta$sample_type <- as.character(miR_meta$sample_type)

#Metastatic Tumor
metastatic_key <- miR_meta[which(miR_meta$sample_type == "Metastatic"),]
miR_meta <- miR_meta[!miR_meta[,2] == metastatic_key[,2],]
miR_df <- miR_df[, -grep(paste0(metastatic_key[,1]), colnames(miR_df))]

#Renaming Conditions
miR_meta$sample_type <- gsub("Primary solid Tumor", "Tumor", miR_meta$sample_type)
miR_meta$sample_type <- gsub("Solid Tissue Normal", "Normal", miR_meta$sample_type)
miR_meta$sample_type <- as.factor(miR_meta$sample_type)
colnames(miR_meta) <- c("cases", "Condition")

## tidy vars
rm(mirna_query)
rm(subset)
rm(number_cols)
rm(metastatic_key)

#Differential Expression Analysis
miR_dds <- DESeqDataSetFromMatrix(miR_df, colData = miR_meta, design = ~ Condition)
miR_dds$Condition <- relevel(miR_dds$Condition, ref = "Normal")
miR_dds <- DESeq(miR_dds)
resultsNames(miR_dds)

#DESeq2
miR_res <- results(miR_dds, alpha = 0.05, name = "Condition_Primary.Tumor_vs_Normal")
miR_res_df <- as.data.frame(miR_res)
summary(miR_res)

#Finding Upreg and Downreg genes
miR_downreg <- get_downregulated(miR_res)
miR_upreg <- get_upregulated(miR_res)

miR_upreg$miRNA_id <- rownames(miR_upreg)
miR_downreg$miRNA_id <- rownames(miR_downreg)

#Volcano Plot
EnhancedVolcano(miR_res_df, x="log2FoldChange", y="padj", lab = miR_res_df$miRNA_id, pCutoff = 0.05, FCcutoff = 2)



