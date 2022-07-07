#Intial Setup

if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)
library(SummarizedExperiment)
library(ggplot2)
BiocManager::install("DESeq2")
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

#Downloading and Data Wrangling

query <- GDCquery(project = "TCGA-PRAD", 
                  data.category = "Gene expression",
                  data.type = "Gene expression quantification",
                  experimental.strategy = "RNA-Seq",
                  platform = "Illumina HiSeq",
                  file.type = "results",
                  legacy = TRUE)
GDCdownload(query, method = "api", files.per.chunk = 100,
            directory = "/content/Data")
mrna_df <- GDCprepare(query, directory = "/content/Data")

# Remove columns we dont need, keep counts

mrna_meta <- mrna_df$sample
mrna_meta <- cbind(mrna_meta, mrna_df$definition)
mrna_df <- assay(mrna_df)
delim_fn = function(x, n, i){
    do.call(c, lapply(x, function(X)
        paste(unlist(strsplit(X, "-"))[(n+1):(i)], collapse = "-")))
}
colnames(mrna_df) <- delim_fn(x = colnames(mrna_df), n = 0, i = 4)
mrna_meta <- as.data.frame(mrna_meta)
mrna_df <- as.data.frame(mrna_df)

# Remove metastatic sample

metastatic_key <- mrna_meta[which(mrna_meta[,2] == "Metastatic"),]
mrna_meta <- mrna_meta[!mrna_meta[,2] == metastatic_key[,2],]
mrna_df <- mrna_df[, -grep(paste0(metastatic_key[,1]), colnames(mrna_df))]
mrna_meta[,2] <- as.character(mrna_meta[,2])
mrna_meta[,2] <- gsub("Primary solid Tumor", "Tumor", mrna_meta[,2])
mrna_meta[,2] <- gsub("Solid Tissue Normal", "Normal", mrna_meta[,2])
mrna_meta[,2] <- as.factor(mrna_meta[,2])
levels(mrna_meta[,2])
colnames(mrna_meta) <- c("cases", "Condition")

# DESeq2 Analysis

mrna_dds <- DESeqDataSetFromMatrix(round(mrna_df), colData = mrna_meta, design = ~ Condition)
mrna_dds$Condition <- relevel(mrna_dds$Condition, ref = "Normal")
mrna_dds <- DESeq(mrna_dds)
vsd <- varianceStabilizingTransformation(mrna_dds, blind=FALSE)

# Dispersions Plot

plotDispEsts(mrna_dds, main="Dispersion plot")
mrna_res <- results(mrna_dds, name = "Condition_Tumor_vs_Normal")

# MA Plot

plotMA(mrna_res)

