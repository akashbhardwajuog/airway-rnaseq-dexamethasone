install.packages(c( 
  "tidyverse", 
  "janitor", 
  "skimr", 
  "here", 
  "pheatmap", 
  "RColorBrewer", 
  "BiocManager" 
)) 

if (!requireNamespace("BiocManager", quietly = TRUE)) { 
  install.packages("BiocManager") 
} 

BiocManager::install(c( 
  "DESeq2", 
  "EnhancedVolcano", 
  "clusterProfiler", 
  "org.Hs.eg.db" 
)) 

BiocManager::install("TCGAbiolinks")
