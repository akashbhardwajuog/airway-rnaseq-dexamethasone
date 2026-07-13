# Install a small public RNA-seq teaching dataset 

if (!requireNamespace("BiocManager", quietly = TRUE)) { 
  install.packages("BiocManager") 
} 

BiocManager::install("airway") 
