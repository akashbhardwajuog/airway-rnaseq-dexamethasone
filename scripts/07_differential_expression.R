# Project: Bulk RNA-seq Analysis of a Public Human Airway Cell Dataset 
# Author: Akash Bhardwaj 
# Purpose: Identify genes differentially expressed after dexamethasone treatment 

library(airway) 
library(SummarizedExperiment) 
library(DESeq2) 
library(tidyverse) 
library(here) 

# Load dataset 
data("airway") 

# Extract counts and metadata 
count_matrix <- assay(airway) 

sample_metadata <- colData(airway) %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("sample_id") %>% 
  mutate( 
    dex = factor(dex, levels = c("untrt", "trt")), 
    cell = factor(cell) 
  ) 

rownames(sample_metadata) <- sample_metadata$sample_id 

# Filter low-expression genes 
keep_genes <- rowSums(count_matrix >= 10) >= 4 
filtered_counts <- count_matrix[keep_genes, ] 

# DESeq2 design controls for donor cell line 
dds <- DESeqDataSetFromMatrix( 
  countData = filtered_counts, 
  colData = sample_metadata, 
  design = ~ cell + dex 
) 

dds <- DESeq(dds) 

# Compare dexamethasone-treated samples against untreated samples 
results_dex <- results( 
  dds, 
  contrast = c("dex", "trt", "untrt"), 
  alpha = 0.05 
) 

# Convert results to table 
results_table <- as.data.frame(results_dex) %>% 
  tibble::rownames_to_column("gene_id") %>% 
  as_tibble() 

# Add gene annotations from the airway dataset 
# Add gene annotations from the airway dataset
gene_annotations <- rowData(airway) %>%
  as.data.frame() %>%
  select(gene_id, gene_name, symbol)

results_table <- results_table %>% 
  left_join(gene_annotations, by = "gene_id") %>% 
  arrange(padj) 

# Define significantly differentially expressed genes 
significant_genes <- results_table %>% 
  filter( 
    !is.na(padj), 
    padj < 0.05, 
    abs(log2FoldChange) >= 1 
  ) 

upregulated_genes <- significant_genes %>% 
  filter(log2FoldChange > 0) 

downregulated_genes <- significant_genes %>% 
  filter(log2FoldChange < 0) 

# Print summary 
cat("Total genes tested:", nrow(results_table), "\n") 
cat("Significant genes (FDR < 0.05 and |log2FC| >= 1):", nrow(significant_genes), "\n") 
cat("Upregulated after treatment:", nrow(upregulated_genes), "\n") 
cat("Downregulated after treatment:", nrow(downregulated_genes), "\n") 

# Save result tables 
write_csv( 
  results_table, 
  here("outputs", "tables", "airway_deseq2_all_results.csv") 
) 

write_csv( 
  significant_genes, 
  here("outputs", "tables", "airway_deseq2_significant_genes.csv") 
) 

write_csv( 
  head(significant_genes, 20), 
  here("outputs", "tables", "airway_top_20_de_genes.csv") 
) 
