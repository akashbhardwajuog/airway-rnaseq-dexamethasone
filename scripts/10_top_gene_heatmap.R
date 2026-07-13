# Project: Bulk RNA-seq Analysis of a Public Human Airway Cell Dataset 
# Author: Akash Bhardwaj 
# Purpose: Create a heatmap of the top differentially expressed genes 

library(airway) 
library(SummarizedExperiment) 
library(DESeq2) 
library(pheatmap) 
library(RColorBrewer) 
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

# Recreate DESeq2 object and variance-stabilised data 
dds <- DESeqDataSetFromMatrix( 
  countData = filtered_counts, 
  colData = sample_metadata, 
  design = ~ cell + dex 
) 

dds <- DESeq(dds) 
vsd <- vst(dds, blind = FALSE) 

# Load DE results and select top 30 genes by adjusted p-value 
results_table <- read_csv( 
  here("outputs", "tables", "airway_deseq2_all_results.csv"), 
  show_col_types = FALSE 
) 

top_genes <- results_table %>% 
  filter(!is.na(padj)) %>% 
  arrange(padj) %>% 
  slice_head(n = 30) 

# Use Ensembl IDs to extract expression values 
top_gene_ids <- top_genes$gene_id 

heatmap_matrix <- assay(vsd)[top_gene_ids, ] 

# Use gene symbols for labels where available 
gene_labels <- ifelse( 
  !is.na(top_genes$symbol) & top_genes$symbol != "", 
  top_genes$symbol, 
  top_genes$gene_id 
) 

rownames(heatmap_matrix) <- make.unique(gene_labels) 

# Row-scale expression values for visual comparison across samples 
heatmap_scaled <- t(scale(t(heatmap_matrix))) 

# Sample annotation 
annotation_col <- sample_metadata %>%  
  dplyr::select(dex, cell)

# Define annotation colours 
annotation_colors <- list( 
  dex = c( 
    untrt = "#F8766D", 
    trt = "#00BFC4" 
  ), 
  cell = setNames( 
    brewer.pal(4, "Set2"), 
    levels(sample_metadata$cell) 
  ) 
) 

# Save heatmap as PNG 
png( 
  filename = here("outputs", "figures", "airway_top_30_de_genes_heatmap.png"), 
  width = 1800, 
  height = 2200, 
  res = 220 
) 

pheatmap( 
  heatmap_scaled, 
  annotation_col = annotation_col, 
  annotation_colors = annotation_colors, 
  cluster_rows = TRUE, 
  cluster_cols = TRUE, 
  show_colnames = TRUE, 
  show_rownames = TRUE, 
  fontsize_row = 8, 
  fontsize_col = 10, 
  main = "Top 30 Differentially Expressed Genes" 
) 

dev.off() 
