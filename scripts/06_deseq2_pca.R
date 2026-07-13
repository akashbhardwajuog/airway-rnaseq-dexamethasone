# Project: Bulk RNA-seq Analysis of a Public Human Airway Cell Dataset 
# Author: Akash Bhardwaj 
# Purpose: Perform DESeq2 normalisation and PCA for RNA-seq quality control 

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

# Confirm sample order 
stopifnot(identical(colnames(count_matrix), rownames(sample_metadata))) 

# Filter low-expression genes 
keep_genes <- rowSums(count_matrix >= 10) >= 4 
filtered_counts <- count_matrix[keep_genes, ] 

# Create DESeq2 object. 
# The design includes cell line so treatment effects are assessed while accounting for donor cell line. 
dds <- DESeqDataSetFromMatrix( 
  countData = filtered_counts, 
  colData = sample_metadata, 
  design = ~ cell + dex 
) 

# Run DESeq2 normalisation and model fitting 
dds <- DESeq(dds) 

# Variance-stabilising transformation for PCA and visualisation 
vsd <- vst(dds, blind = FALSE) 

# Extract transformed values 
vsd_matrix <- assay(vsd) 

# Run PCA on transformed gene expression values 
pca_result <- prcomp(t(vsd_matrix), scale. = FALSE) 

# Calculate variance explained 
pca_variance <- tibble( 
  principal_component = paste0("PC", 1:length(pca_result$sdev)), 
  variance_explained_percent = round( 
    (pca_result$sdev^2 / sum(pca_result$sdev^2)) * 100, 
    1 
  ) 
) 

print(pca_variance) 

# Create PCA plotting data 
pca_scores <- as.data.frame(pca_result$x) %>% 
  tibble::rownames_to_column("sample_id") %>% 
  left_join( 
    sample_metadata %>% 
      select(sample_id, dex, cell), 
    by = "sample_id" 
  ) 


# PCA plot: colour shows treatment and shape shows donor cell line 
pca_plot <- ggplot( 
  pca_scores, 
  aes(x = PC1, y = PC2, colour = dex, shape = cell) 
) + 
  geom_point(size = 4, alpha = 0.9) + 
  labs( 
    title = "PCA of Variance-Stabilised RNA-seq Data", 
    subtitle = "Human airway smooth-muscle cells: dexamethasone treatment and donor cell line", 
    x = paste0("PC1 (", pca_variance$variance_explained_percent[1], "%)"),
    y = paste0("PC2 (", pca_variance$variance_explained_percent[2], "%)"),
    colour = "Treatment", 
    shape = "Cell line" 
  ) + 
  theme_minimal(base_size = 12) 

print(pca_plot) 

# Save outputs 
write_csv( 
  pca_variance, 
  here("outputs", "tables", "airway_pca_variance_explained.csv") 
) 

write_csv( 
  pca_scores, 
  here("outputs", "tables", "airway_pca_scores.csv") 
) 

ggsave( 
  here("outputs", "figures", "airway_pca_plot.png"), 
  plot = pca_plot, 
  width = 8, 
  height = 6, 
  dpi = 300 
) 

saveRDS( 
  dds, 
  here("data_processed", "airway_deseq2_object.rds") 
) 
