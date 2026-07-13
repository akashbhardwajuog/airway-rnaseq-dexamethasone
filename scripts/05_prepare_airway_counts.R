# Project: Bulk RNA-seq Analysis of a Public Human Airway Cell Dataset 
# Author: Akash Bhardwaj 
# Purpose: Prepare count data and perform initial RNA-seq quality-control checks 

library(airway) 
library(SummarizedExperiment) 
library(tidyverse) 
library(here) 

# Load dataset 
data("airway") 

# Extract raw counts and sample metadata 
count_matrix <- assay(airway) 

sample_metadata <- colData(airway) %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("sample_id") %>% 
  mutate( 
    dex = factor(dex, levels = c("untrt", "trt")), 
    cell = factor(cell) 
  ) 

# Confirm sample order matches count-matrix columns 
stopifnot(identical(colnames(count_matrix), sample_metadata$sample_id)) 

# Remove low-expression genes: 
# Keep genes with at least 10 counts in at least 4 samples 
keep_genes <- rowSums(count_matrix >= 10) >= 4 

filtered_counts <- count_matrix[keep_genes, ] 

cat("Genes before filtering:", nrow(count_matrix), "\n") 
cat("Genes after filtering:", nrow(filtered_counts), "\n") 

# Calculate library sizes before filtering 
library_sizes <- tibble( 
  sample_id = colnames(count_matrix), 
  total_counts = colSums(count_matrix) 
) %>% 
  left_join( 
    sample_metadata %>% select(sample_id, dex, cell), 
    by = "sample_id" 
  ) 

print(library_sizes) 

# Create library-size plot 
library_size_plot <- ggplot( 
  library_sizes, 
  aes(x = reorder(sample_id, total_counts), y = total_counts, fill = dex) 
) + 
  geom_col() + 
  coord_flip() + 
  labs( 
    title = "RNA-seq Library Sizes by Sample", 
    x = "Sample", 
    y = "Total Raw Counts", 
    fill = "Treatment" 
  ) + 
  theme_minimal(base_size = 12) 

print(library_size_plot) 

# Save outputs 
write_csv( 
  library_sizes, 
  here("outputs", "tables", "airway_library_sizes.csv") 
) 

ggsave( 
  here("outputs", "figures", "airway_library_sizes.png"), 
  plot = library_size_plot, 
  width = 8, 
  height = 6, 
  dpi = 300 
) 

saveRDS( 
  filtered_counts, 
  here("data_processed", "airway_filtered_counts.rds") 
) 

write_csv( 
  sample_metadata, 
  here("data_processed", "airway_sample_metadata.csv") 
) 
