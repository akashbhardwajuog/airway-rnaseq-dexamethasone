# Project: Bulk RNA-seq Analysis of a Public Human Airway Cell Dataset
# Author: Akash Bhardwaj 
# Purpose: Inspect the Bioconductor airway RNA-seq dataset 

library(airway) 
library(SummarizedExperiment) 
library(tidyverse) 
library(here) 

# Load the dataset 
data("airway") 

# Inspect the full object 
airway 

# Extract raw count matrix 
count_matrix <- assay(airway) 

# Extract sample metadata 
sample_metadata <- colData(airway) %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("sample_id") 

# Check dimensions 
cat("Number of genes:", nrow(count_matrix), "\n") 
cat("Number of samples:", ncol(count_matrix), "\n") 

# Inspect sample metadata 
print(sample_metadata) 

# Check treatment groups 
print(table(sample_metadata$dex)) 

# Save counts and metadata 
write.csv( 
  count_matrix, 
  here("data_raw", "airway_raw_counts.csv"), 
  row.names = TRUE 
) 

write_csv( 
  sample_metadata, 
  here("data_raw", "airway_sample_metadata.csv") 
) 
