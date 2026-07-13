#Project: Bulk RNA-seq Analysis of Dexamethasone Response in Human Airway Smooth-Muscle Cells
# Author: Akash Bhardwaj 
# Purpose: Verify project folders and load core packages 

library(here) 

required_folders <- c( 
  "data_raw", 
  "data_processed", 
  "scripts", 
  "outputs", 
  "outputs/figures", 
  "outputs/tables", 
  "reports", 
  "references" 
) 

folder_check <- data.frame( 
  folder = required_folders, 
  exists = dir.exists(here(required_folders)) 
) 

print(folder_check) 
