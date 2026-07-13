# Project: Bulk RNA-seq Analysis of a Public Human Airway Cell Dataset 
# Author: Akash Bhardwaj 
# Purpose: Perform Gene Ontology biological-process enrichment analysis 

library(tidyverse) 
library(clusterProfiler) 
library(org.Hs.eg.db) 
library(here) 

# Load significant DESeq2 results 
significant_genes <- read_csv( 
  here("outputs", "tables", "airway_deseq2_significant_genes.csv"), 
  show_col_types = FALSE 
) 

# Keep valid Ensembl gene IDs and remove version suffixes if present 
ensembl_ids <- significant_genes %>% 
  filter(!is.na(gene_id)) %>% 
  mutate(gene_id = sub("\\..*$", "", gene_id)) %>% 
  pull(gene_id) %>% 
  unique() 

# Convert Ensembl IDs to Entrez IDs for Gene Ontology analysis 
gene_conversion <- bitr( 
  ensembl_ids, 
  fromType = "ENSEMBL", 
  toType = "ENTREZID", 
  OrgDb = org.Hs.eg.db 
) 

# Keep unique Entrez IDs 
entrez_ids <- unique(gene_conversion$ENTREZID) 

cat("Input Ensembl IDs:", length(ensembl_ids), "\n") 
cat("Mapped Entrez IDs:", length(entrez_ids), "\n") 

# Perform Gene Ontology Biological Process enrichment 
go_enrichment <- enrichGO( 
  gene = entrez_ids, 
  OrgDb = org.Hs.eg.db, 
  keyType = "ENTREZID", 
  ont = "BP", 
  pAdjustMethod = "BH", 
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.05, 
  readable = TRUE 
) 

# Convert results to a table 
go_results <- as.data.frame(go_enrichment) 

# Save full enrichment results 
write_csv( 
  go_results, 
  here("outputs", "tables", "airway_go_bp_enrichment.csv") 
) 

# Display top results 
print(head(go_results, 10)) 

# Create dot plot if results are available 
if (nrow(go_results) > 0) { 
  
  go_plot <- dotplot( 
    go_enrichment, 
    showCategory = 15, 
    title = "Top Enriched Gene Ontology Biological Processes" 
  ) 
  
  print(go_plot) 
  
  ggsave( 
    here("outputs", "figures", "airway_go_bp_dotplot.png"), 
    plot = go_plot, 
    width = 10, 
    height = 8, 
    dpi = 300 
  ) 
} else { 
  message("No significantly enriched GO biological-process terms were identified.") 
} 
