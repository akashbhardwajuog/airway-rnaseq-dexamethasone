# Project: Bulk RNA-seq Analysis of a Public Human Airway Cell Dataset 
# Author: Akash Bhardwaj 
# Purpose: Create volcano and MA plots for DESeq2 results 

library(tidyverse) 
library(EnhancedVolcano) 
library(here) 

# Load DESeq2 results 
results_table <- read_csv( 
  here("outputs", "tables", "airway_deseq2_all_results.csv"),  show_col_types = FALSE 
) 

# Replace missing adjusted p-values for plotting only 
plot_data <- results_table %>% 
  mutate( 
    padj_for_plot = ifelse(is.na(padj), 1, padj), 
    significance = case_when( 
      padj_for_plot < 0.05 & log2FoldChange >= 1 ~ "Upregulated", 
      padj_for_plot < 0.05 & log2FoldChange <= -1 ~ "Downregulated", 
      TRUE ~ "Not significant" 
    ), 
    gene_label = ifelse( 
      !is.na(symbol) & symbol != "", 
      symbol, 
      gene_id 
    ) 
  ) 

# Create volcano plot 
volcano_plot <- ggplot( 
  plot_data, 
  aes( 
    x = log2FoldChange, 
    y = -log10(padj_for_plot), 
    colour = significance 
  ) 
) + 
  geom_point(alpha = 0.6, size = 1.5) + 
  scale_colour_manual( 
    values = c( 
      "Upregulated" = "#D55E00", 
      "Downregulated" = "#0072B2", 
      "Not significant" = "grey70" 
    ) 
  ) + 
  geom_vline( 
    xintercept = c(-1, 1), 
    linetype = "dashed", 
    colour = "black" 
  ) + 
  geom_hline( 
    yintercept = -log10(0.05), 
    linetype = "dashed", 
    colour = "black" 
  ) + 
  labs( 
    title = "Volcano Plot: Dexamethasone-Treated vs Untreated Cells", 
    x = "Log2 Fold Change", 
    y = "-Log10 Adjusted P-value", 
    colour = "Classification" 
  ) + 
  theme_minimal(base_size = 12) 

print(volcano_plot) 

# Create MA plot 
ma_plot <- ggplot( 
  plot_data, 
  aes( 
    x = baseMean, 
    y = log2FoldChange, 
    colour = significance 
  ) 
) + 
  geom_point(alpha = 0.5, size = 1.2) + 
  scale_x_log10() + 
  scale_colour_manual( 
    values = c( 
      "Upregulated" = "#D55E00", 
      "Downregulated" = "#0072B2", 
      "Not significant" = "grey70" 
    ) 
  ) + 
  geom_hline( 
    yintercept = 0, 
    linetype = "dashed", 
    colour = "black" 
  ) + 
  labs( 
    title = "MA Plot: Dexamethasone-Treated vs Untreated Cells", 
    x = "Mean Normalised Expression (log10 scale)", 
    y = "Log2 Fold Change", 
    colour = "Classification" 
  ) + 
  theme_minimal(base_size = 12) 

print(ma_plot) 

# Save plots 
ggsave( 
  here("outputs", "figures", "airway_volcano_plot.png"), 
  plot = volcano_plot, 
  width = 9, 
  height = 7, 
  dpi = 300 
) 

ggsave( 
  here("outputs", "figures", "airway_ma_plot.png"), 
  plot = ma_plot, 
  width = 9, 
  height = 7, 
  dpi = 300 
) 
