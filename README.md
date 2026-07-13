# Bulk RNA-seq Analysis of Dexamethasone Response in Human Airway Smooth-Muscle Cells 
 
## Project objective 
 
This educational project analyses public bulk RNA-seq count data from primary human 
airway smooth-muscle cell lines treated with dexamethasone or left untreated. 
 
The analysis demonstrates an end-to-end RNA-seq workflow: data inspection, gene 
filtering, library-size quality control, DESeq2 normalisation, principal component 
analysis, differential-expression analysis, Gene Ontology enrichment and heatmap 
visualisation. 

## PCA overview ![PCA of variance-stabilised RNA-seq data]
(outputs/figures/airway_pca_plot.png)
 
## Dataset 
 
The analysis uses the `airway` Bioconductor experiment package. The underlying public 
RNA-seq study is available through NCBI Gene Expression Omnibus (GEO), accession 
GSE52778. 
 
- Organism: human 
- Cell type: primary airway smooth-muscle cells 
- Treatment: dexamethasone versus untreated 
- Samples: 8 total; 4 untreated and 4 treated 
- Donor cell lines: 4, with one untreated and one treated sample per donor 
 
## Analysis workflow 
 
1. Loaded raw gene-level count data and sample metadata. 
2. Filtered low-expression genes, retaining genes with at least 10 counts in at least 4 samples. 
3. Checked raw library sizes across samples. 
4. Used DESeq2 with the design `~ cell + dex` to account for donor cell line while testing treatment-associated differences. 
5. Applied variance-stabilising transformation for PCA and visualisation. 
6. Identified differentially expressed genes using adjusted p-value < 0.05 and absolute log2 fold change >= 1. 
7. Performed Gene Ontology Biological Process enrichment analysis. 
8. Visualised results using PCA, volcano, MA and heatmap plots. 
 
## Key results 
 
- Genes before filtering: 63,677 
- Genes after filtering: 16,139 
- PCA variance explained: PC1 = 42.9%; PC2 = 22.3% 
- Differentially expressed genes: 951 
 - Upregulated after dexamethasone: 490 
 - Downregulated after dexamethasone: 461 
- Gene Ontology enrichment highlighted extracellular-matrix organisation, muscle-system processes, actin-filament regulation and signalling-related processes. 
 
## Key outputs 
 
### Figures 
 
- `outputs/figures/airway_library_sizes.png` 
- `outputs/figures/airway_pca_plot.png` 
- `outputs/figures/airway_volcano_plot.png` 
- `outputs/figures/airway_ma_plot.png` 
- `outputs/figures/airway_go_bp_dotplot.png` 
- `outputs/figures/airway_top_30_de_genes_heatmap.png` 
 
### Tables 
 
- `outputs/tables/airway_library_sizes.csv` 
- `outputs/tables/airway_pca_variance_explained.csv` 
- `outputs/tables/airway_pca_scores.csv` 
- `outputs/tables/airway_deseq2_all_results.csv` 
- `outputs/tables/airway_deseq2_significant_genes.csv` 
- `outputs/tables/airway_top_20_de_genes.csv` 
- `outputs/tables/airway_go_bp_enrichment.csv` 
 
## Interpretation and limitations 
 
This is an exploratory educational analysis of a small public dataset. The study includes 
only eight samples from four donor cell lines, so results should not be interpreted as 
clinical evidence or proof of causal mechanisms. Gene Ontology terms are annotation-based 
and may reflect overlapping genes and broad biological processes. 
 
## Technical skills demonstrated 
 
- R/RStudio 
- Bioconductor 
- SummarizedExperiment 
- DESeq2 
- RNA-seq quality control 
- Gene filtering 
- Variance-stabilising transformation 
- Principal component analysis 
- Differential-expression analysis 
- Gene Ontology enrichment using clusterProfiler 
- Data visualisation with ggplot2 and pheatmap 
- Reproducible project organisation with Git/GitHub 
 
## Author 
 
Akash Bhardwaj   
MSc Precision Medicine candidate, University of Glasgow 
 