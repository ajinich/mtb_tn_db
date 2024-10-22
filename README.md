# Project Overview
This code is relevant to the manuscript titled 'Comprehensive genetic interaction maps reveal the landscape of conditional genetic interactions in Mycobacterium tuberculosis,' available at bioRxiv (https://www.biorxiv.org/content/10.1101/2021.03.05.434127v2).

This repository contains code and data for analyzing TnSeq data, including pre-processing, dimensionality reduction, unsupervised and supervised machine learning, and statistical analysis. The goal is to elucidate gene essentiality and interactions in *Mycobacterium tuberculosis* using standardized data from various sources.

## Data Standardization and Pre-processing

The initial standardization of TnSeq data is based on analysis performed using the TRANSIT tool by M. deJesus. The raw output files from TRANSIT can be found in Zenodo at Jinich, A., & DeJesus, M. (2024). Mycobacterium tuberculosis transposon sequencing database (MtbTnSeq) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.13937394. These have been processed into single q-value and log2 fold change tables:

- **Starting Data**: `result_logfc_matrix_2021_05_18.csv`, generated and provided by M. deJesus.
- **Additional Data Stages**:
  - CC-panel data added using `Clare_CC_panel_dataset.ipynb`.
  - BxD panel data added using `BxD_data_prep.ipynb`.

Final processed output files include:
- `result_qval_matrix_2023_02_20_CC_BxD.csv`
- `result_lfc_matrix_2023_02_20_CC_BxD.csv`

## Filtering and Binarization

Data is filtered and binarized using log2 fold change and q-value thresholds:
- **Notebook**: `filter_and_binarize_lfc_qval.ipynb`
- **Output**:
  - `result_qval_matrix_2023_02_20_CC_BxD_processed.csv`
  - `result_lfc_matrix_2023_02_20_CC_BxD_processed.csv`
  - `result_bin_matrix_2023_02_20_CC_BxD_processed.csv`

These processed files serve as inputs for the subsequent analysis.

## Dimensionality Reduction

- **Notebook**: `umap_by_screens.ipynb`
- **Input Files**:
  - `result_lfc_matrix_2023_02_20_CC_BxD_processed.csv`
  - `result_logfc_matrix_2023_02_20_CC_w_mbio_BxD.csv`
- **Output**: UMAP projections by screen type (in vivo or in vitro) and publication.

## Data Preparation for Dash App

- **Notebook**: `Data_for_dash.ipynb`
- **Input Files**: `column_descriptors_standardized_021023.xlsx` and processed q-val/lfc files.
- **Output**: Metadata files for use in a Dash app.

## Unsupervised Learning: Gene Clustering

- **Notebook**: `umap_by_genes.ipynb`
- **Input Files**: Processed q-val, lfc, and binarized data merged with functional annotations.
- **Output Files**: `df_umap_lfc_08_08_24.csv`, `df_umap_10.csv`
- **Results**: UMAP visualizations and k-means clustering of gene data.

## Supervised Machine Learning

- **Notebook**: `make_data.ipynb`
- **Purpose**: Prepares filtered datasets: `lfc_mb_filt_07_22_24.csv`, `qval_mb_filt_07_22_24.csv`, `bin_mb_filt.csv`.
- **Methods**:
  - `Ensemble.ipynb`: Combines logistic regression with L1 regularization and Random Forest.
  - `Tree-based_std.ipynb`: Random Forest classifier.
  - `Log_Reg.ipynb`: XGBoost classifier (used in the associated publication).

## Numerical Analysis

- **counting_screen_types.ipynb**: Classifies screens (in vivo, in vitro, etc.).
- **genomic_distance_corr.ipynb**: Analyzes whether neighboring genes in the genome share similar TnSeq profiles based on UMAP positions.
- **sanity_checks.ipynb**: Compares essentiality calls between standardized and SI datasets.
- **Tn_library_stats.ipynb**: Performs statistical analysis of gene coessentiality counts.

## SI Data Processing

- **Notebook**: `SI_qvals_log2FC.ipynb`
- **Purpose**: Processes SI data tables with TnSeq screens from the literature.

## Visualization

- **Notebook**: `bubble_plots_v2.ipynb`
- **Purpose**: Generates bubble plots of genes that are both essential and unknown in a given screen. These plots are used in the Dash app.

## Repository Structure

- **data/**: Contains input and processed data files.
- **notebooks/**: Jupyter notebooks for analysis and visualization.
- **figures/**: Output figures from various analyses.
- **dash_app/**: Code and resources for the Dash app interface.

## License

This project is licensed under the [MIT License](LICENSE).



