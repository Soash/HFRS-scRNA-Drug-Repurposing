
# In Silico Drug Repurposing for Hantaan Virus Infection Using scRNA-seq Data 🧬💊

![R](https://img.shields.io/badge/R-276DC3?style=for-the-badge&logo=r&logoColor=white)
![Seurat](https://img.shields.io/badge/Seurat-v5.4.0-blue?style=for-the-badge)
![ASGARD](https://img.shields.io/badge/ASGARD-v1.0.0-orange?style=for-the-badge)
![Data](https://img.shields.io/badge/GEO-GSE161354-green?style=for-the-badge)

## 📌 Project Overview
This repository contains the complete computational pipeline and R scripts used for the study **"In Silico Drug Repurposing for Hantaan Virus Infection Using scRNA-seq Data."** Hemorrhagic Fever with Renal Syndrome (HFRS) is a severe viral disease caused by the Hantaan virus (HTNV), characterized by systemic hyper-inflammation and the abnormal peripheral presence of Nucleated Red Blood Cells (NRBCs). Because there are currently no targeted FDA-approved therapeutics for HFRS, this project leverages high-resolution single-cell RNA sequencing (scRNA-seq) to capture the disease's transcriptomic signature and utilizes the **ASGARD** pipeline to identify robust, repurposable drug candidates (notably Vorinostat and Tretinoin) from the LINCS L1000 pharmacogenomic database.

## 📂 Repository Structure

The workflow is divided into four sequential phases. The scripts should be executed in numerical order.

| Script | Description |
| :--- | :--- |
| `1_qc.R` | **Quality Control & Filtering:** Processes raw 10X Genomics data. Applies static thresholds (`nFeature_RNA > 200`, `percent.mt < 20`) and systematically purges mitochondrial genes to prevent downstream clustering bias. |
| `2_integration.R` | **Data Integration & Clustering:** Merges 8 patient and control samples. Performs log-normalization, Canonical Correlation Analysis (CCA) integration via Seurat v5, PCA, and unsupervised graph-based clustering (resolution = 0.5) visualized via UMAP. |
| `3_differential_expression.R` | **DEG Extraction:** Computes cluster-specific Differential Expression Genes (Disease vs. Normal) using the Wilcoxon Rank-Sum test (`lnFC > 0.1`, `min.pct > 0.1`) tailored for ASGARD input. |
| `4_asgard_hfrs.R` | **Drug Repurposing Pipeline:** Maps the specialized single-cell DEGs against the native hematopoietic and lymphoid tissue profiles in the LINCS L1000 database. Calculates the multi-cluster therapeutic scores. |
| `/results/` | *(Directory)* Contains generated CSV files, including the top FDA-approved drug candidates and cluster marker lists. |
| `/plots/` | *(Directory)* Contains UMAP projections, FeaturePlots for hallmark genes (*HBA1*, *HBA2*, *HBB*), and QC violin plots. |

## 📊 Data Availability
The raw single-cell RNA sequencing data analyzed in this study must be downloaded prior to running the pipeline. It is publicly available in the NCBI Gene Expression Omnibus (GEO) database under the accession number:
* **[GSE161354](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161354)**

The LINCS L1000 perturbation profiles required for the ASGARD pipeline must be acquired locally:
* `GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx`
* `GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx`

## 🛠️ Prerequisites & Installation

To run this pipeline, you will need **R (>= 4.1.0)** and the following core packages installed:

```R
# Install core CRAN packages
install.packages(c("dplyr", "ggplot2", "future", "qs2"))

# Install Seurat v5
install.packages("Seurat")

# Install ASGARD (A Single-cell Guided Pipeline to Aid Repurposing of Drugs)
# Refer to the original ASGARD documentation for specific installation instructions via devtools/BiocManager.

```

## 🚀 Usage

1. Clone this repository to your local machine or high-performance cluster.
2. Download the raw patient matrices from GSE161354 and place them in `GSE161354_RAW/disease` and `GSE161354_RAW/normal` directories as configured in the scripts.
3. Download the LINCS L1000 `.gctx` files and place them in the `database/L1000/` directory.
4. Execute the R scripts sequentially from `1_qc.R` through `4_asgard_hfrs.R`.
*(Note: Parallel processing `plan("sequential")` and garbage collection are built into the scripts to manage memory footprint, but a minimum of 32GB RAM is recommended for the CCA integration step).*

## 🔬 Key Findings

Our multi-cluster scoring matrix successfully bypassed the limitations of bulk-tissue RNA sequencing, ensuring rare but clinically vital populations (like NRBCs) were evaluated. The pipeline identified 13 statistically significant FDA-approved drug candidates.

* **Vorinostat:** A potent HDAC inhibitor theoretically capable of suppressing hyper-active immune transcription to mitigate the lethal cytokine storm and vascular collapse.
* **Tretinoin (ATRA):** A differentiation agent theoretically capable of bypassing the HTNV-induced maturation blockade, forcing the terminal differentiation and enucleation of stalled NRBCs.

## ✒️ Citation

If you utilize this pipeline, code, or adapt the methodology in your own work, please cite the upcoming manuscript:

> **Md Abu Sadat Soash** (2026). *In Silico Drug Repurposing for Hantaan Virus Infection Using scRNA-seq Data.* Department of Genetic Engineering and Biotechnology, Shahjalal University of Science and Technology.

## ✉️ Contact

For questions regarding the code, data access, or the computational methodology, please reach out via GitHub Issues or contact:
**Md Abu Sadat Soash** - [mdabu13@student.sust.edu](mailto:mdabu13@student.sust.edu)

