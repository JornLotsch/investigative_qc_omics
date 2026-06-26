# Investigative QC Omics

## A computational investigative quality control framework for omics detecting procedural errors in the data-to-analysis workflow

This repository provides R code for an **investigative QC framework** that detects data structure inconsistent with the study hypothesis, including anomalies arising from procedural errors such as data handover mistakes. The method combines UMAP projection, Ward hierarchical clustering, Voronoi visualization, and supervised random forest classification to test whether data structure reflects biological signals or unexpected technical factors.

### Workflow

1. **Unsupervised structure detection**: UMAP projection, Ward clustering, and Voronoi visualization identify anomalously projected samples and potential subgroups.

2. **Supervised hypothesis testing**: Random forest classification tests whether data structure supports the biological hypothesis or alternative technical groupings (batch identifiers, processing dates, etc.).

## Installation

Download this repository. All required R packages are automatically installed and loaded by the example scripts.

## Code Organization

Example scripts (`*_run.R`) → Core functions (`umap_ward_misclassification_analysis*.R`) → Utility functions (`prepare_dataset.R`, `perform_*.R`, `plot_*.R`)

### Core Functions

- **`umap_ward_misclassification_analysis()`**: Main framework function combining UMAP projection, Ward clustering, and supervised classification. Returns anomaly metrics and flagged samples.
- **`perform_supervised_classification()`**: Tests alternative hypotheses using random forest and SVM with hyperparameter tuning and variable importance extraction.
- **`perform_umap_projection()`, `perform_ward_clustering()`, `plot_umap_with_voronoi()`, `plot_misclassification_heatmap()`**: Individual analysis steps.

All functions automatically load their required libraries.

## Quick Start

```r
# Load data
lipid_profiles <- read.csv("lipid_profiles.csv")
sample_metadata <- read.csv("sample_metadata.csv")

# Run analysis
results <- umap_ward_misclassification_analysis(
  data = lipid_profiles,
  target = sample_metadata$SampleType,
  output_dir = "qc_results",
  file_prefix = "analysis"
)
```

See `example_run.R` for a complete working example.

## Example Scripts

- **`lipid_case_data_run.R`**: Demonstrates detection of a data handover error (uncorrected raw data provided instead of batch-corrected data) in real lipidomics data that passed standard laboratory QC.
- **`lipid_validation_data_run.R`**: Validates framework sensitivity to batch effects of varying magnitude.
- **`example_run.R`**: Simplified introductory example (recommended starting point).

## Dependencies

All required packages are automatically installed and loaded by the code.

**Key packages**: `umap`, `ggplot2`, `ggrepel`, `deldir`, `caret`, `randomForest`, `parallel`, `pbmcapply`, `cABCanalysis`

## License

CC-BY 4.0

## Citation

If you use this code, concept, or framework in your work, please cite:

Lötsch J, Geisslinger G, Himmelspach A, and Kringel D. 
Investigative data science addresses a structural omics quality control gap across the sample-to-dataset research chain. Brief Bioinform 2026 (in revision).


