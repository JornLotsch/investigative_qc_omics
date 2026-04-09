# Investigative QC Omics

## A computational investigative quality control framework for omics detecting hidden batch artifacts beyond the analytical laboratory

This repository provides R code used to develop and validate an 
**investigative quality control (QC) framework for omics data** 
designed to detect hidden laboratory artifacts and batch effects that may generate data 
structure indistinguishable from genuine biological signals. 

The method combines UMAP 
projection, Ward's hierarchical clustering, Voronoi cell visualization, and supervised 
random forest classification to systematically identify anomalously projected samples and 
test whether observed data structure reflects the study hypothesis or technical workflow 
metadata.

### Core Framework

The quality control workflow proceeds in two sequential steps:

1. **Unsupervised structure detection**: UMAP projection combined with Ward clustering 
and Voronoi visualization to identify potential subgroups and anomalously projected 
samples whose cluster assignment contradicts their known class label.

2. **Supervised artifact detection**: Random forest classification across multiple 
alternative class structures, including batch identifiers, processing dates, and 
projection-derived clusters, to test whether data structure supports the biological 
hypothesis or reveals technical confounds. When supervised classifiers achieve higher 
balanced accuracy for workflow-related groupings than for the study hypothesis, this 
signals that technical artifacts may dominate the biological signal, requiring 
investigation before biological interpretation.

### Intended Use

This framework complements existing omics QC pipelines by providing **investigative 
quality control** when standard QC, statistical analysis, and machine learning all appear 
to validate a biological hypothesis. It is particularly useful for:

- Identifying batch effects that escape standard QC procedures, including those of 
  pre-analytical origin outside the analytical laboratory's direct control
- Flagging samples that project anomalously relative to their known class label
- Testing multiple alternative technical explanations for observed data structure
- Validating whether apparent biological signals are robust or driven by technical 
  confounds

**Important**: This framework is a complement to standard data preprocessing and quality 
control, not a replacement for it. Method choice should remain flexible and guided by 
the researcher's knowledge of the data structure under investigation. The framework 
prescribes a principle rather than a fixed pipeline: test the data against alternative 
hypotheses and treat anomalously projected samples as diagnostic signals rather than 
inconveniences.

## Installation

Download this repository to the local hard drive. All required R packages will be 
automatically installed and loaded when running the example scripts.

## Package Architecture

The code is organized in a modular hierarchy:
```
Example Scripts (*_run.R)
    ↓
Core Analysis Functions (umap_ward_misclassification_analysis*.R)
    ↓
Utility Functions (prepare_dataset.R, perform_*.R, plot_*.R)
```

### Utility Functions and Libraries

The core computational functions automatically load their required libraries:

- **`prepare_dataset.R`**: Data standardization and preprocessing
- **`perform_umap_projection.R`**: UMAP dimensionality reduction
  (requires: `umap`)
- **`perform_ward_clustering.R`**: Hierarchical clustering with optimization
  (requires: `stats::hclust`, `clue::solve_LSAP`)
- **`plot_umap_with_voronoi.R`**: Voronoi tessellation visualization
  (requires: `ggplot2`, `ggrepel`, `deldir`)
- **`plot_misclassification_heatmap.R`**: Anomalously projected sample heatmap
  (requires: `ggplot2`, `tidyr`, `scales`)
- **`perform_supervised_classification.R`**: Supervised classification with random forest
  and SVM, including hyperparameter tuning and variable importance extraction
  (requires: `caret`, `randomForest`, `parallel`, `pbmcapply`)
- **`plot_classification_stability_heatmap.R`**: Classification stability
  visualization (requires: `ggplot2`, `tidyr`, `scales`, `grid`)

### Core Method Functions

#### `umap_ward_misclassification_analysis()` — **Main QC Framework Function**

**Description**: Combines UMAP projection, Ward hierarchical clustering, Voronoi
visualization, and anomaly detection to identify samples whose cluster assignment
contradicts their known class label and data structure inconsistent with the study
hypothesis.

**Key Parameters**:
- `data`: Input data (samples × features)
- `target`: Target class labels for comparison
- `output_dir`: Directory for output files
- `file_format`: "svg" or "png"
- `n_neighbors`: UMAP parameter (default: 15)

**Output Files**:
- `{file_prefix}_voronoi.{format}`: UMAP projection with Voronoi visualization
- `{file_prefix}_heatmap.{format}`: Anomalously projected sample heatmap
- `{file_prefix}_combined.{format}`: Combined visualization
- `{file_prefix}_anomalous_samples.csv`: Anomalously projected samples, if any

#### `perform_supervised_classification()` — **Supervised Artifact Detection**

**Description**: Tests multiple alternative hypotheses using random forest and SVM
classifiers with hyperparameter tuning and variable importance extraction, comparing
balanced accuracy across biological targets versus technical groupings including batch
identifiers, processing dates, and projection-derived clusters. SVM may fail to
converge for datasets with very small class sizes.

**Key Parameters**:
- `X`: Predictor features
- `Y`: Target variables to test (alternative hypotheses)
- `n_iter`: Resampling iterations
- `training_size`: Train/test ratio (default: 0.67)

**Returns**: Median balanced accuracy with 95% confidence intervals and variable
importance scores for each hypothesis tested.

**Feature Importance**: Random forest permutation importance scores are extracted across
resampling iterations. Recursive cABC analysis is applied to median importance values to
identify the A-set of most informative features for each classification target.

### Supplementary Functions

**`umap_ward_misclassification_analysis_multi()`**: Multi-target extension for
analyzing multiple metadata variables simultaneously.
**`plot_classification_stability_heatmap()`**: Visualization of classification
stability across resampling iterations.

## Usage Examples

### Basic Usage with Sample Data
```r
# Load the data
lipid_profiles <- read.csv("lipid_profiles.csv")
sample_metadata <- read.csv("sample_metadata.csv")
sample_types <- sample_metadata$SampleType

# Run the analysis
results <- umap_ward_misclassification_analysis(
  data = lipid_profiles,
  target = sample_types,
  labels = sample_metadata$SampleID,
  output_dir = "qc_results",
  file_prefix = "lipid_analysis",
  file_format = "png",
  width = 14,
  height = 10
)

# Check anomalously projected samples
cat("Anomalously projected samples:", 
    sprintf("%.2f%%", results$anomaly_rate * 100), "\n")

# View anomalous samples
if (nrow(results$anomalous_samples) > 0) {
  print(results$anomalous_samples)
}
```

## Example Scripts

This repository includes example R scripts implementing the analyses described 
in the paper.

**`lipid_case_data_run.R`**: Case study on real lipidomics data (psoriatic arthritis
patients versus controls) demonstrating detection of a hidden batch effect that survived
all standard quality control procedures. Performs supervised classification across three
sample stratifications (all samples, controls only, patients only) using both real and
permuted data (6 combinations total), testing whether technical artifacts are consistent
across subgroups or specific to the full cohort. Reproduces the primary case study
described in the paper.

**`lipid_validation_data_run.R`**: Validation pipeline testing framework sensitivity
to batch effects of graded magnitude using simulated technical batches on a
batch-effect-free reference lipidomics dataset.

**`example_run.R`**: Simplified introductory example showing the basic workflow.
Recommended starting point before applying to real data.

## Dependencies

### Core Analysis
All required packages are automatically installed and loaded:
- **Visualization**: `ggplot2`, `ggrepel`, `deldir` (Voronoi), `gridExtra`
- **Analysis**: `umap`, `clue` (Hungarian algorithm)
- **Utilities**: `tidyr`, `scales`, `grid`

### Supervised Classification
- `caret`, `randomForest`, `parallel`, `pbmcapply`, `caTools`, `dplyr`
- Hyperparameter tuning with resampling-based grids (mtry, ntree)

### Demonstrations
- `ComplexHeatmap`, `cowplot`, `plot3D`, `car`, `ggplotify`

### Variable Importance and Feature Selection
- `cABCanalysis`, `twosamples`, `matrixStats` (recursive ABC analysis)

## System Requirements

- **R version**: 4.0 or higher (4.3 or higher recommended)
- **Operating system**: Linux or Unix (macOS supported); not tested on Windows
- **RAM**: Minimum 4 GB (8 GB or more recommended for large datasets)
- **Processing**: Parallel processing recommended; 4 to 8 CPU cores utilized

## License

CC-BY 4.0

## Citation

If you use this code, concept, or framework in your work, please cite:

Lötsch J, Hahnefeld L, ..., Geisslinger G, Himmelspach A, and Kringel D. 
Metadata-driven machine learning exposes preanalytical batch arti-facts that mimic biological discovery in lipidomics. 2026 (in preparation).


