# Investigative QC Omics

## Computational investigative quality control detects hidden technical artifacts 
## that standard omics batch correction misses

This repository provides an **investigative quality control (QC) framework for omics data** 
designed to detect hidden laboratory artifacts and batch effects that may generate data 
structure indistinguishable from genuine biological signals. The method combines UMAP 
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
Core Analysis Functions (umap_ward_investigative_qc*.R)
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
- **`plot_anomaly_heatmap.R`**: Anomalously projected sample heatmap 
  (requires: `ggplot2`, `tidyr`, `scales`)
- **`plot_classification_stability_heatmap.R`**: Classification stability 
  visualization (requires: `ggplot2`, `tidyr`, `scales`, `grid`)

### Core Method Functions

#### `umap_ward_investigative_qc()` — **Main QC Framework Function**

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

#### `perform_supervised_classification()` — **Artifact Detection via Random Forest**

**Description**: Tests multiple alternative hypotheses using random forest classifiers 
with hyperparameter tuning, comparing balanced accuracy across biological targets versus 
technical groupings including batch identifiers, processing dates, and projection-derived 
clusters.

**Key Parameters**:
- `X`: Predictor features
- `Y`: Target variables to test (alternative hypotheses)
- `n_iter`: Resampling iterations (100 in the paper)
- `training_size`: Train/test ratio (default: 0.67)

**Returns**: Median balanced accuracy with 95% confidence intervals for each 
hypothesis tested.

**Batch validation**: Tests batch detection sensitivity across graded technical 
variation (10% to 100% CV).

### Supplementary Functions

**`umap_ward_investigative_qc_multi()`**: Optional multi-target extension.  
**`plot_classification_stability_heatmap()`**: Optional visualization of 
classification stability across resampling iterations.

## Usage Examples

### Basic Usage with Sample Data
```r
# Load the data
lipid_profiles <- read.csv("lipid_profiles.csv")
sample_metadata <- read.csv("sample_metadata.csv")
sample_types <- sample_metadata$SampleType

# Run the analysis
results <- umap_ward_investigative_qc(
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

### Publication Case Study

Case study on real lipidomics data (psoriatic arthritis patients versus controls) 
demonstrating detection of a hidden batch effect of extreme magnitude that survived 
all standard quality control procedures. Shows:

- UMAP projection with Voronoi tessellation visualization
- Ward hierarchical clustering
- Identification of anomalously projected samples
- Random forest supervised classification testing alternative hypotheses including 
  study groups, projection-derived clusters, batch identifiers, and processing dates
- How supervised classifiers revealed that batch effects rather than biological 
  differences dominated the data structure
- Variable importance analysis with recursive cABC selection
- Hyperparameter tuning for random forest classifiers (mtry, ntree)

This script reproduces the primary case study described in the paper.

### Method Validation

**`lipid_validation_data_run.R`**: Full validation pipeline on a batch-effect-free 
reference dataset of 403 healthy controls across 39 serum lipid mediators from a 
multiple sclerosis biomarker study. Tests framework sensitivity across five batch 
effect magnitudes (original, weak 10% CV, medium 20% CV, strong 30% CV, extreme 
100% CV) using simulated technical batches combining batch-day drift and run-order 
signal suppression. Demonstrates:

- Unsupervised UMAP and Ward clustering batch detection index (40.2% to 24.1%)
- Supervised random forest balanced accuracy scaling from 0.51 to 0.86
- Realistic batch simulation matching the range of pre-normalization technical 
  errors reported in large-scale lipidomics studies (17.5 to 34.1% RSD; 
  Fan et al. 2019)

### Demonstration Examples

**`golfball_dataset_run.R`**: Artificial dataset with concentric spherical structure 
demonstrating the visualization approach on synthetic data with fully known structure. 
Useful for understanding framework behavior before applying to real data. Output 
includes a combined visualization showing:

- Left panel: UMAP projection with Voronoi tessellation colored by cluster assignments 
  (for Voronoi tessellation plot, cite Lotsch J and Kringel D (2026). Voronoi 
  tessellation as a complement or replacement for confidence ellipses in the 
  visualization of data projection and clustering results. PLoS One (in revision))
- Right panel: Anomalously projected sample heatmap comparing prior class assignments 
  with Ward-assigned clusters

![Golfball UMAP Analysis Results](golfballs_combined_plot.svg)

**`example_run.R`**: Simplified introductory example showing the basic workflow with 
minimal configuration. Recommended starting point before applying to real data.

**`create_sample_files.R`**: Generates synthetic lipidomics data for testing and 
learning purposes.

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

Lotsch J, Hahnefeld L, Geisslinger G, Himmelspach A, and Kringel D. 
Computational investigative quality control detects hidden technical artifacts 
that standard omics batch correction misses. 2026 (in preparation).


