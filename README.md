# CNVExpression

**Integration of Copy Number Variations with Gene Expression Analysis**

[![R-CMD-check](https://github.com/brunorsci/CNVExpression/workflows/R-CMD-check/badge.svg)](https://github.com/brunorsci/CNVExpression/actions)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-CNVExpression-blue)](https://bioconductor.org/packages/CNVExpression)

## Overview

CNVExpression is an R/Bioconductor package for integrating Copy Number Variation (CNV) data with gene expression (RNA-seq) data. The package enables identification of genes whose expression levels are significantly affected by copy number alterations.

### Key Features

- Complete pipeline for CNV-expression integration analysis
- Support for TCGA and Xena data formats
- Integration with CNVRanger for population-level CNV analysis
- Multiple correlation methods (Spearman, Pearson, Kendall)
- Comprehensive gene annotation (Ensembl, Entrez, gene symbols)
- Publication-ready visualizations

## Installation

### From Bioconductor (when available)

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("CNVExpression")
```

### From GitHub (development version)

```r
devtools::install_github("brunorsci/CNVExpression")
```

## Quick Start

```r
library(CNVExpression)

# Load data
cnv_data <- load_cnv_data("TCGA-LAML.masked_cnv_DNAcopy.tsv.gz")
expr_data <- load_expression_data("TCGA-LAML.star_counts.tsv.gz")

# Run complete pipeline
results <- cnv_expression_pipeline(
    cnv_data = cnv_data,
    expression_data = expr_data,
    genome = "GRCh38",
    logfc_cutoff = 1,
    pvalue_cutoff = 0.05
)

# View summary
print_summary(results)

# Access filtered results
significant_genes <- results$results_filtered
```

## Detailed Workflow

### Step 1: Load and Preprocess CNV Data

```r
# Load CNV data
cnv_raw <- load_cnv_data("TCGA-LAML.masked_cnv_DNAcopy.tsv.gz")

# Preprocess: convert segment mean to copy number, filter diploid, remove sex chromosomes
cnv_processed <- prepare_cnv_data(cnv_raw)
```

### Step 2: Create GRanges Objects

```r
# Create GRangesList for CNV data
cnv_grl <- create_cnv_grangeslist(cnv_processed)

# Calculate population CNV regions
cnvrs <- calculate_population_ranges(cnv_grl)
```

### Step 3: Prepare Expression Data

```r
# Load and prepare expression data
counts <- load_expression_data("TCGA-LAML.star_counts.tsv.gz")
genes_gr <- create_gene_granges(rownames(counts))
rse <- create_expression_se(counts, genes_gr)
```

### Step 4: Run CNV-Expression Analysis

```r
# Run analysis
results <- run_cnv_expression_analysis(cnvrs, cnv_grl, rse)

# Filter significant results
filtered <- filter_cnv_expression_results(results, logfc_cutoff = 1, pvalue_cutoff = 0.05)
```

### Step 5: Annotate and Visualize

```r
# Add gene symbols
annotated <- add_gene_symbols(filtered)

# Generate plots
plot_volcano(filtered)
plot_cnv_by_chromosome(cnv_processed)
```

## Main Functions

| Function | Description |
|----------|-------------|
| `cnv_expression_pipeline()` | Complete analysis pipeline |
| `prepare_cnv_data()` | Preprocess CNV data |
| `create_cnv_grangeslist()` | Create GRangesList from CNV data |
| `run_cnv_expression_analysis()` | Run CNV-expression eQTL analysis |
| `filter_cnv_expression_results()` | Filter significant results |
| `calculate_cnv_expression_correlation()` | Correlation analysis |
| `add_gene_symbols()` | Gene annotation |

## Data Requirements

### CNV Data
- Tab-delimited file with columns: sample, chromosome, start, end, segment_mean
- Supported formats: TCGA DNACopy, Xena CNV files

### Expression Data
- Tab-delimited file with genes as rows, samples as columns
- Gene identifiers: Ensembl IDs or gene symbols
- Supported: Raw counts or normalized expression

## Output

The main output is a data frame with:
- Gene location and identifiers
- Log fold changes for each copy number state (CN0, CN1, CN3, CN4)
- P-values (raw and adjusted)
- Gene annotations (symbols, Entrez IDs)

## Citation

If you use CNVExpression in your research, please cite:

> Assuncao BR (2025). CNVExpression: Integration of Copy Number Variations with
> Gene Expression Analysis. R package version 0.99.0.

## Related Packages

- [CNVRanger](https://bioconductor.org/packages/CNVRanger): Association analysis of CNVs
- [GenomicRanges](https://bioconductor.org/packages/GenomicRanges): Genomic data structures
- [SummarizedExperiment](https://bioconductor.org/packages/SummarizedExperiment): Expression data containers

## License

This package is released under the Artistic-2.0 License.

## Contact

- Issues: https://github.com/brunorsci/CNVExpression/issues
- Author: Bruno Rodrigo Assuncao
