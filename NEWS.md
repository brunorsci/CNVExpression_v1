# CNVExpression News

## Version 0.99.0 (2025-01-29)

### New Features

* Initial submission to Bioconductor
* Complete pipeline for CNV-expression integration analysis via `cnv_expression_pipeline()`
* Support for TCGA and Xena data formats
* Integration with CNVRanger for population-level CNV analysis

### Data Processing Functions

* `load_cnv_data()` - Load CNV data from tab-delimited files
* `load_expression_data()` - Load RNA-seq count data
* `standardize_sample_ids()` - Standardize TCGA sample identifiers
* `match_samples()` - Match samples between CNV and expression datasets

### CNV Analysis Functions

* `convert_segment_mean_to_cn()` - Convert segment mean to integer copy numbers
* `prepare_cnv_data()` - Complete CNV preprocessing pipeline
* `create_cnv_grangeslist()` - Create GRangesList from CNV data
* `calculate_population_ranges()` - Identify recurrent CNV regions
* `remove_sex_chromosomes()` - Filter sex chromosomes
* `filter_diploid()` - Remove diploid (CN=2) entries

### Integration Functions

* `run_cnv_expression_analysis()` - Run CNV-expression eQTL analysis
* `filter_cnv_expression_results()` - Filter by significance and effect size
* `process_integration_results()` - Merge results with CNV annotations
* `create_gene_granges()` - Create gene GRanges using biomaRt
* `create_expression_se()` - Create SummarizedExperiment

### Correlation Functions

* `calculate_cnv_expression_correlation()` - Spearman/Pearson/Kendall correlation
* `prepare_correlation_data()` - Prepare data for correlation analysis
* `filter_significant_correlations()` - Filter significant correlations
* `classify_correlation_direction()` - Classify positive/negative correlations

### Annotation Functions

* `annotate_cnv_with_genes()` - Annotate CNVs with overlapping genes
* `add_gene_symbols()` - Add gene symbol annotations
* `entrez_to_symbol()` - Convert Entrez IDs to gene symbols
* `ensembl_to_symbol()` - Convert Ensembl IDs to gene symbols
* `annotate_with_biomart()` - Comprehensive annotation via biomaRt

### Visualization Functions

* `plot_cnv_types()` - Bar plot of CNV type distribution
* `plot_cnv_by_chromosome()` - CNV distribution by chromosome
* `plot_cnv_length()` - CNV length distribution
* `plot_volcano()` - Volcano plot for CNV-expression results
* `plot_correlation_scatter()` - Scatter plot for gene correlation
* `plot_cnv_per_sample()` - CNV counts per sample
* `generate_report_plots()` - Generate standard analysis plots

### Documentation

* Comprehensive vignette with TCGA-LAML example workflow
* roxygen2 documentation for all exported functions
* README with installation and quick start guide
