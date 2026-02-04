#' @title Utility Functions for CNVExpression
#' @description Helper functions for CNVExpression package
#' @name utils
#' @importFrom rlang .data
#' @importFrom stats cor.test median p.adjust
#' @importFrom methods is
NULL

#' Export Results to TSV
#'
#' Export analysis results to a tab-separated file.
#'
#' @param data Data frame to export
#' @param file Output file path
#' @param ... Additional arguments passed to write.table
#'
#' @return Invisibly returns the file path
#' @export
#'
#' @examples
#' \dontrun{
#' export_results(results, "output/cnv_expression_results.tsv")
#' }
export_results <- function(data, file, ...) {

    # Create directory if needed
    dir_path <- dirname(file)
    if (!dir.exists(dir_path) && dir_path != ".") {
        dir.create(dir_path, recursive = TRUE)
    }

    utils::write.table(
        data,
        file = file,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE,
        ...
    )

    message("Results exported to: ", file)
    invisible(file)
}

#' Check Package Dependencies
#'
#' Check if required packages are installed.
#'
#' @param packages Character vector of package names
#' @param stop_on_missing Logical, stop execution if packages are missing
#'
#' @return Logical indicating if all packages are available
#' @export
#'
#' @examples
#' check_dependencies(c("GenomicRanges", "dplyr"))
check_dependencies <- function(packages, stop_on_missing = FALSE) {

    missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]

    if (length(missing) > 0) {
        msg <- paste0(
            "Missing packages: ", paste(missing, collapse = ", "), "\n",
            "Install with: BiocManager::install(c('",
            paste(missing, collapse = "', '"), "'))"
        )

        if (stop_on_missing) {
            stop(msg)
        } else {
            warning(msg)
            return(FALSE)
        }
    }

    return(TRUE)
}

#' Validate Input Data
#'
#' Validate input data frames for required columns and data types.
#'
#' @param data Data frame to validate
#' @param required_cols Character vector of required column names
#' @param numeric_cols Character vector of columns that should be numeric
#' @param data_name Name of the data for error messages
#'
#' @return Logical TRUE if valid, stops with error otherwise
#' @export
#'
#' @examples
#' \dontrun{
#' validate_input_data(cnv_data, c("sample", "chrom", "start", "end"))
#' }
validate_input_data <- function(data,
                                 required_cols = NULL,
                                 numeric_cols = NULL,
                                 data_name = "data") {

    if (!is.data.frame(data)) {
        stop(data_name, " must be a data frame")
    }

    if (nrow(data) == 0) {
        stop(data_name, " has no rows")
    }

    # Check required columns
    if (!is.null(required_cols)) {
        missing <- setdiff(required_cols, colnames(data))
        if (length(missing) > 0) {
            stop(data_name, " is missing columns: ", paste(missing, collapse = ", "))
        }
    }

    # Check numeric columns
    if (!is.null(numeric_cols)) {
        for (col in numeric_cols) {
            if (col %in% colnames(data) && !is.numeric(data[[col]])) {
                stop("Column '", col, "' in ", data_name, " must be numeric")
            }
        }
    }

    return(TRUE)
}

#' Print Summary Statistics
#'
#' Print formatted summary statistics for analysis results.
#'
#' @param results List of results from analysis pipeline
#'
#' @return Invisibly returns NULL
#' @export
#'
#' @examples
#' \dontrun{
#' print_summary(analysis_results)
#' }
print_summary <- function(results) {

    cat("\n========================================\n")
    cat("CNVExpression Analysis Summary\n")
    cat("========================================\n\n")

    # CNV data summary
    if ("cnv_processed" %in% names(results)) {
        cat("CNV Data:\n")
        cat("  - Total CNV calls: ", nrow(results$cnv_processed), "\n")
        cat("  - Samples: ", length(unique(results$cnv_processed$sample)), "\n")
        if ("annotation" %in% colnames(results$cnv_processed)) {
            cat("  - By type:\n")
            type_table <- table(results$cnv_processed$annotation)
            for (i in seq_along(type_table)) {
                cat("      ", names(type_table)[i], ": ", type_table[i], "\n")
            }
        }
        cat("\n")
    }

    # CNV regions
    if ("cnv_regions" %in% names(results)) {
        cat("Population CNV Regions: ", length(results$cnv_regions), "\n\n")
    }

    # Expression data
    if ("expression_se" %in% names(results)) {
        cat("Expression Data:\n")
        cat("  - Genes: ", nrow(results$expression_se), "\n")
        cat("  - Samples: ", ncol(results$expression_se), "\n\n")
    }

    # Results
    if ("results_raw" %in% names(results)) {
        cat("CNV-Expression Analysis:\n")
        cat("  - Total associations: ", nrow(results$results_raw), "\n")
    }

    if ("results_filtered" %in% names(results)) {
        cat("  - Significant (after filtering): ", nrow(results$results_filtered), "\n")
    }

    if ("results_merged" %in% names(results)) {
        cat("  - Final merged results: ", nrow(results$results_merged), "\n")
    }

    cat("\n========================================\n")

    invisible(NULL)
}

#' Create Sample-Gene Identifier
#'
#' Create a unique identifier combining sample and gene information.
#'
#' @param sample_id Sample identifier
#' @param gene_id Gene identifier
#' @param sep Separator character (default: "_")
#'
#' @return Character vector of combined identifiers
#' @export
#'
#' @examples
#' ids <- create_sample_gene_id(c("S1", "S2"), c("TP53", "BRCA1"))
create_sample_gene_id <- function(sample_id, gene_id, sep = "_") {
    paste(sample_id, gene_id, sep = sep)
}

#' Get Package Version
#'
#' Get the current version of CNVExpression package.
#'
#' @return Character string with package version
#' @export
#'
#' @examples
#' get_cnvexpression_version()
get_cnvexpression_version <- function() {
    as.character(utils::packageVersion("CNVExpression"))
}

#' CNVExpression Package
#'
#' CNVExpression: Integration of Copy Number Variations with Gene Expression
#' Analysis
#'
#' This package provides a comprehensive workflow for integrating Copy Number
#' Variation (CNV) data with gene expression data. It enables identification
#' of genes whose expression levels are significantly affected by copy number
#' alterations.
#'
#' @section Main functions:
#' \itemize{
#'   \item \code{\link{cnv_expression_pipeline}}: Complete analysis pipeline
#'   \item \code{\link{prepare_cnv_data}}: Preprocess CNV data
#'   \item \code{\link{run_cnv_expression_analysis}}: Run CNV-expression analysis
#'   \item \code{\link{filter_cnv_expression_results}}: Filter results
#'   \item \code{\link{calculate_cnv_expression_correlation}}: Correlation analysis
#' }
#'
#' @docType package
#' @name CNVExpression-package
#' @aliases CNVExpression
#'
#' @importFrom dplyr %>% filter select mutate group_by summarise case_when n
#' @importFrom tidyr separate
#' @importFrom stringr str_replace_all
#' @importFrom GenomicRanges GRanges makeGRangesListFromDataFrame findOverlaps
#'   subsetByOverlaps mergeByOverlaps sort
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors mcols mcols<- queryHits subjectHits isEmpty
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom CNVRanger populationRanges cnvEQTL
#' @importFrom BiocGenerics start end width strand
NULL
