#' @title Data Processing Functions for CNVExpression
#' @description Functions for loading and preprocessing CNV and expression data
#' @name data_processing
NULL

#' Load CNV Data from File
#'
#' Load Copy Number Variation data from a tab-delimited file (e.g., from TCGA/Xena).
#'
#' @param file Path to the CNV file (can be gzipped)
#' @param sample_col Name of the column containing sample identifiers
#' @param chrom_col Name of the column containing chromosome information
#' @param start_col Name of the column containing start position
#' @param end_col Name of the column containing end position
#' @param value_col Name of the column containing segment mean values
#'
#' @return A data.frame with standardized column names
#' @export
#'
#' @examples
#' \dontrun{
#' cnv_data <- load_cnv_data("TCGA-LAML.masked_cnv_DNAcopy.tsv.gz")
#' }
load_cnv_data <- function(file,
                          sample_col = "sample",
                          chrom_col = "Chrom",
                          start_col = "Start",
                          end_col = "End",
                          value_col = "value") {

    if (!file.exists(file)) {
        stop("File not found: ", file)
    }

    cnv <- utils::read.delim(file, stringsAsFactors = FALSE)

    # Standardize column names
    col_mapping <- c(
        sample = sample_col,
        chrom = chrom_col,
        start = start_col,
        end = end_col,
        segment_mean = value_col
    )

    # Check if all columns exist
    missing_cols <- setdiff(col_mapping, colnames(cnv))
    if (length(missing_cols) > 0) {
        stop("Missing columns in CNV file: ", paste(missing_cols, collapse = ", "))
    }

    # Select and rename columns
    cnv <- cnv[, col_mapping]
    colnames(cnv) <- names(col_mapping)

    return(cnv)
}

#' Load Expression Data from File
#'
#' Load RNA-seq count data from a tab-delimited file.
#'
#' @param file Path to the expression file (can be gzipped)
#' @param gene_col Name of the column containing gene identifiers (e.g., Ensembl ID)
#' @param remove_par_y Logical, whether to remove PAR_Y genes
#' @param remove_version Logical, whether to remove version suffix from gene IDs
#'
#' @return A data.frame with genes as rows and samples as columns
#' @export
#'
#' @examples
#' \dontrun{
#' expr_data <- load_expression_data("TCGA-LAML.star_counts.tsv.gz")
#' }
load_expression_data <- function(file,
                                  gene_col = "Ensembl_ID",
                                  remove_par_y = TRUE,
                                  remove_version = TRUE) {

    if (!file.exists(file)) {
        stop("File not found: ", file)
    }

    counts <- utils::read.delim(file, stringsAsFactors = FALSE)

    # Check if gene column exists
    if (!gene_col %in% colnames(counts)) {
        stop("Gene column '", gene_col, "' not found in expression file")
    }

    # Set rownames
    rownames(counts) <- counts[[gene_col]]

    # Remove PAR_Y genes if requested
    if (remove_par_y) {
        counts <- counts[!grepl("PAR_Y", rownames(counts)), ]
    }

    # Remove version suffix from Ensembl IDs
    if (remove_version) {
        rownames(counts) <- gsub("\\..*", "", rownames(counts))
    }

    # Remove gene column
    counts <- counts[, !colnames(counts) %in% gene_col, drop = FALSE]

    return(counts)
}

#' Standardize Sample Identifiers
#'
#' Standardize TCGA sample identifiers by replacing hyphens with underscores
#' and optionally truncating to a specific number of parts.
#'
#' @param samples Character vector of sample identifiers
#' @param sep Separator character (default: "-")
#' @param n_parts Number of parts to keep (default: 4, e.g., "TCGA-AB-1234-01")
#' @param output_sep Separator for output (default: "_")
#'
#' @return Character vector with standardized sample identifiers
#' @export
#'
#' @examples
#' samples <- c("TCGA-AB-1234-01A-11R-A12B-07", "TCGA-AB-5678-01A-11R-A12B-07")
#' standardize_sample_ids(samples)
standardize_sample_ids <- function(samples,
                                    sep = "-",
                                    n_parts = 4,
                                    output_sep = "_") {

    standardized <- vapply(samples, function(x) {
        parts <- unlist(strsplit(x, sep, fixed = TRUE))
        if (length(parts) >= n_parts) {
            paste(parts[seq_len(n_parts)], collapse = output_sep)
        } else {
            gsub(sep, output_sep, x, fixed = TRUE)
        }
    }, character(1), USE.NAMES = FALSE)

    # Clean up common suffixes
    standardized <- gsub("_03A$", "_03", standardized)
    standardized <- gsub("_11A$", "_11", standardized)

    return(standardized)
}

#' Match Samples Between CNV and Expression Data
#'
#' Find common samples between CNV and expression datasets.
#'
#' @param cnv_samples Character vector of CNV sample identifiers
#' @param expr_samples Character vector of expression sample identifiers
#' @param standardize Logical, whether to standardize sample IDs before matching
#'
#' @return List with matched sample identifiers and summary
#' @export
#'
#' @examples
#' \dontrun{
#' matched <- match_samples(cnv$sample, colnames(expr_data))
#' }
match_samples <- function(cnv_samples, expr_samples, standardize = TRUE) {

    if (standardize) {
        cnv_samples_std <- standardize_sample_ids(cnv_samples)
        expr_samples_std <- standardize_sample_ids(expr_samples)
    } else {
        cnv_samples_std <- cnv_samples
        expr_samples_std <- expr_samples
    }

    common_samples <- intersect(unique(cnv_samples_std), unique(expr_samples_std))

    list(
        common_samples = common_samples,
        n_common = length(common_samples),
        n_cnv_only = length(setdiff(unique(cnv_samples_std), common_samples)),
        n_expr_only = length(setdiff(unique(expr_samples_std), common_samples))
    )
}
