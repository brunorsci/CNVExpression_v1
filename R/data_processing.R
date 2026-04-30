#' @title Data Processing Functions for CNVExpression
#' @description Functions for loading and preprocessing CNV and expression data
#' @name data_processing
NULL

#' Load CNV Segment Data from File
#'
#' Load Copy Number Variation segment data from a tab-delimited file.
#' The file must be in segment format: one row per genomic segment, with
#' columns for sample ID, chromosome, start/end positions, and segment mean
#' (log2 ratio). Column names can be customised to match any source
#' (TCGA/GDC, UCSC Xena, GATK, DNAcopy, Sequenza, etc.).
#'
#' @param file Path to the CNV segment file (plain text or gzipped).
#' @param sample_col Name of the column containing sample identifiers.
#'   Default \code{"sample"}.
#' @param chrom_col Name of the column containing chromosome labels
#'   (e.g. \code{"chr1"} or \code{"1"} are both accepted downstream).
#'   Default \code{"Chrom"}.
#' @param start_col Name of the column with segment start coordinates (1-based).
#'   Default \code{"Start"}.
#' @param end_col Name of the column with segment end coordinates (inclusive).
#'   Default \code{"End"}.
#' @param value_col Name of the column with segment mean log2-ratio values.
#'   Default \code{"value"}.
#'
#' @return A data.frame with five standardised columns:
#'   \code{sample}, \code{chrom}, \code{start}, \code{end}, \code{segment_mean}.
#' @export
#'
#' @examples
#' \dontrun{
#' # UCSC Xena / DNAcopy format (chr, start, end, value)
#' cnv <- load_cnv_data("SNP6_genomicSegment.tsv.gz",
#'                      chrom_col = "chr", start_col = "start", end_col = "end")
#'
#' # GATK ModelSegments format (CONTIG, START, END, LOG2_COPY_RATIO_POSTERIOR_50)
#' cnv <- load_cnv_data("tumor.cr.seg",
#'                      chrom_col = "CONTIG",
#'                      start_col = "START",
#'                      end_col   = "END",
#'                      value_col = "LOG2_COPY_RATIO_POSTERIOR_50")
#'
#' # Sequenza / custom format
#' cnv <- load_cnv_data("sequenza_segments.txt",
#'                      sample_col = "ID",
#'                      chrom_col  = "chromosome",
#'                      start_col  = "start.pos",
#'                      end_col    = "end.pos",
#'                      value_col  = "CNt")
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
#' Load RNA-seq expression data from a tab-delimited file. The file must have
#' genes as rows and samples as columns, with one leading column containing
#' gene identifiers (Ensembl IDs, gene symbols, or any other identifier).
#'
#' @param file Path to the expression file (plain text or gzipped).
#' @param gene_col Name of the column containing gene identifiers.
#'   Default \code{"Ensembl_ID"}. For UCSC Xena HiSeqV2 files use
#'   \code{"sample"}; for GDC STAR counts use \code{"Ensembl_ID"}.
#' @param remove_par_y Logical. Remove PAR_Y-suffixed Ensembl IDs (GDC-specific
#'   artefact). Set \code{FALSE} when gene IDs are not Ensembl or data do not
#'   come from GDC. Default \code{FALSE}.
#' @param remove_version Logical. Strip the version suffix from Ensembl IDs
#'   (e.g. \code{ENSG00000141510.18} → \code{ENSG00000141510}). Set
#'   \code{FALSE} when using gene symbols or non-versioned IDs. Default
#'   \code{FALSE}.
#'
#' @return A data.frame with gene identifiers as row names and samples as
#'   columns.
#' @export
#'
#' @examples
#' \dontrun{
#' # GDC STAR counts (versioned Ensembl IDs, PAR_Y entries present)
#' expr <- load_expression_data("TCGA-LAML.star_counts.tsv.gz",
#'                              gene_col       = "Ensembl_ID",
#'                              remove_par_y   = TRUE,
#'                              remove_version = TRUE)
#'
#' # UCSC Xena HiSeqV2 (gene symbols, no version suffix)
#' expr <- load_expression_data("HiSeqV2.gz", gene_col = "sample")
#'
#' # Generic count matrix (first column = gene symbol)
#' expr <- load_expression_data("counts.tsv", gene_col = "GeneSymbol")
#' }
load_expression_data <- function(file,
                                  gene_col = "Ensembl_ID",
                                  remove_par_y = FALSE,
                                  remove_version = FALSE) {

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
#' Normalize sample identifiers so that CNV and expression datasets can be
#' matched regardless of how the identifiers were originally formatted.
#' Both hyphens (\code{-}) and dots (\code{.}) are treated as separators,
#' which handles the common case where R silently converts \code{TCGA-AB-1234-01}
#' to \code{TCGA.AB.1234.01} in data frame column names.
#'
#' @param samples Character vector of sample identifiers.
#' @param n_parts Number of identifier parts to keep. Default \code{4}
#'   (e.g. keeps \code{TCGA-AB-1234-01} from a longer barcode). Set
#'   \code{Inf} to keep the full identifier.
#' @param output_sep Separator used in the output. Default \code{"_"}.
#'
#' @return Character vector of standardized sample identifiers.
#' @export
#'
#' @examples
#' # TCGA barcodes with hyphens
#' standardize_sample_ids(c("TCGA-AB-1234-01A-11R-A12B-07"))
#' # R-mangled dots treated identically
#' standardize_sample_ids(c("TCGA.AB.1234.01"))
#' # Non-TCGA IDs are returned unchanged (no hyphen/dot)
#' standardize_sample_ids(c("Sample1", "Sample2"))
standardize_sample_ids <- function(samples,
                                    n_parts = 4,
                                    output_sep = "_") {

    standardized <- vapply(samples, function(x) {
        # Split on either hyphen or dot
        parts <- unlist(strsplit(x, "[-.]"))
        if (length(parts) >= n_parts && is.finite(n_parts)) {
            paste(parts[seq_len(n_parts)], collapse = output_sep)
        } else if (length(parts) > 1) {
            paste(parts, collapse = output_sep)
        } else {
            x  # no separator found — return as-is
        }
    }, character(1), USE.NAMES = FALSE)

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
