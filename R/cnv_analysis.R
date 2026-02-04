#' @title CNV Analysis Functions
#' @description Functions for processing and analyzing Copy Number Variation data
#' @name cnv_analysis
NULL

#' Convert Segment Mean to Copy Number
#'
#' Convert segment mean values to integer copy numbers using the formula:
#' CN = round(2 * 2^segment_mean). Also adds annotation based on copy number.
#'
#' @param df Data frame containing CNV data
#' @param segment_mean_col Name of the column with segment mean values
#' @param cn_col Name for the new copy number column (default: "copy_number")
#' @param annotation_col Name for the new annotation column (default: "annotation")
#'
#' @return Data frame with added copy_number and annotation columns
#' @export
#'
#' @examples
#' cnv <- data.frame(
#'     sample = "sample1",
#'     chrom = "1",
#'     start = 1000,
#'     end = 2000,
#'     segment_mean = c(-2, -0.5, 0, 0.58)
#' )
#' cnv <- convert_segment_mean_to_cn(cnv)
convert_segment_mean_to_cn <- function(df,
                                        segment_mean_col = "segment_mean",
                                        cn_col = "copy_number",
                                        annotation_col = "annotation") {

    if (!segment_mean_col %in% colnames(df)) {
        stop("Column '", segment_mean_col, "' not found in data frame")
    }

    # Calculate copy number: CN = 2 * 2^SM, rounded to nearest integer
    df[[cn_col]] <- round(2 * (2 ^ df[[segment_mean_col]]))

    # Add annotation based on copy number
    df[[annotation_col]] <- dplyr::case_when(
        df[[cn_col]] == 0 ~ "homozygous_deletion",
        df[[cn_col]] == 1 ~ "heterozygous_deletion",
        df[[cn_col]] == 2 ~ "neutral",
        df[[cn_col]] >= 3 ~ "amplification",
        TRUE ~ "unknown"
    )

    return(df)
}

#' Remove Sex Chromosomes from Data
#'
#' Remove entries from sex chromosomes (X, Y) from CNV data.
#'
#' @param df Data frame containing CNV data
#' @param chrom_col Name of the chromosome column
#' @param remove_X Logical, whether to remove chromosome X (default: TRUE)
#' @param remove_Y Logical, whether to remove chromosome Y (default: TRUE)
#'
#' @return Data frame with sex chromosomes removed
#' @export
#'
#' @examples
#' cnv <- data.frame(chrom = c("1", "2", "X", "Y"))
#' cnv <- remove_sex_chromosomes(cnv)
remove_sex_chromosomes <- function(df,
                                    chrom_col = "chrom",
                                    remove_X = TRUE,
                                    remove_Y = TRUE) {

    if (!chrom_col %in% colnames(df)) {
        stop("Column '", chrom_col, "' not found in data frame")
    }

    # Define sex chromosomes to remove
    sex_chroms <- c()
    if (remove_X) sex_chroms <- c(sex_chroms, "X", "chrX")
    if (remove_Y) sex_chroms <- c(sex_chroms, "Y", "chrY")

    df <- df[!df[[chrom_col]] %in% sex_chroms, , drop = FALSE]

    return(df)
}

#' Filter Diploid CNVs
#'
#' Remove diploid (copy number = 2) entries from CNV data.
#'
#' @param df Data frame containing CNV data with copy number information
#' @param cn_col Name of the copy number column
#'
#' @return Data frame with diploid entries removed
#' @export
#'
#' @examples
#' cnv <- data.frame(copy_number = c(0, 1, 2, 3, 4))
#' cnv <- filter_diploid(cnv)
filter_diploid <- function(df, cn_col = "copy_number") {

    if (!cn_col %in% colnames(df)) {
        stop("Column '", cn_col, "' not found in data frame")
    }

    df <- df[df[[cn_col]] != 2, , drop = FALSE]

    return(df)
}

#' Add CNV State Column
#'
#' Add a state column for CNVRanger compatibility (0-4 scale).
#' Values > 4 are capped at 4.
#'
#' @param df Data frame containing CNV data with copy number information
#' @param cn_col Name of the copy number column
#' @param state_col Name for the new state column (default: "state")
#'
#' @return Data frame with added state column
#' @export
#'
#' @examples
#' cnv <- data.frame(copy_number = c(0, 1, 2, 3, 5, 8))
#' cnv <- add_state_column(cnv)
add_state_column <- function(df, cn_col = "copy_number", state_col = "state") {

    if (!cn_col %in% colnames(df)) {
        stop("Column '", cn_col, "' not found in data frame")
    }

    df[[state_col]] <- ifelse(df[[cn_col]] >= 4, 4, df[[cn_col]])

    return(df)
}

#' Standardize Chromosome Names
#'
#' Standardize chromosome names to UCSC format (chr1, chr2, etc.).
#'
#' @param df Data frame containing CNV data
#' @param chrom_col Name of the chromosome column
#'
#' @return Data frame with standardized chromosome names
#' @export
#'
#' @examples
#' cnv <- data.frame(chrom = c("1", "2", "chr3"))
#' cnv <- standardize_chromosomes(cnv)
standardize_chromosomes <- function(df, chrom_col = "chrom") {

    if (!chrom_col %in% colnames(df)) {
        stop("Column '", chrom_col, "' not found in data frame")
    }

    # Add 'chr' prefix if not present
    df[[chrom_col]] <- ifelse(
        grepl("^chr", df[[chrom_col]]),
        df[[chrom_col]],
        paste0("chr", df[[chrom_col]])
    )

    return(df)
}

#' Prepare CNV Data for Analysis
#'
#' Complete preprocessing pipeline for CNV data: convert segment mean to copy
#' number, add annotations, remove sex chromosomes and diploid entries, and
#' standardize chromosome names.
#'
#' @param df Data frame containing raw CNV data
#' @param segment_mean_col Name of the segment mean column
#' @param chrom_col Name of the chromosome column
#' @param remove_diploid Logical, whether to remove diploid entries (default: TRUE)
#' @param remove_sex Logical, whether to remove sex chromosomes (default: TRUE)
#'
#' @return Preprocessed data frame ready for CNVRanger analysis
#' @export
#'
#' @examples
#' \dontrun{
#' cnv_processed <- prepare_cnv_data(cnv_raw)
#' }
prepare_cnv_data <- function(df,
                              segment_mean_col = "segment_mean",
                              chrom_col = "chrom",
                              remove_diploid = TRUE,
                              remove_sex = TRUE) {

    # Convert segment mean to copy number and add annotation
    df <- convert_segment_mean_to_cn(df, segment_mean_col)

    # Add state column for CNVRanger
    df <- add_state_column(df)

    # Standardize chromosome names
    df <- standardize_chromosomes(df, chrom_col)

    # Remove sex chromosomes if requested
    if (remove_sex) {
        df <- remove_sex_chromosomes(df, chrom_col)
    }

    # Remove diploid entries if requested
    if (remove_diploid) {
        df <- filter_diploid(df)
    }

    return(df)
}

#' Create CNV GRangesList
#'
#' Convert CNV data frame to a GRangesList object suitable for CNVRanger.
#'
#' @param df Preprocessed CNV data frame
#' @param sample_col Name of the sample column
#' @param chrom_col Name of the chromosome column
#' @param start_col Name of the start position column
#' @param end_col Name of the end position column
#' @param state_col Name of the state column
#'
#' @return GRangesList object with one GRanges per sample
#' @export
#' @importFrom GenomicRanges makeGRangesListFromDataFrame sort
#' @importFrom GenomeInfoDb seqlevelsStyle<-
#'
#' @examples
#' \dontrun{
#' cnv_grl <- create_cnv_grangeslist(cnv_processed)
#' }
create_cnv_grangeslist <- function(df,
                                    sample_col = "sample",
                                    chrom_col = "chrom",
                                    start_col = "start",
                                    end_col = "end",
                                    state_col = "state") {

    # Select required columns
    required_cols <- c(sample_col, chrom_col, start_col, end_col, state_col)
    missing_cols <- setdiff(required_cols, colnames(df))
    if (length(missing_cols) > 0) {
        stop("Missing columns: ", paste(missing_cols, collapse = ", "))
    }

    # Rename columns for GRanges creation
    df_subset <- df[, required_cols]
    colnames(df_subset) <- c("sample", "seqnames", "start", "end", "state")

    # Create GRangesList
    grl <- GenomicRanges::makeGRangesListFromDataFrame(
        df_subset,
        split.field = "sample",
        keep.extra.columns = TRUE
    )

    # Set seqlevels style to UCSC
    GenomeInfoDb::`seqlevelsStyle<-`(grl, "UCSC")

    # Sort the GRangesList
    grl <- GenomicRanges::sort(grl)

    return(grl)
}

#' Calculate Population CNV Regions
#'
#' Identify recurrent CNV regions across the population using CNVRanger.
#'
#' @param grl GRangesList of CNV calls per sample
#' @param mode Mode for population range identification:
#'   "density" (default), "RO" (reciprocal overlap)
#' @param density Minimum fraction of samples with CNV (for density mode)
#' @param ro_thresh Reciprocal overlap threshold (for RO mode)
#' @param est_recurrence Logical, estimate recurrence significance
#' @param verbose Logical, show progress messages
#'
#' @return GRanges object with population-level CNV regions
#' @export
#' @importFrom CNVRanger populationRanges
#'
#' @examples
#' \dontrun{
#' cnvrs <- calculate_population_ranges(cnv_grl)
#' }
calculate_population_ranges <- function(grl,
                                         mode = c("density", "RO"),
                                         density = 0.1,
                                         ro_thresh = 0.51,
                                         est_recurrence = FALSE,
                                         verbose = TRUE) {

    mode <- match.arg(mode)

    if (mode == "density") {
        cnvrs <- CNVRanger::populationRanges(
            grl,
            density = density,
            est.recur = est_recurrence,
            verbose = verbose
        )
    } else {
        cnvrs <- CNVRanger::populationRanges(
            grl,
            mode = "RO",
            ro.thresh = ro_thresh,
            verbose = verbose
        )
    }

    return(cnvrs)
}

#' Summarize CNV Data
#'
#' Generate summary statistics for CNV data.
#'
#' @param df Preprocessed CNV data frame
#' @param cn_col Name of the copy number column
#' @param annotation_col Name of the annotation column
#' @param sample_col Name of the sample column
#'
#' @return List with summary statistics
#' @export
#'
#' @examples
#' \dontrun{
#' summary <- summarize_cnv_data(cnv_processed)
#' }
summarize_cnv_data <- function(df,
                                cn_col = "copy_number",
                                annotation_col = "annotation",
                                sample_col = "sample") {

    # CNV counts by type
    type_counts <- table(df[[annotation_col]])

    # CNV counts per sample
    sample_counts <- table(df[[sample_col]])

    # CNV size statistics
    if (all(c("start", "end") %in% colnames(df))) {
        sizes <- df$end - df$start
        size_stats <- list(
            mean = mean(sizes),
            median = stats::median(sizes),
            min = min(sizes),
            max = max(sizes)
        )
    } else {
        size_stats <- NULL
    }

    # Summary by copy number and annotation
    summary_by_type <- df |>
        dplyr::group_by(.data[[cn_col]], .data[[annotation_col]]) |>
        dplyr::summarise(
            count = dplyr::n(),
            .groups = "drop"
        )

    list(
        total_cnvs = nrow(df),
        n_samples = length(unique(df[[sample_col]])),
        type_counts = type_counts,
        sample_counts = sample_counts,
        size_stats = size_stats,
        summary_by_type = summary_by_type
    )
}
