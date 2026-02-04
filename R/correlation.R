#' @title CNV-Expression Correlation Functions
#' @description Functions for calculating correlation between CNV and gene expression
#' @name correlation
NULL

#' Calculate CNV-Expression Correlation
#'
#' Calculate correlation between CNV segment mean values and gene expression
#' levels for each gene.
#'
#' @param merged_data Data frame with merged CNV and expression data
#' @param cnv_col Column name for CNV values (segment mean)
#' @param expr_col Column name for expression values
#' @param gene_col Column name for gene identifiers
#' @param method Correlation method: "spearman" (default), "pearson", or "kendall"
#' @param min_samples Minimum number of samples for correlation (default: 5)
#'
#' @return Data frame with correlation results per gene
#' @export
#'
#' @examples
#' \dontrun{
#' cor_results <- calculate_cnv_expression_correlation(merged_data)
#' }
calculate_cnv_expression_correlation <- function(merged_data,
                                                  cnv_col = "segment_mean",
                                                  expr_col = "expression",
                                                  gene_col = "gene_id",
                                                  method = c("spearman", "pearson", "kendall"),
                                                  min_samples = 5) {

    method <- match.arg(method)

    # Check required columns
    required_cols <- c(cnv_col, expr_col, gene_col)
    missing_cols <- setdiff(required_cols, colnames(merged_data))
    if (length(missing_cols) > 0) {
        stop("Missing columns: ", paste(missing_cols, collapse = ", "))
    }

    # Get unique genes
    genes <- unique(merged_data[[gene_col]])

    # Calculate correlation for each gene
    results <- lapply(genes, function(gene) {
        gene_data <- merged_data[merged_data[[gene_col]] == gene, ]

        # Check minimum samples
        if (nrow(gene_data) < min_samples) {
            return(NULL)
        }

        # Calculate correlation
        tryCatch({
            cor_test <- stats::cor.test(
                x = gene_data[[cnv_col]],
                y = gene_data[[expr_col]],
                method = method
            )

            data.frame(
                gene = gene,
                correlation = cor_test$estimate,
                statistic = cor_test$statistic,
                p_value = cor_test$p.value,
                method = method,
                n_samples = nrow(gene_data),
                stringsAsFactors = FALSE
            )
        }, error = function(e) {
            NULL
        })
    })

    # Combine results
    results_df <- do.call(rbind, results[!vapply(results, is.null, logical(1))])

    if (is.null(results_df) || nrow(results_df) == 0) {
        warning("No valid correlations calculated")
        return(data.frame())
    }

    # Add adjusted p-values
    results_df$adj_p_value <- stats::p.adjust(results_df$p_value, method = "BH")

    # Order by p-value
    results_df <- results_df[order(results_df$p_value), ]
    rownames(results_df) <- NULL

    return(results_df)
}

#' Prepare Data for Correlation Analysis
#'
#' Prepare merged CNV and expression data for correlation analysis.
#'
#' @param cnv_data Data frame with CNV data (must have gene annotation)
#' @param expression_data Matrix or data frame with expression counts
#' @param sample_col Sample column in CNV data
#' @param gene_col Gene column in CNV data
#' @param cnv_value_col Column with CNV values (segment mean)
#'
#' @return Data frame ready for correlation analysis
#' @export
#'
#' @examples
#' \dontrun{
#' prep_data <- prepare_correlation_data(cnv_annotated, expr_counts)
#' }
prepare_correlation_data <- function(cnv_data,
                                      expression_data,
                                      sample_col = "sample",
                                      gene_col = "gene_id",
                                      cnv_value_col = "segment_mean") {

    # Check required columns in CNV data
    required_cols <- c(sample_col, gene_col, cnv_value_col)
    missing_cols <- setdiff(required_cols, colnames(cnv_data))
    if (length(missing_cols) > 0) {
        stop("Missing columns in CNV data: ", paste(missing_cols, collapse = ", "))
    }

    # Create sample-gene identifiers
    cnv_data$sample_gene <- paste(cnv_data[[sample_col]], cnv_data[[gene_col]], sep = "_")

    # Get genes present in both datasets
    common_genes <- intersect(cnv_data[[gene_col]], rownames(expression_data))
    if (length(common_genes) == 0) {
        stop("No common genes found between CNV and expression data")
    }

    # Get samples present in both datasets
    cnv_samples <- unique(cnv_data[[sample_col]])
    expr_samples <- colnames(expression_data)
    common_samples <- intersect(cnv_samples, expr_samples)
    if (length(common_samples) == 0) {
        stop("No common samples found between CNV and expression data")
    }

    message("Found ", length(common_genes), " common genes and ",
            length(common_samples), " common samples")

    # Filter CNV data to common genes and samples
    cnv_filtered <- cnv_data[
        cnv_data[[gene_col]] %in% common_genes &
        cnv_data[[sample_col]] %in% common_samples,
    ]

    # Add expression values
    cnv_filtered$expression <- vapply(seq_len(nrow(cnv_filtered)), function(i) {
        gene <- cnv_filtered[[gene_col]][i]
        sample <- cnv_filtered[[sample_col]][i]
        if (gene %in% rownames(expression_data) && sample %in% colnames(expression_data)) {
            as.numeric(expression_data[gene, sample])
        } else {
            NA_real_
        }
    }, numeric(1))

    # Remove rows with NA expression
    cnv_filtered <- cnv_filtered[!is.na(cnv_filtered$expression), ]

    return(cnv_filtered)
}

#' Filter Significant Correlations
#'
#' Filter correlation results to retain only significant associations.
#'
#' @param cor_results Data frame from calculate_cnv_expression_correlation
#' @param p_cutoff P-value cutoff (default: 0.05)
#' @param cor_cutoff Minimum absolute correlation coefficient (default: 0)
#' @param use_adjusted Use adjusted p-values (default: TRUE)
#'
#' @return Filtered data frame
#' @export
#'
#' @examples
#' \dontrun{
#' sig_cors <- filter_significant_correlations(cor_results, p_cutoff = 0.05)
#' }
filter_significant_correlations <- function(cor_results,
                                             p_cutoff = 0.05,
                                             cor_cutoff = 0,
                                             use_adjusted = TRUE) {

    p_col <- if (use_adjusted) "adj_p_value" else "p_value"

    if (!p_col %in% colnames(cor_results)) {
        warning("Column '", p_col, "' not found, using p_value")
        p_col <- "p_value"
    }

    filtered <- cor_results[
        cor_results[[p_col]] <= p_cutoff &
        abs(cor_results$correlation) >= cor_cutoff,
    ]

    message("Filtered from ", nrow(cor_results), " to ", nrow(filtered), " correlations")

    return(filtered)
}

#' Classify Correlation Direction
#'
#' Add classification of correlation direction (positive/negative/neutral).
#'
#' @param cor_results Data frame with correlation results
#' @param threshold Correlation threshold for classification (default: 0.2)
#'
#' @return Data frame with added direction column
#' @export
#'
#' @examples
#' \dontrun{
#' classified <- classify_correlation_direction(cor_results)
#' }
classify_correlation_direction <- function(cor_results, threshold = 0.2) {

    cor_results$direction <- dplyr::case_when(
        cor_results$correlation >= threshold ~ "positive",
        cor_results$correlation <= -threshold ~ "negative",
        TRUE ~ "neutral"
    )

    # Summary
    summary_table <- table(cor_results$direction)
    message("Correlation direction summary: ",
            paste(names(summary_table), "=", summary_table, collapse = ", "))

    return(cor_results)
}

#' Summarize Correlation Results
#'
#' Generate summary statistics for correlation results.
#'
#' @param cor_results Data frame with correlation results
#' @param by_chromosome Logical, summarize by chromosome if available
#'
#' @return List with summary statistics
#' @export
#'
#' @examples
#' \dontrun{
#' summary <- summarize_correlation_results(cor_results)
#' }
summarize_correlation_results <- function(cor_results, by_chromosome = FALSE) {

    summary <- list(
        total_genes = nrow(cor_results),
        significant_genes = sum(cor_results$adj_p_value <= 0.05, na.rm = TRUE),
        mean_correlation = mean(cor_results$correlation, na.rm = TRUE),
        median_correlation = stats::median(cor_results$correlation, na.rm = TRUE),
        positive_correlations = sum(cor_results$correlation > 0, na.rm = TRUE),
        negative_correlations = sum(cor_results$correlation < 0, na.rm = TRUE)
    )

    # Strong correlations
    summary$strong_positive <- sum(cor_results$correlation >= 0.5, na.rm = TRUE)
    summary$strong_negative <- sum(cor_results$correlation <= -0.5, na.rm = TRUE)

    # By direction if available
    if ("direction" %in% colnames(cor_results)) {
        summary$by_direction <- table(cor_results$direction)
    }

    return(summary)
}
