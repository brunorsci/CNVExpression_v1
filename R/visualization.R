#' @title Visualization Functions for CNVExpression
#' @description Functions for visualizing CNV-expression analysis results
#' @name visualization
NULL

#' Plot CNV Type Distribution
#'
#' Create a bar plot showing the distribution of CNV types.
#'
#' @param cnv_data Data frame with CNV data
#' @param type_col Column name for CNV type/annotation
#' @param colors Named vector of colors for each type (optional)
#' @param title Plot title
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' plot_cnv_types(cnv_data)
#' }
plot_cnv_types <- function(cnv_data,
                           type_col = "annotation",
                           colors = NULL,
                           title = "Distribution of CNV Types") {

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' required. Install with: install.packages('ggplot2')")
    }

    if (!type_col %in% colnames(cnv_data)) {
        stop("Column '", type_col, "' not found in data")
    }

    # Default colors
    if (is.null(colors)) {
        colors <- c(
            "homozygous_deletion" = "#d73027",
            "heterozygous_deletion" = "#fc8d59",
            "neutral" = "#ffffbf",
            "amplification" = "#4575b4"
        )
    }

    p <- ggplot2::ggplot(cnv_data, ggplot2::aes(x = .data[[type_col]])) +
        ggplot2::geom_bar(fill = "steelblue") +
        ggplot2::labs(
            title = title,
            x = "CNV Type",
            y = "Count"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
        )

    return(p)
}

#' Plot CNV Distribution by Chromosome
#'
#' Create a stacked bar plot showing CNV distribution across chromosomes.
#'
#' @param cnv_data Data frame with CNV data
#' @param chrom_col Column name for chromosome
#' @param type_col Column name for CNV type/annotation
#' @param colors Named vector of colors for each type (optional)
#' @param title Plot title
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' plot_cnv_by_chromosome(cnv_data)
#' }
plot_cnv_by_chromosome <- function(cnv_data,
                                    chrom_col = "chrom",
                                    type_col = "annotation",
                                    colors = NULL,
                                    title = "CNV Distribution by Chromosome") {

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' required. Install with: install.packages('ggplot2')")
    }

    # Default colors
    if (is.null(colors)) {
        colors <- c(
            "homozygous_deletion" = "#d73027",
            "heterozygous_deletion" = "#fc8d59",
            "neutral" = "#ffffbf",
            "amplification" = "#4575b4"
        )
    }

    # Order chromosomes
    cnv_data[[chrom_col]] <- factor(
        cnv_data[[chrom_col]],
        levels = c(paste0("chr", 1:22), "chrX", "chrY",
                   as.character(1:22), "X", "Y")
    )

    p <- ggplot2::ggplot(cnv_data, ggplot2::aes(
        x = .data[[chrom_col]],
        fill = .data[[type_col]]
    )) +
        ggplot2::geom_bar() +
        ggplot2::scale_fill_manual(values = colors, na.value = "gray50") +
        ggplot2::labs(
            title = title,
            x = "Chromosome",
            y = "Count",
            fill = "CNV Type"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
        )

    return(p)
}

#' Plot CNV Length Distribution
#'
#' Create a boxplot or histogram showing CNV length distribution.
#'
#' @param cnv_data Data frame with CNV data
#' @param start_col Column name for start position
#' @param end_col Column name for end position
#' @param type_col Column name for CNV type (for grouping)
#' @param plot_type Type of plot: "boxplot" or "histogram"
#' @param log_scale Use log10 scale for length
#' @param title Plot title
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' plot_cnv_length(cnv_data)
#' }
plot_cnv_length <- function(cnv_data,
                            start_col = "start",
                            end_col = "end",
                            type_col = "annotation",
                            plot_type = c("boxplot", "histogram"),
                            log_scale = TRUE,
                            title = "CNV Length Distribution") {

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' required. Install with: install.packages('ggplot2')")
    }

    plot_type <- match.arg(plot_type)

    # Calculate length
    cnv_data$length <- cnv_data[[end_col]] - cnv_data[[start_col]]

    if (plot_type == "boxplot" && type_col %in% colnames(cnv_data)) {
        p <- ggplot2::ggplot(cnv_data, ggplot2::aes(
            x = .data[[type_col]],
            y = .data$length
        )) +
            ggplot2::geom_boxplot(fill = "lightblue") +
            ggplot2::labs(
                title = title,
                x = "CNV Type",
                y = "Length (bp)"
            )
    } else {
        p <- ggplot2::ggplot(cnv_data, ggplot2::aes(x = .data$length)) +
            ggplot2::geom_histogram(bins = 50, fill = "steelblue", color = "white") +
            ggplot2::labs(
                title = title,
                x = "Length (bp)",
                y = "Count"
            )
    }

    if (log_scale) {
        p <- p + ggplot2::scale_y_continuous(trans = "log10")
    }

    p <- p +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
        )

    return(p)
}

#' Plot Volcano Plot for CNV-Expression Results
#'
#' Create a volcano plot showing log fold change vs -log10(p-value).
#'
#' @param results Data frame with CNV-expression results
#' @param logfc_col Column name for log fold change
#' @param pvalue_col Column name for p-value (or adjusted p-value)
#' @param logfc_cutoff Log fold change threshold for significance
#' @param pvalue_cutoff P-value threshold for significance
#' @param highlight_genes Character vector of genes to highlight
#' @param title Plot title
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' plot_volcano(results, logfc_col = "logFC.CN1")
#' }
plot_volcano <- function(results,
                         logfc_col = "logFC.CN1",
                         pvalue_col = "AdjPValue",
                         logfc_cutoff = 1,
                         pvalue_cutoff = 0.05,
                         highlight_genes = NULL,
                         title = "Volcano Plot") {

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' required. Install with: install.packages('ggplot2')")
    }

    # Prepare data
    plot_data <- results[, c(logfc_col, pvalue_col), drop = FALSE]
    colnames(plot_data) <- c("logFC", "pvalue")
    plot_data$neg_log_p <- -log10(plot_data$pvalue)

    # Classify points
    plot_data$significance <- "Not significant"
    plot_data$significance[abs(plot_data$logFC) >= logfc_cutoff &
                              plot_data$pvalue <= pvalue_cutoff] <- "Significant"

    p <- ggplot2::ggplot(plot_data, ggplot2::aes(
        x = .data$logFC,
        y = .data$neg_log_p,
        color = .data$significance
    )) +
        ggplot2::geom_point(alpha = 0.6, size = 1.5) +
        ggplot2::scale_color_manual(
            values = c("Not significant" = "gray50", "Significant" = "red")
        ) +
        ggplot2::geom_hline(
            yintercept = -log10(pvalue_cutoff),
            linetype = "dashed",
            color = "blue"
        ) +
        ggplot2::geom_vline(
            xintercept = c(-logfc_cutoff, logfc_cutoff),
            linetype = "dashed",
            color = "blue"
        ) +
        ggplot2::labs(
            title = title,
            x = paste("Log2 Fold Change (", logfc_col, ")", sep = ""),
            y = "-Log10(P-value)",
            color = "Significance"
        ) +
        ggplot2::theme_minimal()

    return(p)
}

#' Plot Correlation Scatter
#'
#' Create a scatter plot showing CNV-expression correlation for a gene.
#'
#' @param data Data frame with merged CNV and expression data
#' @param gene Gene identifier to plot
#' @param cnv_col Column name for CNV values
#' @param expr_col Column name for expression values
#' @param gene_col Column name for gene identifier
#' @param title Plot title (auto-generated if NULL)
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' plot_correlation_scatter(merged_data, gene = "TP53")
#' }
plot_correlation_scatter <- function(data,
                                      gene,
                                      cnv_col = "segment_mean",
                                      expr_col = "expression",
                                      gene_col = "gene_id",
                                      title = NULL) {

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' required. Install with: install.packages('ggplot2')")
    }

    # Filter to specific gene
    gene_data <- data[data[[gene_col]] == gene, ]

    if (nrow(gene_data) == 0) {
        stop("Gene '", gene, "' not found in data")
    }

    # Calculate correlation
    cor_test <- stats::cor.test(gene_data[[cnv_col]], gene_data[[expr_col]],
                                method = "spearman")

    # Auto-generate title
    if (is.null(title)) {
        title <- sprintf("%s: r = %.3f, p = %.2e",
                         gene, cor_test$estimate, cor_test$p.value)
    }

    p <- ggplot2::ggplot(gene_data, ggplot2::aes(
        x = .data[[cnv_col]],
        y = .data[[expr_col]]
    )) +
        ggplot2::geom_point(alpha = 0.6) +
        ggplot2::geom_smooth(method = "lm", se = TRUE, color = "red") +
        ggplot2::labs(
            title = title,
            x = "CNV Segment Mean",
            y = "Expression"
        ) +
        ggplot2::theme_minimal()

    return(p)
}

#' Plot CNV Summary per Sample
#'
#' Create a bar plot showing CNV counts per sample.
#'
#' @param cnv_data Data frame with CNV data
#' @param sample_col Column name for sample identifier
#' @param type_col Column name for CNV type (optional, for stacking)
#' @param top_n Show only top N samples (default: NULL for all)
#' @param title Plot title
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' plot_cnv_per_sample(cnv_data, top_n = 20)
#' }
plot_cnv_per_sample <- function(cnv_data,
                                 sample_col = "sample",
                                 type_col = "annotation",
                                 top_n = NULL,
                                 title = "CNV Counts per Sample") {

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' required. Install with: install.packages('ggplot2')")
    }

    # Count CNVs per sample
    sample_counts <- cnv_data |>
        dplyr::group_by(.data[[sample_col]]) |>
        dplyr::summarise(count = dplyr::n(), .groups = "drop") |>
        dplyr::arrange(dplyr::desc(.data$count))

    # Select top N if specified
    if (!is.null(top_n) && top_n < nrow(sample_counts)) {
        top_samples <- sample_counts[[sample_col]][seq_len(top_n)]
        cnv_data <- cnv_data[cnv_data[[sample_col]] %in% top_samples, ]
    }

    # Order samples by count
    cnv_data[[sample_col]] <- factor(
        cnv_data[[sample_col]],
        levels = sample_counts[[sample_col]]
    )

    if (type_col %in% colnames(cnv_data)) {
        p <- ggplot2::ggplot(cnv_data, ggplot2::aes(
            x = .data[[sample_col]],
            fill = .data[[type_col]]
        )) +
            ggplot2::geom_bar()
    } else {
        p <- ggplot2::ggplot(cnv_data, ggplot2::aes(x = .data[[sample_col]])) +
            ggplot2::geom_bar(fill = "steelblue")
    }

    p <- p +
        ggplot2::labs(
            title = title,
            x = "Sample",
            y = "CNV Count"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 6)
        )

    return(p)
}

#' Generate Analysis Report Plots
#'
#' Generate a set of standard plots for CNV-expression analysis.
#'
#' @param results List of results from cnv_expression_pipeline
#' @param output_dir Directory to save plots (NULL for no saving)
#' @param format Output format: "png", "pdf", or "both"
#'
#' @return List of ggplot objects
#' @export
#'
#' @examples
#' \dontrun{
#' plots <- generate_report_plots(analysis_results, output_dir = "plots")
#' }
generate_report_plots <- function(results,
                                   output_dir = NULL,
                                   format = c("png", "pdf", "both")) {

    format <- match.arg(format)
    plots <- list()

    # CNV type distribution
    if ("cnv_processed" %in% names(results)) {
        plots$cnv_types <- plot_cnv_types(results$cnv_processed)
        plots$cnv_by_chrom <- plot_cnv_by_chromosome(results$cnv_processed)
        plots$cnv_length <- plot_cnv_length(results$cnv_processed)
    }

    # Volcano plot
    if ("results_filtered" %in% names(results) && nrow(results$results_filtered) > 0) {
        logfc_cols <- grep("logFC", colnames(results$results_filtered), value = TRUE)
        if (length(logfc_cols) > 0) {
            plots$volcano <- plot_volcano(results$results_filtered, logfc_col = logfc_cols[1])
        }
    }

    # Save plots if output directory specified
    if (!is.null(output_dir)) {
        if (!dir.exists(output_dir)) {
            dir.create(output_dir, recursive = TRUE)
        }

        for (name in names(plots)) {
            if (format %in% c("png", "both")) {
                ggplot2::ggsave(
                    filename = file.path(output_dir, paste0(name, ".png")),
                    plot = plots[[name]],
                    width = 8,
                    height = 6,
                    dpi = 300
                )
            }
            if (format %in% c("pdf", "both")) {
                ggplot2::ggsave(
                    filename = file.path(output_dir, paste0(name, ".pdf")),
                    plot = plots[[name]],
                    width = 8,
                    height = 6
                )
            }
        }
        message("Plots saved to: ", output_dir)
    }

    return(plots)
}
