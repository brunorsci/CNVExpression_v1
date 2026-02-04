#' @title CNV-Expression Integration Functions
#' @description Functions for integrating CNV and gene expression data
#' @name integration
NULL

#' Create Gene GRanges from Expression Data
#'
#' Create a GRanges object for genes based on expression data using biomaRt.
#'
#' @param gene_ids Character vector of gene identifiers (e.g., Ensembl IDs)
#' @param id_type Type of gene ID: "ensembl" (default) or "symbol"
#' @param genome Reference genome: "GRCh38" (default) or "GRCh37"
#'
#' @return GRanges object with gene coordinates
#' @export
#'
#' @examples
#' \dontrun{
#' genes_gr <- create_gene_granges(rownames(expr_data))
#' }
create_gene_granges <- function(gene_ids,
                                 id_type = c("ensembl", "symbol"),
                                 genome = c("GRCh38", "GRCh37")) {

    id_type <- match.arg(id_type)
    genome <- match.arg(genome)

    if (!requireNamespace("biomaRt", quietly = TRUE)) {
        stop("Package 'biomaRt' is required. Install with: ",
             "BiocManager::install('biomaRt')")
    }

    # Clean gene IDs (remove version for Ensembl)
    if (id_type == "ensembl") {
        gene_ids_clean <- gsub("\\..*", "", gene_ids)
        filter_name <- "ensembl_gene_id"
    } else {
        gene_ids_clean <- gene_ids
        filter_name <- "external_gene_name"
    }

    # Connect to Ensembl
    if (genome == "GRCh37") {
        ensembl <- biomaRt::useEnsembl(
            biomart = "ensembl",
            dataset = "hsapiens_gene_ensembl",
            GRCh = 37
        )
    } else {
        ensembl <- biomaRt::useEnsembl(
            biomart = "ensembl",
            dataset = "hsapiens_gene_ensembl"
        )
    }

    # Get gene information
    gene_info <- biomaRt::getBM(
        attributes = c("ensembl_gene_id", "chromosome_name",
                       "start_position", "end_position", "strand"),
        filters = filter_name,
        values = gene_ids_clean,
        mart = ensembl
    )

    # Filter to standard chromosomes
    standard_chroms <- c(as.character(1:22), "X", "Y")
    gene_info <- gene_info[gene_info$chromosome_name %in% standard_chroms, ]

    # Create GRanges object
    genes_gr <- GenomicRanges::GRanges(
        seqnames = paste0("chr", gene_info$chromosome_name),
        ranges = IRanges::IRanges(
            start = gene_info$start_position,
            end = gene_info$end_position
        ),
        strand = "*",
        gene_id = gene_info$ensembl_gene_id
    )

    names(genes_gr) <- genes_gr$gene_id

    return(genes_gr)
}

#' Create SummarizedExperiment from Expression Data
#'
#' Create a SummarizedExperiment object from expression counts and gene
#' coordinates.
#'
#' @param counts Matrix or data.frame of expression counts (genes x samples)
#' @param genes_gr GRanges object with gene coordinates
#' @param match_genes Logical, match genes between counts and GRanges
#'
#' @return SummarizedExperiment object
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
#'
#' @examples
#' \dontrun{
#' rse <- create_expression_se(counts, genes_gr)
#' }
create_expression_se <- function(counts, genes_gr, match_genes = TRUE) {

    if (match_genes) {
        # Find common genes
        common_genes <- intersect(rownames(counts), names(genes_gr))

        if (length(common_genes) == 0) {
            stop("No common genes found between counts and gene GRanges")
        }

        message("Found ", length(common_genes), " common genes out of ",
                nrow(counts), " in expression data")

        # Subset to common genes
        counts <- counts[common_genes, , drop = FALSE]
        genes_gr <- genes_gr[common_genes]
    }

    # Ensure counts is a matrix
    if (!is.matrix(counts)) {
        counts <- as.matrix(counts)
    }

    # Create SummarizedExperiment
    rse <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = counts),
        rowRanges = genes_gr
    )

    return(rse)
}

#' Run CNV-Expression Integration Analysis
#'
#' Perform CNV-expression quantitative trait loci (eQTL) analysis using
#' CNVRanger to identify genes whose expression is associated with copy
#' number changes.
#'
#' @param cnvrs GRanges object with CNV regions (from populationRanges)
#' @param cnv_grl GRangesList of individual CNV calls
#' @param rse SummarizedExperiment with expression data
#' @param window Window size for CNV-gene association (default: "1Mbp")
#' @param min_samples Minimum number of samples with CNV (default: 3)
#' @param min_state_freq Minimum frequency per CNV state (default: 3)
#' @param filter_by_expr Logical, filter lowly expressed genes (default: TRUE)
#' @param verbose Logical, show progress messages (default: TRUE)
#'
#' @return Data frame with CNV-expression association results
#' @export
#' @importFrom CNVRanger cnvEQTL
#'
#' @examples
#' \dontrun{
#' results <- run_cnv_expression_analysis(cnvrs, cnv_grl, rse)
#' }
run_cnv_expression_analysis <- function(cnvrs,
                                         cnv_grl,
                                         rse,
                                         window = "1Mbp",
                                         min_samples = 3,
                                         min_state_freq = 3,
                                         filter_by_expr = TRUE,
                                         verbose = TRUE) {

    # Run cnvEQTL analysis
    results <- CNVRanger::cnvEQTL(
        cnvrs = cnvrs,
        grl = cnv_grl,
        rse = rse,
        window = window,
        min.samples = min_samples,
        min.state.freq = min_state_freq,
        filter.by.expr = filter_by_expr,
        verbose = verbose
    )

    # Convert to data frame
    results_df <- as.data.frame(results)

    return(results_df)
}

#' Filter CNV-Expression Results
#'
#' Filter CNV-expression association results by significance and effect size.
#'
#' @param results Data frame from run_cnv_expression_analysis
#' @param logfc_cutoff Minimum absolute log fold change (default: 1)
#' @param pvalue_cutoff Maximum adjusted p-value (default: 0.05)
#' @param logfc_cols Columns containing logFC values
#'
#' @return Filtered data frame
#' @export
#'
#' @examples
#' \dontrun{
#' filtered <- filter_cnv_expression_results(results, logfc_cutoff = 1)
#' }
filter_cnv_expression_results <- function(results,
                                           logfc_cutoff = 1,
                                           pvalue_cutoff = 0.05,
                                           logfc_cols = c("logFC.CN0", "logFC.CN1",
                                                          "logFC.CN3", "logFC.CN4")) {

    # Find available logFC columns
    available_logfc <- intersect(logfc_cols, colnames(results))

    if (length(available_logfc) == 0) {
        warning("No logFC columns found in results")
        # Filter by p-value only
        filtered <- results[results$AdjPValue <= pvalue_cutoff, ]
    } else {
        # Create logFC filter: at least one logFC exceeds threshold
        logfc_matrix <- as.matrix(results[, available_logfc, drop = FALSE])
        logfc_filter <- apply(logfc_matrix, 1, function(x) {
            any(abs(x) >= logfc_cutoff, na.rm = TRUE)
        })

        # Apply both filters
        filtered <- results[results$AdjPValue <= pvalue_cutoff & logfc_filter, ]
    }

    message("Filtered from ", nrow(results), " to ", nrow(filtered), " results")

    return(filtered)
}

#' Process CNV-Expression Integration Results
#'
#' Process and merge CNV-expression results with original CNV data to create
#' a comprehensive output.
#'
#' @param results Filtered results from filter_cnv_expression_results
#' @param cnv_data Original preprocessed CNV data frame
#' @param sample_col Sample column in CNV data
#' @param chrom_col Chromosome column in CNV data
#'
#' @return Merged data frame with CNV and expression information
#' @export
#'
#' @examples
#' \dontrun{
#' merged <- process_integration_results(filtered_results, cnv_data)
#' }
process_integration_results <- function(results,
                                         cnv_data,
                                         sample_col = "sample",
                                         chrom_col = "chrom") {

    # Parse group_name to extract gene coordinates
    if ("group_name" %in% colnames(results)) {
        results <- results |>
            tidyr::separate(
                group_name,
                into = c("gene_chrom", "gene_positions"),
                sep = ":",
                remove = FALSE
            ) |>
            tidyr::separate(
                gene_positions,
                into = c("gene_start", "gene_end"),
                sep = "-",
                convert = TRUE
            )
    }

    # Create GRanges for genes
    if (all(c("gene_chrom", "gene_start", "gene_end") %in% colnames(results))) {
        gene_gr <- GenomicRanges::GRanges(
            seqnames = results$gene_chrom,
            ranges = IRanges::IRanges(
                start = results$gene_start,
                end = results$gene_end
            )
        )

        # Standardize chromosome names in CNV data
        cnv_data_std <- cnv_data
        if (!any(grepl("^chr", cnv_data_std[[chrom_col]]))) {
            cnv_data_std[[chrom_col]] <- paste0("chr", cnv_data_std[[chrom_col]])
        }

        # Create GRanges for CNVs
        cnv_gr <- GenomicRanges::GRanges(
            seqnames = cnv_data_std[[chrom_col]],
            ranges = IRanges::IRanges(
                start = cnv_data_std$start,
                end = cnv_data_std$end
            ),
            sample = cnv_data_std[[sample_col]],
            copy_number = cnv_data_std$copy_number,
            annotation = cnv_data_std$annotation
        )

        # Find overlaps
        overlaps <- GenomicRanges::findOverlaps(gene_gr, cnv_gr)

        if (length(overlaps) > 0) {
            # Merge data
            merged <- data.frame(
                results[S4Vectors::queryHits(overlaps), ],
                cnv_sample = cnv_data_std[[sample_col]][S4Vectors::subjectHits(overlaps)],
                cnv_chrom = cnv_data_std[[chrom_col]][S4Vectors::subjectHits(overlaps)],
                cnv_start = cnv_data_std$start[S4Vectors::subjectHits(overlaps)],
                cnv_end = cnv_data_std$end[S4Vectors::subjectHits(overlaps)],
                cnv_copy_number = cnv_data_std$copy_number[S4Vectors::subjectHits(overlaps)],
                cnv_annotation = cnv_data_std$annotation[S4Vectors::subjectHits(overlaps)]
            )

            return(merged)
        }
    }

    warning("Could not process integration results. Returning original results.")
    return(results)
}

#' Complete CNV-Expression Analysis Pipeline
#'
#' Run the complete CNV-expression integration pipeline from raw data to
#' filtered results.
#'
#' @param cnv_data Data frame with raw CNV data
#' @param expression_data Matrix or data frame with expression counts
#' @param segment_mean_col Name of segment mean column in CNV data
#' @param sample_col Name of sample column in CNV data
#' @param genome Reference genome ("GRCh38" or "GRCh37")
#' @param window Window for CNV-gene association
#' @param logfc_cutoff Minimum logFC for filtering
#' @param pvalue_cutoff Maximum p-value for filtering
#' @param verbose Show progress messages
#'
#' @return List with analysis results and intermediate objects
#' @export
#'
#' @examples
#' \dontrun{
#' results <- cnv_expression_pipeline(cnv_raw, expr_counts)
#' }
cnv_expression_pipeline <- function(cnv_data,
                                     expression_data,
                                     segment_mean_col = "segment_mean",
                                     sample_col = "sample",
                                     genome = "GRCh38",
                                     window = "1Mbp",
                                     logfc_cutoff = 1,
                                     pvalue_cutoff = 0.05,
                                     verbose = TRUE) {

    if (verbose) message("Step 1/6: Preparing CNV data...")
    cnv_processed <- prepare_cnv_data(cnv_data, segment_mean_col)

    # Standardize sample IDs in CNV data
    cnv_processed[[sample_col]] <- standardize_sample_ids(cnv_processed[[sample_col]])

    # Standardize sample IDs in expression data
    colnames(expression_data) <- standardize_sample_ids(colnames(expression_data))

    if (verbose) message("Step 2/6: Creating CNV GRangesList...")
    cnv_grl <- create_cnv_grangeslist(cnv_processed, sample_col = sample_col)

    if (verbose) message("Step 3/6: Calculating population CNV regions...")
    cnvrs <- calculate_population_ranges(cnv_grl, verbose = verbose)

    if (verbose) message("Step 4/6: Creating gene GRanges...")
    genes_gr <- create_gene_granges(rownames(expression_data), genome = genome)

    if (verbose) message("Step 5/6: Creating SummarizedExperiment...")
    rse <- create_expression_se(expression_data, genes_gr)

    if (verbose) message("Step 6/6: Running CNV-expression analysis...")
    results <- run_cnv_expression_analysis(
        cnvrs, cnv_grl, rse,
        window = window,
        verbose = verbose
    )

    # Filter results
    filtered <- filter_cnv_expression_results(
        results,
        logfc_cutoff = logfc_cutoff,
        pvalue_cutoff = pvalue_cutoff
    )

    # Process and merge with CNV data
    merged <- process_integration_results(filtered, cnv_processed, sample_col)

    list(
        cnv_processed = cnv_processed,
        cnv_grangeslist = cnv_grl,
        cnv_regions = cnvrs,
        genes_granges = genes_gr,
        expression_se = rse,
        results_raw = results,
        results_filtered = filtered,
        results_merged = merged
    )
}
