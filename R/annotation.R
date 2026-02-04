#' @title Gene Annotation Functions
#' @description Functions for annotating CNV regions with gene information
#' @name annotation
NULL

#' Annotate CNV Regions with Genes
#'
#' Annotate CNV regions with overlapping genes using a reference genome.
#'
#' @param cnv_gr GRanges object with CNV regions
#' @param genome Reference genome: "hg38" (default) or "hg19"
#'
#' @return GRanges object with gene annotations added
#' @export
#'
#' @examples
#' \dontrun{
#' cnv_annotated <- annotate_cnv_with_genes(cnv_gr, genome = "hg38")
#' }
annotate_cnv_with_genes <- function(cnv_gr, genome = c("hg38", "hg19")) {

    genome <- match.arg(genome)

    # Load appropriate TxDb
    if (genome == "hg38") {
        if (!requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)) {
            stop("Package 'TxDb.Hsapiens.UCSC.hg38.knownGene' required. ",
                 "Install with: BiocManager::install('TxDb.Hsapiens.UCSC.hg38.knownGene')")
        }
        txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
    } else {
        if (!requireNamespace("TxDb.Hsapiens.UCSC.hg19.knownGene", quietly = TRUE)) {
            stop("Package 'TxDb.Hsapiens.UCSC.hg19.knownGene' required. ",
                 "Install with: BiocManager::install('TxDb.Hsapiens.UCSC.hg19.knownGene')")
        }
        txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
    }

    # Get genes from TxDb
    genes_ref <- GenomicFeatures::genes(txdb)

    # Find overlaps
    overlaps <- GenomicRanges::findOverlaps(cnv_gr, genes_ref)

    if (length(overlaps) == 0) {
        warning("No overlapping genes found")
        S4Vectors::mcols(cnv_gr)$gene_id <- NA_character_
        return(cnv_gr)
    }

    # Add gene IDs to CNV regions
    gene_ids <- split(
        names(genes_ref)[S4Vectors::subjectHits(overlaps)],
        S4Vectors::queryHits(overlaps)
    )

    # Initialize with NA
    S4Vectors::mcols(cnv_gr)$entrez_id <- NA_character_

    # Add gene IDs
    for (i in names(gene_ids)) {
        S4Vectors::mcols(cnv_gr)$entrez_id[as.integer(i)] <-
            paste(gene_ids[[i]], collapse = ";")
    }

    return(cnv_gr)
}

#' Convert Entrez IDs to Gene Symbols
#'
#' Convert Entrez gene IDs to gene symbols using org.Hs.eg.db.
#'
#' @param entrez_ids Character vector of Entrez IDs
#'
#' @return Named character vector with gene symbols
#' @export
#'
#' @examples
#' \dontrun{
#' symbols <- entrez_to_symbol(c("1", "2", "3"))
#' }
entrez_to_symbol <- function(entrez_ids) {

    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
        stop("Package 'org.Hs.eg.db' required. ",
             "Install with: BiocManager::install('org.Hs.eg.db')")
    }

    # Remove NAs and split if multiple IDs are concatenated
    entrez_ids_clean <- unique(unlist(strsplit(entrez_ids, ";")))
    entrez_ids_clean <- entrez_ids_clean[!is.na(entrez_ids_clean) & entrez_ids_clean != ""]

    if (length(entrez_ids_clean) == 0) {
        return(character(0))
    }

    symbols <- AnnotationDbi::mapIds(
        org.Hs.eg.db::org.Hs.eg.db,
        keys = entrez_ids_clean,
        column = "SYMBOL",
        keytype = "ENTREZID",
        multiVals = "first"
    )

    return(symbols)
}

#' Convert Ensembl IDs to Gene Symbols
#'
#' Convert Ensembl gene IDs to gene symbols using org.Hs.eg.db.
#'
#' @param ensembl_ids Character vector of Ensembl IDs
#' @param remove_version Logical, remove version suffix (default: TRUE)
#'
#' @return Named character vector with gene symbols
#' @export
#'
#' @examples
#' \dontrun{
#' symbols <- ensembl_to_symbol(c("ENSG00000141510", "ENSG00000134982"))
#' }
ensembl_to_symbol <- function(ensembl_ids, remove_version = TRUE) {

    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
        stop("Package 'org.Hs.eg.db' required. ",
             "Install with: BiocManager::install('org.Hs.eg.db')")
    }

    # Remove version suffix if requested
    if (remove_version) {
        ensembl_ids <- gsub("\\..*", "", ensembl_ids)
    }

    # Remove NAs and duplicates
    ensembl_ids_clean <- unique(ensembl_ids[!is.na(ensembl_ids)])

    if (length(ensembl_ids_clean) == 0) {
        return(character(0))
    }

    symbols <- AnnotationDbi::mapIds(
        org.Hs.eg.db::org.Hs.eg.db,
        keys = ensembl_ids_clean,
        column = "SYMBOL",
        keytype = "ENSEMBL",
        multiVals = "first"
    )

    return(symbols)
}

#' Add Gene Symbols to Results
#'
#' Add gene symbol annotations to CNV-expression results.
#'
#' @param results Data frame with CNV-expression results
#' @param gene_col Column name containing gene identifiers
#' @param id_type Type of gene ID: "entrez", "ensembl", or "auto" (detect)
#'
#' @return Data frame with added gene_symbol column
#' @export
#'
#' @examples
#' \dontrun{
#' results_annotated <- add_gene_symbols(results)
#' }
add_gene_symbols <- function(results,
                              gene_col = "gene_id",
                              id_type = c("auto", "entrez", "ensembl")) {

    id_type <- match.arg(id_type)

    if (!gene_col %in% colnames(results)) {
        warning("Gene column '", gene_col, "' not found")
        return(results)
    }

    gene_ids <- results[[gene_col]]

    # Auto-detect ID type
    if (id_type == "auto") {
        if (any(grepl("^ENSG", gene_ids, ignore.case = TRUE))) {
            id_type <- "ensembl"
        } else {
            id_type <- "entrez"
        }
        message("Detected gene ID type: ", id_type)
    }

    # Convert to symbols
    if (id_type == "ensembl") {
        symbols <- ensembl_to_symbol(gene_ids)
    } else {
        symbols <- entrez_to_symbol(as.character(gene_ids))
    }

    # Add to results
    results$gene_symbol <- symbols[gene_ids]

    return(results)
}

#' Annotate Results Using biomaRt
#'
#' Annotate CNV-expression results with comprehensive gene information
#' using biomaRt.
#'
#' @param results Data frame with CNV-expression results
#' @param chrom_col Column with chromosome information
#' @param start_col Column with start position
#' @param end_col Column with end position
#' @param genome Reference genome: "GRCh38" or "GRCh37"
#'
#' @return Data frame with added gene annotations
#' @export
#'
#' @examples
#' \dontrun{
#' annotated <- annotate_with_biomart(results)
#' }
annotate_with_biomart <- function(results,
                                   chrom_col = "gene_chrom",
                                   start_col = "gene_start",
                                   end_col = "gene_end",
                                   genome = c("GRCh38", "GRCh37")) {

    genome <- match.arg(genome)

    if (!requireNamespace("biomaRt", quietly = TRUE)) {
        stop("Package 'biomaRt' required. Install with: BiocManager::install('biomaRt')")
    }

    # Check columns
    required_cols <- c(chrom_col, start_col, end_col)
    missing_cols <- setdiff(required_cols, colnames(results))
    if (length(missing_cols) > 0) {
        stop("Missing columns: ", paste(missing_cols, collapse = ", "))
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

    # Query for each region
    results$ensembl_id <- NA_character_
    results$gene_symbol <- NA_character_
    results$gene_description <- NA_character_

    for (i in seq_len(nrow(results))) {
        chrom <- gsub("chr", "", results[[chrom_col]][i])
        start <- results[[start_col]][i]
        end <- results[[end_col]][i]

        tryCatch({
            gene_info <- biomaRt::getBM(
                attributes = c("ensembl_gene_id", "external_gene_name", "description"),
                filters = c("chromosome_name", "start", "end"),
                values = list(chrom, start, end),
                mart = ensembl
            )

            if (nrow(gene_info) > 0) {
                results$ensembl_id[i] <- paste(gene_info$ensembl_gene_id, collapse = ";")
                results$gene_symbol[i] <- paste(gene_info$external_gene_name, collapse = ";")
                results$gene_description[i] <- gene_info$description[1]
            }
        }, error = function(e) {
            # Continue on error
        })
    }

    return(results)
}

#' Get Gene Annotation Summary
#'
#' Generate summary statistics for gene annotations.
#'
#' @param results Annotated results data frame
#' @param symbol_col Column with gene symbols
#'
#' @return List with annotation summary
#' @export
#'
#' @examples
#' \dontrun{
#' summary <- get_annotation_summary(annotated_results)
#' }
get_annotation_summary <- function(results, symbol_col = "gene_symbol") {

    total <- nrow(results)

    if (!symbol_col %in% colnames(results)) {
        return(list(
            total_entries = total,
            annotated = NA,
            not_annotated = NA,
            unique_genes = NA
        ))
    }

    annotated <- sum(!is.na(results[[symbol_col]]) & results[[symbol_col]] != "")
    unique_genes <- length(unique(results[[symbol_col]][!is.na(results[[symbol_col]])]))

    list(
        total_entries = total,
        annotated = annotated,
        not_annotated = total - annotated,
        annotation_rate = annotated / total,
        unique_genes = unique_genes
    )
}
