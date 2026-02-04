# Tests for data processing functions

test_that("standardize_sample_ids works correctly", {
    samples <- c(
        "TCGA-AB-1234-01A-11R-A12B-07",
        "TCGA-AB-5678-01A-11R-A12B-07"
    )

    result <- standardize_sample_ids(samples)

    expect_equal(result, c("TCGA_AB_1234_01", "TCGA_AB_5678_01"))
})

test_that("standardize_sample_ids handles already standardized IDs", {
    samples <- c("TCGA_AB_1234_01", "TCGA_AB_5678_01")

    result <- standardize_sample_ids(samples)

    expect_equal(result, samples)
})

test_that("match_samples finds common samples", {
    cnv_samples <- c("S1", "S2", "S3", "S4")
    expr_samples <- c("S2", "S3", "S5", "S6")

    result <- match_samples(cnv_samples, expr_samples, standardize = FALSE)

    expect_equal(result$n_common, 2)
    expect_equal(sort(result$common_samples), c("S2", "S3"))
    expect_equal(result$n_cnv_only, 2)
    expect_equal(result$n_expr_only, 2)
})

test_that("create_sample_gene_id creates correct identifiers", {
    result <- create_sample_gene_id(c("S1", "S2"), c("TP53", "BRCA1"))

    expect_equal(result, c("S1_TP53", "S2_BRCA1"))
})

test_that("create_sample_gene_id respects custom separator", {
    result <- create_sample_gene_id(c("S1", "S2"), c("TP53", "BRCA1"), sep = "-")

    expect_equal(result, c("S1-TP53", "S2-BRCA1"))
})
