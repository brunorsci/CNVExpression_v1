# Tests for integration functions

test_that("filter_cnv_expression_results filters by p-value", {
    results <- data.frame(
        gene = c("A", "B", "C", "D"),
        logFC.CN1 = c(1.5, 0.5, -2.0, 0.3),
        logFC.CN3 = c(-0.5, 1.2, 0.8, -1.5),
        AdjPValue = c(0.01, 0.1, 0.03, 0.2)
    )

    filtered <- filter_cnv_expression_results(
        results,
        logfc_cutoff = 1,
        pvalue_cutoff = 0.05
    )

    # Should keep A (logFC.CN1 = 1.5, p = 0.01) and C (logFC.CN1 = -2.0, p = 0.03)
    expect_equal(nrow(filtered), 2)
    expect_true(all(filtered$AdjPValue <= 0.05))
})

test_that("filter_cnv_expression_results handles missing logFC columns", {
    results <- data.frame(
        gene = c("A", "B"),
        other_col = c(1, 2),
        AdjPValue = c(0.01, 0.1)
    )

    expect_warning(
        filter_cnv_expression_results(results, logfc_cutoff = 1, pvalue_cutoff = 0.05),
        "No logFC columns found"
    )
})
