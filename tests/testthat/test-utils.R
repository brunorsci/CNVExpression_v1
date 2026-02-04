# Tests for utility functions

test_that("create_sample_gene_id creates correct identifiers", {
    result <- create_sample_gene_id(c("S1", "S2"), c("TP53", "BRCA1"))
    expect_equal(result, c("S1_TP53", "S2_BRCA1"))
})

test_that("validate_input_data checks for data frame", {
    expect_error(validate_input_data("not a data frame"))
})

test_that("validate_input_data checks for empty data frame", {
    empty_df <- data.frame()
    expect_error(validate_input_data(empty_df))
})

test_that("validate_input_data checks required columns", {
    df <- data.frame(a = 1, b = 2)
    expect_error(validate_input_data(df, required_cols = c("a", "c")))
    expect_true(validate_input_data(df, required_cols = c("a", "b")))
})

test_that("check_dependencies returns TRUE for installed packages", {
    result <- check_dependencies(c("stats", "methods"))
    expect_true(result)
})

test_that("check_dependencies warns for missing packages", {
    expect_warning(
        check_dependencies(c("nonexistent_package_12345")),
        "Missing packages"
    )
})
