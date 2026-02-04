# Tests for CNV analysis functions

test_that("convert_segment_mean_to_cn calculates correct copy numbers", {
    # Create test data
    df <- data.frame(
        sample = c("S1", "S1", "S1", "S1", "S1"),
        chrom = c("1", "2", "3", "4", "5"),
        start = c(1000, 2000, 3000, 4000, 5000),
        end = c(1500, 2500, 3500, 4500, 5500),
        segment_mean = c(-3, -1, 0, 0.58, 1)
    )

    result <- convert_segment_mean_to_cn(df)

    # Check that columns were added
    expect_true("copy_number" %in% colnames(result))
    expect_true("annotation" %in% colnames(result))

    # Check copy number calculation: CN = round(2 * 2^SM)
    # SM = -3: CN = round(2 * 2^-3) = round(0.25) = 0
    # SM = -1: CN = round(2 * 2^-1) = round(1) = 1
    # SM = 0: CN = round(2 * 2^0) = round(2) = 2
    # SM = 0.58: CN = round(2 * 2^0.58) = round(3) = 3
    # SM = 1: CN = round(2 * 2^1) = round(4) = 4
    expect_equal(result$copy_number, c(0, 1, 2, 3, 4))

    # Check annotations
    expect_equal(result$annotation[1], "homozygous_deletion")
    expect_equal(result$annotation[2], "heterozygous_deletion")
    expect_equal(result$annotation[3], "neutral")
    expect_equal(result$annotation[4], "amplification")
    expect_equal(result$annotation[5], "amplification")
})

test_that("remove_sex_chromosomes filters correctly", {
    df <- data.frame(
        chrom = c("1", "2", "X", "Y", "chr1", "chrX", "chrY")
    )

    result <- remove_sex_chromosomes(df)

    expect_equal(nrow(result), 3)
    expect_false(any(result$chrom %in% c("X", "Y", "chrX", "chrY")))
})

test_that("filter_diploid removes CN=2 entries", {
    df <- data.frame(
        copy_number = c(0, 1, 2, 2, 3, 4)
    )

    result <- filter_diploid(df)

    expect_equal(nrow(result), 4)
    expect_false(any(result$copy_number == 2))
})

test_that("add_state_column caps values at 4", {
    df <- data.frame(
        copy_number = c(0, 1, 2, 3, 4, 5, 6, 10)
    )

    result <- add_state_column(df)

    expect_equal(result$state, c(0, 1, 2, 3, 4, 4, 4, 4))
})

test_that("standardize_chromosomes adds chr prefix", {
    df <- data.frame(
        chrom = c("1", "2", "chr3", "X")
    )

    result <- standardize_chromosomes(df)

    expect_equal(result$chrom, c("chr1", "chr2", "chr3", "chrX"))
})

test_that("prepare_cnv_data performs complete preprocessing", {
    df <- data.frame(
        sample = rep("S1", 6),
        chrom = c("1", "2", "3", "X", "Y", "4"),
        start = 1:6 * 1000,
        end = 1:6 * 1000 + 500,
        segment_mean = c(-2, -0.5, 0, 0.3, 0.6, 1)
    )

    result <- prepare_cnv_data(df, remove_diploid = TRUE, remove_sex = TRUE)

    # Should have removed sex chromosomes
    expect_false(any(result$chrom %in% c("chrX", "chrY")))

    # Should have removed diploid
    expect_false(any(result$copy_number == 2))

    # Should have copy_number, annotation, and state columns
    expect_true(all(c("copy_number", "annotation", "state") %in% colnames(result)))
})
