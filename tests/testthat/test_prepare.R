test_that("loadPeakFile can produce GRange with feature_id column", {
    peakfile <- system.file("extdata", "ChIP.bed.gz", package = "FindIT2")
    peak_GR <- loadPeakFile(peakfile)
    expect_s4_class(peak_GR, "GRanges")
    expect_match(paste(colnames(mcols(peak_GR)), collapse = " "), "feature_id")
})
