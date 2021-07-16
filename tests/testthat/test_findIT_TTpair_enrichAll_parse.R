data("TF_target_database")
data("test_geneSet")
data("test_featureSet")
peak_path <- system.file("extdata", "ATAC.bed.gz", package = "FindIT2")
peak_GR <- loadPeakFile(peak_path)
ChIP_peak_path <- system.file("extdata", "ChIP.bed.gz", package = "FindIT2")
ChIP_peak_GR <- loadPeakFile(ChIP_peak_path)
ChIP_peak_GR$TF_id <- "AT1G28300"

TF_GR_database_path <- system.file("extdata",
                                   "TF_GR_database.bed.gz",
                                   package = "FindIT2")
TF_GR_database <- loadPeakFile(TF_GR_database_path)
TF_GR_database$TF_id <- TF_GR_database$feature_id

test_that("jaccard findIT TTpair test",{
    result_findIT_TTPair <- findIT_TTPair(
        input_genes = test_geneSet,
        TF_target_database = TF_target_database
    )

    jaccard_data <- jaccard_findIT_TTpair(
        input_genes = test_geneSet,
        TF_target_database = TF_target_database,
        input_TF_id = result_findIT_TTPair$TF_id[1:3]
    )

    expect_true(is.matrix(jaccard_data))

    expect_equal(jaccard_data[1, 2], 0.1363636, tolerance = 1e-5)
    expect_equal(jaccard_data[1, 1], 0)


})

test_that("jaccard findIT enrich in All test",{

    expect_true(is.matrix(
        jaccard_findIT_enrichFisher(
            input_feature_id = test_featureSet,
            peak_GR = peak_GR,
            TF_GR_database = ChIP_peak_GR,
            input_TF_id = "AT1G28300"
        )
    ))

    jaccard_result <- jaccard_findIT_enrichFisher(
        input_feature_id = test_featureSet,
        peak_GR = peak_GR,
        TF_GR_database = TF_GR_database,
        input_TF_id = c("AT2G36270", "AT3G59060", "AT5G24110")
    )

    expect_equal(jaccard_result[1, 1], 0)
    expect_equal(jaccard_result[1, 2], 0.2857143, tolerance = 1e-5)

})

