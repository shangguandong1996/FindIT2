suppressPackageStartupMessages(library(TxDb.Athaliana.BioMart.plantsmart28))

data("ATAC_normCount")
data("test_geneSet")
data("test_featureSet")
data("TF_target_database")

Txdb <- TxDb.Athaliana.BioMart.plantsmart28
seqlevels(Txdb) <- paste0("Chr", c(1:5, "M", "C"))

peak_path <- system.file("extdata", "ATAC.bed.gz", package = "FindIT2")
peak_GR <- loadPeakFile(peak_path)

ChIP_peak_path <- system.file("extdata", "ChIP.bed.gz", package = "FindIT2")
ChIP_peak_GR <- loadPeakFile(ChIP_peak_path)
ChIP_peak_GR$TF_id <- "AT1G28300"

TF_GR_database_path <- system.file("extdata", "TF_GR_database.bed.gz", package = "FindIT2")
TF_GR_database <- loadPeakFile(TF_GR_database_path)
colnames(mcols(TF_GR_database))[1] <- "TF_id"
TF_GR_database <- c(TF_GR_database, ChIP_peak_GR)

test_that("findIT regionRP test",{
    mmAnno <- mm_geneScan(peak_GR, Txdb)

    regionRP <- calcRP_region(
        mmAnno = mmAnno,
        peakScoreMt = ATAC_normCount,
        Txdb = Txdb,
        Chrs_included = "Chr5"
    )

    set.seed(20160806)

    result_findIT_regionRP <- findIT_regionRP(
        regionRP = regionRP,
        Txdb = Txdb,
        TF_GR_database = ChIP_peak_GR,
        input_genes = test_geneSet,
        background_number = 3000
    )

    expect_equal(assay(result_findIT_regionRP)[1], 0.47925271)
    expect_equal(
        log10(assays(result_findIT_regionRP)$TF_pvalue[1]),
        -15.85757
    )

})


test_that("findIT TTpair test", {
    result_findIT_TTPair <- findIT_TTPair(
        input_genes = test_geneSet,
        TF_target_database = TF_target_database
    )

    expect_equal(result_findIT_TTPair$qvalue[1], 0.030707091)
})

test_that("findIT TFHit test", {
    set.seed(20160806)
    result_findIT_TFHit <- findIT_TFHit(
        input_genes = test_geneSet,
        Txdb = Txdb,
        TF_GR_database = ChIP_peak_GR,
        scan_dist = 2e4,
        decay_dist = 1e3
    )

    expect_equal(
        log10(result_findIT_TFHit$pvalue[1]),
        -62.56263,
        tolerance = 0.0001)
})


test_that("findIT enrich In shuffle test", {
    set.seed(20160806)

    result_findIT_enrichInShuffle <- findIT_enrichInShuffle(
        input_feature_id = test_featureSet,
        peak_GR = peak_GR,
        TF_GR_database = ChIP_peak_GR,
        shuffleN = 10
    )

    expect_equal(result_findIT_enrichInShuffle$shuffleMean_hitN[1], 41.1)
})


test_that("findIT enrich In All test", {
    result_findIT_enrichInAll <- findIT_enrichInAll(
        input_feature_id = test_featureSet,
        peak_GR = peak_GR,
        TF_GR_database = ChIP_peak_GR
    )

    expect_equal(result_findIT_enrichInAll$inputRatio[1], "64/77")
    expect_equal(is.na(result_findIT_enrichInAll$qvalue[1]), TRUE)
})


test_that("findIT MARA test", {

    set.seed(20160806)
    result_findIT_MARA <- findIT_MARA(
        input_feature_id = test_featureSet,
        peak_GR = peak_GR,
        peakScoreMt = ATAC_normCount,
        TF_GR_database = TF_GR_database
    )

    expect_equal(result_findIT_MARA$E5_0h_R1[1], 1.1503494, tolerance = 0.0001)
})
