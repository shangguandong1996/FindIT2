suppressPackageStartupMessages(library(TxDb.Athaliana.BioMart.plantsmart28))
suppressPackageStartupMessages(library(SummarizedExperiment))
Txdb <- TxDb.Athaliana.BioMart.plantsmart28
seqlevels(Txdb) <- paste0("Chr", c(1:5, "M", "C"))

bwFile <- system.file("extdata", "E50h_sampleChr5.bw", package = "FindIT2")

if (.Platform$OS.type != "windows") {
    test_that("calcRP_coverage test",{
        quiet(RP_df <- calcRP_coverage(
            bwFile = bwFile,
            Txdb = Txdb,
            Chrs_included = "Chr5"
        ))

        expect_equal(RP_df$sumRP[2], 0.65769128)
        expect_equal(RP_df$gene_id[2], "AT5G01060")
    })
}

test_that("calcRP_region test",{
    data("ATAC_normCount")
    peak_path <- system.file("extdata", "ATAC.bed.gz", package = "FindIT2")
    peak_GR <- loadPeakFile(peak_path)
    quiet(mmAnno <- mm_geneScan(peak_GR, Txdb))
    quiet(calcRP_region(
        mmAnno = mmAnno,
        peakScoreMt = ATAC_normCount,
        Txdb = Txdb,
        Chrs_included = "Chr5"
    ) -> regionRP)

    expect_s4_class(regionRP, "MultiAssayExperiment")
    expect_equal(assays(regionRP)$sumRP[5,5], 235.44406)
    expect_equal(assays(regionRP)$fullRP[5,5], 138.69528)

    quiet(calcRP_region(
        mmAnno = mmAnno,
        peakScoreMt = ATAC_normCount,
        Txdb = Txdb,
        Chrs_included = "Chr5",
        log_transform = TRUE
    ) -> regionRP)

    expect_equal(assays(regionRP)$sumRP[5,5], 1.12383564)

})


test_that("calcRP_TFHit", {
    peak_path <- system.file("extdata", "ChIP.bed.gz", package = "FindIT2")
    peak_GR <- loadPeakFile(peak_path)
    quiet(mmAnno <- mm_geneScan(peak_GR, Txdb))

    quiet(fullRP_hit <- calcRP_TFHit(
        mmAnno = mmAnno,
        Txdb = Txdb,
        report_fullInfo = TRUE
    ))

    expect_equal(fullRP_hit$RP[100], 0.65474271)

    quiet(geneRP_df <- calcRP_TFHit(
        mmAnno = mmAnno,
        Txdb = Txdb,
        report_fullInfo = FALSE
    ))

    expect_equal(geneRP_df$withPeakN[100], 4)

})
