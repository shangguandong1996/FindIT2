suppressPackageStartupMessages(library(TxDb.Athaliana.BioMart.plantsmart28))
Txdb <- TxDb.Athaliana.BioMart.plantsmart28
seqlevels(Txdb) <- paste0("Chr", c(1:5, "M", "C"))

peak_path <- system.file("extdata", "ChIP.bed.gz", package = "FindIT2")
peak_GR <- loadPeakFile(peak_path)[1:100]
peakAnno <- mm_nearestGene(peak_GR, Txdb)
peakAnno_scan <- mm_geneScan(peak_GR, Txdb)

test_that("getAssocPairNumber test", {
    expect_equal(getAssocPairNumber(peakAnno)$peakNumber[1], 1)

    expect_equal(
        getAssocPairNumber(peakAnno, output_type = "feature_id")$geneNumber[1], 1
        )

    expect_equal(
        getAssocPairNumber(peakAnno_scan, output_summary = TRUE)$gene_freq[6], 2
        )

})

# how to test ggplot2 ?
test_that("plot_annoDistance test", {
    expect_error(
        plot_annoDistance(peakAnno_scan),
        "sorry, it only accept mmAnno from mm_nearestGene"
    )


})

