suppressPackageStartupMessages(library(TxDb.Athaliana.BioMart.plantsmart28))
Txdb <- TxDb.Athaliana.BioMart.plantsmart28
seqlevels(Txdb) <- paste0("Chr", c(1:5, "M", "C"))

peak_path <- system.file("extdata", "ChIP.bed.gz", package = "FindIT2")
peak_GR <- loadPeakFile(peak_path)[1:100]
peakAnno <- mm_nearestGene(peak_GR, Txdb, verbose = FALSE)
peakAnno_scan <- mm_geneScan(peak_GR, Txdb, verbose = FALSE)

test_that("getAssocPairNumber test", {
    expect_equal(getAssocPairNumber(peakAnno)$peakNumber[1], 1)

    expect_equal(
        getAssocPairNumber(peakAnno, output_type = "feature_id")$geneNumber[1], 1
        )

    expect_equal(
        getAssocPairNumber(peakAnno_scan, output_summary = TRUE)$gene_freq[6], 2
        )

})

test_that("plot_annoDistance test", {
    expect_error(
        plot_annoDistance(peakAnno_scan),
        "sorry, it only accept mmAnno from mm_nearestGene"
    )

    p <- plot_annoDistance(peakAnno)
    pb_1 <- ggplot2::ggplot_build(p[[1]])
    pb_2 <- ggplot2::ggplot_build(p[[2]])


    expect_equal(pb_1$data[[1]]$density[1], 3.372373, tolerance = 1e-5)
    expect_match(pb_1$data[[2]]$label, "\nMedian :  200.5")

    expect_equal(pb_2$data[[1]]$density[250], 0.0001044554, tolerance = 1e-5)
    expect_match(pb_2$data[[2]]$label[2], "\nMin. :  -3165")

})

