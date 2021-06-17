suppressPackageStartupMessages(library(TxDb.Athaliana.BioMart.plantsmart28))
Txdb <- TxDb.Athaliana.BioMart.plantsmart28
seqlevels(Txdb) <- paste0("Chr", c(1:5, "M", "C"))

peak_path <- system.file("extdata", "ChIP.bed.gz", package = "FindIT2")
peak_GR <- loadPeakFile(peak_path)[1:100]

test_that("mmAnno nearest mode test", {
    seqlevels(peak_GR) <- "5"
    expect_error(
        quiet(mm_nearestGene(peak_GR = peak_GR,Txdb = Txdb)),
        "have no sequence levels in common"
    )

    seqlevels(peak_GR) <- "Chr5"
    peak_GR_add <- suppressWarnings(c(peak_GR, GRanges(seqnames = "Chr6",
                                                       ranges = "1-5",
                                                       feature_id = "peak_0",
                                                       score = 1)))


    expect_warning(
        quiet(mm_nearestGene(peak_GR = peak_GR_add,Txdb = Txdb)),
        "some peak's Chr is nor in your Txdb, for example: Chr6"
    )

    mmAnno <- suppressMessages(quiet(
        mm_nearestGene(peak_GR = peak_GR,Txdb = Txdb)
        ))

    expect_equal(mmAnno$distanceToTSS[1], -344)
    expect_equal(mmAnno$gene_id[5], "AT5G01040")

})


test_that("mmAnno geneScan mode test", {

    mmAnno <- suppressMessages(quiet(
        mm_geneScan(peak_GR = peak_GR, Txdb = Txdb)
    ))

    expect_equal(mmAnno$distanceToTSS[1], -1174)
    expect_equal(mmAnno$gene_id[1], "AT5G01010")

})


test_that("mmAnno gene Bound mode test", {
    expect_message(
        quiet(mm_geneBound(peak_GR, Txdb, c("AT5G01015"))),
        "all your input gene have been annotated by nearestGene mode")

    expect_error(
        quiet(mm_geneBound(peak_GR, Txdb, c("AT2G17950"))),
        "sorry, all of your input genes all not in your peak_GR chrosome")

    expect_message(
        quiet(mm_geneBound(peak_GR, Txdb, c("AT5G01015", "AT5G67570"))),
        "It seems that there 1 genes")

    expect_warning(
        quiet(mm_geneBound(peak_GR, Txdb, c("AT5G01015", "AT5G67570", "AT2G17950"))),
        "some of your input gene are not in your peak_GR chrosome")

    mmAnno <- mm_geneBound(peak_GR = peak_GR,
                           Txdb = Txdb,
                           input_genes = c("AT5G01015", "AT5G67570"))

    expect_equal(mmAnno$distanceToTSS_abs[2], 26516130)


})
