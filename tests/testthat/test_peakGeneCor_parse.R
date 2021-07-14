suppressPackageStartupMessages(library(TxDb.Athaliana.BioMart.plantsmart28))
Txdb <- TxDb.Athaliana.BioMart.plantsmart28
seqlevels(Txdb) <- paste0("Chr", c(1:5, "M", "C"))

data("RNA_normCount")
data("ATAC_normCount")
peak_path <- system.file("extdata", "ATAC.bed.gz", package = "FindIT2")
peak_GR <- loadPeakFile(peak_path)[1:100]
mmAnno <- mm_geneScan(peak_GR, Txdb)
ATAC_colData <- data.frame(
    row.names = colnames(ATAC_normCount),
    type = gsub("_R[0-9]", "", colnames(ATAC_normCount))
)
integrate_replicates(ATAC_normCount, ATAC_colData) -> ATAC_normCount_merge
RNA_colData <- data.frame(
    row.names = colnames(RNA_normCount),
    type = gsub("_R[0-9]", "", colnames(RNA_normCount))
)
integrate_replicates(RNA_normCount, RNA_colData) -> RNA_normCount_merge

mmAnnoCor <- suppressWarnings(peakGeneCor(
    mmAnno = mmAnno,
    peakScoreMt = ATAC_normCount_merge,
    geneScoreMt = RNA_normCount_merge,
    parallel = FALSE
))


mm_ePLink <- enhancerPromoterCor(
    peak_GR = peak_GR,
    Txdb = Txdb,
    peakScoreMt = ATAC_normCount,
    parallel = FALSE)

test_that("plot_peakGeneCor test",{
    p <- plot_peakGeneCor(mmAnnoCor, select_gene = "AT5G01010")
    pb <- ggplot2::ggplot_build(p)
    expect_s3_class(p, "ggplot")
    expect_equal(pb$data[[1]]$x[1], 91.89544, tolerance = 1e-5)
    # there 4 facet
    expect_equal(nrow(pb$data[[3]]), 4)

    p <- plot_peakGeneCor(mm_ePLink, select_gene = "AT5G01010")
    pb <- ggplot2::ggplot_build(p)
    expect_s3_class(p, "ggplot")
    expect_equal(pb$data[[1]]$x[1], 96.44832, tolerance = 1e-5)
    # there 8 facet
    expect_equal(nrow(pb$data[[3]]), 8)
})
