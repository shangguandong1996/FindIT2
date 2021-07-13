suppressPackageStartupMessages(library(TxDb.Athaliana.BioMart.plantsmart28))
Txdb <- TxDb.Athaliana.BioMart.plantsmart28
seqlevels(Txdb) <- paste0("Chr", c(1:5, "M", "C"))

data("RNADiff_LEC2_GR")

peak_path <- system.file("extdata", "ChIP.bed.gz", package = "FindIT2")
peak_GR <- loadPeakFile(peak_path)


test_that("integrate ChIP RNA test", {
    mmAnno <- mm_geneScan(peak_GR, Txdb)

    calcRP_TFHit(
        mmAnno = mmAnno,
        Txdb = Txdb
    ) -> result_geneRP

    integrate_ChIP_RNA(
        result_geneRP = result_geneRP,
        result_geneDiff = RNADiff_LEC2_GR
    ) -> merge_data

    pb <- ggplot2::ggplot_build(merge_data)

    expect_s3_class(merge_data, "ggplot")
    expect_equal(merge_data$data$rankProduct[1], 248)

    # test Ks test whether correct
    expect_match(pb$data[[2]]$label,
                 "1.86e-12")

    expect_warning(
        integrate_ChIP_RNA(
            result_geneRP = result_geneRP,
            result_geneDiff = RNADiff_LEC2_GR,
            lfc_threshold = 2
        ),
        "no significant down genes, just returing rank product result"
    )

})

set.seed(19960203)
mt <- matrix(runif(100, 0, 100), nrow = 10)
colnames(mt) <- paste0(paste0("type", 1:5), "_", rep(1:2, 5))
rownames(mt) <- paste0("TF", 1:10)
colData <- data.frame(
    type = gsub("_[0-9]", "", colnames(mt)),
    row.names = colnames(mt)
)

test_that("integrate replicates test",{
    expect_equal(
        integrate_replicates(mt, colData, type = "value")[1, 2],
        33.36415,
        tolerance = 1e-5
        )

    expect_equal(
        integrate_replicates(mt, colData, type = "rank")[1, 2],
        4,
        tolerance = 1e-5
    )

    mt <- apply(mt, 2, INT)
    expect_equal(
        integrate_replicates(mt, colData, type = "rank_zscore")[1, 2],
        -0.8217251,
        tolerance = 1e-5
    )

})
