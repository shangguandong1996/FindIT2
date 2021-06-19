suppressPackageStartupMessages(library(TxDb.Athaliana.BioMart.plantsmart28))
Txdb <- TxDb.Athaliana.BioMart.plantsmart28
seqlevels(Txdb) <- paste0("Chr", c(1:5, "M", "C"))

data("RNA_normCount")
data("ATAC_normCount")
peak_path <- system.file("extdata", "ATAC.bed.gz", package = "FindIT2")
peak_GR <- loadPeakFile(peak_path)[1:100]
mmAnno <- quiet(mm_geneScan(peak_GR, Txdb))
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

suppressWarnings(quiet(peakGeneCor(
    mmAnno = mmAnno,
    peakScoreMt = ATAC_normCount_merge,
    geneScoreMt = RNA_normCount_merge,
    parallel = FALSE
))) -> mmAnnoCor

quiet(
    enhancerPromoterCor(
        peak_GR = peak_GR,
        Txdb = Txdb,
        peakScoreMt = ATAC_normCount,
        parallel = FALSE)
    ) -> mm_ePLink


test_that("peakGeneCor test",{
    expect_warning(
        peakGeneCor(
            mmAnno = mmAnno,
            peakScoreMt = ATAC_normCount_merge,
            geneScoreMt = RNA_normCount_merge,
            parallel = FALSE
        ),

        "final cor and pvalue of these gene_id or feature_id pair will be NA"
    )
    expect_equal(mmAnnoCor$feature_id[100], "WFX_peak_19682")
    expect_equal(mmAnnoCor$cor[1], -0.24446226)
    expect_equal(mmAnnoCor$qvalue[1], 0.7149576)
})


test_that("enhancerPromoter test",{
    expect_equal(mm_ePLink$feature_id[20], "WFX_peak_19632")
    expect_equal(mm_ePLink$cor[1], 0.62563266)
    expect_equal(mm_ePLink$qvalue[1], 0.045637884)
})

