#' peakGeneCor
#'
#' @import BiocParallel
#' @importFrom stats cor.test p.adjust
#'
#' @param mmAnno the annotated GRange object from mm_geneScan or mm_nearestGene
#' @param peakScoreMt peak count matrix. The rownames are feature_id in mmAnno,
#' while the colnames are sample names.
#' @param geneScoreMt gene count matirx. The rownames are gene_id in mmAnno,
#' while the colnames are sample names.
#' @param parallel whehter you want to using bplapply to speed up calculation
#'
#' @return mmAnno with Cor, pvalue,padj,qvalue column
#' @export
#'
#' @examples
#' data("RNA_normCount")
#' data("ATAC_normCount")
#'
#' library(TxDb.Athaliana.BioMart.plantsmart28)
#' Txdb <- TxDb.Athaliana.BioMart.plantsmart28
#' seqlevels(Txdb) <- c(paste0("Chr", 1:5), "M", "C")
#' peak_path <- system.file("extdata", "ATAC.bed.gz", package = "FindIT2")
#' peak_GR <- loadPeakFile(peak_path)[1:100]
#' mmAnno <- mm_geneScan(peak_GR,Txdb)
#'
#' ATAC_colData <- data.frame(row.names = colnames(ATAC_normCount),
#'                            type = gsub("_R[0-9]", "", colnames(ATAC_normCount))
#'                            )
#'
#' integrate_replicates(ATAC_normCount, ATAC_colData) -> ATAC_normCount_merge
#' RNA_colData <- data.frame(row.names = colnames(RNA_normCount),
#'                           type = gsub("_R[0-9]", "", colnames(RNA_normCount))
#'                           )
#' integrate_replicates(RNA_normCount, RNA_colData) -> RNA_normCount_merge
#' peakGeneCor(mmAnno = mmAnno,
#'             peakScoreMt = ATAC_normCount_merge,
#'             geneScoreMt = RNA_normCount_merge,
#'             parallel = FALSE) -> mmAnnoCor
#'
#' mmAnnoCor
#'
peakGeneCor <- function(mmAnno,
                        peakScoreMt,
                        geneScoreMt,
                        parallel = FALSE) {

    # check colnames match
    if (all(colnames(peakScoreMt) == colnames(geneScoreMt))) {
        message("Good, your two matrix colnames matchs")
    } else {
        stop("Sorry, your two matrix colnames do not match",
             call. = FALSE
        )
    }

    # if someone using gff as Txdb, then will appear miRNA,lncRNA gene
    # which will not appear in geneScoreMt
    withr::local_options(list(warn = 1))
    if (mean(mmAnno$gene_id %in% rownames(geneScoreMt)) < 1 |
        mean(mmAnno$feature_id %in% rownames(peakScoreMt)) < 1) {

        left_index <- mmAnno$gene_id %in% rownames(geneScoreMt) &
            mmAnno$feature_id %in% rownames(peakScoreMt)

        mmAnno_left <- mmAnno[left_index]
        mmAnno_drop <- mmAnno[!left_index]
        msg <- paste("some gene_id or feature_id in your mmAnno is not in your",
                     "geneScoreMt or peakScore Mt,",
                     "final cor and pvalue of these gene_id or feature_id pair will be NA")

        warning(msg, call. = FALSE)
    } else {
        mmAnno_left <- mmAnno
        mmAnno_drop <- GRanges()
    }


    cat(
        ">> calculating cor and pvalue, which may be time consuming...\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n"
    )

    mt_gene <- geneScoreMt[mmAnno_left$gene_id, ]
    mt_peak <- peakScoreMt[mmAnno_left$feature_id, ]

    if (parallel) {
        # print_msgParallel()
        # though cor(t(mt_gene["gene",]), t(mt_peak["peak,]))
        # may quickyly produce the whole result
        # but it will produce a very very huge matrix
        suppressWarnings(bplapply(seq_len(dim(mt_gene)[1]), function(index) {
            cor <- suppressWarnings(cor.test(
                mt_gene[index, ],
                mt_peak[index, ]
            ))

            return(c(cor$estimate, cor$p.value))
        })) -> cor_result
    } else {
        suppressWarnings(lapply(seq_len(dim(mt_gene)[1]), function(index) {
            cor <- suppressWarnings(cor.test(
                mt_gene[index, ],
                mt_peak[index, ]
            ))

            return(c(cor$estimate, cor$p.value))
        })) -> cor_result
    }

    cor_result <- unlist(cor_result)
    cor_mt <- matrix(cor_result, ncol = 2, byrow = TRUE)

    cat(
        ">> merging all info together...\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n"
    )

    mmAnno_left$cor <- cor_mt[, 1]
    mmAnno_left$pvalue <- cor_mt[, 2]
    mmAnno_left$p_adj <- p.adjust(mmAnno_left$pvalue)

    # compare with BH, qvalue is more soft
    qvalue_result <- tryCatch(qvalue::qvalue(
        p = mmAnno_left$pvalue,
        fdr.level = 0.05,
        pi0.method = "bootstrap"
    ),
    error = function(e) NULL
    )

    if (class(qvalue_result) == "qvalue") {
        qvalues <- qvalue_result$qvalues
    } else {
        qvalues <- NA
    }
    mmAnno_left$qvalue <- qvalues

    mmAnno <- sort(c(mmAnno_drop, mmAnno_left))
    metadata(mmAnno)$peakScoreMt <- peakScoreMt
    metadata(mmAnno)$geneScoreMt <- geneScoreMt

    return(mmAnno)
}
