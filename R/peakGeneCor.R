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


#' enhancerPromoterCor
#'
#' @param peak_GR peak GRange with a column named feature_id representing you peak name
#' @param Txdb Txdb
#' @param up_scanPromoter the scan distance which is used to scan nearest promoter
#' @param down_scanPromoter the scan distance which is used to scan nearest promoter
#' @param up_scanEnhancer the scan distance which is used to scan feature
#' @param down_scanEnhacner the scan distance which is used to scan feature
#' @param peakScoreMt peak count matrix. The rownames are feature_id in peak_GR
#'
#' @return mmAnno with Cor, pvalue,padj,qvalue column
#' @export
#'
#' @examples
enhancerPromoterCor <- function(peak_GR,
                                Txdb,
                                up_scanPromoter = 500,
                                down_scanPromoter = 500,
                                up_scanEnhancer = 2e4,
                                down_scanEnhacner = 2e4,
                                peakScoreMt,
                                parallel){

    mm_geneScan(peak_GR = peak_GR,
                Txdb = Txdb,
                upstream = up_scanPromoter,
                downstream = down_scanPromoter) -> mm_promoter

    mm_promoter %>%
        data.frame() %>%
        dplyr::select(feature_id, gene_id, distanceToTSS) %>%
        dplyr::mutate(distanceToTSS_abs = abs(distanceToTSS)) %>%
        dplyr::group_by(gene_id) %>%
        dplyr::mutate(
            promoter_feature = feature_id[which.min(distanceToTSS_abs)],
            peak_gene_pair = paste0(feature_id, ":", gene_id)
            ) %>%
        dplyr::ungroup() -> mm_promoter_pair

    pair_info <- mm_promoter_pair %>%
        dplyr::select(gene_id, promoter_feature) %>%
        dplyr::distinct() %>%
        dplyr::mutate(
            peak_gene_pair = paste0(promoter_feature, ":", gene_id)
        )

    mm_promoter_pair %>%
        dplyr::filter(peak_gene_pair %in% pair_info$peak_gene_pair) %>%
        dplyr::select(gene_id, distanceToTSS, promoter_feature) -> mm_promoter_tidy


    mm_geneScan(peak_GR = peak_GR,
                Txdb = Txdb,
                upstream = up_scanEnhancer,
                downstream = down_scanEnhacner) -> mm_scan

    # some gene may not have promoter(no open peak) in here
    mm_scan <- subset(mm_scan, gene_id %in% mm_promoter_tidy$gene_id)

    mcols(mm_scan) %>%
        data.frame() %>%
        dplyr::left_join(mm_promoter_tidy[, c(1, 3)]) -> enhancer_peak_pair


    mt_promoter <- peakScoreMt[enhancer_peak_pair$promoter_feature, ]
    mt_enhancer <- peakScoreMt[enhancer_peak_pair$feature_id, ]


    if (parallel) {
        # print_msgParallel()
        # though cor(t(mt_gene["gene",]), t(mt_peak["peak,]))
        # may quickyly produce the whole result
        # but it will produce a very very huge matrix
        suppressWarnings(bplapply(seq_len(dim(mt_promoter)[1]), function(index) {
            cor <- suppressWarnings(cor.test(
                mt_promoter[index, ],
                mt_enhancer[index, ]
            ))

            return(c(cor$estimate, cor$p.value))
        })) -> cor_result
    } else {
        suppressWarnings(lapply(seq_len(dim(mt_promoter)[1]), function(index) {
            cor <- suppressWarnings(cor.test(
                mt_promoter[index, ],
                mt_enhancer[index, ]
            ))

            return(c(cor$estimate, cor$p.value))
        })) -> cor_result
    }

    cor_result <- unlist(cor_result)
    cor_mt <- matrix(cor_result, ncol = 2, byrow = TRUE)

    mm_scan$promoter_feature <- enhancer_peak_pair$promoter_feature
    mm_scan$cor <- cor_mt[, 1]
    mm_scan$pvalue <- cor_mt[, 2]
    # we do not need result about feature_id's own correlation
    mm_scan <- subset(mm_scan, feature_id != promoter_feature)

    mm_scan$p_adj <- p.adjust(mm_scan$pvalue)

    # compare with BH, qvalue is more soft
    qvalue_result <- tryCatch(qvalue::qvalue(
        p = mm_scan$pvalue,
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
    mm_scan$qvalue <- qvalues

    metadata(mm_scan)$mm_promoter_tidy <- mm_promoter_tidy

    return(mm_scan)


}
