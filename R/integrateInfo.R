utils::globalVariables(c("padj", "RP_rank", "diff_rank", "rankProduct",
                         "rankOf_rankProduct", "gene_category"))

#' integrate_ChIP_RNA
#'
#' integrate ChIP-Seq and RNA-Seq data to find TF target genes
#'
#' @importFrom stats ks.test
#'
#' @param result_geneRP the simplify result from calcRP_TFHit(report_fullInfo = FALSE)
#' @param result_geneDiff the result from RNA diff result with three column gene_id,
#' log2FoldChange, padj
#' @param lfc_threshold the threshold which decide significant genes
#' @param padj_threshold the threshold which decide significant genes
#'
#' @return a ggplot object if having significant genes in your result. If not, it
#' will report a data.frame with integrated info.
#' @export
#'
#' @examples
#' data("RNADiff_LEC2_GR")
#' library(TxDb.Athaliana.BioMart.plantsmart28)
#' Txdb <- TxDb.Athaliana.BioMart.plantsmart28
#' seqlevels(Txdb) <- c(paste0("Chr", 1:5), "M", "C")
#' peak_path <- system.file("extdata", "ChIP.bed.gz", package = "FindIT2")
#' peak_GR <- loadPeakFile(peak_path)
#' mmAnno <- mm_geneScan(peak_GR,Txdb)
#'
#' calcRP_TFHit(mmAnno = mmAnno,
#'              Txdb = Txdb) -> result_geneRP
#' # output a plot
#' integrate_ChIP_RNA(result_geneRP = result_geneRP,
#'                   result_geneDiff = RNADiff_LEC2_GR) -> merge_data
#' # if you want to extract merge target data
#' target_data <- merge_data$data
integrate_ChIP_RNA <- function(result_geneRP,
                               result_geneDiff,
                               lfc_threshold = 1,
                               padj_threshold = 0.05) {
    check_colnames(
        colnames = c("gene_id", "log2FoldChange", "padj"),
        result_geneDiff
    )

    check_colnames(
        colnames = "gene_id",
        result_geneRP
    )

    # add GOSemSim ? so fun(1/3)
    # fun <- function(exp){
    #   function(x) {
    #     x ^ exp
    #   }
    # }
    # fun_sqrt <- fun(1/2)

    # allGenes_N <- as.double(length(union(result_geneRP$gene_id,
    #                                      result_geneDiff$gene_id)))
    # dplyr::full_join(result_geneRP,
    #                  result_geneDiff) %>%
    #   dplyr::mutate(RP_rank = dplyr::case_when(
    #     is.na(RP_rank) ~ allGenes_N,
    #     TRUE ~ RP_rank)) %>%
    #   dplyr::mutate(diff_rank = rank(padj, na.last = "keep"),
    #                 diff_rank = dplyr::case_when(
    #                   is.na(diff_rank) ~ allGenes_N,
    #                   TRUE ~ diff_rank),
    #                 rankProduct = sqrt(RP_rank * diff_rank),
    #                 rankOf_rankProduct = rank(rankProduct)
    #                 ) -> merge_result



    dplyr::left_join(result_geneRP,
                     result_geneDiff,
                     by = "gene_id"
    ) -> merge_result
    allGenes_N <- as.double(nrow(merge_result))

    merge_result %>%
        dplyr::mutate(
            diff_rank = rank(padj, na.last = "keep"),
            diff_rank = dplyr::case_when(
                is.na(diff_rank) ~ allGenes_N,
                TRUE ~ diff_rank
            ),
            rankProduct = RP_rank * diff_rank,
            rankOf_rankProduct = rank(rankProduct)
        ) %>%
        dplyr::arrange(rankOf_rankProduct) %>%
        dplyr::mutate(
            gene_category = dplyr::case_when(
                log2FoldChange > lfc_threshold & padj < padj_threshold ~ "up",
                log2FoldChange < -lfc_threshold & padj < padj_threshold ~ "down",
                TRUE ~ "static"
            ),
            gene_category = factor(gene_category,
                                   levels = c("up", "down", "static")
            )
        ) -> merge_result


    # in case of no up or down genes
    upGenes_rank <- filter(merge_result, gene_category == "up")$RP_rank
    downGenes_rank <- filter(merge_result, gene_category == "down")$RP_rank
    staticGenes_rank <- filter(merge_result, gene_category == "static")$RP_rank

    if (length(upGenes_rank) == 0 & length(downGenes_rank) == 0) {
        warning("no significant genes, just returing rank product result",
                call. = FALSE
        )
        return(merge_result)
    } else if (length(upGenes_rank) == 0) {
        warning("no significant up genes, just returing rank product result",
                call. = FALSE
        )
        return(merge_result)
    } else if (length(downGenes_rank) == 0) {
        warning("no significant down genes, just returing rank product result",
                call. = FALSE
        )
        return(merge_result)
    }

    suppressWarnings(ks.test(upGenes_rank,
                             staticGenes_rank,
                             alternative = "greater"
    )$p.value) -> up_static_pvalue

    suppressWarnings(ks.test(downGenes_rank,
                             staticGenes_rank,
                             alternative = "greater"
    )$p.value) -> down_static_pvalue


    paste0(
        "Kolmogorov-Smirnov Tests ",
        "\npvalue of up vs static: ",
        format(up_static_pvalue, digits = 3, scientific = TRUE),
        "\npvalue of down vs static: ",
        format(down_static_pvalue, digits = 3, scientific = TRUE)
    ) -> ks_test



    annotate_df <- data.frame(
        xpos = -Inf,
        ypos = Inf,
        annotateText = ks_test,
        hjustvar = 0,
        vjustvar = 1
    )

    merge_result %>%
        ggplot2::ggplot(aes(x = RP_rank)) +
        ggplot2::stat_ecdf(aes(color = gene_category), geom = "line") +
        ggplot2::geom_text(data = annotate_df, aes(
            x = xpos,
            y = ypos,
            hjust = hjustvar,
            vjust = vjustvar,
            label = annotateText
        )) +
        ggplot2::scale_x_continuous(expand = c(0, 0)) -> p

    return(p)
}






#' integrate_replicates
#'
#' integrate value from replicates
#'
#' @param mt value matrix
#' @param colData a data.frame with a single column. Rows of colData
#' correspond to columns of mt.
#' @param fun the function you want to use. If set NULL, program will decide integrate
#' method according to your 'type' parameter.
#' @param type one of 'value', 'rank', 'pvalue'. value will use mean to integrate replicates,
#' rank will use product, and pvalue will use CCT(Cauchy distribution)
#'
#'
#' @return matrix
#' @export
#'
#' @examples
#' mt <- matrix(runif(100, 0, 1), nrow = 10)
#' colnames(mt) <- paste0(paste0("type", 1:5), "_", rep(1:2, 5))
#' rownames(mt) <- paste0("TF", 1:10)
#'
#' colData <- data.frame(
#' type = gsub("_[0-9]","", colnames(mt)),
#' row.names = colnames(mt))
#'
#'integrate_replicates(mt, colData, type = "value")
#'
integrate_replicates <- function(mt,
                                 colData,
                                 fun = NULL,
                                 type = "value") {
    if (!all(colnames(mt) == rownames(colData))) {
        stop("Please make sure order between colnames(mt) and rownames(colData) are consistent",
             call. = FALSE
        )
    } else if (colnames(colData) > 1 & colnames(colData) != "type") {
        stop("only support one column colData and its column name must be 'type'",
             call. = FALSE
        )
    }

    if (is.null(fun)) {
        fun <- switch(type,
                      value = mean,
                      rank = prod,
                      pvalue = CCT,
                      stop("Invalid type value")
        )
    } else {
        fun <- fun
    }

    replicates <- colData$type
    sample <- unique(colData$type)

    sapply(sample, function(x) {
        apply(mt[, x == replicates, drop = FALSE], 1, fun)
    }) -> result

    if (length(result) == 1){
        return(result)
    }

    if (type == "rank") {
        result <- apply(result, 2, function(x) rank(-x))
        return(result)
    } else {
        return(result)
    }
}

#' @importFrom stats pcauchy
CCT <- function(pvals, weights = NULL) {
    # This code is from
    # https://github.com/xihaoli/STAAR/blob/dc4f7e509f4fa2fb8594de48662bbd06a163108c/R/CCT.R

    # This code is used to combinde replicates or different source result
    # The idea is from
    # Qin, Q. et al. (2020). Lisa: inferring transcriptional regulators through integrative modeling of public chromatin accessibility and ChIP-seq data. Genome Biol 21, 32.
    # Combined statistics method for TR ranking


    #### check if there is NA
    if (sum(is.na(pvals)) > 0) {
        stop("Cannot have NAs in the p-values!", call. = FALSE)
    }

    #### check if all p-values are between 0 and 1
    if ((sum(pvals < 0) + sum(pvals > 1)) > 0) {
        stop("All p-values must be between 0 and 1!", call. = FALSE)
    }

    #### check if there are p-values that are either exactly 0 or 1.
    is.zero <- (sum(pvals == 0) >= 1)
    is.one <- (sum(pvals == 1) >= 1)
    if (is.zero && is.one) {
        stop("Cannot have both 0 and 1 p-values!", call. = FALSE)
    }
    if (is.zero) {
        return(0)
    }
    if (is.one) {
        warning("There are p-values that are exactly 1!",
                call. = FALSE
        )
        return(1)
    }

    #### check the validity of weights (default: equal weights) and standardize them.
    if (is.null(weights)) {
        weights <- rep(1 / length(pvals), length(pvals))
    } else if (length(weights) != length(pvals)) {
        stop("The length of weights should be the same as that of the p-values!", call. = FALSE)
    } else if (sum(weights < 0) > 0) {
        stop("All the weights must be positive!", call. = FALSE)
    } else {
        weights <- weights / sum(weights)
    }

    #### check if there are very small non-zero p-values
    is.small <- (pvals < 1e-16)
    if (sum(is.small) == 0) {
        cct.stat <- sum(weights * tan((0.5 - pvals) * pi))
    } else {
        cct.stat <- sum((weights[is.small] / pvals[is.small]) / pi)
        cct.stat <- cct.stat + sum(weights[!is.small] * tan((0.5 - pvals[!is.small]) * pi))
    }

    #### check if the test statistic is very large.
    if (cct.stat > 1e+15) {
        pval <- (1 / cct.stat) / pi
    } else {
        pval <- 1 - pcauchy(cct.stat)
    }
    return(pval)
}
