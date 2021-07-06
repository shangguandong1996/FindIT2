utils::globalVariables(c("N", "freqNumber"))
utils::globalVariables(c(
    "abs_dist", "distanceToTSS", "xpos", "ypos",
    "hjustvar", "vjustvar", "annotateText"
))
utils::globalVariables(c(
    "N", "freq", "xpos", "annotateText",
    "ypos", "xend", "yend"
))
#' getAssocPairNumber
#'
#' get associated peak number of gene and vice verse.
#'
#' @import rlang
#' @importFrom magrittr %>%
#'
#' @param mmAnno the annotated GRange object from mm_geneScan or mm_nearestGene
#' @param output_type one of 'feature_id' or 'gene_id'
#' @param output_summary whether you want to detailed info
#'
#' @return data.frame
#' @export
#'
#' @examples
#' if (require(TxDb.Athaliana.BioMart.plantsmart28)) {
#'     Txdb <- TxDb.Athaliana.BioMart.plantsmart28
#'     seqlevels(Txdb) <- paste0("Chr", c(1:5, "M", "C"))
#'
#'     peak_path <- system.file("extdata", "ChIP.bed.gz", package = "FindIT2")
#'     peak_GR <- loadPeakFile(peak_path)
#'     peakAnno <- mm_nearestGene(peak_GR, Txdb)
#'
#'     getAssocPairNumber(peakAnno)
#'
#' }
getAssocPairNumber <- function(mmAnno,
                               output_type = "gene_id",
                               output_summary = FALSE) {

    type <- match.arg(output_type, c("feature_id", "gene_id"))

    if (type == "gene_id") {
        colName <- "peakNumber"
    } else if (type == "feature_id") {
        colName <- "geneNumber"
    }

    # In R 3.6 GenomicGRanges, if duplicted names in GRanges
    # as.data.frame will report error(R 4.10 not)
    # So I set mmAnno into NULL to avoid this error
    names(mmAnno) <- NULL
    mmAnno %>%
        as.data.frame(stringsAsFactors = FALSE) %>%
        dplyr::select(feature_id, gene_id) -> mmAnno_df

    if ("mmCor_mode" %in% names(metadata(mmAnno))){
        if (metadata(mmAnno)$mmCor_mode == "enhancerPromoterCor") {
            mm_promoter_tidy <- metadata(mmAnno)$mm_promoter_tidy %>%
                dplyr::select(promoter_feature, gene_id) %>%
                dplyr::rename(feature_id = promoter_feature)

            mmAnno_df <- rbind(mmAnno_df,
                               mm_promoter_tidy) %>%
                dplyr::distinct()

        }
    }

    mmAnno_df %>%
        dplyr::group_by(!!sym(type)) %>%
        dplyr::summarise(!!sym(colName) := dplyr::n()) %>%
        dplyr::mutate(!!sym(type) := stringr::str_sort(!!sym(type),
            numeric = TRUE
        )) -> pairNumber


    freqName <- gsub("_id", "_freq", type)

    pairNumber %>%
        dplyr::group_by(!!sym(colName)) %>%
        dplyr::summarise(freqNumber = dplyr::n()) %>%
        dplyr::mutate(N = dplyr::case_when(
            !!sym(colName) >= 8 ~ ">=8",
            TRUE ~ as.character(!!sym(colName))
        )) %>%
        dplyr::group_by(N) %>%
        dplyr::summarise(!!sym(freqName) := sum(freqNumber)) %>%
        dplyr::mutate(N = factor(N, levels = c(">=8", 7:1))) %>%
        dplyr::mutate(N = droplevels(N)) -> pairNumber_summary
    # when nearest mode, one peak-one gene
    # which means N only has 1


    if (output_summary) {
        return(pairNumber_summary)
    } else {
        return(pairNumber)
    }
}


#' plot_annoDistance
#'
#' plot the distance distribution of mmAnno from mm_nearestGene, which helps you
#' decide whehter your TF is promoter or enhancer dominant
#'
#' @importFrom ggplot2 aes
#' @importFrom patchwork plot_layout
#'
#' @param mmAnno the annotated GRange object from mm_nearestGene
#' @param quantile the quantile of distanceToTSS you want to show
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples
#' if (require(TxDb.Athaliana.BioMart.plantsmart28)) {
#'     Txdb <- TxDb.Athaliana.BioMart.plantsmart28
#'     seqlevels(Txdb) <- paste0("Chr", c(1:5, "M", "C"))
#'
#'     peak_path <- system.file("extdata", "ChIP.bed.gz", package = "FindIT2")
#'     peak_GR <- loadPeakFile(peak_path)
#'     peakAnno <- mm_nearestGene(peak_GR, Txdb)
#'     plot_annoDistance(peakAnno)
#'
#' }
plot_annoDistance <- function(mmAnno,
                              quantile = c(0.01, 0.99)) {

    if (metadata(mmAnno)$mmAnno_mode != "nearestGene") {
        stop("sorry, it only accept mmAnno from mm_nearestGene",
            call. = FALSE
        )
    }

    # TODO: maybe I can take suggestion from
    # https://stackoverflow.com/questions/58785930/r-find-maximum-of-density-plot
    # p1
    mmAnno$abs_dist <- abs(mmAnno$distanceToTSS) + 1
    summary_value <- summary(mmAnno$abs_dist) - 1
    paste(names(summary_value), ": ",
        round(summary_value, digits = 2), "  ",
        collapse = "\n"
    ) -> annotate_text

    annotate_text <- paste0(
        "\nsummary of abs(distanceToTSS)  \n",
        annotate_text
    )

    ggplot2::ggplot(data.frame(mmAnno), aes(x = abs_dist)) +
        ggplot2::geom_density() +
        ggplot2::scale_x_log10() +
        ggplot2::theme_bw() +
        ggplot2::xlab("abs(distanceToTSS) + 1") +
        ggplot2::annotate("text", x = Inf, y = Inf, label = annotate_text,
                          vjust = 1, hjust = 1) -> p1


    quantile_value <- quantile(mmAnno$distanceToTSS, quantile)
    lim_text <- paste0(
        "\n  The plot is zoomed based on",
        "\n  quantile_", quantile[1], ": ", round(quantile_value[1], digits = 2),
        "\n  quantile_", quantile[2], ": ", round(quantile_value[2], digits = 2)
    )

    summary_value <- summary(mmAnno$distanceToTSS)
    paste(names(summary_value), ": ",
        round(summary_value, digits = 2), "  ",
        collapse = "\n"
    ) -> annotate_text

    annotate_text <- paste0(
        "\nsummary of distanceToTSS  \n",
        annotate_text
    )

    # https://stackoverflow.com/questions/32123288/position-ggplot-text-in-each-corner
    annotate_df <- data.frame(
        xpos = c(-Inf, Inf),
        ypos = c(Inf, Inf),
        annotateText = c(
            lim_text,
            annotate_text
        ),
        hjustvar = c(0, 1),
        vjustvar = c(1, 1)
    )

    ggplot2::ggplot(data.frame(mmAnno), aes(x = distanceToTSS)) +
        ggplot2::geom_density() +
        ggplot2::coord_cartesian(xlim = quantile_value) +
        ggplot2::theme_bw() +
        ggplot2::geom_text(data = annotate_df, aes(
            x = xpos,
            y = ypos,
            hjust = hjustvar,
            vjust = vjustvar,
            label = annotateText
        )) -> p2

    # https://github.com/thomasp85/patchwork/issues/246
    # for the patchwork operator
    p1 / p2
}






#' plot_peakGeneAlias_summary
#'
#' @importFrom magrittr %>%
#' @importFrom ggplot2 aes
#'
#' @param mmAnno the annotated GRange object from mm_geneScan or mm_nearestGene
#' @param mmAnno_corFilter the filter mmAnno object according to p-value or cor, defalut is NULL
#' @param output_type one of 'feature_id' or 'gene_id'
#' @param fillColor the bar plot color
#'
#' @return a ggplot object
#' @export
#'
#' @examples
#' if (require(TxDb.Athaliana.BioMart.plantsmart28)) {
#'     Txdb <- TxDb.Athaliana.BioMart.plantsmart28
#'     seqlevels(Txdb) <- paste0("Chr", c(1:5, "M", "C"))
#'
#'     peak_path <- system.file("extdata", "ChIP.bed.gz", package = "FindIT2")
#'     peak_GR <- loadPeakFile(peak_path)
#'     peakAnno <- mm_nearestGene(peak_GR, Txdb)
#'
#'     plot_peakGeneAlias_summary(peakAnno)
#'
#' }
plot_peakGeneAlias_summary <- function(mmAnno,
                                       mmAnno_corFilter = NULL,
                                       output_type = "gene_id",
                                       fillColor = "#ca6b67") {

    if (output_type == "gene_id") {
        ylabs <- paste0("The freq of each gene has N peak")
    } else if (output_type == "feature_id") {
        ylabs <- paste0("The freq of each peak has N gene")
    }

    getAssocPairNumber(
        mmAnno = mmAnno,
        output_type = output_type,
        output_summary = TRUE
    ) -> summary
    colnames(summary)[2] <- "freq"
    # because I want to plot
    # bidirectional bar chart with positive labels on both sides ggplot2
    # so I have to use aes instead of aes_string, so I change the colnames

    if (!is.null(mmAnno_corFilter)) {
        summary[, 2] <- -summary[, 2]

        getAssocPairNumber(
            mmAnno = mmAnno_corFilter,
            output_type = output_type,
            output_summary = TRUE
        ) -> summary_filter
        colnames(summary_filter)[2] <- "freq"
    }


    if (is.null(mmAnno_corFilter)) {
        ggplot2::ggplot(summary, aes(x = N, y = freq)) +
            ggplot2::geom_bar(
                stat = "identity", color = "black",
                fill = fillColor
            ) +
            ggplot2::geom_text(aes(label = abs(freq), y = freq)) +
            ggplot2::ylab(ylabs) +
            ggplot2::coord_flip() -> p

        return(p)
    } else {
        ggplot2::ggplot(summary, aes(x = N, y = freq)) +
            ggplot2::geom_bar(stat = "identity", color = "black") +
            ggplot2::geom_text(aes(label = abs(freq), y = freq),
                hjust = 1.1
            ) +
            ggplot2::geom_bar(
                data = summary_filter,
                stat = "identity", fill = fillColor, color = "black"
            ) +
            ggplot2::geom_text(aes(label = freq, y = freq),
                data = summary_filter, color = "black",
                hjust = -0.1
            ) +
            ggplot2::scale_y_continuous(
                labels = abs,
                limits = c(
                    min(summary$freq, -summary_filter$freq) * 1.1,
                    max(-summary$freq, summary_filter$freq) * 1.1
                )
            ) +
            ggplot2::ylab(ylabs) -> p1


        annotate_df <- data.frame(
            xpos = c(levels(summary$N)[1], levels(summary$N)[1]),
            ypos = c(max(abs(summary$freq)), -max(abs(summary$freq))),
            xend = c(levels(summary$N)[3], levels(summary$N)[3]),
            yend = c(max(abs(summary$freq) / 1.2), -max(abs(summary$freq) / 1.2)),
            annotateText = c(
                "The result of filtered cor",
                "The result of origin cor"
            ),
            curve = c(0.3, -0.3)
        )

        # deal with nearest mode + feature_id
        annotate_df[is.na(annotate_df)] <- 1


        p1 +
            ggplot2::geom_text(data = annotate_df, aes(
                x = xpos,
                y = yend,
                label = annotateText
            ), ) +
            ggplot2::geom_curve(
                data = annotate_df[1, ], aes(x = xpos, y = ypos, xend = xend, yend = yend),
                arrow = ggplot2::arrow(length = ggplot2::unit(0.07, "inch")), size = 0.5,
                color = "gray20",
                curvature = annotate_df[1, 6]
            ) +
            ggplot2::geom_curve(
                data = annotate_df[2, ], aes(x = xpos, y = ypos, xend = xend, yend = yend),
                arrow = ggplot2::arrow(length = ggplot2::unit(0.07, "inch")), size = 0.5,
                color = "gray20",
                curvature = annotate_df[2, 6]
            ) +
            ggplot2::coord_flip() +
            ggplot2::theme_bw() -> p

        return(p)
    }
}
