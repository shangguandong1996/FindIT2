utils::globalVariables(c("hit"))

#' jaccard_findIT_TTpair
#'
#' @param input_feature_id a character vector which represent peaks set
#' which you want to find influential TF for (same as your find_IT_TTpair parameter)
#' @param TF_target_database TF_target pair data
#' @param input_TF_id TF_id which you want to calculate jaccard index for
#'
#' @return jaccard similarity matrix
#' @export
#'
#' @examples
#'
#' data("TF_target_database")
#' data("test_geneSet")
#' findIT_TTPair(
#'     input_genes = test_geneSet,
#'     TF_target_database = TF_target_database
#' ) -> result_findIT_TTPair
#'
#' jaccard_findIT_TTpair(
#'     input_genes = test_geneSet,
#'     TF_target_database = TF_target_database,
#'     input_TF_id = result_findIT_TTPair$TF_id[1:3]
#' )
jaccard_findIT_TTpair <- function(input_genes,
                                  TF_target_database,
                                  input_TF_id) {

    input_genes <- as.character(unique(input_genes))
    input_TF_id <- as.character(unique(input_TF_id))

    fill_df <- tibble::tibble(
        TF_id = rep(input_TF_id, each = length(input_genes)),
        target_gene = rep(input_genes, length(input_TF_id))
    )

    TF_target_select <- TF_target_database %>%
        dplyr::filter(
            TF_id %in% input_TF_id,
            target_gene %in% input_genes
        ) %>%
        dplyr::select(TF_id, target_gene) %>%
        dplyr::mutate(hit = 1)

    jaccard_similarity <- calc_jaccard(
        fill_df = fill_df,
        TF_target_select = TF_target_select
    )

    return(jaccard_similarity)
}


#' jaccard_findIT_enrichInAll
#'
#' @param input_feature_id a character vector which represent peaks set
#' which you want to find influential TF for (same as your find_IT_enrichInAll parameter)
#' @param peak_GR a GRange object represent your whole feature location with a
#' column named feature_id, which your input_feature_id should a part of it.
#' @param TF_GR_database TF peak GRange with a column named TF_id representing you TF name
#' @param input_TF_id TF_id which you want to calculate jaccard index for
#'
#' @return jaccard similarity matrix
#' @export
#'
#' @examples
#' data("test_featureSet")
#' peak_path <- system.file("extdata", "ATAC.bed.gz", package = "FindIT2")
#' peak_GR <- loadPeakFile(peak_path)
#'
#' ChIP_peak_path <- system.file("extdata", "ChIP.bed.gz", package = "FindIT2")
#' ChIP_peak_GR <- loadPeakFile(ChIP_peak_path)
#' ChIP_peak_GR$TF_id <- "AT1G28300"
#' findIT_enrichInAll(
#'     input_feature_id = test_featureSet,
#'     peak_GR = peak_GR,
#'     TF_GR_database = ChIP_peak_GR
#' ) -> result_findIT_enrichInAll
#'
#' jaccard_findIT_enrichInAll(
#'     input_feature_id = test_featureSet,
#'     peak_GR = peak_GR,
#'     TF_GR_database = ChIP_peak_GR,
#'     input_TF_id = result_findIT_enrichInAll$TF_id[1]
#' )
jaccard_findIT_enrichInAll <- function(input_feature_id,
                                       peak_GR,
                                       TF_GR_database,
                                       input_TF_id) {

    input_feature_id <- as.character(unique(input_feature_id))
    input_TF_id <- as.character(unique(input_TF_id))

    fill_df <- tibble::tibble(
        TF_id = rep(input_TF_id,
            each = length(input_feature_id)
        ),
        feature_id = rep(
            input_feature_id,
            length(input_TF_id)
        )
    )

    hits <- GenomicRanges::findOverlaps(peak_GR, TF_GR_database)
    peakTF_pair <- data.frame(
        TF_id = TF_GR_database$TF_id[subjectHits(hits)],
        feature_id = peak_GR$feature_id[queryHits(hits)]
    ) %>%
        dplyr::distinct(feature_id, TF_id, .keep_all = TRUE)

    TF_target_select <- peakTF_pair %>%
        dplyr::filter(
            TF_id %in% input_TF_id,
            feature_id %in% input_feature_id
        ) %>%
        dplyr::mutate(hit = 1)

    jaccard_similarity <- calc_jaccard(
        fill_df = fill_df,
        TF_target_select = TF_target_select
    )


    return(jaccard_similarity)
}


#' @importFrom stats dist
calc_jaccard <- function(fill_df,
                         TF_target_select) {

    suppressMessages(TF_hit_wide <- fill_df %>%
                         dplyr::left_join(TF_target_select) %>%
                         replace(is.na(.), 0) %>%
                         tidyr::pivot_wider(
                             names_from = TF_id,
                             values_from = hit
                         ))

    TF_hit_mt <- as.matrix(TF_hit_wide[, -1])
    rownames(TF_hit_mt) <- TF_hit_wide[[1]]

    jaccard_dist <- dist(t(TF_hit_mt), method = "binary")
    jaccard_dist_matrix <- as.matrix(jaccard_dist)

    jaccard_similarity <- 1 - jaccard_dist_matrix

    # the diag of similarity matrix is equal to 0
    # which will be helpful for subsequent heatmap plot
    diag(jaccard_similarity) <- 0

    return(jaccard_similarity)
}
