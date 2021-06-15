utils::globalVariables(c("hit"))

#' jaccard_findIT_TTpair
#'
#' @param input_genes genes which you want to find influentail TF for
#' (same as your find_IT_TTpair parameter)
#' @param TF_target_database TF_target pair data
#' @param input_TF_id TF_id which you want to calculate jaccard index for
#'
#' @return jaccard similarity matrix
#' @export
#'
#' @examples
jaccard_findIT_TTpair <- function(input_genes,
                                  TF_target_database,
                                  input_TF_id){

    input_genes <- as.character(unique(input_genes))
    input_TF_id <- as.character(unique(input_TF_id))

    tibble::tibble(TF_id = rep(input_TF_id, each = length(input_genes)),
                   target_gene = rep(input_genes, length(input_TF_id))) -> fill_df

    TF_target_database %>%
        dplyr::filter(TF_id %in% input_TF_id,
                      target_gene %in% input_genes) %>%
        dplyr::select(TF_id, target_gene) %>%
        dplyr::mutate(hit = 1) -> TF_target_select

    calc_jaccard(fill_df = fill_df,
                 TF_target_select = TF_target_select) -> jaccard_similarity

    return(jaccard_similarity)

}


#' jaccard_findIT_enrichInAll
#'
#' @param input_feature_id peaks which you want to find influentail TF for
#' (same as your find_IT_enrichInAll parameter)
#' @param peak_GR all peak sets.Your input_feature_id is a part of it.
#' @param TF_GR_database TF peak GRange with a column named TF_id representing you TF name
#' @param input_TF_id TF_id which you want to calculate jaccard index for
#'
#' @return
#' @export
#'
#' @examples
jaccard_findIT_enrichInAll <- function(input_feature_id,
                                       peak_GR,
                                       TF_GR_database,
                                       input_TF_id) {

    input_feature_id <- as.character(unique(input_feature_id))
    input_TF_id <- as.character(unique(input_TF_id))

    tibble::tibble(TF_id = rep(input_TF_id,
                               each = length(input_feature_id)),

                   feature_id = rep(input_feature_id,
                                     length(input_TF_id))) -> fill_df

    hits <- GenomicRanges::findOverlaps(peak_GR, TF_GR_database)
    peakTF_pair <- data.frame(
        TF_id = TF_GR_database$TF_id[subjectHits(hits)],
        feature_id = peak_GR$feature_id[queryHits(hits)]
    ) %>%
        dplyr::distinct(feature_id, TF_id, .keep_all = TRUE)

    peakTF_pair %>%
        dplyr::filter(TF_id %in% input_TF_id,
                      feature_id %in% input_feature_id) %>%
        dplyr::mutate(hit = 1)-> TF_target_select

    calc_jaccard(fill_df = fill_df,
                 TF_target_select = TF_target_select) -> jaccard_similarity


    return(jaccard_similarity)




}


#' @importFrom magrittr %>%
#' @importFrom stats dist
calc_jaccard <- function(fill_df,
                         TF_target_select){
    fill_df %>%
        dplyr::left_join(TF_target_select) %>%
        replace(is.na(.), 0) %>%
        tidyr::pivot_wider(names_from = TF_id,
                           values_from = hit) -> TF_hit_wide

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
