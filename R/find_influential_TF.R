utils::globalVariables(c("gene_type"))
utils::globalVariables(c("TF_hit"))
utils::globalVariables(c("TF", "target_gene", "num_TFHit_bg",
                         "num_TFHit_inputGene", "num_topLeft", "num_bottomLeft",
                         "num_topRight", "num_TFHit_input", "inputRatio",
                         "bgRatio", "odds_ratio"))
utils::globalVariables(c("TF_id", "num_TFHit_bg", "num_topLeft", "num_bottomLeft",
                         "num_topRight", "num_TFHit_inputFeature"))
utils::globalVariables(c("hit_N"))


#' findI(nfluential)T(F)_regionRP
#'
#' find Influential TF of your input gene set based on regulatory potential data and TF ChIP-Seq or
#' motif data
#'
#' @importFrom magrittr %>%
#' @importFrom stats t.test
#'
#' @param regionRP the MultiAssayExperiment object from calcRP_region
#' @param Txdb Txdb
#' @param TF_GR_database TF peak GRange with a column named TF_id representing you TF name
#' @param input_genes genes which you want to find influentail TF for
#' @param background_genes If you do not assign background gene , program will sample
#' background_number genes as background genes from all gene sets.
#' @param background_number background genes number
#'
#'
#' @return a MultiAssayExperiment object containg detailed TF-percent and TF-pvalue
#' @export
#'
#' @examples
#' #see in vignettes
#'
findIT_regionRP <- function(regionRP,
                            Txdb,
                            TF_GR_database,
                            input_genes,
                            background_genes = NULL,
                            background_number = 3000) {

    check_colnames("TF_id", TF_GR_database)

    Chrs_included <- S4Vectors::metadata(regionRP)$Chrs_included
    input_genes <- as.character(unique(input_genes))

    # warning will appear firstly
    withr::local_options(list(warn = 1))
    all_geneSets <- GenomicFeatures::genes(Txdb)
    all_geneSets <- names(all_geneSets[GenomeInfoDb::seqnames(all_geneSets) %in% Chrs_included])

    if (mean(input_genes %in% all_geneSets) < 1) {
        dropN <- sum(!input_genes %in% all_geneSets)
        msg <- paste0(dropN, " input genes are not in your Txdb, I will filter these genes")
        warning(msg, call. = FALSE)
        input_genes <- input_genes[input_genes %in% all_geneSets]
    }

    if (is.null(background_genes)) {
        background_pools <- all_geneSets[!all_geneSets %in% input_genes]
        # # This seed is anniversary of me and Daisy :)
        # set.seed(seed)
        background_genes <- sample(background_pools,
                                   size = background_number
        )
    } else if (mean(background_genes %in% all_geneSets) < 1) {
        dropN <- sum(!background_genes %in% all_geneSets)
        msg <- paste0(dropN,
                      " background genes are not in your Txdb,
                      I will filter these genes")
        warning(msg,call. = FALSE)
    }

    if (any(input_genes %in% background_genes)) {
        warning("your input and background have overlapping genes", call. = FALSE)
    }

    fullRP_mt <- SummarizedExperiment::assays(regionRP)$fullRP
    fullRP_GR <- SummarizedExperiment::rowRanges(MultiAssayExperiment::experiments(regionRP)[[1]])
    # pre-fill gene without peak with percent 0
    genes_twoSets <- c(input_genes, background_genes)
    data.frame(
        gene_id = genes_twoSets[!genes_twoSets %in% unique(fullRP_GR$gene_id)],
        percent = 0,
        stringsAsFactors = FALSE
    ) -> TF_percent_fill

    data.frame(fullRP_GR, stringsAsFactors = FALSE) %>%
        dplyr::select(c("seqnames", "start", "end", "width", "strand", "feature_id")) %>%
        dplyr::distinct(.keep_all = TRUE) %>%
        GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE) -> peakGR_origin

    hits <- GenomicRanges::findOverlaps(TF_GR_database, peakGR_origin)

    # although some motif may hits more than one times in one peak
    # matain these info is not useful for find IT based on percent

    # these info can be useful for ridge analysis, but not here

    data.frame(
        TF_id = TF_GR_database$TF_id[queryHits(hits)],
        feature_id = peakGR_origin$feature_id[subjectHits(hits)],
        stringsAsFactors = FALSE) %>%
        # dplyr::group_by(TF_id, feature_id) %>%
        # dplyr::summarise(hit_N = dplyr::n()) %>%
        # dplyr::mutate(hitN_norm = hit_N - sum(hit_N) / length(peakGR_origin)) %>%
        # dplyr::ungroup() -> hits_df
        dplyr::distinct(TF_id, feature_id, .keep_all = TRUE) %>%
        dplyr::as_tibble() -> hits_df

    # select related feature_id to speed up below calculation
    input_related_feature <- subset(fullRP_GR, gene_id %in% genes_twoSets)$feature_id
    hits_df <- subset(hits_df, feature_id %in% input_related_feature)

    # select related fullRP_GR to speed up below calculation
    fullRP_GR_select <- subset(fullRP_GR, gene_id %in% genes_twoSets)
    mmAnno_peakGenePair <- S4Vectors::mcols(fullRP_GR_select)[, c("feature_id", "gene_id")]
    fullRP_mt_select <- fullRP_mt[names(fullRP_GR_select), ]

    all_TFs <- unique(TF_GR_database$TF_id)

    allSample_pvalue_list <- list()
    allSample_mean_list <- list()
    allSample_percent_list <- list()

    for (i in seq_len(dim(fullRP_mt)[2])) {
        sample_name <- colnames(fullRP_mt_select)[i]

        fullRP_mt_sample <- data.frame(mmAnno_peakGenePair,
                                       RP = fullRP_mt_select[, sample_name],
                                       stringsAsFactors = FALSE
        )

        cat(
            ">> dealing with", sample_name, "...\t\t",
            format(Sys.time(), "%Y-%m-%d %X"), "\n"
        )

        pvalue_list <- list()
        mean_list <- list()
        percent_list <- list()
        pb <- progress::progress_bar$new(total = length(all_TFs))

        for (TF in all_TFs) {
            pb$tick()
            suppressMessages(fullRP_mt_sample %>%
                                 dplyr::left_join(subset(hits_df, TF_id == TF)) %>%
                                 dplyr::mutate(TF_hit = dplyr::case_when(
                                     is.na(TF_id) ~ 0,
                                     TRUE ~ 1
                                 )) %>%
                                 dplyr::group_by(gene_id) %>%
                                 dplyr::summarise(percent = sum(TF_hit * RP) / sum(RP)) %>%
                                 dplyr::bind_rows(TF_percent_fill) -> TF_percent)

            TF_percent_input <- subset(TF_percent, gene_id %in% input_genes)
            # TF_percent_input$gene_type <- "input_genes"
            TF_percent_background <- subset(TF_percent, gene_id %in% background_genes)
            # TF_percent_background$gene_type <- "background_genes"

            # percent_list[[TF]] <- dplyr::bind_rows(
            #     TF_percent_input,
            #     TF_percent_background
            # )

            percent_list[[TF]] <- TF_percent_input

            input_percent <- TF_percent_input$percent
            background_percent <- TF_percent_background$percent

            mean_list[[TF]] <- mean(input_percent) - mean(background_percent)

            # maybe we should not use wilcox.test
            # https://github.com/AllenWLynch/lisa/issues/3
            # http://archive.bio.ed.ac.uk/jdeacon/statistics/tress4.html#Transformation%20of%20data
            # https://stackoverflow.com/questions/55984490/arcsine-transformation-of-percentage-data
            # pvalue_list[[TF]] <- wilcox.test(input_percent,
            #                                  background_percent,
            #                                  alternative = "greater")$p.value

            # one side ? or two side ?
            pvalue_list[[TF]] <- t.test(
                asin(sqrt(input_percent)),
                asin(sqrt(background_percent))
            )$p.value
        }

        lapply(names(percent_list), function(x) {
            data <- percent_list[[x]]
            data$TF_id <- x
            return(data)
        }) -> percent_df

        do.call(rbind, percent_df) -> allSample_percent_list[[sample_name]]

        allSample_mean_list[[sample_name]] <- unlist(mean_list)
        allSample_pvalue_list[[sample_name]] <- unlist(pvalue_list)
    }

    allSample_pvalue_se <- SummarizedExperiment::SummarizedExperiment(do.call(cbind, allSample_pvalue_list))
    allSample_mean_se <- SummarizedExperiment::SummarizedExperiment(do.call(cbind, allSample_mean_list))

    lapply(names(allSample_percent_list), function(x) {
        data <- allSample_percent_list[[x]]
        data$sample <- x
        return(data)
    }) -> allSample_percent_df
    do.call(rbind, allSample_percent_df) -> allSample_percent_df

    final_result <- MultiAssayExperiment::MultiAssayExperiment(
        list(
            "TF_percentMean" = allSample_mean_se,
            "TF_pvalue" = allSample_pvalue_se
        )
    )

    S4Vectors::metadata(final_result)$percent_df <- allSample_percent_df
    S4Vectors::metadata(final_result)$hits_df <- hits_df


    cat(
        ">> done\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n\n"
    )
    return(final_result)
}



#' findI(nfluential)T(F)_T(F)T(arget)Pair
#'
#' find influential TF of your input gene set based on public TF-Target data
#'
#' @importFrom magrittr %>%
#' @importFrom stats fisher.test
#'
#' @param input_genes genes which you want to find influentail TF for
#' @param TF_target_database TF_target pair data
#' @param gene_background If you do not assign background gene, program will consider
#' all target gene as background
#' @param TFHit_min minimal size of target gene regulated by TF
#' @param TFHit_max maximal size of target gene regulated by TF
#'
#' @return data.frame
#' @export
#'
#' @examples
#' #see in vignettes
#'
findIT_TTPair <- function(input_genes,
                          TF_target_database,
                          gene_background = NULL,
                          TFHit_min = 5,
                          TFHit_max = 10000) {

    # maybe we can include TF-pair score ?

    # input_genes and TF-target pair may have some overlap
    # some human gene maybe look like number, like 12340?
    input_genes <- as.character(unique(input_genes))

    if (is.null(gene_background)) {
        gene_background <- unique(TF_target_database$target_gene)
    }

    # drop genes may not in your TF database
    # just like GO enrichment
    input_genes <- input_genes[input_genes %in% gene_background]

    num_gene_bg <- length(gene_background)
    num_gene_input <- length(input_genes)

    # some human gene maybe look like number, like 12340?
    TF_target_database %>%
        dplyr::distinct(TF_id, target_gene,
                        .keep_all = TRUE
        ) %>%
        dplyr::mutate(
            TF_id = as.character(TF_id),
            target_gene = as.character(target_gene)
        ) %>%
        dplyr::group_by(TF_id) %>%
        dplyr::mutate(num_TFHit_bg = dplyr::n()) %>%
        dplyr::filter(
            num_TFHit_bg > TFHit_min,
            num_TFHit_bg < TFHit_max
        ) %>%
        # this step can drop TF whose target is not in inputGenes
        dplyr::filter(target_gene %in% input_genes) -> TF_target_filter

    if (nrow(TF_target_filter) == 0) {
        stop("sorry, your input genes is not in your target genes",
             call. = FALSE
        )
    }

    TF_target_filter %>%
        dplyr::mutate(num_TFHit_inputGene = dplyr::n()) %>%
        dplyr::ungroup() %>%
        dplyr::select(TF_id, num_TFHit_bg, num_TFHit_inputGene) %>%
        dplyr::distinct(.keep_all = TRUE) %>%
        dplyr::rename(num_topLeft = num_TFHit_inputGene) %>%
        # N of TF Hit in inputGene
        dplyr::mutate(
            num_bottomLeft = num_TFHit_bg - num_topLeft, # N of TF Hit in leftGene
            num_topRight = num_gene_input - num_topLeft, # N of TF not Hit in inputGene
            # N of TF not Hit in leftGene
            num_bottomRight = num_gene_bg - num_topLeft - num_bottomLeft - num_topRight
        ) -> TFHit_summary


    # # https://stackoverflow.com/questions/40093595/dplyr-group-by-and-convert-multiple-columns-to-a-vector
    # This TF_split list will be large object, for 1.7M TF-target, it will be 122.4MB
    # TF_target_database %>%
    # split(.$TF) %>%
    #   lapply(function(x){x$target_gene}) -> TF_split



    TF_fihserMt <- TFHit_summary[, -c(1:2)]
    tibble::tibble(
        TF_id = TFHit_summary$TF_id,
        pvalue = apply(TF_fihserMt, 1, function(x) fisher.test(matrix(x, nrow = 2),
                                                               alternative = "greater")$p.value),
        odds_ratio = apply(TF_fihserMt, 1, function(x) fisher.test(matrix(x, nrow = 2),
                                                                   alternative = "greater")$estimate)
    ) %>%
        dplyr::mutate(pvalue = dplyr::case_when(
            pvalue > 1 ~ 1, # https://github.com/StoreyLab/qvalue/issues/7
            TRUE ~ pvalue
        )) -> TF_fisher_result

    # If qvalue error, maybe
    # https://github.com/YuLab-SMU/clusterProfiler/blob/97ddd55fb4aa131956e3323fda9e59cb43a6a8ec/R/enrichDAVID.R#L131-L132

    dplyr::inner_join(TFHit_summary, TF_fisher_result) %>%
        dplyr::mutate(
            inputRatio = paste0(num_topLeft, "/", num_gene_input),
            bgRatio = paste0(num_TFHit_bg, "/", num_gene_bg),
            num_TFHit_input = num_topLeft
        ) %>%
        dplyr::select(TF_id, num_TFHit_input, inputRatio, bgRatio, pvalue, odds_ratio) %>%
        dplyr::mutate(p.adj = p.adjust(pvalue)) -> final_result

    qvalue_result <- tryCatch(qvalue::qvalue(
        p = final_result$pvalue,
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

    final_result %>%
        dplyr::mutate(
            qvalue = qvalues,
            rank = rank(pvalue)
        ) %>%
        dplyr::arrange(pvalue) -> final_result

    return(final_result)
}


#' findI(nfluential)T(F)_TFHit
#'
#' find influential TF of your input gene set based on TF ChIP-Seq or motif data
#'
#' @importFrom stats wilcox.test
#'
#' @param input_genes genes which you want to find influentail TF for
#' @param Txdb Txdb
#' @param TF_GR_database TF peak GRange with a column named TF_id representing you TF name
#' @param scan_dist scan distance
#' @param decay_dist decay distance
#' @param Chrs_included chromsomes which you want to sample background genes from
#' @param background_genes If you do not assign background gene , program will sample
#' background_number genes as background genes from all gene sets.
#' @param background_number background genes number
#'
#' @return data.frame
#' @export
#'
#' @examples
#' #see in vignettes
#'
findIT_TFHit <- function(input_genes,
                         Txdb,
                         TF_GR_database,
                         scan_dist = 1e3,
                         decay_dist = 2e4,
                         Chrs_included,
                         background_genes = NULL,
                         background_number = 3000) {
    input_genes <- as.character(unique(input_genes))

    # TODO: add TF_id report

    if (missing(Chrs_included)) {
        Chrs_included <- seqlevels(Txdb)
    }

    # warning will appear firstly
    withr::local_options(list(warn = 1))
    gene_location <- GenomicFeatures::genes(Txdb)
    all_geneSets <- names(gene_location[GenomeInfoDb::seqnames(gene_location) %in% Chrs_included])

    if (mean(input_genes %in% all_geneSets) < 1) {
        dropN <- sum(!input_genes %in% all_geneSets)
        msg <- paste0(dropN,
                      " input genes are not in your Txdb, I will filter these genes")
        warning(msg, call. = FALSE)
        input_genes <- input_genes[input_genes %in% all_geneSets]
    }

    if (is.null(background_genes)) {
        background_pools <- all_geneSets[!all_geneSets %in% input_genes]
        # This seed is anniversary of me and Daisy :)
        # set.seed(seed)
        background_genes <- sample(background_pools,
                                   size = background_number
        )
    } else if (mean(background_genes %in% all_geneSets) < 1) {
        dropN <- sum(!background_genes %in% all_geneSets)
        msg <- paste0(dropN,
                      " background genes are not in your Txdb,
                      I will filter these genes")
        warning(msg, call. = FALSE)
    }

    if (any(input_genes %in% background_genes)) {
        warning("your input and background have overlapping genes", call. = FALSE)
    }

    genes_twoSets <- c(input_genes, background_genes)
    # 我只要看input gene那里的peak就可以了，不需要所有的peak
    gene_start <- GenomicRanges::resize(gene_location[genes_twoSets],
                                        width = 1, fix = "start"
    )

    gene_promoter <- suppressWarnings(
        GenomicRanges::promoters(gene_location[genes_twoSets],
                                 upstream = scan_dist,
                                 downstream = scan_dist
        )
    )

    if (all(is.na(seqinfo(Txdb)@seqlengths))) {
        warning("your chr length is not set, gene_promoter may cross bound",
                call. = FALSE
        )
    } else {
        cat(
            ">> some scan range may cross Chr bound, trimming...\t\t",
            format(Sys.time(), "%Y-%m-%d %X"), "\n"
        )
        gene_promoter <- GenomicRanges::trim(gene_promoter)
    }

    TF_GR_list <- split(TF_GR_database, TF_GR_database$TF_id)

    purrr::map(names(TF_GR_list), function(x) {
        TF_GR <- TF_GR_list[[x]]
        overlapHits <- GenomicRanges::findOverlaps(gene_promoter, TF_GR)

        hit_genes <- gene_promoter[S4Vectors::queryHits(overlapHits)]$gene_id

        RP_fill <- data.frame(
            gene_id = genes_twoSets[!genes_twoSets %in% unique(hit_genes)],
            sumRP = 0
        )

        hit_features_location <- TF_GR[S4Vectors::subjectHits(overlapHits)]

        GenomicRanges::distance(
            gene_start[hit_genes],
            GenomicRanges::resize(hit_features_location, fix = "center", width = 1)
        ) -> centerToTSS

        RP <- 2^(-centerToTSS / decay_dist)

        data.frame(
            gene_id = hit_genes,
            RP = RP
        ) %>%
            dplyr::group_by(gene_id) %>%
            dplyr::summarise(sumRP = sum(RP)) %>%
            dplyr::bind_rows(RP_fill) %>%
            dplyr::mutate(gene_type = dplyr::case_when(
                gene_id %in% input_genes ~ "input_genes",
                TRUE ~ "bg_genes"
            )) %>%
            dplyr::group_split(gene_type) -> RP_list

        wilcox.test(RP_list[[2]]$sumRP,
                    RP_list[[1]]$sumRP,
                    alternative = "greater"
        )$p.value -> pvalue

        mean(RP_list[[2]]$sumRP) - mean(RP_list[[1]]$sumRP) -> mean_value

        return(c(mean_value, pvalue, length(TF_GR)))
    }) -> result

    final_result <- data.frame(matrix(unlist(result), ncol = 3, byrow = TRUE))
    colnames(final_result) <- c("mean_value", "pvalue", "TFPeak_number")

    final_result$rank <- rank(final_result$pvalue)
    final_result$TF_id <- names(TF_GR_list)
    final_result <- dplyr::arrange(final_result, pvalue)
}




#' findI(nfluential)T(F)_enrichInShuffle
#'
#' find influential TF of your input peak set compared with your whole peak sets
#' based on TF ChIP-Seq or motif data.
#'
#' @import BiocParallel
#' @importFrom stats ecdf
#'
#' @param input_feature_id peaks which you want to find influentail TF for
#' @param peak_GR all peak sets.Your input_feature_id is a part of it.
#' @param TF_GR_database TF peak GRange with a column named TF_id representing you TF name
#' @param shuffleN shuffle number
#'
#' @return data.frame
#' @export
#'
#' @examples
#'
#' #see in vignettes
#'
findIT_enrichInShuffle <- function(input_feature_id,
                                   peak_GR,
                                   TF_GR_database,
                                   shuffleN = 1000) {
    input_feature_id <- unique(input_feature_id)
    all_feature_id <- unique(peak_GR$feature_id)
    names(peak_GR) <- peak_GR$feature_id
    input_GR <- peak_GR[input_feature_id]

    # # This seed is anniversary of me and Daisy :)
    # set.seed(20160806)
    purrr::map(seq_len(shuffleN), function(x) {
        sample(all_feature_id, size = length(input_feature_id))
    }) -> shuffle_feature
    TF_GR_list <- split(TF_GR_database, TF_GR_database$TF_id)

    # 要用jaccard index，而不是单纯的overlap吧？
    # 不用，因为我们使用随机的shuffle来做的
    cat(
        ">> done\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n"
    )
    bplapply(names(TF_GR_list), function(x) {
        TF_GR <- TF_GR_list[[x]]

        # Because we use the shuffle
        # we can consider more than one Tf hit in one peak

        # true_number <- length(unique(
        #   S4Vectors::queryHits(GenomicRanges::findOverlaps(input_GR, TF_GR))
        #   ))

        true_number <- length(queryHits(GenomicRanges::findOverlaps(input_GR, TF_GR)))

        purrr::map_dbl(shuffle_feature, function(x) {
            length(queryHits(GenomicRanges::findOverlaps(peak_GR[x], TF_GR)))
        }) -> background_number

        # https://stats.stackexchange.com/questions/50080/estimate-quantile-of-value-in-a-vector
        pvalue <- 1 - ecdf(c(true_number, background_number))(true_number)
        data.frame(
            TF_id = x,
            pvalue = pvalue,
            hitInput_N = true_number,
            shuffle_meanN = mean(background_number),
            TFPeak_number = length(TF_GR)
        ) -> df
        return(df)
    }) -> pvalue_result

    final_result <- do.call(rbind, pvalue_result)
    final_result$rank <- rank(final_result$pvalue)
    final_result <- dplyr::arrange(final_result, pvalue)

    cat(
        ">> done\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n"
    )

    return(final_result)
}



#' findI(nfluential)T(F)_enrichInAll
#'
#' find influential TF of your input peak set compared with your whole peak sets
#' based on TF ChIP-Seq or motif data.
#'
#' @importFrom stats fisher.test
#'
#' @param input_feature_id peaks which you want to find influentail TF for
#' @param peak_GR all peak sets.Your input_feature_id is a part of it.
#' @param TF_GR_database TF peak GRange with a column named TF_id representing you TF name
#'
#' @return data.frame
#' @export
#'
#' @examples
#' #see in vignettes
#'
findIT_enrichInAll <- function(input_feature_id,
                               peak_GR,
                               TF_GR_database) {
    num_feature_input <- length(input_feature_id)
    num_feature_bg <- length(peak_GR$feature_id)

    hits <- GenomicRanges::findOverlaps(peak_GR, TF_GR_database)
    peakTF_pair <- data.frame(
        feature_id = peak_GR$feature_id[queryHits(hits)],
        TF_id = TF_GR_database$TF_id[subjectHits(hits)]
    ) %>%
        dplyr::distinct(feature_id, TF_id, .keep_all = TRUE)


    TFHit_summary <- peakTF_pair %>%
        dplyr::group_by(TF_id) %>%
        dplyr::mutate(num_TFHit_bg = dplyr::n()) %>%
        dplyr::filter(feature_id %in% input_feature_id) %>%
        dplyr::mutate(num_topLeft = dplyr::n()) %>%
        dplyr::ungroup() %>%
        dplyr::select(TF_id, num_TFHit_bg, num_topLeft) %>%
        dplyr::distinct(.keep_all = TRUE) %>%
        dplyr::mutate(
            num_bottomLeft = num_TFHit_bg - num_topLeft,
            num_topRight = num_feature_input - num_topLeft,
            num_bottomRight = num_feature_bg - num_topLeft - num_bottomLeft - num_topRight
        )


    TF_fihserMt <- as.matrix(TFHit_summary[, -c(1:2)])

    data.frame(
        TF_id = TFHit_summary$TF_id,
        pvalue = apply(TF_fihserMt, 1, function(x) fisher.test(matrix(x, nrow = 2),
                                                               alternative = "greater")$p.value),
        odds_ratio = apply(TF_fihserMt, 1, function(x) fisher.test(matrix(x, nrow = 2),
                                                                   alternative = "greater")$estimate)
    ) %>%
        dplyr::mutate(pvalue = dplyr::case_when(
            pvalue > 1 ~ 1, # https://github.com/StoreyLab/qvalue/issues/7
            TRUE ~ pvalue
        )) -> TF_fisher_result


    final_result <- dplyr::inner_join(TFHit_summary, TF_fisher_result) %>%
        dplyr::mutate(
            inputRatio = paste0(num_topLeft, "/", num_feature_input),
            bgRatio = paste0(num_TFHit_bg, "/", num_feature_bg),
            num_TFHit_inputFeature = num_topLeft
        ) %>%
        dplyr::select(TF_id, num_TFHit_inputFeature, inputRatio, bgRatio, pvalue, odds_ratio) %>%
        dplyr::mutate(p.adj = p.adjust(pvalue))

    qvalue_result <- tryCatch(qvalue::qvalue(
        p = final_result$pvalue,
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

    final_result %>%
        dplyr::mutate(
            qvalue = qvalues,
            rank = rank(pvalue)
        ) %>%
        dplyr::arrange(pvalue) -> final_result

    return(final_result)
}

#' findIT_MARA
#'
#' @importFrom glmnet cv.glmnet glmnet
#' @importFrom stats model.matrix
#' @importFrom stats coef
#'
#' @param input_feature_id peaks which you want to find influentail TF for
#' @param peak_GR all peak sets.Your input_feature_id is a part of it.
#' @param peakScoreMt peak count matrix.
#' @param TF_GR_database TF peak GRange with a column named TF_id representing you TF name
#' @param log whether you want to log your peakScoreMt
#' @param meanScale whether you want to mean-centered per row
#'
#' @return a data.frame
#' @export
#'
#' @examples
#'
#'
findIT_MARA <- function(input_feature_id,
                        peak_GR,
                        peakScoreMt,
                        TF_GR_database,
                        log = TRUE,
                        meanScale = TRUE) {
    input_feature_id <- unique(input_feature_id)

    hits <- GenomicRanges::findOverlaps(TF_GR_database, peak_GR)

    all_TF <- unique(TF_GR_database$TF_id)

    data.frame(
        TF_id = TF_GR_database$TF_id[queryHits(hits)],
        feature_id = peak_GR$feature_id[subjectHits(hits)],
        stringsAsFactors = FALSE) -> hits_df

    hits_df %>%
        dplyr::group_by(TF_id) %>%
        dplyr::summarise(norm = dplyr::n() / length(peak_GR)) -> TF_normFactor

    hits_df %>%
        dplyr::group_by(TF_id, feature_id) %>%
        dplyr::summarise(hit_N = dplyr::n()) -> TF_hitN

    # some feature_id may not hit by your TF
    # which means feature in all TF will be 0
    fill_hitN <- data.frame(TF_id = rep(all_TF, each = length(input_feature_id)),
                            feature_id = rep(input_feature_id, length(all_TF))
                            )

    fill_hitN %>%
        dplyr::left_join(TF_hitN) %>%
        replace(is.na(.), 0) %>%
        tidyr::pivot_wider(names_from = TF_id,
                           values_from = hit_N) -> hitN_select_wide

    hitN_mt <- as.data.frame(hitN_select_wide[, -1])
    rownames(hitN_mt) <- hitN_select_wide$feature_id
    hitN_mt <- hitN_mt[, TF_normFactor$TF_id, drop = FALSE]
    hitN_mt <- sweep(hitN_mt, 2, TF_normFactor$norm, "-")

    peakScoreMt_select <- peakScoreMt[input_feature_id, ]

    if (log) {
        peakScoreMt_select <- log(peakScoreMt_select + 1)
    }

    if (meanScale){
        peakScoreMt_select <- t(apply(peakScoreMt_select, 1,
                                      function(x) x - mean(x))
                                )
    }

    if (mean(rownames(peakScoreMt_select) == rownames(hitN_mt)) < 1){
        stop("something wrong")
    }

    suppressMessages(hitN_mt %>%
        dplyr::as_tibble(rownames = "feature_id") %>%
        dplyr::inner_join(
            dplyr::as_tibble(peakScoreMt_select,
                             rownames = "feature_id")
            )) -> merge_df

    purrr::map(colnames(peakScoreMt_select), function(sample){

        train_data <- merge_df[, c(sample, all_TF)]
        colnames(train_data)[1] <- "value"

        if (length(all_TF) == 1) {
            tibble::tibble(TF_id = all_TF,
                           cor = cor(train_data)[1, -1]) -> activity_df

            colnames(activity_df)[2] <- sample

            return(activity_df)

        }


        x = model.matrix(value ~ ., train_data)[, -1, drop = FALSE]
        y = train_data$value

        cv.out = cv.glmnet(x, y, alpha = 0)

        model <- glmnet(x, y, alpha = 0, lambda = cv.out$lambda.min)
        data.frame(coef(model)[-1, 1]) %>%
            dplyr::as_tibble(rownames = "TF_id") -> activity_df
        activity_df[, 2] <- INT(activity_df[, 2])
        colnames(activity_df)[2] <- sample

        # should we need cor info ?
        # cor_data <- cor(train_data)[1, -1]
        #
        # suppressMessages(
        #     activity_df %>%
        #         dplyr::inner_join(
        #             data.frame(TF_id = names(cor_data),
        #                        cor = cor_data)
        #         )
        # ) -> activity_df
        #
        # colnames(activity_df)[c(2,3)] <- paste0(sample, c("", "_cor"))


        return(activity_df)
        #activity_df[, 2] <-

    }) -> TF_activity_list

    suppressMessages(
        purrr::reduce(TF_activity_list, inner_join)
        ) -> TF_activity_df

    return(TF_activity_df)
}
