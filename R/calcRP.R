utils::globalVariables(c("RP", "sumRP"))
utils::globalVariables(c("bw", "centerToTSS", "RP", "sumRP", "logRP", "normRP"))

#' calcRP_coverage
#'
#' calculate regulatory potential using big wig files, which is useful for ATAC
#' or H3K27ac histone modification data.
#'
#' @importFrom IRanges Views viewApply
#'
#' @param bwFile bw file
#' @param Txdb Txdb
#' @param gene_included genes which you want to calcluate RP for
#' @param Chrs_included chromosome where you want to calcluate gene RP in
#' @param decay_dist decay distance
#' @param scan_dist scan distance
#'
#' @return data.frame
#' @details
#' Please note that because of rtracklayer::import has some issue on 32 bit R of windows, so the
#' calcRP_coverage can not work on this system. But if your R is 64 bit,
#' which now be applied on the most windows R, this function still work.
#'
#' @export
#' @examples
#' if (.Platform$OS.type != "windows" & require(TxDb.Athaliana.BioMart.plantsmart28)) {
#'     Txdb <- TxDb.Athaliana.BioMart.plantsmart28
#'     seqlevels(Txdb) <- paste0("Chr", c(1:5, "M", "C"))
#'     bwFile <- system.file("extdata", "E50h_sampleChr5.bw", package = "FindIT2")
#'
#'     RP_df <- calcRP_coverage(
#'         bwFile = bwFile,
#'         Txdb = Txdb,
#'         Chrs_included = "Chr5"
#'     )
#'
#' }
calcRP_coverage <- function(bwFile,
                            Txdb,
                            gene_included,
                            Chrs_included,
                            decay_dist = 1e3,
                            scan_dist = 2e4) {
    cat(
        ">> preparing gene features information...\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n"
    )
    if (missing(gene_included)) {
        genes_GR <- genes(Txdb)
    } else {
        genes_GR <- genes(Txdb)[gene_included]
    }

    if (missing(Chrs_included)) {
        Chrs_included <- seqlevels(Txdb)
    }

    if (all(is.na(seqinfo(Txdb)@seqlengths))) {
        stop("your chr length is not set, scan region may cross bound. please set your chr length in Txdb",
            call. = FALSE
        )
    }

    scan_region <- suppressWarnings(promoters(genes_GR,
        upstream = scan_dist,
        downstream = scan_dist + 1
    ))

    TSS_location <- resize(genes_GR,
        width = 1,
        fix = "start"
    )

    cat(
        ">> some scan range may cross Chr bound, trimming...\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n"
    )
    scan_region <- trim(scan_region)

    cat(
        ">> preparing weight info...\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n"
    )
    # https://github.com/qinqian/lisa/blob/9fef7cb682264bdef5cca6847d09acf6c92a08f2/lisa/utils.py#L79-L98
    alpha <- -log(1 / 3) * (scan_dist / decay_dist)
    d <- -scan_dist:scan_dist
    weight <- (2 * exp(-alpha * abs(d) / scan_dist)) / (1 + exp(-alpha * abs(d) / scan_dist))

    cat(
        ">> loading", basename(bwFile), "info...\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n"
    )
    bwInfo <- rtracklayer::import(bwFile,
        format = "Bigwig",
        as = "RleList"
    )

    signal_list <- list()

    for (Chr in Chrs_included) {
        cat("------------\n")
        cat(
            ">> extracting and calcluating", Chr, "signal...\t\t",
            format(Sys.time(), "%Y-%m-%d %X"), "\n"
        )

        scan_region_included <- scan_region[seqnames(scan_region) == Chr &
            width(scan_region) == (2 * scan_dist + 1)]
        view <- Views(bwInfo[Chr][[1]], ranges(scan_region_included))
        signal_Part1 <- viewApply(view, function(x) sum(as.numeric(x) * weight))

        # # a test
        # value <- vector(mode = "numeric", length = 1000)
        #
        # for (i in 1:1000) {
        #   sum(as.numeric(view[i][[1]]) * weight) -> value[i]
        # }
        #
        # mean(value == viewApply(view, function(x) sum(as.numeric(x) * weight))[1:1000])

        scan_region_left <- scan_region[seqnames(scan_region) == Chr &
            width(scan_region) < (2 * scan_dist + 1)]

        if (length(scan_region_left) > 0) {
            cat(
                ">> dealing with", Chr, "left gene signal...\t\t",
                format(Sys.time(), "%Y-%m-%d %X"), "\n"
            )
            gene_left <- names(scan_region_left)

            signal_Part2 <- vector(mode = "numeric", length = length(scan_region_left))

            for (i in seq_along(gene_left)) {
                gene_IR <- ranges(scan_region_left[i])
                value <- as.numeric(Views(bwInfo[Chr][[1]], gene_IR)[[1]])

                # deal with cross boundary
                if (start(gene_IR) == 1) {
                    weight_index <- (scan_dist * 2 - width(gene_IR) + 2):(2 * scan_dist + 1)
                } else {
                    weight_index <- 1:width(gene_IR)
                }

                signal_Part2[i] <- sum(value * weight[weight_index])
            }


            names(signal_Part2) <- gene_left

            signal_merge <- c(signal_Part1, signal_Part2)
        } else {
            signal_merge <- signal_Part1
        }

        # normalized per Chr
        # different Chrs has own RP profile
        cat(
            ">> norming", Chr, "RP accoring to the whole Chr RP...\t\t",
            format(Sys.time(), "%Y-%m-%d %X"), "\n"
        )
        signal_norm <- log(signal_merge + 1) - mean(log(signal_merge + 1))

        data.frame(
            gene_id = names(signal_norm),
            sumRP = signal_norm,
            stringsAsFactors = FALSE
        ) -> signal_list[[Chr]]
    }

    cat(
        ">> merging all Chr RP together...\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n"
    )
    final_RP <- do.call(rbind, signal_list)
    rownames(final_RP) <- NULL

    cat(
        ">> done\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n\n"
    )
    return(final_RP)
}


#' calcRP_region
#'
#' calculate regulatory potential based on mm_geneScan result and peakCount matrix,
#' which is useful for ATAC or H3K27ac histone modification data.
#'
#' @importFrom magrittr %>%
#' @import rlang
#'
#' @param mmAnno the annotated GRange object from mm_geneScan
#' @param peakScoreMt peak count matrix. The rownames are feature_id in mmAnno,
#' while the colnames are sample names
#' @param Txdb Txdb
#' @param Chrs_included chromosome where you want to calculate gene RP in.
#' If Chromosome is not be set, it will calculate gene RP in all chromosomes in Txdb.
#' @param decay_dist decay distance
#' @param log_transform whether you want to log and norm your RP
#'
#' @return a MultiAssayExperiment object containg detailed peak-RP-gene relationship
#' and sumRP info
#' @export
#'
#' @examples
#'
#' if (require(TxDb.Athaliana.BioMart.plantsmart28)) {
#'     data("ATAC_normCount")
#'     library(SummarizedExperiment)
#'     Txdb <- TxDb.Athaliana.BioMart.plantsmart28
#'     seqlevels(Txdb) <- paste0("Chr", c(1:5, "M", "C"))
#'
#'     peak_path <- system.file("extdata", "ATAC.bed.gz", package = "FindIT2")
#'     peak_GR <- loadPeakFile(peak_path)
#'     mmAnno <- mm_geneScan(peak_GR, Txdb)
#'
#'     calcRP_region(
#'         mmAnno = mmAnno,
#'         peakScoreMt = ATAC_normCount,
#'         Txdb = Txdb,
#'         Chrs_included = "Chr5"
#'     ) -> regionRP
#'
#'     sumRP <- assays(regionRP)$sumRP
#'     fullRP <- assays(regionRP)$fullRP
#' }
#'
calcRP_region <- function(mmAnno,
                          peakScoreMt,
                          Txdb,
                          Chrs_included,
                          decay_dist = 1e3,
                          log_transform = FALSE) {

    check_parameter_length(mmAnno, decay_dist)

    cat(
        "\n>> calculating peakCenter to TSS using peak-gene pair...\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n"
    )
    # names(peak_GR) <- peak_GR$feature_id
    # peak_scaned <- peak_GR[scan_result$feature_id]
    peak_scaned_center <- GenomicRanges::resize(mmAnno,
        width = 1,
        fix = "center"
    )

    if (missing(Chrs_included)) {
        Chrs_included <- GenomeInfoDb::seqlevels(Txdb)
    }
    all_gene_location <- GenomicFeatures::genes(Txdb)
    all_gene_location <- all_gene_location[GenomeInfoDb::seqnames(all_gene_location) %in% Chrs_included]


    gene_scaned <- all_gene_location[mmAnno$gene_id]
    gene_scaned_TSS <- GenomicRanges::resize(gene_scaned,
        width = 1,
        fix = "start"
    )

    mmAnno$centerToTSS <- GenomicRanges::distance(peak_scaned_center, gene_scaned_TSS)


    # In R 3.6 GenomicGRanges, if duplicted names in GRanges
    # as.data.frame will report error(R 4.10 not)
    # So I set mmAnno into NULL to avoid this error, so that I can test on server
    names(mmAnno) <- NULL

    scan_result <- as.data.frame(mmAnno)
    scan_result <- scan_result[, c(
        "seqnames", "feature_id",
        "gene_id", "centerToTSS"
    )]


    geneAll <- all_gene_location$gene_id
    noRPGene <- geneAll[!geneAll %in% unique(mmAnno$gene_id)]
    cat(
        ">> pre-filling", length(noRPGene), "noAssoc peak gene's RP with 0...\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n"
    )
    noRPGene_df <- data.frame(
        seqnames = as.character(GenomeInfoDb::seqnames(all_gene_location[noRPGene])),
        gene_id = noRPGene,
        sumRP = 0,
        stringsAsFactors = FALSE
    )

    RP_list <- list()

    fullRP_mt <- matrix(
        data = NA_integer_,
        nrow = nrow(scan_result),
        ncol = dim(peakScoreMt)[2]
    )
    colnames(fullRP_mt) <- colnames(peakScoreMt)

    cat(
        ">> calculating RP using centerToTSS and peak score",
        format(Sys.time(), "%Y-%m-%d %X"), "\n"
    )


    for (i in seq_len(dim(peakScoreMt)[2])) {
        sample_name <- colnames(peakScoreMt)[i]

        sample_peakScore <- data.frame(
            feature_id = rownames(peakScoreMt),
            score = peakScoreMt[, i],
            stringsAsFactors = FALSE
        )
        # colnames(sample_peakScore)[2] <- sample_name

        # The gene in Txdb which do not have assoc peak should be fill 0 at first !!
        # This fill can be useful when comapred with geneExpr data
        # some interesting gene may no RP

        suppressMessages(dplyr::left_join(scan_result, sample_peakScore)) %>%
            dplyr::mutate(RP = score * 2^(-centerToTSS / decay_dist)) -> peakRP_gene
        # mutate(RP = !!sym(sample_name) * 2^(-centerToTSS/decay_dist)) -> peakRP_gene

        fullRP_mt[, sample_name] <- peakRP_gene$RP

        peakRP_gene %>%
            dplyr::group_by(gene_id) %>%
            dplyr::mutate(sumRP = sum(RP)) %>%
            dplyr::select(c("seqnames", "gene_id", "sumRP")) %>%
            dplyr::distinct(gene_id, .keep_all = TRUE) %>%
            dplyr::bind_rows(noRPGene_df) -> calc_result

        if (log_transform) {
            # this log and per chrom norm idea is from
            # Qin, Q. et al. (2020). Lisa: inferring transcriptional regulators through integrative modeling of public chromatin accessibility and ChIP-seq data. Genome Biol 21, 32.
            # Methos:Preprocessing of chromatin profiles

            # log norm can help you find high RP gene in several samples or

            calc_result %>%
                dplyr::mutate(logRP = log(sumRP + 1)) %>%
                # which makes normRP has meaning
                dplyr::ungroup() %>%
                dplyr::group_by(seqnames) %>%
                dplyr::mutate(normRP = logRP - mean(logRP)) %>%
                # the min logRP is 0
                dplyr::ungroup() %>%
                # a gene has negative RP means this gene RP is smaller than mean
                dplyr::select(gene_id, normRP) %>%
                dplyr::rename(!!sym(sample_name) := normRP) -> RP_list[[i]]
        } else {
            # no log norm may be suitable for finding tissue-specific RP gene
            # when ploting heatmap

            # you can use sampleN/median(sample_all) to plot heatmap
            # this idea is from
            # Wang, S. et al. (2016). Modeling cis-regulation with a compendium of genome-wide histone H3K27ac profiles. Genome Res. 26, 1417–1429.
            # Fig.2B

            # no log norm RP is suitable for latter analysis, like findIT

            calc_result %>%
                dplyr::ungroup() %>%
                dplyr::select(gene_id, sumRP) %>%
                dplyr::rename(!!sym(sample_name) := sumRP) -> RP_list[[i]]
        }
    }

    cat(
        ">> merging all info together\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n"
    )
    sum_RP <- suppressMessages(purrr::reduce(RP_list, inner_join))
    rownames_gene <- sum_RP$gene_id
    sum_RP <- as.matrix(sum_RP[, -1])
    rownames(sum_RP) <- rownames_gene

    rownames(fullRP_mt) <- paste0(mmAnno$feature_id, ":", mmAnno$gene_id)
    names(mmAnno) <- paste0(mmAnno$feature_id, ":", mmAnno$gene_id)


    data_fullRP <- SummarizedExperiment::SummarizedExperiment(
        assays = fullRP_mt,
        rowRanges = mmAnno
    )

    final_result <- MultiAssayExperiment::MultiAssayExperiment(list(
        "fullRP" = data_fullRP,
        "sumRP" = sum_RP
    ))

    S4Vectors::metadata(final_result)$Chrs_included <- Chrs_included
    S4Vectors::metadata(final_result)$decay_dist <- decay_dist
    S4Vectors::metadata(final_result)$log_transform <- log_transform

    cat(
        ">> done\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n\n"
    )
    return(final_result)
}

#' calcRP_TFHit
#'
#' calculate regulatory potential based on ChIP-Seq peak data, which is useful
#' for TF ChIP-seq data.
#'
#' @importFrom magrittr %>%
#'
#' @param mmAnno the annotated GRange object from mm_geneScan
#' @param Txdb Txdb
#' @param decay_dist decay distance
#' @param report_fullInfo whether you want to report full peak-RP-gene info
#'
#' @return if report_fullInfo is TRUE, it will output GRanges with detailed info.
#' While FALSE, it will output data frame
#' @export
#'
#' @details
#' If your origin peak_GR of mmAnno have column named feature_score, calcRP_TFHit
#' will consider this column when calculating sumRP. Otherwise, it will consider
#' all peak Hit feature_score is 1.
#'
#' @examples
#' if (require(TxDb.Athaliana.BioMart.plantsmart28)){
#'     Txdb <- TxDb.Athaliana.BioMart.plantsmart28
#'     seqlevels(Txdb) <- paste0("Chr", c(1:5, "M", "C"))
#'     peak_path <- system.file("extdata", "ChIP.bed.gz", package = "FindIT2")
#'     peak_GR <- loadPeakFile(peak_path)
#'     mmAnno <- mm_geneScan(peak_GR, Txdb)
#'
#'     # if you just want to get RP_df, you can set report_fullInfo FALSE
#'     fullRP_hit <- calcRP_TFHit(
#'         mmAnno = mmAnno,
#'         Txdb = Txdb,
#'         report_fullInfo = TRUE
#'     )
#'
#'     RP_df <- metadata(fullRP_hit)$peakRP_gene
#'
#' }
calcRP_TFHit <- function(mmAnno,
                         Txdb,
                         decay_dist = 1000,
                         report_fullInfo = FALSE) {
    check_parameter_length(mmAnno, decay_dist)

    cat(
        "\n>> calculating peakCenter to TSS using peak-gene pair...\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n"
    )

    peak_scaned_center <- GenomicRanges::resize(mmAnno,
        width = 1,
        fix = "center"
    )

    all_gene_location <- GenomicFeatures::genes(Txdb)

    gene_scaned <- all_gene_location[mmAnno$gene_id]
    gene_scaned_TSS <- GenomicRanges::resize(gene_scaned,
        width = 1,
        fix = "start"
    )

    mmAnno$centerToTSS <- GenomicRanges::distance(peak_scaned_center, gene_scaned_TSS)

    scan_result <- as.data.frame(mmAnno)


    if ("feature_score" %in% colnames(mcols(mmAnno))) {
        cat(
            ">> calculating RP using centerToTSS and feature_score\t\t",
            format(Sys.time(), "%Y-%m-%d %X"), "\n"
        )
        mmAnno$RP <- mmAnno$feature_score * 2^(-mmAnno$centerToTSS / decay_dist)
    } else {
        cat(
            ">> calculating RP using centerToTSS and TF hit\t\t",
            format(Sys.time(), "%Y-%m-%d %X"), "\n"
        )
        mmAnno$RP <- 2^(-mmAnno$centerToTSS / decay_dist)
    }

    cat(
        ">> merging all info together\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n"
    )
    mcols(mmAnno) %>%
        data.frame(stringsAsFactors = FALSE) %>%
        dplyr::group_by(gene_id) %>%
        dplyr::summarise(
            withPeakN = dplyr::n(),
            sumRP = sum(RP)
        ) %>%
        dplyr::arrange(-sumRP) %>%
        dplyr::mutate(RP_rank = rank(-sumRP)) -> peakRP_gene

    cat(
        ">> done\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n\n"
    )

    if (report_fullInfo) {
        metadata(mmAnno)$peakRP_gene <- peakRP_gene
        return(mmAnno)
    } else {
        return(peakRP_gene)
    }

}
