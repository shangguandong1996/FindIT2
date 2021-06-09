utils::globalVariables(c("feature_id", "gene_id"))

#' mm_nearestGene
#'
#' Annotate peaks using nearest gene mode, which means every peak only have one
#' related gene.
#'
#' @import GenomeInfoDb GenomicFeatures S4Vectors BiocGenerics GenomicRanges
#'
#' @param peak_GR peak GRange with a column named feature_id representing you
#' peak name
#' @param Txdb Txdb
#' @param reportGeneInfo whether you want to report full gene info
#' @param ... additional arguments in distanceToNearest
#'
#' @return Granges object with annotated info
#' @export
#'
#' @examples
#' library(TxDb.Athaliana.BioMart.plantsmart28)
#' Txdb <- TxDb.Athaliana.BioMart.plantsmart28
#' seqlevels(Txdb) <- c(paste0("Chr", 1:5), "M", "C")
#'
#' peak_path <- system.file("extdata", "ChIP.bed.gz", package = "FindIT2")
#' peak_GR <- loadPeakFile(peak_path)
#' peakAnno <- mm_nearestGene(peak_GR,Txdb)
#' peakAnno
#'
mm_nearestGene <- function(peak_GR,
                           Txdb,
                           reportGeneInfo = FALSE,
                           ...) {
    if ("mmAnno_mode" %in% names(metadata(peak_GR))) {
        warning("it seems you have done mmAnno before, please not do again",
                call. = FALSE
        )
        return(peak_GR)
    }

    if (!"feature_id" %in% colnames(mcols(peak_GR))) {
        addfeatureIdSuggsesion()
    }

    # warning will appear firstly
    withr::local_options(list(warn = 1))
    check_seqlevel(peak_GR, Txdb)
    peak_GR <- peak_GR[seqnames(peak_GR) %in% seqlevels(Txdb)]

    cat("------------\n")
    cat("annotatePeak using nearest gene mode begins\n")
    cat(
        ">> preparing gene features information...\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n"
    )
    gene_location <- genes(Txdb)
    gene_start <- resize(gene_location, width = 1, fix = "start")

    cat(
        ">> finding nearest gene and calculating distance...\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n"
    )

    # It seems that there is two distanceToNearest method,
    # one is from GenomicGRanges, another is from IRanges
    # both support ignore.strand while the latter do not show in help
    # if peak is unstrand(for most cases), the ignore.strand doesn't matter
    nearestHit <- suppressWarnings(distanceToNearest(
        peak_GR,
        gene_start, ...
    ))

    if (length(nearestHit) == 0) {
        stop("Sorry, there is no hits.",
             call. = FALSE
        )
    }

    cat(
        ">> dealing with gene strand ...\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n"
    )
    data.frame(
        start_dis = start(peak_GR[queryHits(nearestHit)]) -
            start(gene_start[subjectHits(nearestHit)]),
        end_dis = end(peak_GR[queryHits(nearestHit)]) -
            start(gene_start[subjectHits(nearestHit)])
    ) -> a

    # if start and end are all in TSS left, sum(x > 0)) it will be 2.
    # if start is in left while end is in right, sum(x > 0)) it will be 1
    # if start and end are all in TSS right, sum(x > 0)) it will be 0
    sign_1 <- ifelse(apply(a, 1, function(x) sum(x > 0)) == 2, 1, -1)

    # sign_1 and sign_2 will decide distance's sign
    sign_2 <- ifelse(strand(gene_start[subjectHits(nearestHit)]) == "+", 1, -1)

    gene_id <- names(gene_start)[subjectHits(nearestHit)]

    cat(
        ">> merging all info together ...\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n"
    )

    if (reportGeneInfo) {
        # someone's peakN may be too big, in this case,
        # they may do not want full gene info
        # so they can save disk or memory
        mcols(peak_GR) <- cbind(mcols(peak_GR),
                                report_geneInfo(gene_location[gene_id]))
        mcols(peak_GR)$distanceToTSS <- mcols(nearestHit)$distance * sign_1 * sign_2
    } else {
        mcols(peak_GR)$gene_id <- gene_id
        mcols(peak_GR)$distanceToTSS <- mcols(nearestHit)$distance * sign_1 * sign_2
    }

    metadata(peak_GR)$mmAnno_mode <- "nearestGene"

    cat(
        ">> done\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n"
    )
    return(peak_GR)
}


#' mm_geneScan
#'
#' Annotate peaks using geneScan mode, which means every peak have more than one
#' related genes.
#'
#' @import GenomeInfoDb GenomicFeatures S4Vectors BiocGenerics GenomicRanges
#'
#' @param peak_GR peak GRange with a column named feature_id representing you peak name
#' @param Txdb Txdb
#' @param upstream distance to start site(upstream)
#' @param downstream distance to start site(downstream)
#' @param reportGeneInfo whether you want to add gene info
#' @param ... additional arguments in findOverlaps
#'
#' @return Granges object with annotated info
#' @export
#'
#' @examples
#' library(TxDb.Athaliana.BioMart.plantsmart28)
#' Txdb <- TxDb.Athaliana.BioMart.plantsmart28
#' seqlevels(Txdb) <- c(paste0("Chr", 1:5), "M", "C")
#' peak_path <- system.file("extdata", "ChIP.bed.gz", package = "FindIT2")
#' peak_GR <- loadPeakFile(peak_path)
#' peakAnno <- mm_geneScan(peak_GR,Txdb)
#' peakAnno
#'
mm_geneScan <- function(peak_GR,
                        Txdb,
                        upstream = 3000,
                        downstream = 3000,
                        reportGeneInfo = FALSE,
                        ...) {
    if ("mmAnno_mode" %in% names(peak_GR@metadata)) {
        warning("it seems you have done mmAnno before, please not do again",
                call. = FALSE
        )
        return(peak_GR)
    }

    if (!"feature_id" %in% colnames(mcols(peak_GR))) {
        addfeatureIdSuggsesion()
    }

    # warning will appear firstly
    withr::local_options(list(warn = 1))
    check_seqlevel(peak_GR, Txdb)
    peak_GR <- peak_GR[seqnames(peak_GR) %in% seqlevels(Txdb)]

    cat("------------\n")
    cat("annotatePeak using geneScan mode begins\n")
    cat(
        ">> preparing gene features information and scan region...\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n"
    )
    gene_location <- GenomicFeatures::genes(Txdb)
    gene_start <- resize(gene_location, width = 1, fix = "start")

    gene_promoter <- suppressWarnings(promoters(gene_location,
                                                               upstream = upstream,
                                                               downstream = downstream
    ))

    if (all(is.na(seqinfo(Txdb)@seqlengths))) {
        warning("your chr length is not set, gene_promoter may cross bound",
                call. = FALSE
        )
    } else {
        cat(
            ">> some scan range may cross Chr bound, trimming...\t\t",
            format(Sys.time(), "%Y-%m-%d %X"), "\n"
        )
        gene_promoter <- trim(gene_promoter)
    }

    cat(
        ">> finding overlap peak in gene scan region...\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n"
    )
    overlapHits <- suppressWarnings(findOverlaps(
        gene_promoter,
        peak_GR,
        ...
    ))

    cat(
        ">> dealing with left peak not your gene scan region...\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n"
    )
    # cat(">> these peak will be annotated by nearest gene\n")
    peak_GR_left <- subset(
        peak_GR,
        !feature_id %in% unique(peak_GR[subjectHits(overlapHits)]$feature_id)
    )

    # if I use nearest, it will report 0 when no peak_GR left
    # so I use distanceToNearest
    nearest_hit <- suppressWarnings(distanceToNearest(
        peak_GR_left,
        gene_start
    ))
    gene_TSS_left <- gene_start[subjectHits(nearest_hit)]

    cat(
        ">> merging two set peaks...\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n"
    )
    final_gene_GR <- c(gene_start[queryHits(overlapHits)], gene_TSS_left)
    final_peak_GR <- c(peak_GR[subjectHits(overlapHits)], peak_GR_left)

    cat(
        ">> calculating distance and dealing with gene strand...\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n"
    )
    dist <- suppressWarnings(distance(
        final_peak_GR,
        final_gene_GR
    ))

    data.frame(
        start_dis = start(final_peak_GR) - start(final_gene_GR),
        end_dis = end(final_peak_GR) - start(final_gene_GR)
    ) -> a

    sign_1 <- ifelse(apply(a, 1, function(x) sum(x > 0)) == 2, 1, -1)
    sign_2 <- ifelse(strand(final_gene_GR) == "+", 1, -1)


    gene_id <- names(final_gene_GR)

    cat(
        ">> merging all info together ...\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n"
    )
    if (reportGeneInfo) {
        mcols(final_peak_GR) <- cbind(
            mcols(final_peak_GR),
            report_geneInfo(gene_location[gene_id])
        )
        mcols(final_peak_GR)$distanceToTSS <- dist * sign_1 * sign_2
    } else {
        mcols(final_peak_GR)$gene_id <- gene_id
        mcols(final_peak_GR)$distanceToTSS <- dist * sign_1 * sign_2
    }

    metadata(final_peak_GR)$mmAnno_mode <- "geneScan"
    metadata(final_peak_GR)$upstream <- upstream
    metadata(final_peak_GR)$downstream <- downstream

    cat(
        ">> done\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n"
    )
    return(final_peak_GR)
}


#' mm_geneBound
#'
#' find related peaks of your input genes, which is useful when you want to plot
#' volcano plot or heatmap of peaks.
#'
#' @import GenomeInfoDb GenomicFeatures S4Vectors BiocGenerics GenomicRanges
#' @importFrom utils capture.output
#'
#' @param peak_GR peak GRange with a column named feature_id representing you peak name
#' @param Txdb Txdb
#' @param input_genes genes which you want to find related peak for
#' @param ... additional arguments in distanceToNearest
#'
#' @return data.frame with three column: related peak id, your input gene id,
#' and distance
#' @export
#'
#' @examples
#' library(TxDb.Athaliana.BioMart.plantsmart28)
#' Txdb <- TxDb.Athaliana.BioMart.plantsmart28
#' seqlevels(Txdb) <- c(paste0("Chr", 1:5), "M", "C")
#' peak_path <- system.file("extdata", "ChIP.bed.gz", package = "FindIT2")
#' peak_GR <- loadPeakFile(peak_path)
#' peak_pair <- mm_geneBound(peak_GR, Txdb, c("AT5G01015", "AT5G67570"))
#' peak_pair
#'
mm_geneBound <- function(peak_GR,
                         Txdb,
                         input_genes,
                         ...) {
    # some genes maybe 12345……
    input_genes <- as.character(unique(input_genes))

    if ("gene_id" %in% colnames(mcols(peak_GR))) {
        mmAnno <- peak_GR
    } else {
        cat(
            ">> using mm_nearestGene to annotate Peak...\t\t",
            format(Sys.time(), "%Y-%m-%d %X"), "\n"
        )
        # https://stackoverflow.com/questions/34208564/how-to-hide-or-disable-in-function-printed-message
        invisible(capture.output(mm_nearestGene(
            peak_GR = peak_GR,
            Txdb = Txdb,
            reportGeneInfo = FALSE,
            ...
        ) -> mmAnno))
    }


    noPeakBind_genes <- input_genes[!input_genes %in% mmAnno$gene_id]

    gene_location <- GenomicFeatures::genes(Txdb)
    all_genes <- gene_location$gene_id

    mmAnno$distanceToTSS_abs <- abs(mmAnno$distanceToTSS)
    mmAnno_inputGenes <- mcols(subset(mmAnno, gene_id %in% input_genes))[, c("feature_id", "gene_id", "distanceToTSS_abs")]
    mmAnno_inputGenes <- data.frame(mmAnno_inputGenes, stringsAsFactors = FALSE)

    if (length(noPeakBind_genes) == 0) {
        message("It seems that all your input gene have been annotated",
             call. = FALSE
        )
        return(mmAnno_inputGenes)


    } else if (mean(noPeakBind_genes %in% all_genes) < 1) {
        warning("some of your input genes is not in your Txdb\n",
                call. = FALSE
        )
        noPeakBind_genes <- noPeakBind_genes[noPeakBind_genes %in% all_genes]
    }

    cat(
        ">> It seems that there", length(noPeakBind_genes), "genes do not have realted peak...\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n"
    )

    cat(
        ">> using distanceToNearest to find nearest peak of these genes...\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n"
    )
    gene_start <- resize(gene_location, width = 1, fix = "start")
    gene_GR <- subset(gene_start, gene_id %in% noPeakBind_genes)



    nearestHit <- suppressWarnings(distanceToNearest(
        gene_GR,
        peak_GR, ...
    ))


    data.frame(
        feature_id = peak_GR[subjectHits(nearestHit)]$feature_id,
        gene_id = gene_GR[queryHits(nearestHit)]$gene_id,
        distanceToTSS_abs = mcols(nearestHit)$distance,
        stringsAsFactors = FALSE
    ) -> left_anno

    cat(
        ">> merging all anno...\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n"
    )
    final_anno <- rbind(
        mmAnno_inputGenes,
        left_anno
    )

    cat(
        ">> done\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n"
    )
    return(final_anno)
}
