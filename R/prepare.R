#' loadPeakFile
#'
#' read peak file and transform it into GRanges object
#'
#' @param filePath peak Path
#' @param TFBS_database whether your peak file is a TFBS database file. If you
#' want the final GRanges have a column named "TF_id", you should set TFBS_database
#' TRUE. The GRanges with TF_id can be applied in "TF_GR_database" parameter of
#' findIT_TFHit, findIT_enrichFisher, findIT_enrichWilcox, findIT_regionRP. If
#' FALSE, the GRanges will have a column named "feature_id", which always be the input
#' of "peak_GR" parameter.
#'
#'
#' @return GRanges object with a column named feature_id or TF_id
#' @export
#'
#' @examples
#' peakfile <- system.file("extdata", "ChIP.bed.gz", package = "FindIT2")
#' loadPeakFile(peakfile)
#'
#' @details
#' The GRanges with TF_id always be the input of "TF_GR_database" parameter. It
#' represents the TFBS database like motif scan result, public database ChIP-seq
#' site and so on.
#'
#' The GRanges with feature_id always be the input of "peak_GR" parameter.
#'
loadPeakFile <- function(filePath,
                         TFBS_database = FALSE) {

    peak_GR <- rtracklayer::import(con = filePath)
    if (TFBS_database){
        colnames(S4Vectors::mcols(peak_GR))[1] <- "TF_id"
    } else {
        colnames(S4Vectors::mcols(peak_GR))[1] <- "feature_id"
    }
    return(peak_GR)
}
