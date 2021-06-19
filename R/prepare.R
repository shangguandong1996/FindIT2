#' loadPeakFile
#'
#' read peak file and transform it into GRanges object
#'
#' @param filePath peak Path
#' @param fromMACS2 whether peak from MACS2.If it set TRUE, extra column in bed
#' file will be named according to the content in MACS2 manual
#' @param narrowPeak whether peak is narrowPeak or broadPeak.
#'
#' @return GRanges object with a column named feature_id, which is important
#' for later analysis
#' @export
#'
#' @examples
#' peakfile <- system.file("extdata", "ChIP.bed.gz", package = "FindIT2")
#' loadPeakFile(peakfile)
loadPeakFile <- function(filePath,
    fromMACS2 = FALSE,
    narrowPeak = TRUE) {
    if (fromMACS2) {
        if (narrowPeak) {
            # This idea is from
            # https://charlesjb.github.io/How_to_import_narrowPeak/
            extraCols <- c(
                Fold = "numeric", pValue = "numeric",
                qValue = "numeric", relative_summit = "integer"
            )
            peak_GR <- rtracklayer::import(
                con = filePath,
                format = "BED",
                extraCols = extraCols
            )
        } else {
            extraCols <- c(
                Fold = "numeric", pValue = "numeric",
                qValue = "numeric"
            )
            peak_GR <- rtracklayer::import(
                con = filePath,
                format = "BED",
                extraCols = extraCols
            )
        }
    } else {
        peak_GR <- rtracklayer::import(con = filePath)
    }
    colnames(S4Vectors::mcols(peak_GR))[1] <- "feature_id"
    return(peak_GR)
}
