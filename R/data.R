# https://stackoverflow.com/questions/24434438/r-package-build-undocumented-code-objects/42094220

#' RNA normCount of E50h-72h in Chr5
#'
#' @usage data(RNA_normCount)
#'
#' @format A matrix
#'
#' @source \url{https://doi.org/10.1016/j.devcel.2020.07.003}
#'
"RNA_normCount"


#' ATAC normCount of E50h-72h in Chr5
#'
#' @usage data(ATAC_normCount)
#'
#' @format A matrix
#'
#' @source \url{https://doi.org/10.1016/j.devcel.2020.07.003}
#'
"ATAC_normCount"

#' RNA diff result from LEC2_GR VS LEC2_DMSO
#'
#' @usage data(RNADiff_LEC2_GR)
#'
#' @format a data frame
#'
#' @source \url{https://doi.org/10.1016/j.devcel.2020.07.003}
#'
"RNADiff_LEC2_GR"


#' TF-target database
#'
#' @usage data(TF_target_database)
#'
#' @format a data frame
#'
#' @source \url{http://bioinformatics.psb.ugent.be/webtools/iGRN/pages/download}
"TF_target_database"


#' input_genes
#'
#' @usage data(input_genes)
#'
#' @format gene vector
#'
#' @examples
#' \dontrun{
#' library(FindIT2)
#' library(TxDb.Athaliana.BioMart.plantsmart28)
#' Txdb <- TxDb.Athaliana.BioMart.plantsmart28
#' seqlevels(Txdb) <- c(paste0("Chr", 1:5), "M", "C")
#' ChIP_peak_path <- system.file("extdata", "ChIP.bed.gz", package = "FindIT2")
#' ChIP_peak_GR <- loadPeakFile(ChIP_peak_path)
#' ATAC_peak_path <- system.file("extdata", "ATAC.bed.gz", package = "FindIT2")
#' ATAC_peak_GR <- loadPeakFile(ATAC_peak_path)
#'
#' mmAnno_geneScan <- mm_geneScan(peak_GR = ChIP_peak_GR,
#'                                Txdb = Txdb,
#'                                upstream = 2e4,
#'                                downstream = 2e4)
#'
#' peakRP_gene <- calcRP_TFHit(mmAnno = mmAnno_geneScan,
#'                             Txdb = Txdb,
#'                             report_fullInfo = FALSE)
#'
#' data("RNADiff_LEC2_GR")
#' integrate_ChIP_RNA(result_geneRP = peakRP_gene,
#' result_geneDiff = RNADiff_LEC2_GR) -> merge_result
#'
#' target_result <- merge_result$data
#' input_genes <- target_result$gene_id[1:50]
#'
#' related_peaks <- mm_geneBound(peak_GR = ATAC_peak_GR,
#'                               Txdb = Txdb,
#'                               input_genes = input_genes)
#' input_feature_id <- related_peaks$feature_id
#' #save(input_genes, file = "data/input_genes.rda", version = 2)
#' #save(input_feature_id, file = "data/input_feature_id.rda", version = 2)
#' }
#'
#'
"input_genes"


#' input_feature_id
#'
#' @usage data(input_feature_id)
#'
#' @format feature_id vector
#'
#' @examples
#' \dontrun{
#' library(FindIT2)
#' library(TxDb.Athaliana.BioMart.plantsmart28)
#' Txdb <- TxDb.Athaliana.BioMart.plantsmart28
#' seqlevels(Txdb) <- c(paste0("Chr", 1:5), "M", "C")
#' ChIP_peak_path <- system.file("extdata", "ChIP.bed.gz", package = "FindIT2")
#' ChIP_peak_GR <- loadPeakFile(ChIP_peak_path)
#' ATAC_peak_path <- system.file("extdata", "ATAC.bed.gz", package = "FindIT2")
#' ATAC_peak_GR <- loadPeakFile(ATAC_peak_path)
#'
#' mmAnno_geneScan <- mm_geneScan(peak_GR = ChIP_peak_GR,
#'                                Txdb = Txdb,
#'                                upstream = 2e4,
#'                                downstream = 2e4)
#'
#' peakRP_gene <- calcRP_TFHit(mmAnno = mmAnno_geneScan,
#'                             Txdb = Txdb,
#'                             report_fullInfo = FALSE)
#'
#' data("RNADiff_LEC2_GR")
#' integrate_ChIP_RNA(result_geneRP = peakRP_gene,
#' result_geneDiff = RNADiff_LEC2_GR) -> merge_result
#'
#' target_result <- merge_result$data
#' input_genes <- target_result$gene_id[1:50]
#'
#' related_peaks <- mm_geneBound(peak_GR = ATAC_peak_GR,
#'                               Txdb = Txdb,
#'                               input_genes = input_genes)
#' input_feature_id <- related_peaks$feature_id
#' #save(input_genes, file = "data/input_genes.rda", version = 2)
#' #save(input_feature_id, file = "data/input_feature_id.rda", version = 2)
#' }
#'
#'
"input_feature_id"
