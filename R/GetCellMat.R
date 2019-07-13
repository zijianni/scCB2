#' Extract real cell matrix from `CBFindCell` output
#' 
#' Simply function to extract real cell matrix from `CBFindCell` output. It also provides
#' the option to filter out broken cells based on proportion of mitochondrial gene expressions.
#'
#' @param CBout Output object from `CBFindCell`.
#' 
#' @param MTfilter Numeric value between 0 to 1. Mitochondrial gene expression proportion filtering threshold for broken cells. The proportion of mitochondrial gene
#' expressions is calculated from the scaled sum of all genes starting with "MT-" (human) or "mt-" (mouse).
#'
#' @param MTgene Character vector. User may use customized mitochondrial gene IDs to perform the filtering. 
#'
#' @example
#' data(mbrainSub)
#' CBOut <- CBFindCell(mbrainSub, FDR_threshold = 0.01, pooling_threshold = 100, Ncores = 4)
#' RealCell <- GetCellMat(CBOut, MTfilter = 0.05)
#' 
#' @importFrom 
#' 
#' 
#' @export


GetCellMat <- function(CBout,
                       MTfilter = 1,
                       MTgene = NULL) {
  if (!(is.numeric(MTfilter) && MTfilter >= 0 && MTfilter <= 1)) {
    stop("\"MTfilter\" should be a numeric value between 0 and 1.")
  }
  xmat <- cbind(CBout$cluster_matrix, CBout$cell_matrix)
  if (is.null(MTgene)) {
    MTgene <- c(grep(pattern = "\\<MT-", CBout = rownames(dat_norm)),
                grep(pattern = "\\<mt-", CBout = rownames(dat_norm)))
  } else if (!all(is.character(MTgene))) {
    stop("\"MTgene\" should be a character vector of mitochondrial gene names.")
  }
  if (length(MTgene) > 0) {
    MTprop <- colSums(CBout[MTgene, ]) / colSums(CBout)
    brokenCell <- MTprop > MTfilter
    xmat <- xmat[, !brokenCell]
  }
  return(xmat)
}