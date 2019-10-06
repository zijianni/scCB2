#' Extract real cell matrix from \code{CB2FindCell} output
#'
#' Handy function to extract real cell matrix from \code{CB2FindCell} output. 
#' It also provides the option to filter out broken cells based on proportion 
#' of mitochondrial gene expressions.
#'
#' @param CBout Output object from \code{CB2FindCell}.
#'
#' @param MTfilter Numeric value between 0 to 1. Default: \code{1} 
#' (No filtering). Mitochondrial gene expression proportion filtering 
#' threshold for broken cells. The proportion of mitochondrial gene
#' expressions is calculated from the scaled sum of all genes starting 
#' with "MT-" (human) or "mt-" (mouse).
#'
#' @param MTgene Character vector. User may specify customized mitochondrial
#' gene IDs to perform the filtering.
#'
#' @return A \code{dgCMatrix} count matrix of real cells.
#'
#' @examples
#' data(mbrainSub)
#' CBOut <- CB2FindCell(mbrainSub, FDR_threshold = 0.01,
#'     background_threshold = 100, Ncores = 2)
#' RealCell <- GetCellMat(CBOut, MTfilter = 0.05)
#' str(RealCell)
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
        MTgene <- c(grep(pattern = "\\<MT-", x = rownames(xmat)),
                    grep(pattern = "\\<mt-", x = rownames(xmat)))
    } else if (!all(is.character(MTgene))) {
        stop("\"MTgene\" should be a character vector of gene names.")
    }
    if (length(MTgene) > 0) {
        MTprop <- colSums(xmat[MTgene,]) / colSums(xmat)
        brokenCell <- MTprop > MTfilter
        xmat <- xmat[,!brokenCell]
    }
    return(xmat)
}