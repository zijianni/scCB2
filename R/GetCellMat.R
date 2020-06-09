#' Extract real cell matrix from \code{CB2FindCell} output and optionally
#' filter out low-quality cells
#'
#' Handy function to extract real cell matrix from \code{CB2FindCell} output. 
#' It provides the option to filter out broken cells based on proportion 
#' of mitochondrial gene expressions. The input can also be a sparse matrix
#' only for cell filtering.
#'
#' @param CBout Output object from \code{CB2FindCell}, or a sparse matrix 
#' (for example, from \code{QuickCB2}).
#'
#' @param MTfilter Numeric value between 0 and 1. Default: \code{1} 
#' (No filtering). For each barcode, if the proportion of mitochondrial 
#' gene expression exceeds \code{MTfilter}, this barcode will be filtered out.
#' No barcode exceeds 100\% mitochondrial gene expression, thus the default
#' (100\%) corresponds to no filtering. The proportion of mitochondrial 
#' gene expressions is usually a criterion for evaluating cell quality, 
#' and is calculated using the scaled sum of all genes starting 
#' with "MT-" (human) or "mt-" (mouse) if row names are gene symbols, 
#' or customized mitochondrial genes specified by \code{MTgene}.
#'
#' @param MTgene Character vector. User may specify customized mitochondrial
#' gene IDs to perform the filtering. This should correspond to a subset 
#' of row names in raw data.
#'
#' @return A \code{dgCMatrix} count matrix of real cells.
#'
#' @examples
#' 
#' # Please also refer to the example in function CB2FindCell.
#' 
#' # Simulate CB2FindCell output object.
#' data(mbrainSub)
#' mbrainReal <- mbrainSub[,Matrix::colSums(mbrainSub)>500]
#' CBOut <- list(
#'  cluster_matrix = mbrainReal[,
#'  sample(ncol(mbrainReal), 100, replace = TRUE)],
#'  cell_matrix = mbrainReal[,
#'  sample(ncol(mbrainReal), 100, replace = TRUE)])
#'                  
#' # Get cell matrix, filtering out barcodes with 
#' # more than 10% of counts from mitochondrial genes.     
#' 
#' RealCell <- GetCellMat(CBOut, MTfilter = 0.1)
#' str(RealCell)             
#'
#' @export


GetCellMat <- function(CBout,
                        MTfilter = 1,
                        MTgene = NULL) {
    if (!(is.numeric(MTfilter) && MTfilter >= 0 && MTfilter <= 1)) {
        stop("\"MTfilter\" should be a numeric value between 0 and 1.")
    }
    
    if(is.list(CBout)){
        xmat <- cbind(CBout$cluster_matrix, CBout$cell_matrix)
    }else{
        xmat <- CBout
    }

    if (is.null(MTgene)) {
        MTgene <- c(grep(pattern = "\\<MT-", x = rownames(xmat)),
                    grep(pattern = "\\<mt-", x = rownames(xmat)))
    } else if (!all(is.character(MTgene))) {
        stop("\"MTgene\" should be a character vector of row names.")
    }
    if (length(MTgene) > 0) {
        MTprop <- colSums(xmat[MTgene,]) / colSums(xmat)
        brokenCell <- MTprop > MTfilter
        xmat <- xmat[,!brokenCell]
    }
    return(xmat)
}