#' All-in-one function from raw data to filtered cell matrix
#'
#' All-in-one function for \code{scCB2}. Take 10x output raw data as input 
#' and return either a matrix of real cells identified by CB2 or 
#' a Seurat object containing this matrix, which can be incorporated 
#' with downstream analysis using Seurat pipeline.  
#'
#' @param dir The directory of 10x output data. For Cell Ranger version <3,
#' directory should include three files: barcodes.tsv, genes.tsv, matrix.mtx.
#' For Cell Ranger version >=3, directory should include three
#' files: barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz.
#' 
#' @param h5file The path of 10x output HDF5 file (ended with .h5).
#' 
#' @param FDR_threshold Numeric between 0 and 1. Default: 0.01. 
#' The False Discovery Rate (FDR) to be controlled for multiple testing.
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
#' @param AsSeurat Logical. Default: \code{FALSE}. Decides whether a Seurat
#' object is returned instead of cell matrix. Set to \code{TRUE} so that 
#' users can directly apply Seurat pipeline for downstream analyses.
#' 
#' @param Ncores Positive integer. Default: \code{detectCores() - 2}. 
#' Number of cores for parallel computation.
#' 
#' @param ... Additional arguments to be passed to \code{CB2FindCell}.
#'
#' @return Either a sparse matrix of real cells identified by CB2 or a
#' Seurat object containing real cell matrix.
#' 
#' @details 
#' 
#' \code{QuickCB2} is a quick function to apply CB2 on 10x Cell Ranger 
#' raw data by combining \code{Read10xRaw}, \code{Read10xRawH5}, 
#' \code{CB2FindCell} and \code{GetCellMat} into one simple function 
#' under default parameters. 
#' 
#' 
#' @examples
#' 
#' # Since this function is a combination of others, to avoid duplication
#' # please refer to the help pages of other functions 
#' # for runnable examples. 
#' 
#' \dontrun{
#' 
#' # Run CB2 on 10x raw data and get cell matrix.
#' # Control FDR at 0.1%. Use 4-core parallel computation.
#' 
#' RealCell <- QuickCB2(dir = "path/to/raw/data/directory", 
#'                      FDR_threshold = 0.001,
#'                      Ncores = 4)
#' str(RealCell)
#' }
#' 
#' @import BiocStyle
#' @importFrom Seurat CreateSeuratObject
#' @importFrom parallel detectCores
#' @export

QuickCB2 <- function(dir = NULL, 
                    h5file = NULL, 
                    FDR_threshold = 0.01,
                    MTfilter = 1,
                    MTgene = NULL,
                    AsSeurat = FALSE, 
                    Ncores = detectCores() - 2,
                    ...){

    message("Loading data...")
    if(is.null(dir)){
        RawDat <- Read10xRawH5(h5file = h5file)    
    }else if(is.null(h5file)){
        RawDat <- Read10xRaw(dir = dir)
    }else{
        stop("dir and h5file can't be specified together.")
    }
    message("Done.\n")
    
    message("Running CB2...")
    CBOut <- CB2FindCell(RawDat = RawDat,
                        FDR_threshold = FDR_threshold,
                        Ncores = Ncores,
                        ...)
    message("Done.\n")
    
    message("Generating final cell matrix...")
    CBMat <- GetCellMat(CBOut, 
                        MTfilter = MTfilter, 
                        MTgene = MTgene)
    message("Done.\n")
    
    if(AsSeurat){
        return(CreateSeuratObject(counts = CBMat, 
                                project = "CB2output"))
    }else{
        return(CBMat)
    }
    
}


