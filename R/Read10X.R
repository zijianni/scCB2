#' Read 10X output data
#'
#' One-line simple function to read 10X output data. Works under both old 
#' (<3) and new (>=3) Cell Ranger version.
#' 
#' @param dir The directory of 10X output data. For Cell Ranger version <3, 
#' directory should include three files: barcodes.tsv, genes.tsv, matrix.mtx. 
#' For Cell Ranger version >=3, directory should include three 
#' files: barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz.
#' @param row.name Specify either using gene names (\code{row.name = "name"}) or 
#' gene Ensembl IDs (\code{row.name = "id"}) as row names of the count matrix. 
#' Default is \code{row.name = "name"}.
#' @param meta Logical. If \code{TRUE}, returns a list containing both the count matrix
#' and metadata of genes (features). Metadata includes feature names, IDs and other 
#' additional information depending on Cell Ranger output. If \code{FALSE} (default),
#' only returns the count matrix.
#' 
#' @return If \code{meta = T}, returns a list of two elements: a "dgCMatrix" 
#' sparse matrix containing expression counts and a data frame containing metadata of 
#' genes (features). For the count matrix, each row is a gene (feature) and 
#' each column is a barcode.  If \code{meta = F}, only returns the count matrix.
#' 
#' @example 
#' \dontrun{
#' data_dir <- 'path/to/data/directory'
#' 
#' #count matrix
#' data_mat <- Read10X(data.dir = data_dir)
#' str(data_mat)
#' 
#' #a list including count matrix and feature metadata
#' data_list <- Read10X(data.dir = data_dir, meta = T)
#' str(data_list)
#' }
#' 
#' @importFrom utils read.delim
#' @importFrom Matrix readMM
#' @importFrom Matrix as
#' 
#' @export

Read10X <- function(dir,
                    row.name = "name",
                    meta = F) {
  if (!row.name %in% c("name", "id"))
    stop("row.name should be either \"name\" or \"id\".")
  if (!dir.exists(dir))
    stop("Directory does not exist.")
  dir <- gsub("/$", "", dir)
  fname <- list.files(dir)
  #check if Cell Ranger version >=3
  V3 <- "features.tsv.gz" %in% fname
  if (V3) {
    Barcode <- file.path(dir, "barcodes.tsv.gz")
    Gene <- file.path(dir, "features.tsv.gz")
    CountMat <- file.path(dir, "matrix.mtx.gz")
  } else{
    Barcode <- file.path(dir, "barcodes.tsv")
    Gene <- file.path(dir, "genes.tsv")
    CountMat <- file.path(dir, "matrix.mtx")
  }
  barcode <- readLines(Barcode)
  gene.meta <- read.delim(Gene, header = F, colClasses = "character")
  if (V3) {
    colnames(gene.meta) <- c("id", "name", "type")
  } else{
    colnames(gene.meta) <- c("id", "name")
  }
  if (row.name == "name") {
    #gene names as row names of count matrix
    gene <- gene.meta$name
  } else{
    #gene ids as row names of count matrix
    gene <- gene.meta$id
  }
  countmat <- as(readMM(CountMat), "dgCMatrix")
  colnames(countmat) <- barcode
  rownames(countmat) <- make.unique(gene)
  if (meta) {
    #return a list including count matrix and gene metadata
    return(list(CountMatrix = countmat, Metadata = gene.meta))
  } else{
    #only return count matrix
    return(countmat)
  }
}
