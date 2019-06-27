#' Read 10X output data in hdf5 format
#'
#' One-line simple function to read 10X output hdf5 file. Works under both old 
#' (<3) and new (>=3) CellRanger version.
#' 
#' @param h5file The path of 10X output hdf5 file (ended with .h5). 
#' @param row.name Specify either using gene names (\code{row.name = "name"}) or 
#' gene Ensembl IDs (\code{row.name = "id"}) as row names of the count matrix. 
#' Default is \code{row.name = "name"}.
#' @param meta Logical. If \code{TRUE}, returns a list containing both the count matrix
#' and metadata of genes (features). Metadata includes feature names, IDs and other 
#' additional information depending on CellRanger output. If \code{FALSE} (default),
#' only returns the count matrix.
#' 
#' @return If \code{meta = T}, returns a list of two elements: a "dgCMatrix" 
#' sparse matrix containing expression counts and a data frame containing metadata of 
#' genes (features). For the count matrix, each row is a gene (feature) and 
#' each column is a barcode.  If \code{meta = F}, only returns the count matrix.
#' 
#' @example 
#' \dontrun{
#' data_path <- 'path/to/data/directory/data.h5'
#' 
#' #count matrix
#' data_mat <- Read10X(h5file = data_path)
#' str(data_mat)
#' 
#' #a list including count matrix and feature metadata
#' data_list <- Read10X(h5file = data_path, meta = T)
#' str(data_list)
#' }
#' 
#' @importFrom rhdf5 h5ls
#' @importFrom rhdf5 h5read
#' @importFrom Matrix sparseMatrix
#' 
#' @export

Read10Xh5 <- function(h5file,
                      row.name = "name",
                      meta = F) {
  if (!row.name %in% c("name", "id"))
    stop("row.name should be either \"name\" or \"id\".")
  if (!file.exists(h5file))
    stop("File does not exist.")
  fname <- h5ls(h5file)
  #check if CellRanger version >=3
  V3 <- "/matrix" %in% fname$group
  if (V3) {
    data.temp <-  h5read(h5file, "/matrix")
    barcode <- data.temp$barcodes
    gene.meta <- data.temp$features
    if (row.name == "name") {
      #gene names as row names of count matrix
      gene <- gene.meta$name
    } else{
      #gene ids as row names of count matrix
      gene <- gene.meta$id
    }
    countmat <-
      sparseMatrix(
        i = data.temp$indices,
        p = data.temp$indptr,
        x = data.temp$data,
        index1 = F,
        dims = data.temp$shape
      )
    colnames(countmat) <- as.vector(barcode)
    rownames(countmat) <- make.unique(gene)
    
    if (meta) {
      #return a list including count matrix and gene metadata
      return(list(CountMatrix = countmat, Metadata = gene.meta))
    } else{
      #only return count matrix
      return(countmat)
    }
  } else{
    subname <- fname$name[fname$group == "/"]
    #CellRanger V2 data might include multiple count matrices
    CountMatList <- vector("list", length(subname))
    names(CountMatList) <- subname
    for (i in seq_along(subname)) {
      data.temp <- h5read(h5file, paste0("/", subname[i]))
      barcode <- data.temp$barcodes
      gene.meta <-
        data.frame(name = data.temp$gene_names, id = data.temp$genes)
      if (row.name == "name") {
        #gene names as row names of count matrix
        gene <- gene.meta$name
      } else{
        #gene ids as row names of count matrix
        gene <- gene.meta$id
      }
      countmat <-
        sparseMatrix(
          i = data.temp$indices,
          p = data.temp$indptr,
          x = data.temp$data,
          dims = data.temp$shape,
          index1 = F,
          giveCsparse = T
        )
      colnames(countmat) <- as.vector(barcode)
      rownames(countmat) <- make.unique(gene)
      if (meta) {
        #return a list including count matrix and gene metadata
        CountMatList[[i]] <-
          list(CountMatrix = countmat, Metadata = gene.meta)
      } else{
        #only return count matrix
        CountMatList[[i]] <- countmat
      }
    }
    if (length(unique(fname$group)) > 2) {
      warning("Detected multiple datasets. Returning a list with multiple matrices.")
      return(CountMatList)
    } else{
      return(CountMatList[[1]])
    }
  }
}
