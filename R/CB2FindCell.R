#' Distinguish real cells from empty droplets using clustering-based Monte-Carlo test.
#'
#' The main function of \code{scCB2} package. Distinguish real cells from empty droplets 
#' using clustering-based Monte-Carlo test. 
#'
#' @param RawDat Matrix. Supports standard matrix or sparse matrix. This is the raw feature-by-barcode count matrix.
#' 
#' @param FDR_threshold Numeric between 0 and 1. Default: 0.01. The False Discovery Rate (FDR) to be controlled for multiple testing.
#' 
#' @param background_threshold Positive integer. Default: 100. All barcodes whose total count below or equal to this threshold
#' are defined as background empty droplets. They will be used to estimate the background distribution. The remaining barcodes
#' will be test against background distribution. If sequencing depth is deliberately made higher (lower)
#' than usual, this threshold can be leveled up (down) correspondingly to get reasonable number of cells.
#' Recommended sequencing depth for this default threshold: 40,000~80,000 reads per cell.
#' 
#' @param RemoveProtein Logical. Default: \code{TRUE}. For 10X Cell Ranger version >=3, extra features (surface proteins) besides genes 
#' are measured simultaneously. If \code{RemoveProtein = TRUE}, only genes are used for testing. Removing extra features are recommended
#' because the default pooling threshold (100) is chosen only for handling gene expression. Protein expression level is hugely different
#' from gene expression level. If using the default pooling threshold while keeping proteins, the estimated background distribution
#' will be hugely biased and does not reflect the real background distribution of empty droplets.   
#' 
#' @param retain Positive integer. Default: \code{NULL}. This is the retain threshold for large barcodes. All barcodes whose
#' total counts are larger or equal to retain threshold are directly classified as real cells prior to testing. If \code{retain = NULL}, the knee point
#' of the log rank curve of barcodes total counts will serve as the retain threshold, which is calculated using package
#' \code{DropletUtils}'s method. If \code{retain = Inf}, no barcodes will be retained prior to testing. If manually specified,
#' it should be greater than pooling threshold. 
#' 
#' @param Ncores Positive integer. Default: \code{detectCores() - 2}. Number of cores for parallel computation.
#' 
#' @param PrintProg Logical. Default: \code{TRUE}. If \code{PrintProg = TRUE}, progressing messages will be printed.
#'
#' @return A list of (1) real cell barcode matrix distinguished during cluster-level test, (2) real cell barcode matrix
#' distinguished during single-barcode-level test, (3) testing statistics (Pearson correlation to the background) 
#' for all candidate barcode clusters, (4) barcode IDs for all candidate barcode clusters, the name of each cluster is its median 
#' barcode size, (5) testing statistics (log likelihood under background distribution) for remaining single barcodes not clustered,
#' (6) estimated background distribution count vector.
#' 
#' @details 
#' 
#' Input data is a feature-by-barcode matrix. Background barcodes are defined based
#' on \code{background_threshold}. Large barcodes are automatically treated as real cells
#' based on \code{retain}. Remaining barcodes will be first clustered into subgroups, then 
#' tested against background using Monte-Carlo p-values simulated from Multinomial distribution.
#' The rest barcodes will be further tested using EmptyDrops (Aaron T. L. Lun \emph{et. al. 2019}).
#' FDR is controlled based on \code{FDR_threshold}.
#' 
#' This function supports parallel computation. \code{Ncores} is used to specify
#' number of cores. 
#' 
#' Under CellRanger version >=3, extra features other than genes are simultaneously measured 
#' (e.g. surface protein). We recommend filtering them out using \code{RemoveProtein = TRUE}
#' because the measurement of protein abundance is not in the same level as gene expression counts.
#' If using the default pooling threshold while keeping proteins, the estimated background distribution
#' will be hugely biased and does not reflect the real background distribution of empty droplets.
#' The resulting matrix will contain lots of barcodes who have almost zero gene expression
#' and relatively high protein expression, which are usually not useful for RNA-Seq study.
#' 
#' @examples
#' # raw data, all barcodes
#' data(mbrainSub)
#' str(mbrainSub)
#' 
#' # run CB2
#' CBOut <- CB2FindCell(mbrainSub, FDR_threshold = 0.01, 
#'     background_threshold = 100, Ncores = 2)
#' RealCell <- GetCellMat(CBOut, MTfilter = 0.05)
#' 
#' # real cells
#' str(RealCell)
#' 
#' @import doParallel
#' @import foreach
#' @import parallel
#' @import SingleCellExperiment
#' @import Matrix
#' @importFrom DropletUtils barcodeRanks
#' @importFrom DropletUtils emptyDrops
#' @importFrom edgeR goodTuringProportions
#' @importFrom methods as
#' @importFrom methods is
#' @importFrom stats cor
#' @importFrom stats median
#' @importFrom stats p.adjust
#' @importFrom stats rmultinom
#' @export

CB2FindCell <- function(RawDat,
                        FDR_threshold = 0.01,
                        background_threshold = 100,
                        RemoveProtein = TRUE,
                        retain = NULL,
                        Ncores = detectCores() - 2,
                        PrintProg = TRUE) {
    time_begin <- Sys.time()
    
    if (any(duplicated(colnames(RawDat))))
        stop("Detected duplicated barcode ID.")
    if (any(duplicated(rownames(RawDat))))
        stop("Detected duplicated gene ID.")
    
    if (PrintProg) {
        message("FDR threshold: ", FDR_threshold)
        message("Background threshold: ", background_threshold)
        message("Cores allocated: ", Ncores, "\n")
    }
    
    if (Ncores < 1)
        Ncores <- 1
    
    #############filter 0-expressed genes and barcodes
    if (PrintProg)
        message("\n(1/5) Filtering empty barcodes and features in raw data...")
    if (!is(RawDat, "dgCMatrix")) {
        RawDat <- as(RawDat, "dgCMatrix")
    }
    
    dat_filter <- FilterGB(RawDat)
    if (RemoveProtein) {
        protein <- grep(pattern = "TotalSeqB", x = rownames(dat_filter))
        if (length(protein) > 0) {
            dat_filter <- dat_filter[-protein,]
        }
    }
    
    if (is.null(retain)) {
        brank <- barcodeRanks(dat_filter, lower = background_threshold)
        
        knee <-
            as.integer(ifelse(is.null(brank$knee),
                            brank@metadata$knee, brank$knee))
        
        if (is.null(knee)) {
            stop("Failed to calculate knee point. Check input data.")
        }
        retain <- knee
    }
    
    if(retain <= background_threshold){
        stop("Retain threshold should be larger than background threshold.")
    }
    
    bc <- colSums(dat_filter)
    if (length(names(which(bc >= retain))) > 0) {
        retain_cell <- names(bc)[bc >= retain]
        retain_mat <- RawDat[, retain_cell]
    }
    B0_cell <- names(bc)[bc <= background_threshold]
    B1_cell <- names(bc)[(bc > background_threshold & bc < retain)]
    B0 <- dat_filter[, B0_cell]
    dat <- dat_filter[, B1_cell]
    
    nzero_gene_B1 <- which(rowSums(dat) > 0)
    nzero_gene_B0 <- which(rowSums(B0) > 0)
    dat_filter <-
        dat_filter[sort(union(nzero_gene_B0, nzero_gene_B1)),]
    B0 <- dat_filter[, B0_cell]
    
    null_count <- rowSums(B0)
    
    null_prob <- goodTuringProportions(null_count)
    dat <- dat_filter[, B1_cell]
    
    
    #automatically calculate cluster cutoff
    c_threshold <-
        mean(unlist(replicate(1000, cor(
            rmultinom(1, 2 * background_threshold, null_prob),
            rmultinom(1, 2 * background_threshold, null_prob)
        ))))
    
    stat_fun <- function(x) {
        cor(x, null_prob)
    }
    
    if (PrintProg)
        message("Done.\n")
    ##############################estimate test statistic for all barcodes
    
    if (PrintProg)
        message("(2/5) Calculating test statistics for barcodes...")
    
    dat_Cor <-
        Calc_stat(dat,
                size = 1000,
                Ncores = Ncores,
                null_prob = null_prob)
    dat_temp <- data.frame(cbind(dat_Cor, colSums(dat)))
    #Here the test statistic is correlation, but name it as logLH for step 5
    if (PrintProg)
    colnames(dat_temp) <- c("logLH", "count") 
        message("Done.\n")
    
    ###############################construct clusters
    clust_mat <- NULL
    if (PrintProg) {
        message("(3/5) Constructing highly-correlated clusters...")
        message("Cluster threshold: ", round(c_threshold, 3))
    }
    output_cl <- Calc_cluster_Paral(
        dat,
        size_hc = 2000,
        Ncores = Ncores,
        Nclust = 100,
        c_threshold = c_threshold,
        dat_cor = dat_Cor
    )
    if(is.null(output_cl$barcode)){
        if(PrintProg){
            message("Failed to construct any cluster. Skipping to Step 5.")
            message("Does raw data have too few barcodes to construct clusters?\n")
        }
        warning("No cluster constructed. Clustering tests are skipped.")
    } else{
        if (PrintProg)
            message("Done.\n")
        ##############################estimate test statistic for every cluster
        cl_Cor <- output_cl$stat
        
        cl_temp <- data.frame(cbind(cl_Cor, names(cl_Cor)))
        
        colnames(cl_temp) <- c("corr", "count")
        cl_temp$corr <- as.numeric(as.character(cl_temp$corr))
        cl_temp$count <- as.numeric(as.character(cl_temp$count))
        
        clust_b <- unlist(output_cl$barcode)
        ##############################cluster test
        if (PrintProg)
            message("(4/5) Calculating empirical p-value for each cluster...")
        
        cand_count <- sort(unique(cl_temp$count))
        cl <- makeCluster(Ncores)
        registerDoParallel(cl)
        cfun <- function(a, b) {
            sim <- c(a$sim, b$sim)
            pval <- c(a$pval, b$pval)
            return(list(sim = sim, pval = pval))
        }
        
        output_ <-
            foreach(
                i = seq_along(cand_count),
                .combine = "cfun",
                .packages = c("Matrix","edgeR")
            ) %dopar% {
                input_stat = cl_Cor[cl_temp$count == cand_count[i]]
                v_temp <- rmultinom(1000, cand_count[i], null_prob)
                e_temp <- apply(v_temp, 2, stat_fun)
                
                pval <- vapply(input_stat, function(x)
                        (sum(e_temp <= x) + 1) / 1001, 0)
                
                names(pval) <- rep(cand_count[i], length(pval))
                list(sim = c(cand_count[i]), pval = pval)
            }
        stopCluster(cl)
        
        #####################find significant clusters
        for (i in seq_along(cand_count)) {
            cl_temp$pval[cl_temp$count == cand_count[i]] <-
                output_$pval[names(output_$pval) == cand_count[i]]
        }
        
        cl_temp$padj <- p.adjust(cl_temp$pval, method = "BH")
        
        rm(output_)
        
        ##get real clusters
        cell_cluster <- which(cl_temp$padj <= FDR_threshold)
        
        res_b <- c()
        for (j in cell_cluster) {
            res_b <- c(res_b, output_cl$barcode[[j]])
        }
        
        ##check whether significant cluster exists
        if (!is.null(res_b)) {
            clust_mat <- RawDat[, res_b]
            
            dat_Cor <- dat_Cor[setdiff(colnames(dat), res_b)]
            
            dat <- dat[, setdiff(colnames(dat), res_b)]
            dat_temp <- dat_temp[setdiff(colnames(dat), res_b),]
        }
        
        if (PrintProg)
            message("Done.\n")
        
    }
    ##########Permute null distribution of the test statistic
    if (PrintProg)
        message("(5/5) Calculating empirical p-value for each barcode...")
    
    if ((ncol(dat)) == 0) {
        if (PrintProg)
            message("\nNo remaining individual barcodes.\nAll finished.\n")
        if (PrintProg)
            print(Sys.time() - time_begin)
        return(
            list(
                cluster_matrix = clust_mat,
                cell_matrix = retain_mat,
                ClusterStat = cl_temp,
                Cluster = output_cl$barcode,
                background = null_count
            )
        )
        
    } else{
        cand_barcode <- c(colnames(cbind(dat,B0)), colnames(retain_mat))
        ED_out <- emptyDrops(RawDat[, cand_barcode], 
                    lower = background_threshold, retain = retain)
        dat_temp$logLH <- ED_out$LogProb[seq_len(ncol(dat))]
        dat_temp$pval <- ED_out$PValue[seq_len(ncol(dat))]
        dat_temp$padj <- ED_out$FDR[seq_len(ncol(dat))]
        cell_barcode <-
            cand_barcode[ifelse(is.na(ED_out$FDR), FALSE,
                                ED_out$FDR <= FDR_threshold)]
        cell_mat <- RawDat[, cell_barcode]
        if (PrintProg){
            message("Done.\n")
            message("All finished.\n")
            print(Sys.time() - time_begin)
        }
           
        return(
            list(
                cluster_matrix = clust_mat,
                cell_matrix = cell_mat,
                ClusterStat = cl_temp,
                Cluster = output_cl$barcode,
                BarcodeStat = dat_temp,
                background = null_count
            )
        )
    }
}

#' @importFrom iterators iter
Calc_stat <-
    function(dat,
            size = 1000,
            Ncores = detectCores() - 2,
            null_prob) {
        stat_fun <- function(x) {
            cor(x, null_prob)
        }
        #####parallel computation
        apply_ind <-
            split(seq_len(ncol(dat)), ceiling(seq_len(ncol(dat)) / size))
        
        cl <- makeCluster(Ncores)
        registerDoParallel(cl)
        dat_tmp <- list()
        for (i in seq_along(apply_ind)) {
            dat_tmp[[i]] <- dat[, apply_ind[[i]], drop = FALSE]
        }
        
        dat_tmp_it <- iter(dat_tmp)
        x <- NULL #initialize x to prevent note in R CMD check
        Cor <-
            foreach(
                x = dat_tmp_it,
                .combine = 'c',
                .packages = c("Matrix", "edgeR")
            ) %dopar% {
                apply(x, 2, stat_fun)
            }
        stopCluster(cl)
        return(Cor)
    }







#' @importFrom stats as.dist
#' @importFrom stats hclust
#' @importFrom stats cutree

Calc_cluster_Paral <- function(dat,
                                size_hc = 2000,
                                Ncores = detectCores() - 2,
                                Nclust = 100,
                                c_threshold,
                                dat_cor) {
    bc <- colSums(dat)
    dat <- dat[, names(sort(bc))]
    
    #####parallel computation
    apply_ind <-
        split(seq_len(ncol(dat)), ceiling(seq_len(ncol(dat)) / size_hc))
    
    cl <- makeCluster(Ncores)
    registerDoParallel(cl)
    dat_tmp <- list()
    for (i in seq_along(apply_ind)) {
        dat_tmp[[i]] <- dat[, apply_ind[[i]], drop = FALSE]
    }
    
    dat_tmp_it <- iter(dat_tmp)
    
    cfun <- function(a, b) {
        barcode <- c(a$barcode, b$barcode)
        stat <- c(a$stat, b$stat)
        return(list(barcode = barcode, stat = stat))
    }
    
    x <- NULL #initialize x to prevent note in R CMD check
    Output <- foreach(x = dat_tmp_it, .combine = "cfun", 
        .packages = c("Matrix")) %dopar% {
        ave_cor <- function(cor_list) {
            #remove diagonal 1, remove single-point clusters
            mean_vec <- c()
            size_clust <- size_cor(cor_list)
            for (i in seq_along(size_clust)) {
                if (size_clust[i] == 1) {
            #single-point clusters will be filtered out
            #by correlation threshold
                    mean_vec <-
                        c(mean_vec, 0)
                } else{
                    cor_temp <- cor_list[[i]]
                    diag(cor_temp) <- NA
                    mean_vec <- c(mean_vec, mean(cor_temp, na.rm = TRUE))
                }
            }
            return(mean_vec)
        }
        
        #Note: I know it looks stupid to define functions inside functions, but 
        #it seems to be the only way I can make the package `foreach` work
        #correctly. The `.export` parameter in `foreach` function works in R
        #scripts, but does not work inside my package. I would appreciate any
        #suggestions for fixing it.
        
        #return size of each cluster in a list
        size_cor <- function(cor_list) {
            unlist(lapply(cor_list,
                        function(x)
                            ifelse(is.null(dim(x)), 1, dim(x)[1])))
        }
        
        #cluster barcodes based on correlation
        cor_clust <- function(cor_mat, Nclust = 100) {
            dist_mat <- as.dist(1 - cor_mat)
            hc <- hclust(dist_mat)
            hc1 <- cutree(hc, k = Nclust)
            cor_list <- list()
            for (i in seq_len(Nclust)) {
                cor_list[[i]] <- cor_mat[which(hc1 == i), which(hc1 == i)]
                if (length(cor_list[[i]]) == 1) {
                    cor_list[[i]] <- as.matrix(cor_list[[i]])
                    colnames(cor_list[[i]]) <-
                        colnames(cor_mat)[which(hc1 == i)]
                    rownames(cor_list[[i]]) <- colnames(cor_list[[i]])
                }
            }
            return(cor_list)
        }
        
        Calc_cluster <-
            function(dat,
                    Nclust = 100,
                    c_threshold,
                    dat_cor,
                    dat_bc) {
                sparse_cor <- function(x) {
                    #efficient correlation calculation of large sparse matrix
                    n <- nrow(x)
                    m <- ncol(x)
                    # non-empty rows
                    ii <- unique(x@i) + 1
                    Ex <- colMeans(x)
                    #centralize
                    nozero <- as.vector(x[ii, ]) - rep(Ex, each = length(ii))
                    covmat <- (crossprod(matrix(nozero, ncol = m)) +
                                crossprod(t(Ex)) * (n - length(ii))) / (n - 1)
                    sdvec <- sqrt(diag(covmat))
                    return(covmat / crossprod(t(sdvec)))
                }
                
                Output_stat <- c()
                cor_temp <- sparse_cor(dat)
                clust_temp <-
                    cor_clust(cor_temp, Nclust = min(Nclust, ncol(dat)))
                cand_clust <- #highly correlated, more than 1 barcode clusters
                    which(ave_cor(clust_temp) > c_threshold) 
                #Under some extreme distribution, no cluster exceeds threshold.
                #Lower threshold by 10%
                if (length(cand_clust) == 0) { 
                    cand_clust <-
                        which(ave_cor(clust_temp) > c_threshold*0.9)  
                }
                
                #If still no cluster, then keep this result
                if (length(cand_clust) == 0) {
                    return(NULL)
                } else{
                    cand_bclust <- list()
                    Output_stat <- numeric(length(cand_clust))
                    for (i in seq_along(cand_clust)) {
                        cand_bclust[[i]] <- colnames(
                            clust_temp[[cand_clust[i]]])
                    }
                    for (i in seq_along(cand_bclust)) {
                        Output_stat[i] <- median(dat_cor[cand_bclust[[i]]])
                        WM <-
                            which.min(abs(dat_cor[cand_bclust[[i]]] - 
                                        Output_stat[i]))
                        names(Output_stat)[i] <-  dat_bc[cand_bclust[[i]]][WM]
                    }
                    names(cand_bclust) <- names(Output_stat)
                    return(list(barcode = cand_bclust, stat = Output_stat))
                }
            }
        Calc_cluster(x,
                    Nclust,
                    c_threshold,
                    dat_cor = dat_cor,
                    dat_bc = bc)
    }
    
    stopCluster(cl)
    return(Output)
}
