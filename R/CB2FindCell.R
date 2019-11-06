#' Distinguish real cells from empty droplets using clustering-based 
#' Monte-Carlo test.
#'
#' The main function of \code{scCB2} package. Distinguish real cells 
#' from empty droplets using clustering-based Monte-Carlo test. 
#'
#' @param RawDat Matrix. Supports standard matrix or sparse matrix. 
#' This is the raw feature-by-barcode count matrix.
#' 
#' @param FDR_threshold Numeric between 0 and 1. Default: 0.01. 
#' The False Discovery Rate (FDR) to be controlled for multiple testing.
#' 
#' @param lower Positive integer. Default: 100. All barcodes 
#' whose total count below or equal to this threshold are defined as 
#' background empty droplets. They will be used to estimate the background 
#' distribution. The remaining barcodes will be test against background 
#' distribution. If sequencing depth is deliberately made higher (lower)
#' than usual, this threshold can be leveled up (down) correspondingly to 
#' get reasonable number of cells. Recommended sequencing depth for this 
#' default threshold: 40,000~80,000 reads per cell.
#' 
#' @param upper Positive numeric. Default: \code{NULL}. This is the upper 
#' threshold for large barcodes. All barcodes whose total counts are larger 
#' or equal to upper threshold are directly classified as real cells prior 
#' to testing. If \code{upper = NULL}, the knee point of the log rank curve 
#' of barcodes total counts will serve as the upper threshold, which is 
#' calculated using package \code{DropletUtils}'s method. If 
#' \code{upper = Inf}, no barcodes will be retained prior to testing. 
#' If manually specified, it should be greater than pooling threshold. 
#' 
#' @param RemoveProtein Logical. Default: \code{TRUE}. For 10x Cell Ranger 
#' version >=3, extra features (surface proteins) besides genes 
#' are measured simultaneously. If \code{RemoveProtein = TRUE}, only genes 
#' are used for testing. Removing extra features are recommended
#' because the default pooling threshold (100) is chosen only for handling 
#' gene expression. Protein expression level is hugely different
#' from gene expression level. If using the default pooling threshold 
#' while keeping proteins, the estimated background distribution
#' will be hugely biased and does not reflect the real background distribution 
#' of empty droplets.   
#' 
#' @param Ncores Positive integer. Default: \code{detectCores() - 2}. 
#' Number of cores for parallel computation.
#' 
#' @param verbose Logical. Default: \code{TRUE}. If \code{verbose = TRUE}, 
#' progressing messages will be printed.
#'
#' @return A list of (1) real cell barcode matrix distinguished during 
#' cluster-level test, (2) real cell barcode matrix 
#' distinguished during single-barcode-level 
#' test plus large cells who exceed upper threshold, (3) testing statistics 
#' (Pearson correlation to the background) for all candidate barcode clusters, 
#' (4) barcode IDs for all candidate barcode clusters, the name of each cluster 
#' is its median barcode size, (5) testing statistics 
#' (log likelihood under background distribution) for remaining single 
#' barcodes not clustered, (6) background distribution count vector,
#' without Good-Turing correction.
#' 
#' @details 
#' 
#' Input data is a feature-by-barcode matrix. Background barcodes are 
#' defined based on \code{lower}. Large barcodes are 
#' automatically treated as real cells based on \code{upper}. Remaining 
#' barcodes will be first clustered into subgroups, then 
#' tested against background using Monte-Carlo p-values simulated from 
#' Multinomial distribution. The rest barcodes will be further tested 
#' using EmptyDrops (Aaron T. L. Lun \emph{et. al. 2019}).
#' FDR is controlled based on \code{FDR_threshold}.
#' 
#' This function supports parallel computation. \code{Ncores} is used to specify
#' number of cores. 
#' 
#' Under CellRanger version >=3, extra features other than genes are 
#' simultaneously measured (e.g. surface protein). We recommend filtering 
#' them out using \code{RemoveProtein = TRUE} because the measurement of 
#' protein abundance is not in the same level as gene expression counts.
#' If using the default pooling threshold while keeping proteins, the 
#' estimated background distribution will be hugely biased and does not 
#' reflect the real background distribution of empty droplets. The resulting 
#' matrix will contain lots of barcodes who have almost zero gene expression
#' and relatively high protein expression, which are usually not useful for 
#' RNA-Seq study.
#' 
#' @examples
#' # raw data, all barcodes
#' data(mbrainSub)
#' str(mbrainSub)
#' 
#' # run CB2
#' CBOut <- CB2FindCell(mbrainSub, FDR_threshold = 0.01, 
#'     lower = 100, Ncores = 2)
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
                        lower = 100,
                        upper = NULL,
                        RemoveProtein = TRUE,
                        Ncores = detectCores() - 2,
                        verbose = TRUE) {
    time_begin <- Sys.time()
    
    if(!(is(RawDat,"matrix") | is(RawDat,"dgCMatrix")))
        stop("Incorrect raw data format.")
    
    if (any(duplicated(colnames(RawDat))))
        stop("Detected duplicated barcode ID.")
    if (any(duplicated(rownames(RawDat))))
        stop("Detected duplicated gene ID.")
    
    if (verbose) {
        message("FDR threshold: ", FDR_threshold)
        message("Lower threshold: ", lower)
        message("Cores allocated: ", Ncores, "\n")
    }
    
    if (Ncores < 1)
        Ncores <- 1
    
    #############filter proteins, 0-expressed genes and barcodes
    if (verbose)
        message("(1/5) Filtering empty barcodes and features in raw data...")
    if (is(RawDat, "matrix")) {
        RawDat <- as(RawDat, "dgCMatrix")
    }
    
    if (RemoveProtein) {
        protein <- grep(pattern = "TotalSeqB", x = rownames(RawDat))
    }else{
        protein <- integer(0)
    }
    
    if (length(protein) > 0) {
        dat_filter <- FilterGB(RawDat[-protein,],0,0)
    }else{
        dat_filter <- FilterGB(RawDat,0,0) 
    }

    if (is.null(upper)) {
        
        brank <- Calc_upper(dat_filter, lower = lower)
        #check convergence of knee point
        repeat{
            upper_temp <- brank$knee
            brank <- Calc_upper(dat_filter, brank$inflection + 100)
            if(brank$knee==upper_temp) break
        }
        
        if (is.null(upper_temp)) {
            stop("Probably not enough barcodes to calculate knee point.")
        }
        upper <- upper_temp
    }
    
    if (verbose) {
        message("Upper threshold: ", upper)
    }
    if(upper <= lower){
        stop("Upper threshold should be larger than lower threshold.")
    }
    
    #####Large real cells
    bc <- colSums(dat_filter)
    if (any(bc >= upper)) {
        upper_cell <- names(bc)[bc >= upper]
        upper_mat <- RawDat[, upper_cell]
        dat_filter <- FilterGB(dat_filter[,bc < upper], 0, 0)
    }else{
        upper_mat <- NULL
    }
    
    #####remaining barcodes
    B0_bc <- names(bc)[bc <= lower]
    B1_bc <- names(bc)[(bc > lower & bc < upper)]
    B0 <- dat_filter[, B0_bc]
    dat <- dat_filter[, B1_bc]
    
    null_count <- rowSums(B0)
    null_prob <- goodTuringProportions(null_count)[,1]
    c_prob <- null_prob
    
    if (any(bc >= upper)) {
        upper_count <- rowSums(upper_mat[rownames(B0),])
        upper_prob <- goodTuringProportions(upper_count)[,1]
        
        if(c_entropy(null_prob)<=c_entropy(upper_prob)){
            c_prob <- upper_prob
        }
    }
    
    #function for calculating test statistic
    stat_fun <- function(x) {
        cor(x, null_prob)
    }
    
    if (verbose)
        message("Done.\n")
    
    ####################estimate test statistic for all barcodes
    
    if (verbose)
        message("(2/5) Calculating test statistics for barcodes...")
    
    dat_Cor <-
        Calc_stat(dat,
                size = 1000,
                Ncores = Ncores,
                null_prob = null_prob)
    dat_temp <- data.frame(cbind(dat_Cor, colSums(dat)))
    #Here the test statistic is correlation, but 
    #we name it as logLH for step 5 use.
    if (verbose)
    colnames(dat_temp) <- c("logLH", "count") 
        message("Done.\n")
    
    ###############################construct clusters
    clust_mat <- NULL
    if (verbose) {
        message("(3/5) Constructing highly-correlated clusters...")
    }
    
    output_cl_raw <- Calc_cluster_Paral(
        dat,
        size_hc = 1000,
        Ncores = Ncores,
        dat_cor = dat_Cor,
        cor_clust, ave_cor, size_cor, 
        cfun, Calc_cluster, sparse_cor
    )
    
    #automatically calculate cluster cutoff
    c_threshold <- numeric(10)
    for(bg in seq_len(10)){
        c_mat <- cor(rmultinom(100, bg * 100, c_prob))
        diag(c_mat) <- NA
        c_threshold[bg] <- mean(c_mat, na.rm = TRUE)
    }
    if(verbose){
        message("Baseline clustering threshold: ",round(c_threshold[2],3))
    }
    c_size_int <- as.numeric(names(output_cl_raw$stat))%/%100
    c_size_int <- ifelse(c_size_int>10,9,c_size_int)
    c_size_int <- ifelse(c_size_int<3,2,c_size_int)
    
    c_within <- as.numeric(names(output_cl_raw$barcode))

    #Filter out bad clusters. Only keep well-grouped clusters
    good_clust <- which(c_within>c_threshold[c_size_int])
    output_cl <- list(barcode=NULL,stat=NULL)
    for(gc in good_clust){
        output_cl$barcode <- c(output_cl$barcode,
                            list(output_cl_raw$barcode[[gc]]))
        output_cl$stat <- c(output_cl$stat,output_cl_raw$stat[gc])
    }
        
    if(is.null(output_cl$barcode)){
        warning("Failed to construct any cluster. Skipping to Step 5.")
        if(verbose){
            message("Raw data may not have enough barcodes for clustering.\n")
            message("Also check for and remove outlier overexpressed gene.\n")
        }
        cl_temp <- NULL
    } else{
        if (verbose)
            message("Done.\n")
        ##############estimate test statistic for every cluster
        cl_Cor <- output_cl$stat
        
        cl_temp <- data.frame(cbind(cl_Cor, names(cl_Cor)))
        
        colnames(cl_temp) <- c("corr", "count")
        cl_temp$corr <- as.numeric(as.character(cl_temp$corr))
        cl_temp$count <- as.numeric(as.character(cl_temp$count))
        
        clust_b <- unlist(output_cl$barcode)
        ##############################cluster test
        if (verbose)
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
        
        if (verbose)
            message("Done.\n")
        
    }
    ##########Permute null distribution of the test statistic
    if (verbose)
        message("(5/5) Calculating empirical p-value for each barcode...")
    
    if ((ncol(dat)) == 0) {
        if (verbose)
            message("\nNo remaining individual barcodes.
                    \nAll finished.\n")
        if (verbose)
            print(Sys.time() - time_begin)
        return(
            list(
                cluster_matrix = clust_mat,
                cell_matrix = upper_mat,
                ClusterStat = cl_temp,
                Cluster = output_cl$barcode,
                background = null_count
            )
        )

    } else{
        cand_barcode <- c(colnames(cbind(dat,B0)), 
                        colnames(upper_mat))
        ED_out <- emptyDrops(RawDat[, cand_barcode], 
                    lower = lower, retain = upper)
        dat_temp$logLH <- ED_out$LogProb[seq_len(ncol(dat))]
        dat_temp$pval <- ED_out$PValue[seq_len(ncol(dat))]
        dat_temp$padj <- ED_out$FDR[seq_len(ncol(dat))]
        cell_barcode <-
            cand_barcode[ifelse(is.na(ED_out$FDR), FALSE,
                                ED_out$FDR <= FDR_threshold)]
        cell_mat <- RawDat[, cell_barcode]
        if (verbose){
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

#Calculate test statistc for each barcode in parallel 
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

#Calculate test statistc for each cluster in parallel 
Calc_cluster_Paral <- function(dat,
                                size_hc = 1000,
                                Ncores = detectCores() - 2,
                                dat_cor,
                                cor_clust, ave_cor, size_cor, 
                                cfun, Calc_cluster, sparse_cor
                                ) {
    bc <- colSums(dat)
    dat <- dat[, names(sort(bc))]
    
    #####parallel computation
    apply_ind <-
        split(seq_len(ncol(dat)), ceiling(seq_len(ncol(dat)) / size_hc))
    
    #if the last group is too small, combine with the second last.
    if( (length(apply_ind[[length(apply_ind)]]) < size_hc / 2) &&
        (length(apply_ind) > 1) ){
        apply_ind[[length(apply_ind)-1]] <- 
            c(apply_ind[[length(apply_ind)-1]], 
            apply_ind[[length(apply_ind)]])
        apply_ind[[length(apply_ind)]] <- NULL
    }
    
    cl <- makeCluster(Ncores)
    registerDoParallel(cl)
    dat_tmp <- list()
    for (i in seq_along(apply_ind)) {
        dat_tmp[[i]] <- dat[, apply_ind[[i]], drop = FALSE]
    }
    
    dat_tmp_it <- iter(dat_tmp)
    

    
    x <- NULL #initialize x to prevent note in R CMD check
    Output <- foreach(x = dat_tmp_it, .combine = "cfun", 
        .packages = c("Matrix")) %dopar% {

        Calc_cluster(x,
                    dat_cor = dat_cor,
                    dat_bc = bc,
                    sparse_cor,
                    cor_clust,
                    ave_cor,
                    size_cor)
    }
    
    stopCluster(cl)
    return(Output)
}

#Construct clusters and calculate test statistics for barcode subsets
Calc_cluster <-
    function(dat,
            dat_cor,
            dat_bc,
            sparse_cor,
            cor_clust,
            ave_cor,
            size_cor
    ) {
        
        if(ncol(dat)<=20) return(NULL)
        Nclust <- ncol(dat)/20
        Output_stat <- c() #save test statistic for each cluster
        cor_temp <- sparse_cor(dat) #first calculate pairwise correlation
        clust_temp <- #then hierarchical clustering 
            cor_clust(cor_temp, Nclust = min(Nclust, ncol(dat)))
        cand_bclust <- lapply(clust_temp,colnames)
        
        Output_stat <- numeric(Nclust)
        for (i in seq_len(Nclust)) {
            Output_stat[i] <- median(dat_cor[cand_bclust[[i]]])
            WhichMedian <-
                which.min(abs(dat_cor[cand_bclust[[i]]] - 
                                Output_stat[i]))
            names(Output_stat)[i] <-  dat_bc[cand_bclust[[i]]][WhichMedian]
        }
        
        names(cand_bclust) <- ave_cor(clust_temp,size_cor)
        
        return(list(barcode = cand_bclust, stat = Output_stat))
    }

#cluster barcodes based on correlation and return a list
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

#return size of each cluster in a list
size_cor <- function(cor_list) {
    unlist(lapply(cor_list,
            function(x)
                ifelse(is.null(dim(x)), 1, dim(x)[1])))
}

#return averaged correlation of each cluster in a list
ave_cor <- function(cor_list,size_cor) {
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

#efficient correlation calculation of large sparse matrix
sparse_cor <- function(x) {
    
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

#combine results in foreach parallel computation
cfun <- function(a, b) {
    barcode <- c(a$barcode, b$barcode)
    stat <- c(a$stat, b$stat)
    return(list(barcode = barcode, stat = stat))
}

#entropy
c_entropy <- function(prob){
    prob[prob==0] <- 1
    -sum(prob*log(prob),na.rm = TRUE)
}

#' @importFrom stats smooth.spline
#' @importFrom stats predict

#Calculate knee point of a dataset
#Note: This function is a modified version of barcodeRanks() in 
#package DropletUtils. We fixed a minor bug causing returned knee point
#being larger than actual knee point. We also moved the smooth spline fitting
#to the beginning to avoid unstable knee point estimation when lower threshold
#changes. For its original script, see 
#https://github.com/MarioniLab/DropletUtils/blob/master/R/barcodeRanks.R
Calc_upper <- function(dat,lower){
    dat <- FilterGB(dat)
    totals <- unname(colSums(dat))
    o <- order(totals, decreasing=TRUE)
    
    stuff <- rle(totals[o])
    # Get mid-rank of each run.
    run.rank <- cumsum(stuff$lengths) - (stuff$lengths-1)/2 
    run.totals <- stuff$values
    
    keep <- run.totals > lower
    if (sum(keep)<3) { 
        stop("insufficient unique points for computing knee/inflection points")
    }

    x <- log10(run.rank[keep])
    fit <- smooth.spline(log10(run.rank), log10(run.totals), df=20)
    y <- predict(fit)$y[keep]
    # Numerical differentiation to identify bounds for spline fitting.
    # The upper/lower bounds are defined at the 
    # plateau and inflection, respectively.
    d1n <- diff(y)/diff(x)
    right.edge <- which.min(d1n)
    left.edge <- which.max(d1n[seq_len(right.edge)])
    # We restrict to this region, thereby simplifying the shape of the curve.
    # This allows us to get a decent fit with low df for stable differentiation.
    new.keep <- left.edge:right.edge
    # Smoothing to avoid error multiplication upon differentiation.
    # Minimizing the signed curvature and returning 
    # the total for the knee point.
    d1 <- predict(fit, deriv=1)$y[keep][new.keep]
    d2 <- predict(fit, deriv=2)$y[keep][new.keep]
    curvature <- d2/(1 + d1^2)^1.5
    knee <- 10^(y[new.keep][which.min(curvature)])
    inflection <- 10^(y[right.edge])
    return(list(knee=round(knee),inflection=round(inflection)))
}
