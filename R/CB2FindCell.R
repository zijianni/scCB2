CB2FindCell <- function(raw_dat,
                        FDR_threshold = 0.01,
                        pooling_threshold = 100,
                        remove_protein = T,
                        retain = NULL,
                        Ncores = detectCores() - 2,
                        print_prog = T) {
  time_begin <- Sys.time()
  
  if (any(duplicated(colnames(raw_dat))))
    stop("Detected duplicated barcode ID.")
  if (any(duplicated(rownames(raw_dat))))
    stop("Detected duplicated gene ID.")
  
  
  if (print_prog) {
    message("FDR threshold: ", FDR_threshold)
    message("Background threshold: ", pooling_threshold)
    message("Cores allocated: ", Ncores, "\n")
  }
  
  if (Ncores < 1)
    Ncores <- 1
  
  #############filter 0-expressed genes and barcodes
  if (print_prog)
    message("\n(1/5) Filtering empty barcodes and features in raw data...")
  if (!is(raw_dat, "dgCMatrix")) {
    raw_dat <- as(raw_dat, "dgCMatrix")
  }
  
  dat_filter <- Filter_GB(raw_dat)
  if (remove_protein) {
    protein <- grep(pattern = "TotalSeqB", x = rownames(dat_filter))
    if (length(protein) > 0) {
      dat_filter <- dat_filter[-protein,]
    }
  }
  
  if (is.null(retain)) {
    brank <- barcodeRanks(dat_filter, lower = pooling_threshold)
    
    knee <-
      as.integer(ifelse(is.null(brank$knee), brank@metadata$knee, brank$knee))
    
    if (is.null(knee)) {
      stop("Failed to calculate knee point. Check input data.")
    }
    retain <- knee
  }
  
  
  bc <- Matrix::colSums(dat_filter)
  #inflection <- ifelse(is.null(brank$inflection),brank@metadata$inflection,brank$inflection)
  if (length(names(which(bc >= retain))) > 0) {
    retain_cell <- names(bc)[bc >= retain]
    retain_mat <- raw_dat[, retain_cell]
    
    B0_cell <- names(bc)[bc <= pooling_threshold]
    B1_cell <- names(bc)[(bc > pooling_threshold & bc < retain)]
    B0 <- dat_filter[, B0_cell]
    dat <- dat_filter[, B1_cell]
    
    
    
    nzero_gene_B1 <- which(rowSums(dat) > 0)
    nzero_gene_B0 <- which(rowSums(B0) > 0)
    dat_filter <-
      dat_filter[sort(union(nzero_gene_B0, nzero_gene_B1)),]
    B0 <- dat_filter[, B0_cell]
    dat <- dat_filter[, B1_cell]
    
  }
  
  
  null_count <- Matrix::rowSums(B0)
  
  null_prob <- goodTuringProportions(null_count)
  
  
  
  #automatically calculate cluster cutoff
  c_threshold <-
    mean(unlist(replicate(1000, cor(
      rmultinom(1, 2 * pooling_threshold, null_prob),
      rmultinom(1, 2 * pooling_threshold, null_prob)
    ))))
  
  
  
  
  stat_fun <- function(x) {
    cor(x, null_prob)
  }
  
  
  if (print_prog)
    message("Done.\n")
  
  
  ##############################estimate test statistic for all barcodes
  
  if (print_prog)
    message("(2/5) Calculating test statistics for barcodes...")
  
  
  dat_Cor <-
    Calc_stat(dat,
              size = 1000,
              Ncores = Ncores,
              
              null_prob = null_prob)
  
  
  dat_temp <- data.frame(cbind(dat_Cor, Matrix::colSums(dat)))
  colnames(dat_temp) <- c("statistic", "count")
  if (print_prog)
    message("Done.\n")
  
  ###############################construct clusters
  clust_mat <- NULL
  if (print_prog) {
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
  if (print_prog)
    message("Done.\n")
  ##############################estimate test statistic for every cluster
  cl_Cor <- output_cl$stat
  
  cl_temp <- data.frame(cbind(cl_Cor, names(cl_Cor)))
  
  colnames(cl_temp) <- c("statistic", "count")
  cl_temp$statistic <- as.numeric(as.character(cl_temp$statistic))
  cl_temp$count <- as.numeric(as.character(cl_temp$count))
  
  clust_b <- unlist(output_cl$barcode)
  ##############################cluster test
  if (print_prog)
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
      i = 1:length(cand_count),
      .combine = "cfun",
      .packages = c("edgeR")
    ) %dopar% {
      input_stat = cl_Cor[cl_temp$count == cand_count[i]]
      v_temp <-
        rmultinom(1000, cand_count[i], null_prob)  #might improve. Redundant with barcode-wise test
      e_temp <- apply(v_temp, 2, stat_fun)
      
      pval <-
        sapply(input_stat, function(x)
          (sum(e_temp <= x) + 1) / 1001)
      
      names(pval) <- rep(cand_count[i], length(pval))
      list(sim = c(cand_count[i]), pval = pval)
      
    }
  stopCluster(cl)
  
  #####################find significant clusters
  for (i in 1:length(cand_count)) {
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
    clust_mat <- raw_dat[, res_b]
    
    dat_Cor <- dat_Cor[setdiff(colnames(dat), res_b)]
    
    dat <- dat[, setdiff(colnames(dat), res_b)]
    dat_temp <- dat_temp[setdiff(colnames(dat), res_b),]
    
    
  }
  
  if (print_prog)
    message("Done.\n")
  
  
  #########################################Permute null distribution of the test statistic
  if (print_prog)
    message("(5/5) Calculating empirical p-value for each barcode...")
  
  if ((ncol(dat)) == 0) {
    if (print_prog)
      message("\nNo remaining individual barcodes.\nAll finished.\n")
    if (print_prog)
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
    cand_barcode <- c(colnames(cbind(dat, B0)), colnames(retain_mat))
    ED_out <-
      emptyDrops(raw_dat[, cand_barcode], lower = pooling_threshold, retain = knee)
    dat_temp$pval <- ED_out$PValue[1:ncol(dat)]
    dat_temp$padj <- ED_out$FDR[1:ncol(dat)]
    cell_barcode <-
      cand_barcode[ifelse(is.na(ED_out$FDR), FALSE, ED_out$FDR <= FDR_threshold)]
    cell_mat <- raw_dat[, cell_barcode]
    if (print_prog)
      message("Done.\n")
    
    if (print_prog)
      message("All finished.\n")
    if (print_prog)
      print(Sys.time() - time_begin)
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


Filter_GB <- function(dat,g_threshold = 0, b_threshold = 0) {
  #filter barcodes and genes
  bc <- Matrix::colSums(dat)
  dat <- dat[, bc > b_threshold]
  gc <- Matrix::rowSums(dat)
  dat <- dat[gc > g_threshold, ]
}



Calc_stat <- function(dat, size = 1000, Ncores = detectCores() - 2, null_prob) {
  
    stat_fun <- function(x) {
      cor(x, null_prob)
    }

  
  #####parallel computation
  apply_ind <- split(1:ncol(dat), ceiling(1:ncol(dat) / size))  
  
  cl <- makeCluster(Ncores)
  registerDoParallel(cl)
  dat_tmp <- list()
  for (i in seq_along(apply_ind)) {
    dat_tmp[[i]] <- dat[,apply_ind[[i]],drop=F]
  }
  
  dat_tmp_it <- iter(dat_tmp)
  Cor <- foreach(x = dat_tmp_it, .combine = 'c', .packages = c("Matrix", "edgeR")) %dopar% {
    apply(x, 2, stat_fun)
  }
  stopCluster(cl)
  return(Cor)
}

#return size of each cluster in a list
size_cor <- function(cor_list) { 
  unlist(lapply(cor_list, function(x) ifelse(is.null(dim(x)), 1, dim(x)[1])))
}

ave_cor <- function(cor_list) {
  #remove diagonal 1, remove single-point clusters
  mean_vec <- c()
  size_clust <- size_cor(cor_list)
  
  for (i in 1:length(size_clust)) {
    if (size_clust[i] == 1) {
      mean_vec <- c(mean_vec, 0)   #single-point clusters will be filtered out by correlation threshold
    } else{
      cor_temp <- cor_list[[i]]
      diag(cor_temp) <- NA
      mean_vec <- c(mean_vec, mean(cor_temp, na.rm = T))
    }
  }
  return(mean_vec)
}

cor_clust <- function(cor_mat, Nclust = 100) {
  dist_mat <- as.dist(1 - cor_mat)
  hc <- hclust(dist_mat)
  hc1 <- cutree(hc, k = Nclust)
  
  cor_list <- list()
  for (i in 1:Nclust) {
    cor_list[[i]] <- cor_mat[which(hc1 == i), which(hc1 == i)]
    if (length(cor_list[[i]]) == 1) {
      cor_list[[i]] <- as.matrix(cor_list[[i]])
      colnames(cor_list[[i]]) <- colnames(cor_mat)[which(hc1 == i)]
      rownames(cor_list[[i]]) <- colnames(cor_list[[i]])
    }
  }
  
  return(cor_list)
}

sparse_cor <- function(x) {
  n <- nrow(x)
  m <- ncol(x)
  # non-empty rows
  ii <- unique(x@i) + 1 
  
  Ex <- colMeans(x)
  nozero <- as.vector(x[ii, ]) - rep(Ex, each = length(ii))        # colmeans
  
  covmat <- (crossprod(matrix(nozero, ncol = m)) +
               crossprod(t(Ex)) * (n - length(ii))) / (n - 1)
  sdvec <- sqrt(diag(covmat))
  return(covmat / crossprod(t(sdvec)))
}

Calc_cluster <- function(dat, Nclust = 100, c_threshold, dat_cor, dat_bc) {
  Output_stat <- c()
  cor_temp <- sparse_cor(dat)
  clust_temp <- cor_clust(cor_temp, Nclust = min(Nclust, ncol(dat)))
  cand_clust <-
    which(ave_cor(clust_temp) > c_threshold)  #highly correlated, more than 1 barcode
  if (length(cand_clust) == 0) {
    return(NULL)
  } else{
    cand_bclust <- list()
    Output_stat <- numeric(length(cand_clust))
    
    for (i in 1:length(cand_clust)) {
      cand_bclust[[i]] <- colnames(clust_temp[[cand_clust[i]]])
    }
    
    for (i in 1:length(cand_bclust)) {

        Output_stat[i] <- median(dat_cor[cand_bclust[[i]]])
        WM <-
          which.min(abs(dat_cor[cand_bclust[[i]]] - Output_stat[i]))
        names(Output_stat)[i] <-  dat_bc[cand_bclust[[i]]][WM]

      
    }
    names(cand_bclust) <- names(Output_stat)
    return(list(barcode = cand_bclust, stat = Output_stat))
  }
  
}


Calc_cluster_Paral <- function(dat,
                         size_hc = 2000,
                         Ncores = detectCores()-2,
                         Nclust = 100,
                         c_threshold,
                         dat_cor) {
  
  bc <- Matrix::colSums(dat)
  dat <- dat[,names(sort(bc))]
  
  #####parallel computation
  apply_ind <- split(1:ncol(dat), ceiling(1:ncol(dat) / size_hc))  
  
  cl <- makeCluster(Ncores)
  registerDoParallel(cl)
  dat_tmp <- list()
  for (i in seq_along(apply_ind)) {
    dat_tmp[[i]] <- dat[,apply_ind[[i]],drop=F]
  }
  
  dat_tmp_it <- iter(dat_tmp)
  
  cfun <- function(a, b) {
    barcode <- c(a$barcode, b$barcode)
    stat <- c(a$stat, b$stat)
    return(list(barcode = barcode, stat = stat))
  }
  
  Output <-
    foreach(
      aa = (dat_tmp_it),
      .combine = "cfun",
      .packages = c("Matrix"),
      .export = c(
        "sparse_cor",
        "Calc_cluster",
        "size_cor",
        "ave_cor",
        "cor_clust"
      )
    ) %dopar% {
      Calc_cluster(
        aa,
        Nclust,
        c_threshold,
        dat_cor = dat_cor,
        dat_bc = bc
      )
      
    }
  
  stopCluster(cl)
  return(Output)
}

