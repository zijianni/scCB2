#' Check different background cutoffs and recommend an appropriate one
#' 
#' The key parameter of CB2 as well as other similar methods is the 
#' background cutoff, which divides barcodes 
#' into two groups: (1) small barcodes that are most likely to be 
#' background; (2) the rest barcodes that can be either background or cell, and 
#' remain to be tested. Those small barcodes will be used to estimate a 
#' background distribution, which guides the identification of cells from 
#' background. It is crucial to have an unbiased estimation of the background 
#' distribution. 
#' 
#' An appropriate background cutoff should be reasonably large to contain 
#' enough background information, but shouldn't be too large to mistakenly 
#' include real cells. We recommend a background cutoff which 
#' (1) puts more than 90% barcodes into background, or 
#' (2) puts more than 10% UMI counts into background. 
#' The smallest cutoff satisfying either condition is the recommended cutoff.
#' 
#' @param RawDat Matrix. Supports standard matrix or sparse matrix. 
#' This is the raw feature-by-barcode count matrix.
#' 
#' @return A list containing a data frame summarizing background information 
#' under different background cutoffs, and the recommended background cutoff 
#' for the input data. For the data frame, `n_bg_bcs` is the number of 
#' barcodes less or equal to the cutoff, `n_bg_counts` is the number of UMI 
#' counts within the barcodes less or equal to the cutoff, 
#' `prop_bg_bcs` and `prop_bg_counts` are the corresponding proportions.

#' @examples 
#' data(mbrainSub)
#' CheckBackgroundCutoff(mbrainSub)
#' 
#' @importFrom Matrix colSums
#'
#' @export

CheckBackgroundCutoff <- function(RawDat){
    cand_cutoffs <- 2:10*50
    summary_df <- data.frame(cutoff=cand_cutoffs, n_bg_bcs=0, 
                             n_bg_counts=0, prop_bg_bcs=0, prop_bg_counts=0)
    for(idx in seq_along(cand_cutoffs)){
        
        col_sums <- colSums(RawDat)
        bg_bcs <- col_sums<=cand_cutoffs[idx]
        summary_df[idx,2] <- sum(bg_bcs)
        summary_df[idx,3] <- sum(col_sums[bg_bcs])
        summary_df[idx,4] <- mean(bg_bcs)
        summary_df[idx,5] <- summary_df[idx,3]/sum(col_sums)
    }
    recommended_cutoff <- min(summary_df$cutoff[
        summary_df$prop_bg_bcs>=0.9 | summary_df$prop_bg_counts>=0.1
    ])
    return(list(summary_table=summary_df, recommended_cutoff=recommended_cutoff))
    
}
