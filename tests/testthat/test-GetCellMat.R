data(mbrainSub)

mbrainReal <- mbrainSub[, Matrix::colSums(mbrainSub) > 500]
CBOut <- list(cluster_matrix = mbrainReal[,
                    sample(ncol(mbrainReal), 100, replace = TRUE)],
                    cell_matrix = mbrainReal[,
                    sample(ncol(mbrainReal), 100, replace = TRUE)])


test_that("Output matrix format", {
    expect_identical(class(GetCellMat(CBOut)), class(CBOut$cell_matrix))
})

test_that("Mitochondrial filtering", {
    expect_error(GetCellMat(CBOut, MTgene = "1"),
                 "invalid character indexing")
    expect_error(GetCellMat(CBOut, MTfilter = -1), "between 0 and 1")
    
    MTgene <- sample(rownames(CBOut$cell_matrix), 10)
    for (k in c(0.5, 0.1, 0.01, 0.001, 0)) {
        CBMat <- GetCellMat(CBOut, k, MTgene)
        MTprop <-
            Matrix::colSums(CBMat[MTgene, ]) / Matrix::colSums(CBMat)
        expect_true(all(MTprop <= k))
    }
})
