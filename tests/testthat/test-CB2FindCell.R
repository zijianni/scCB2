data(mbrainSub)
CBOut <- CB2FindCell(mbrainSub, FDR_threshold = 0.01,
                    lower = 100, Ncores = 2,PrintProg = FALSE)


test_that("Output CBout format", {
    expect_true(all(c("cluster_matrix","cell_matrix")%in%names(CBOut)))
})

test_that("Background threshold", {
    for(pooling in c(50, 100)){
        CBOut_temp <- CB2FindCell(mbrainSub,
                            lower = pooling, 
                            Ncores = 2, PrintProg = FALSE)
        expect_true(all(Matrix::colSums(CBOut_temp$cluster_matrix)>pooling))
        expect_true(all(Matrix::colSums(CBOut_temp$cell_matrix)>pooling))
    }
})
                                                                                

test_that("Upper threshold", {
    for(Upper in c(500, 1000)){
        CBOut_temp <- CB2FindCell(mbrainSub,
                            upper = Upper,
                            Ncores = 2, PrintProg = FALSE)
        Upper_b <- colnames(mbrainSub)[Matrix::colSums(mbrainSub)>=Upper]
        Cell_b <- colnames(cbind(CBOut_temp$cluster_matrix,
                            CBOut_temp$cell_matrix))
        expect_true(all(Upper_b%in%Cell_b))
    }
})
