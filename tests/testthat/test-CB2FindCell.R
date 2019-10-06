data(mbrainSub)
CBOut <- CB2FindCell(mbrainSub, FDR_threshold = 0.01,
                    background_threshold = 100, Ncores = 2,PrintProg = FALSE)


test_that("Output CBout format", {
    expect_true(all(c("cluster_matrix","cell_matrix")%in%names(CBOut)))
})

test_that("Background threshold", {
    for(pooling in c(50, 100)){
        CBOut_temp <- CB2FindCell(mbrainSub,
                            background_threshold = pooling, 
                            Ncores = 2, PrintProg = FALSE)
        expect_true(all(Matrix::colSums(CBOut_temp$cluster_matrix)>pooling))
        expect_true(all(Matrix::colSums(CBOut_temp$cell_matrix)>pooling))
    }
})
                                                                                

test_that("Retain threshold", {
    for(Retain in c(500, 1000)){
        CBOut_temp <- CB2FindCell(mbrainSub,
                            retain = Retain,
                            Ncores = 2, PrintProg = FALSE)
        Retain_b <- colnames(mbrainSub)[Matrix::colSums(mbrainSub)>=Retain]
        Cell_b <- colnames(cbind(CBOut_temp$cluster_matrix,
                            CBOut_temp$cell_matrix))
        expect_true(all(Retain_b%in%Cell_b))
    }
})
