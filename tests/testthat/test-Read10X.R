hgmmh5 <- "~/Google Drive/Hallu/DATA/10x/testReadFunctions/hgmm100/hgmm_100_raw_gene_bc_matrices_h5.h5"
hg19 <- "~/Google Drive/Hallu/DATA/10x/testReadFunctions/hgmm100/raw_gene_bc_matrices/hg19/"
pbmch5 <- "~/Google Drive/Hallu/DATA/10x/testReadFunctions/PBMC5K/5k_pbmc_v3_filtered_feature_bc_matrix.h5"
pbmc <- "~/Google Drive/Hallu/DATA/10x/testReadFunctions/PBMC5K/filtered_feature_bc_matrix/"

test_that("Directory/files existence", {
    expect_error(Read10X("~/a"),"Directory does not exist")
    expect_error(Read10X("~"),"No 10X output file detected")
})

test_that("Data loading", {
    expect_identical(Read10X(hg19),Read10Xh5(hgmmh5)$hg19)
    expect_identical(Read10X(pbmc),Read10Xh5(pbmch5))
})

test_that("Use gene ID", {
    expect_match( rownames(Read10X(hg19,row.name = "id")),"ENSG")
    expect_match( rownames(Read10X(pbmc,row.name = "id")),"ENSG")
})


test_that("Metadata loading", {
    expect_true(is.list(Read10X(hg19,meta = TRUE)))
    expect_true(is.list(Read10X(pbmc,meta = TRUE)))
})
