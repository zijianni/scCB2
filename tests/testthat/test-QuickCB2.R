hgmmh5 <- "~/Google Drive/Hallu/DATA/10x/testReadFunctions/hgmm100/hgmm_100_raw_gene_bc_matrices_h5.h5"
hg19 <- "~/Google Drive/Hallu/DATA/10x/testReadFunctions/hgmm100/raw_gene_bc_matrices/hg19/"


test_that("Both directory and HDF5 are used", {
    expect_error(QuickCB2(dir=hg19, h5file=hgmmh5),"dir and h5file can't be specified together.")
})

