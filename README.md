# R package 'CB'

'CB' is an R package for pre-processing raw droplet-based single cell RNA-seq data. It helps distinguish real cells from empty droplets using a clustering-based Monte-Carlo test. 

'CB' provides easy-to-use functions `Read10X` and `Read10Xh5` to read 10X Chromium output files and generate feature-by-barcode count matrix. The main function `CBFindCell` takes count matrix as input and return testing information for each barcode. `GetCellMat` extracts real cell matrix from `CBFindCell` output. Since the output format is still matrix, it can be easily integrated into downstream statistical analysis pipelines like 'Seurat'.
