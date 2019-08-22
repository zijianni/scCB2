# scCB2

CB2 is a cluster-based approach for distinguishing true cells from background barcodes in droplet-based single cell RNA-seq experiments (especially for 10X Chromium output), while scCB2 is its corresponding R package. It is based on clustering similar barcodes and calculating Monte-Carlo p-value for each cluster to test against background distribution. This cluster-level test outperforms single-barcode-level tests not only for high count barcodes, but also in dealing with low count barcodes and homogeneous sequencing library, while keeping FDR well controlled.
