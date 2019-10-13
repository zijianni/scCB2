# scCB2

Droplet-based single cell RNA-seq technologies provide a novel insight in transcriptome profiles of individual cells with a much larger cell numbers and cheaper cost compared with microfluidic-based single cell RNA-seq. During library preparation, each cell is expected to be captured by one droplet. The number of droplets are usually much more than the number of cells, thus most droplets do not contain real cells. However, there are always free-floating RNA fragments in the library due to broken cells or library contamination. Empty droplets will capture them and have non-zero expression values after sequencing.


**CB2** is a cluster-based approach for distinguishing true cells from background barcodes in droplet-based single cell RNA-seq experiments (especially for 10X Chromium output), while `scCB2` is its corresponding R package. It is based on clustering similar barcodes and calculating Monte-Carlo p-value for each cluster to test against background distribution. This cluster-level test outperforms single-barcode-level tests not only for high count barcodes, but also in dealing with low count barcodes and homogeneous sequencing library, while keeping FDR well controlled.

Install `scCB2` via Github into R using `devtools`: 

`devtools::install_github("zijianni/scCB2")`
