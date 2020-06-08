# Updates

## scCB2 0.99.0

---------------------

* Created initial github page.
* Created package.
* Created vignettes.
* Checked package... passed R CMD check and BiocCheck.
* Ready to submit to Bioconductor.

## scCB2 0.99.11

---------------------

* Customized knee point calculation function.
* CB2FindCell will print out retain threshold in standard output.
* Minor bug fixes when there is no cluster.

## scCB2 0.99.12

---------------------

* Minor changes on vignette.
* Baseline (2 * background) clustering threshold will be printed out.

## scCB2 0.99.13

---------------------

* Rounded baseline threshold.
* Minor changes on vignette.
* Bug fix when retain threshold is larger than any barcode.

## scCB2 0.99.14

---------------------

* Change parameter names to match those in the paper: retain -> upper. background -> lower.

## scCB2 0.99.15

---------------------

* When dividing barcodes into groups with similar barcode counts, the last group will be combined with second last group if the number of barcodes in the last group is less than half of that in the second last group. 

## scCB2 0.99.16

---------------------

* Change function name `Read10X` to be `Read10xRaw`, `Read10Xh5` to be `Read10xRawH5`.
* Change parameter name `PrintProg` to be `verbose`.
* Minor edits on parameter description of function `GetCellMat`.
* Added a quick function `QuickCB2` by combining all necessary functions into one. Input: directory of raw data. Output: filtered matrix, or a Seurat object containing filtered matrix.

## scCB2 0.99.17

---------------------

* Initial version preparing to submit to Bioconductor.

## scCB2 0.99.18

---------------------

* Minor bug fix.
* Minor edits on package vignettes.

## scCB2 0.99.19

---------------------

* Minor edits on package vignettes.

## scCB2 0.99.20

---------------------

* Updated citation.

