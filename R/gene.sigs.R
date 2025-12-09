#' Data Frame of 5 gene signatures compatible with CCLE/DepMap gene nomenclature
#'
#' This data was imported from the file `gene_signature_set.csv` using the `sigImport` fuction to confirm that the gene nomenclature used matches that in CCLE/DepMap files.
#' @description
#' Example gene signatures for use in other functions from this package. Each successive signature is a "trimmed" version of the one before it. When used with `sig2target`, the `gsva.vs.exp` output highlights that the "bottom" genes of each list contribute least to GSVA.
#' @source Access original file using `system.file("extdata","gene_signature_set.csv",package = "GeneSigEnhancer")`. This file can be directly fed into `sigImport`, which will output the data frame as part of a list.
#' @docType data
#' @keywords datasets
#' @name gene.sigs
#' @usage sig2target(gene.sigs, exp.dat, dep.dat)
NULL
