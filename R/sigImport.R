#' Import gene signature(s), gene expression data, and gene dependency data for downstream use
#'
#' Function is only compatible with .csv files stored in the current working directory. When prompted, user provides full file names and `sigImport` checks that they are compatible with the requirements of functions in this package.
#' @param NULL No parameters passed to function; function prompts user for file names and checks that they are suitable
#' @return A list of 3 objects: 1) a data frame of gene signatures (`gene.sigs`), 2) a data frame of gene expression data (`exp.dat`), and 3) a data frame of gene dependency data (`dep.dat`)
#' @description
#' Guides user in importing key data for use in the `sig2target`, `target2sig` and `sigAccuracy` functions from this package. For data formatting, please pay close attention to each prompt!
#'
#' @export
sigImport <- function(){
library(stringr)
gene.sigs.file <- readline(prompt = "Enter the name of your gene signature file EXACTLY and include the \".csv\" extension.\nNOTE: Make sure each distinct signature/list is in its own column with a title (header) at the top. ")
while(!(gene.sigs.file %in% list.files())){
  gene.sigs.file <- readline(prompt = paste0("No such file found in this directory (",getwd() ,"). Try new name or ESC and reset directory. "))
}
gene.sigs <- read.csv(gene.sigs.file, header = TRUE)
if(any(sapply(gene.sigs, is.character) == FALSE)){
  stop("Gene signature object must be a data frame with ONLY character values.")
}
exp.dat.file <- readline(prompt = "Enter the name of your expression data file EXACTLY and include the \".csv\" extension.\nNOTE: This should be formatted like CCLE/DepMap expression data (columns by gene name, rows by cell line). ")
while(!(exp.dat.file %in% list.files())){
  exp.dat.file <- readline(prompt = paste0("No such file found in this directory (",getwd() ,"). Try new name or ESC and reset directory. "))
}
exp.dat <- read.csv(exp.dat.file, header = TRUE, row.names = 1)
colnames(exp.dat) <- word(colnames(exp.dat), 1, sep=fixed('..'))
if(ncol(exp.dat)==0){
  stop("No data found. Was a single gene signature accidentally imported?")
}
if(any(sapply(exp.dat, is.numeric) == FALSE)){
  stop("Gene expression object must be a data frame with ONLY numeric values.")
}
gene.sigs.mat <- as.matrix(gene.sigs)
gene.sigs.mat <- gene.sigs.mat[gene.sigs.mat!=""]
if(sum(gene.sigs.mat %in% colnames(exp.dat)) == length(gene.sigs.mat)){
  print("All genes in signature(s) match those in expression data set.")
}else{
  unmatched <- unique(gene.sigs.mat[!(gene.sigs.mat %in% colnames(exp.dat))])
  unmatched.trim <- str_sub(unmatched, start = 1, end = 3)
  poss.matches <- lapply(unmatched.trim, function(x){
    grep(paste0("^",x), colnames(exp.dat), value =TRUE)
  })
  names(poss.matches) <- unmatched
  print(poss.matches)
  stop(paste0(length(poss.matches)," signature gene(s) not found in expression data frame.\n\tSee above printout for name(s) and possible matches.\n\tConsider converting to HUGO nomenclature."))
}
dep.dat.file <- readline(prompt = "Enter the name of your gene-dependency data file EXACTLY and include the \".csv\" extension.\nNOTE: This should be formatted like CCLE/DepMap dependency data (columns by gene name, rows by cell line). ")
while(!(dep.dat.file %in% list.files())){
  dep.dat.file <- readline(prompt = paste0("No such file found in this directory (",getwd() ,"). Try new name or ESC and reset directory. "))
}
dep.dat <- read.csv(dep.dat.file, header = TRUE, row.names = 1)
colnames(dep.dat) <- word(colnames(dep.dat), 1, sep=fixed('..'))
if(any(sapply(dep.dat, is.numeric) == FALSE)){
  stop("Gene dependency object must be a data frame with ONLY numeric values.")
}
genes.present <- sum(colnames(exp.dat) %in% colnames(dep.dat))
print(paste0("Of ", length(colnames(exp.dat)), " unique genes in the expression data set there are ", genes.present, " genes matched to the gene dependency data set."))
sig.list <- list(gene.sigs,exp.dat, dep.dat)
names(sig.list) <- c("gene.sigs","exp.dat","dep.dat")
return(sig.list)
}
