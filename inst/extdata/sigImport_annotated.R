### Import gene signatures, gene expression data, and gene dependency data from .csv files and prepare for use in `sig2target`

## Prompt user to import the target files

# Prompt user for name of gene signature file
gene.sigs.file <- readline(prompt = "Enter the name of your gene signature file EXACTLY and include the \".csv\" extension.\nNOTE: Make sure each distinct signature/list is in its own column with a title (header) at the top. ")
# Check to make sure file exits in current working directory
while(!(gene.sigs.file %in% list.files())){
  gene.sigs.file <- readline(prompt = paste0("No such file found in this directory (",getwd() ,"). Try new name or ESC and reset directory. "))
}
# Use file name to import file as a data frame
gene.sigs <- read.csv(gene.sigs.file, header = TRUE)
# Check to make sure the imported data is of the correct type (columns of character data)
if(any(sapply(gene.sigs, is.character) == FALSE)){
  stop("Gene signature object must be a data frame with ONLY character values.")
}
# Prompt user for name of gene expression file
exp.dat.file <- readline(prompt = "Enter the name of your expression data file EXACTLY and include the \".csv\" extension.\nNOTE: This should be formatted like CCLE/DepMap expression data (columns by gene name, rows by cell line). ")
# Check to make sure file exists in current working directory
while(!(exp.dat.file %in% list.files())){
  exp.dat.file <- readline(prompt = paste0("No such file found in this directory (",getwd() ,"). Try new name or ESC and reset directory. "))
}
# Use file name to import file as a data frame and then reformat the column (gene) names
exp.dat <- read.csv(exp.dat.file, header = TRUE, row.names = 1)
colnames(exp.dat) <- word(colnames(exp.dat), 1, sep=fixed('..'))
# Check to make sure there is at lest one column of data (if a single signature is accidentally imported in this way, there will be 0 columns)
if(ncol(exp.dat)==0){
  stop("No data found. Was a single gene signature accidentally imported?")
}
# Check to make sure the imported data is of the correct type (columns of character data)
if(any(sapply(exp.dat, is.numeric) == FALSE)){
  stop("Gene expression object must be a data frame with ONLY numeric values.")
}


## Confirm that all signature genes are in the expression data frame

# Conversion to matrix simplifies comparison to gene names in the expression data frame
gene.sigs.mat <- as.matrix(gene.sigs)
# When signature lengths are unequal, blanks ("") are inserted into the list that must be removed
gene.sigs.mat <- gene.sigs.mat[gene.sigs.mat!=""]
# Check the genes in the signatures against the expression data frame
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


##Import gene dependency data and check overlap with expression data

#Prompt user for name of gene dependency file
dep.dat.file <- readline(prompt = "Enter the name of your gene-dependency data file EXACTLY and include the \".csv\" extension.\nNOTE: This should be formatted like CCLE/DepMap dependency data (columns by gene name, rows by cell line). ")
# Check to make sure file exits in current working directory
while(!(dep.dat.file %in% list.files())){
  dep.dat.file <- readline(prompt = paste0("No such file found in this directory (",getwd() ,"). Try new name or ESC and reset directory. "))
}
# Use file name to import file as a data frame and then reformat the column (gene) names
dep.dat <- read.csv(dep.dat.file, header = TRUE, row.names = 1)
colnames(dep.dat) <- word(colnames(dep.dat), 1, sep=fixed('..'))
# Check the structure of the expression data frame
if(any(sapply(dep.dat, is.numeric) == FALSE)){
  stop("Gene dependency object must be a data frame with ONLY numeric values.")
}


## Report the number of genes from the expression set that are also in the dependency set

genes.present <- sum(colnames(exp.dat) %in% colnames(dep.dat))
print(paste0("Of ", length(colnames(exp.dat)), " unique genes in the expression data set there are ", genes.present, " genes matched to the gene dependency data set."))


## Return a list of the key objects, which can be used in other functions

sig.list <- list(gene.sigs,exp.dat, dep.dat)
names(sig.list) <- c("gene.sigs","exp.dat","dep.dat")
return(sig.list)
