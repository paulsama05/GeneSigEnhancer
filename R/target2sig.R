#' Uses target gene dependency data to create/test new gene signature(s)
#'
#' User inputs a target gene name and data frames of gene dependency and gene expression data (e.g. from `sigImport`) and function both returns a list of key results and exports those results to .csv files starting with "target2sig".
#' @param dep.dat Data frame of gene dependency data organized with genes in columns, samples (e.g. cell lines) in rows. Sample name formatting must match that of expression data (although not all samples must match)
#' @param gene Name of target gene. Must exactly match that in dependency data set.
#' @param exp.dat Data frame of gene expression data organized with genes in columns, samples (e.g. cell lines) in rows. Sample name formatting must match that of dependency data (although not all samples must match)
#' @param cutoff  Positive integer value that sets the maximum size of the output gene signature. Default is 100.
#' @param gene.sig OPTIONAL. A single gene signature vector. Names must match those in dependency/expression data sets.
#' @return A list 3 objects: 1) List of new gene signature(s), 2) matrix of GSVA scores with each row corresponding to new signature(s), 3) data frame of GSVA-vs.-dependency Pearson correlations limited to significantly correlated genes
#' @description
#' Creates new gene signature(s) by finding genes whose expression correlates to target gene dependency. New signature(s) used to run GSVA and then compare back to gene dependency data.
#'
#' @export
target2sig <- function(dep.dat, gene, exp.dat, cutoff = 100, gene.sig = NA){
if(!requireNamespace("BiocManager", quietly = TRUE)){install.packages("BiocManager")}
if(!requireNamespace("GSVA", quietly = TRUE)){BiocManager::install("GSVA")}
library(GSVA)
library(dplyr)
if(any(sapply(dep.dat, is.numeric) == FALSE)){
  stop("First parameter, 'dep.dat', must be a data frame with ONLY numeric values.")
}
avg.dep <- sapply(dep.dat, mean)
if(!(gene %in% colnames(dep.dat))){
  stop(paste0("Provided `gene`, \"",gene,"\", not found in dependency data set, `dep.dat`."))
}
if(any(sapply(exp.dat, is.numeric) == FALSE)){
  stop("Third parameter, 'exp.dat', must be a data frame with ONLY numeric values.")
}
if(!is.numeric(cutoff) | cutoff <1){
  stop("\"Cutoff\" value for maximum output gene signature length must be a positive number. Entries rounded down to integers")
}
cutoff <- as.integer(cutoff)
if(!any(is.na(gene.sig))){
  if(!is.character(gene.sig) | !is.vector(gene.sig)){
    stop("Unless NA, 'gene.sig' parameter must be a single vector of ONLY character values.")
  }
}
gene.dep <- dep.dat[[gene]]
names(gene.dep) <- rownames(dep.dat)
gene.dep.trim <- gene.dep[names(gene.dep) %in% rownames(exp.dat)]
gene.dep.trim <- gene.dep.trim[order(names(gene.dep.trim))]
exp.dat.trim <- exp.dat[rownames(exp.dat) %in% names(gene.dep),]
exp.dat.trim <- exp.dat.trim[order(rownames(exp.dat.trim)),]
dep.exp.cor <- sapply(exp.dat.trim, function(x){
  test.out <- cor.test(gene.dep.trim, x, method = "pearson")
  test.key <- c(test.out$estimate, test.out$p.value)
  return(test.key)
})
dep.exp.cor.df <- as.data.frame(t(dep.exp.cor))
colnames(dep.exp.cor.df) <- c('Pearson.Cor','P.value')
dep.exp.cor.df$FDR <- p.adjust(dep.exp.cor.df$`P.value`)
dep.exp.cor.df <- dep.exp.cor.df %>% filter(FDR<0.05) %>% arrange(FDR)
if(cutoff > nrow(dep.exp.cor.df)){
  cutoff <- nrow(dep.exp.cor.df)
}
if(any(is.na(gene.sig))){
  new.sigs <- as.data.frame(rownames(dep.exp.cor.df)[1:cutoff], ncol = 1)
  names(new.sigs) <- "top.cors"
  write.csv(new.sigs,"target2sig_new.sigs.csv",row.names = FALSE)
}else{
  top.cors <- rownames(dep.exp.cor.df)[1:cutoff]
  full.overlap <- gene.sig[gene.sig %in% rownames(dep.exp.cor.df)]
  cutoff.overlap <- gene.sig[gene.sig %in% top.cors]
  new.sigs <- list(top.cors, gene.sig, full.overlap, cutoff.overlap)
  names(new.sigs) <- c("top.cors","input.sig","full.overlap","top.overlap")
  max.length <- max(lengths(new.sigs))
  new.sigs.equal <- sapply(new.sigs, function(x){
    length(x) <- max.length
    return(x)
  })
  new.sigs.equal[is.na(new.sigs.equal)] <- ""
  write.csv(new.sigs.equal,"target2sig_new.sigs.csv",row.names = FALSE)
}
exp.dat.mat <- as.matrix(t(exp.dat))
new.sigs.list <- lapply(new.sigs, as.list)
gsva.scores <- gsva(exp.dat.mat, new.sigs.list, verbose=TRUE,kcdf="Gaussian")
write.csv(gsva.scores,"target2sig_gsva.scores.csv")
if(nrow(gsva.scores)==1){
  gsva.scores.trim <- gsva.scores[,colnames(gsva.scores) %in% rownames(dep.dat)]
  gsva.scores.final <- as.data.frame(gsva.scores.trim[order(names(gsva.scores.trim))], ncol=1)
  colnames(gsva.scores.final) <- rownames(gsva.scores)
}else{
  gsva.scores.t <- as.data.frame(t(gsva.scores))
  gsva.scores.trim <- gsva.scores.t[rownames(gsva.scores.t) %in% rownames(dep.dat),]
  gsva.scores.final <- gsva.scores.trim[order(rownames(gsva.scores.trim)),]
}
dep.dat.trim <- dep.dat[rownames(dep.dat) %in% colnames(gsva.scores),]
dep.dat.final <- dep.dat.trim[order(rownames(dep.dat.trim)),]
gsva.vs.dep.summary <-as.data.frame(matrix(NA, ncol = 6, nrow = 1))
colnames(gsva.vs.dep.summary) <- c('Pearson.Cor','P.value','FDR','Dependency','Gene','Signature')
for(i in 1:ncol(gsva.scores.final)){
  gsva.dep.cor <- sapply(dep.dat.final, function(x){
    test.out <- cor.test(gsva.scores.final[,i], x, method = "pearson")
    test.key <- c(test.out$estimate, test.out$p.value)
    return(test.key)
  })
  gsva.dep.cor.df <- as.data.frame(t(gsva.dep.cor))
  colnames(gsva.dep.cor.df) <- c('Pearson.Cor','P.value')
  gsva.dep.cor.df$FDR <- p.adjust(gsva.dep.cor.df$`P.value`)
  gsva.dep.cor.df$Dependency <- avg.dep
  gsva.dep.cor.df <- gsva.dep.cor.df %>% filter(Pearson.Cor<0, FDR<0.05) %>% arrange(Pearson.Cor)
  gsva.dep.cor.df <- gsva.dep.cor.df %>% mutate(Gene = rownames(gsva.dep.cor.df), Signature = colnames(gsva.scores.final)[i])
  gsva.vs.dep.summary <- bind_rows(gsva.vs.dep.summary, gsva.dep.cor.df)
}
gsva.vs.dep.summary <- gsva.vs.dep.summary[-1,]
gsva.vs.dep.summary <- gsva.vs.dep.summary %>% select(c(6,5,4,1:3)) %>% arrange(Signature, FDR)
write.csv(gsva.vs.dep.summary,"target2sig_gsva.vs.dep.csv",row.names = FALSE)
summary.list <- list(new.sigs, gsva.scores, gsva.vs.dep.summary)
names(summary.list) <- c("new.sigs","gsva.scores","gsva.vs.dep")
return(summary.list)
}
