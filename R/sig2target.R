#' Uses gene signatures, expression data, and dependency data to identify target genes
#'
#' User inputs data frames of gene signatures, gene expression, and gene dependency data (e.g. from `sigImport`) and function both returns a list of key results and exports those results to .csv files starting with "sig2target".
#' @param gene.sigs Data frame of gene signatures (each column a unique signature with signature/column name). All gene names must match those in expression data frame.
#' @param exp.dat Data frame of gene expression data organized with genes in columns, samples (e.g. cell lines) in rows. Sample name formatting must match that of dependency data (although not all samples must match)
#' @param dep.dat Data frame of gene dependency data organized with genes in columns, samples (e.g. cell lines) in rows. Sample name formatting must match that of expression data (although not all samples must match)
#' @return A list 3 objects: 1) Matrix of GSVA scores with each row corresponding to an input signature, 2) data frame of GSVA-vs.-expression Pearson correlations limited to signature genes, 3) data frame of GSVA-vs.-dependency Pearson correlations limited to significantly correlated genes
#' @description
#' Runs GSVA using gene signatures and gene expression data, then uses GSVA and dependency data to identify candidate targets.
#'
#' @export
sig2target <- function(gene.sigs, exp.dat, dep.dat){
if(!requireNamespace("BiocManager", quietly = TRUE)){install.packages("BiocManager")}
if(!requireNamespace("GSVA", quietly = TRUE)){BiocManager::install("GSVA")}
library(GSVA)
library(dplyr)
if(any(sapply(gene.sigs, is.character) == FALSE)){
  stop("First parameter, 'gene.sigs', must be a data frame with ONLY character values.")
}
if(any(sapply(exp.dat, is.numeric) == FALSE)){
  stop("Second parameter, 'exp.dat', must be a data frame with ONLY numeric values.")
}
if(any(sapply(dep.dat, is.numeric) == FALSE)){
  stop("Third parameter, 'dep.dat', must be a data frame with ONLY numeric values.")
}
gene.sigs.mat <- as.matrix(gene.sigs)
gene.sigs.mat <- gene.sigs.mat[gene.sigs.mat!=""]
if(sum(gene.sigs.mat %in% colnames(exp.dat)) == length(gene.sigs.mat)){
  print("All genes in signature(s) match those in expression data frame.")
}else{
  print(str_c(unique(gene.sigs.mat[!(gene.sigs.mat %in% colnames(exp.dat))]), collapse = ", "))
  stop("The above genes were NOT found in data frame. Please alter signature to match nomenclature in data frame and try again.")
}
exp.dat.mat <- as.matrix(t(exp.dat))
gene.sigs.list <- lapply(gene.sigs, as.list)
gene.sigs.list <- lapply(gene.sigs.list, function(x){x[x!=""]})
gsva.scores <- gsva(exp.dat.mat, gene.sigs.list, verbose=TRUE,kcdf="Gaussian")
write.csv(gsva.scores,"sig2target_gsva.scores.csv")
gsva.vs.exp.summary <-as.data.frame(matrix(NA, ncol = 5, nrow = 1))
colnames(gsva.vs.exp.summary) <- c('Pearson.Cor','P.value','FDR','Gene','Signature')
for(i in 1:nrow(gsva.scores)){
  gsva.exp.cor <- sapply(exp.dat, function(x){
    test.out <- cor.test(gsva.scores[i,], x, method = "pearson")
    test.key <- c(test.out$estimate, test.out$p.value)
    return(test.key)
  })
  gsva.exp.cor.df <- as.data.frame(t(gsva.exp.cor))
  colnames(gsva.exp.cor.df) <- c('Pearson.Cor','P.value')
  gsva.exp.cor.df$FDR <- p.adjust(gsva.exp.cor.df$`P.value`)
  gsva.exp.cor.df <- gsva.exp.cor.df[rownames(gsva.exp.cor.df) %in% gene.sigs[,i],]
  gsva.exp.cor.df <- gsva.exp.cor.df %>% mutate(Gene = rownames(gsva.exp.cor.df), Signature = rownames(gsva.scores)[i])
  gsva.vs.exp.summary <- bind_rows(gsva.vs.exp.summary, gsva.exp.cor.df)
}
gsva.vs.exp.summary <- gsva.vs.exp.summary[-1,]
gsva.vs.exp.summary <- gsva.vs.exp.summary %>% select(c(5,4,1:3)) %>% arrange(Signature, FDR)
write.csv(gsva.vs.exp.summary,"sig2target_gsva.vs.exp.csv",row.names = FALSE)
avg.dep <- sapply(dep.dat, mean)
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
write.csv(gsva.vs.dep.summary,"sig2target_gsva.vs.dep.csv",row.names = FALSE)
summary.list <- list(gsva.scores, gsva.vs.exp.summary, gsva.vs.dep.summary)
names(summary.list) <- c("gsva.scores","gsva.vs.exp","gsva.vs.dep")
return(summary.list)
}
