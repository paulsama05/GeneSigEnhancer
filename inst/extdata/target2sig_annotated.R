### Uses a candidate gene (e.g. from `sig2target` output) to tailor a gene signature and then rerun GSVA and dependency correlations

## Parameters: 1) the dependency data frame, 2) the target gene,  3) the expression data frame, 4) a max cutoff for the new signature size (default 100), and 5) a starting gene signature (optional)

# Example dependency data taken from `sigCheck` output
dep.dat <- sig.list$dep.dat
if(any(sapply(dep.dat, is.numeric) == FALSE)){
  stop("First parameter, 'dep.dat', must be a data frame with ONLY numeric values.")
}
#Calculate the mean dependency of each gene
avg.dep <- sapply(dep.dat, mean)
# Define target gene and confirm it's presence in dependency data set
gene <- "ADAR"
if(!(gene %in% colnames(dep.dat))){
  stop(paste0("Provided `gene`, \"",gene,"\", not found in dependency data set, `dep.dat`."))
}
# Example expression data taken from `sigCheck` output
exp.dat <- sig.list$exp.dat
if(any(sapply(exp.dat, is.numeric) == FALSE)){
  stop("Third parameter, 'exp.dat', must be a data frame with ONLY numeric values.")
}
# Define cutoff
cutoff <- 100
if(!is.numeric(cutoff) | cutoff <1){
  stop("\"Cutoff\" value for maximum output gene signature length must be a positive number. Entries rounded down to integers")
}
cutoff <- as.integer(cutoff)
# Example gene signature taken from `sigCheck` output (first signature)
gene.sig <- sig.list$gene.sigs$A
# Check the structure of the gene signature (single vector of character type)
if(!any(is.na(gene.sig))){
  if(!is.character(gene.sig) | !is.vector(gene.sig)){
    stop("Unless NA, 'gene.sig' parameter must be a single vector of ONLY character values.")
  }
}

## The gene name is used to excise the relevant dependency data (with associated cell line identifiers)
gene.dep <- dep.dat[[gene]]
names(gene.dep) <- rownames(dep.dat)

## The relevant dependency data is then used to run a Pearson correlation against the expression of every gene and identify significantly correlated genes

# Expression and dependency data are "trimmed" to have the same cell line/sample identifiers in the same order
gene.dep.trim <- gene.dep[names(gene.dep) %in% rownames(exp.dat)]
gene.dep.trim <- gene.dep.trim[order(names(gene.dep.trim))]
exp.dat.trim <- exp.dat[rownames(exp.dat) %in% names(gene.dep),]
exp.dat.trim <- exp.dat.trim[order(rownames(exp.dat.trim)),]
# Pearson correlation is run between the target gene's dependency data and the expression data for every gene
dep.exp.cor <- sapply(exp.dat.trim, function(x){
  test.out <- cor.test(gene.dep.trim, x, method = "pearson")
  test.key <- c(test.out$estimate, test.out$p.value)
  return(test.key)
})
# Convert to a data frame where columns are variables
dep.exp.cor.df <- as.data.frame(t(dep.exp.cor))
# Label the initial data frame
colnames(dep.exp.cor.df) <- c('Pearson.Cor','P.value')
# Calculate the FDR (adjusted p-value) for each gene and add it to a new column
dep.exp.cor.df$FDR <- p.adjust(dep.exp.cor.df$`P.value`)
# Trim the data to only include genes with FDR values < 0.05 and reorder by FDR (small to large)
dep.exp.cor.df <- dep.exp.cor.df %>% filter(FDR<0.05) %>% arrange(FDR)


## Build the new signature list

#Confirm that the "cutoff" value supplied isn't larger than the number of significantly correlated genes
if(cutoff > nrow(dep.exp.cor.df)){
  cutoff <- nrow(dep.exp.cor.df)
}
# If no signature supplied, turn top X ("cutoff") hits into a new signature (stored as a data frame with one column)
if(any(is.na(gene.sig))){
  new.sigs <- as.data.frame(rownames(dep.exp.cor.df)[1:cutoff], ncol = 1)
  names(new.sigs) <- "top.cors"
  write.csv(new.sigs,"target2sig_new.sigs.csv",row.names = FALSE)
}else{
# When a comparison signature is provided, create a list of the following signatures:
  # 1) The top X ("cutoff") correlated genes
  # 2) The input gene signature
  # 3) Genes found in both the comparison signature and the COMPLETE list of significantly correlated genes
  # 4) Genes found in both the comparison signature and the top X correlations
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


## Reformat the expression data and the "new" gene signature list then run GSVA

# Transpose the expression data frame and convert to a matrix
exp.dat.mat <- as.matrix(t(exp.dat))
  # Genes by row, cell lines by column
# Convert the gene signature data frame to a list of lists and remove any blanks
new.sigs.list <- lapply(new.sigs, as.list)
# Run the GSVA
gsva.scores <- gsva(exp.dat.mat, new.sigs.list, verbose=TRUE,kcdf="Gaussian")
# Export a summary file that contains cell-line level GSVA scores for each unique signature
write.csv(gsva.scores,"target2sig_gsva.scores.csv")


## Align the cell lines found in the GSVA scores and dependency data frames

# Reduce GSVA output to only cell lines in dependency data and reorder (depends on # of input GSVA signatures)
if(nrow(gsva.scores)==1){
  gsva.scores.trim <- gsva.scores[,colnames(gsva.scores) %in% rownames(dep.dat)]
  gsva.scores.final <- as.data.frame(gsva.scores.trim[order(names(gsva.scores.trim))], ncol=1)
  colnames(gsva.scores.final) <- rownames(gsva.scores)
}else{
  gsva.scores.t <- as.data.frame(t(gsva.scores))
  gsva.scores.trim <- gsva.scores.t[rownames(gsva.scores.t) %in% rownames(dep.dat),]
  gsva.scores.final <- gsva.scores.trim[order(rownames(gsva.scores.trim)),]
}
# Reduce dependency data to only cell lines in GSVA output and reorder
dep.dat.trim <- dep.dat[rownames(dep.dat) %in% colnames(gsva.scores),]
dep.dat.final <- dep.dat.trim[order(rownames(dep.dat.trim)),]


## Use the GSVA signatures to run Pearson correlations and linear regressions against the cell line dependency data, then store the key data in a data frame

# Create a placeholder for the summary data
gsva.vs.dep.summary <-as.data.frame(matrix(NA, ncol = 6, nrow = 1))
colnames(gsva.vs.dep.summary) <- c('Pearson.Cor','P.value','FDR','Dependency','Gene','Signature')
# Use a for loop to call each GSVA signature one by one
for(i in 1:ncol(gsva.scores.final)){
  # And sapply to call the expression data of each gene one by one to run a pearson correlation and store the pearson scores and p-values
  gsva.dep.cor <- sapply(dep.dat.final, function(x){
    test.out <- cor.test(gsva.scores.final[,i], x, method = "pearson")
    test.key <- c(test.out$estimate, test.out$p.value)
    return(test.key)
  })
  # Convert to a data frame where columns are variables
  gsva.dep.cor.df <- as.data.frame(t(gsva.dep.cor))
  # Label the initial data frame
  colnames(gsva.dep.cor.df) <- c('Pearson.Cor','P.value')
  # Calculate the FDR (adjusted p-value) for each gene and add it to a new column
  gsva.dep.cor.df$FDR <- p.adjust(gsva.dep.cor.df$`P.value`)
  # Add in the average dependency scores for each gene
  gsva.dep.cor.df$Dependency <- avg.dep
  # Trim the data to only include genes with:
  # 1) Negative Pearson correlation scores (Looking for positive signature with negative dependency score)
  # 2) FDR values < 0.05
  #And rearrange by Pearson score (low)
  gsva.dep.cor.df <- gsva.dep.cor.df %>% filter(Pearson.Cor<0, FDR<0.05) %>% arrange(Pearson.Cor)
  # Move the gene names to a column and add a column to identify the corresponding GSVA signature
  gsva.dep.cor.df <- gsva.dep.cor.df %>% mutate(Gene = rownames(gsva.dep.cor.df), Signature = colnames(gsva.scores.final)[i])
  # Bind the single-signature data frame to the "master" data frame
  gsva.vs.dep.summary <- bind_rows(gsva.vs.dep.summary, gsva.dep.cor.df)
}
# Trim off the first row (leftover from the placeholder) and reorganize the data
gsva.vs.dep.summary <- gsva.vs.dep.summary[-1,]
gsva.vs.dep.summary <- gsva.vs.dep.summary %>% select(c(6,5,4,1:3)) %>% arrange(Signature, FDR)
# Export the complete summary table
write.csv(gsva.vs.dep.summary,"target2sig_gsva.vs.dep.csv",row.names = FALSE)


## Compile the key summary data (new gene signatures, GSVA scores,dependency correlations) into a list and return
summary.list <- list(new.sigs, gsva.scores, gsva.vs.dep.summary)
names(summary.list) <- c("new.sigs","gsva.scores","gsva.vs.dep")
return(summary.list)
