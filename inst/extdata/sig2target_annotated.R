### Use gene signatures (lists) and gene expression data to generate GSVA signatures then combine with gene dependency data to identify candidate target genes


## Import the relevant data and perform initial quality checks
  ## NOT IN FUNCTION: import steps, which are part of the function `sigCheck` or done manually by the user
  ## Where possible (size), example data objects will be included in the package

# Import gene signature(s) as a data frame
gene.sigs <- read.csv("gene_signature_set.csv", header = T)
# Check the structure of the signature data frame (all columns should of class "character")
if(any(sapply(gene.sigs, is.character) == FALSE)){
  stop("First parameter, 'gene.sigs', must be a data frame with ONLY character values.")
}

# Import expression data frame and isolate gene names
exp.dat <- read.csv("CCLE_expression.csv", header = TRUE, row.names = 1)
colnames(exp.dat) <- word(colnames(exp.dat), 1, sep=fixed('..'))
# Check the structure of the expression data frame
if(any(sapply(exp.dat, is.numeric) == FALSE)){
  stop("Second parameter, 'exp.dat', must be a data frame with ONLY numeric values.")
}

# Confirm that all signature genes are in the expression data frame
gene.sigs.mat <- as.matrix(gene.sigs)
  # Conversion to matrix simplifies comparison to gene names in the expression data frame
gene.sigs.mat <- gene.sigs.mat[gene.sigs.mat!=""]
  # When signature lengths are unequal, blanks ("") are inserted into the list that must be removed
if(sum(gene.sigs.mat %in% colnames(exp.dat)) == length(gene.sigs.mat)){
  print("All genes in signature(s) match those in expression data frame.")
}else{
  print(str_c(unique(gene.sigs.mat[!(gene.sigs.mat %in% colnames(exp.dat))]), collapse = ", "))
  stop("The above genes were NOT found in data frame. Please alter signature to match nomenclature in data frame and try again.")
}


## Convert signature and expression data to GSVA-appropriate formats

# Transpose the expression data frame and convert to a matrix
exp.dat.mat <- as.matrix(t(exp.dat))
  #Genes by row, cell lines by column
# Convert the gene signature data frame to a list of lists and remove any blanks
gene.sigs.list <- lapply(gene.sigs, as.list)
gene.sigs.list <- lapply(gene.sigs.list, function(x){x[x!=""]})


## Run GSVA across all signatures and check signature performance

# Run the GSVA
gsva.scores <- gsva(exp.dat.mat, gene.sigs.list, verbose=TRUE,kcdf="Gaussian")
# Export a summary file that contains cell-line level GSVA scores for each unique signature
write.csv(gsva.scores,"sig2target_gsva.scores.csv")


## Check "performance" of each gene in each signature by running Pearson correlations between each GSVA-score set and the gene expression data

# Create a placeholder for the summary data
gsva.vs.exp.summary <-as.data.frame(matrix(NA, ncol = 5, nrow = 1))
colnames(gsva.vs.exp.summary) <- c('Pearson.Cor','P.value','FDR','Gene','Signature')

# Use a for loop to call each GSVA signature one by one
for(i in 1:nrow(gsva.scores)){
  # And sapply to call the expression data of each gene one by one to run a pearson correlation  and store the pearson scores and p-values
  gsva.exp.cor <- sapply(exp.dat, function(x){
    test.out <- cor.test(gsva.scores[i,], x, method = "pearson")
    test.key <- c(test.out$estimate, test.out$p.value)
    return(test.key)
  })
  #Convert to a data frame where columns are variables
  gsva.exp.cor.df <- as.data.frame(t(gsva.exp.cor))
  #Label the initial data frame
  colnames(gsva.exp.cor.df) <- c('Pearson.Cor','P.value')
  #Calculate the FDR (adjusted p-value) for each gene and add it to a new column
  gsva.exp.cor.df$FDR <- p.adjust(gsva.exp.cor.df$`P.value`)
  #Trim the data to only include genes that were present in the signature
  gsva.exp.cor.df <- gsva.exp.cor.df[rownames(gsva.exp.cor.df) %in% gene.sigs[,i],]
  # Move the gene names to a column and add a column to identify the corresponding GSVA signature
  gsva.exp.cor.df <- gsva.exp.cor.df %>% mutate(Gene = rownames(gsva.exp.cor.df), Signature = rownames(gsva.scores)[i])
  # Bind the single-signature data frame to the "master" data frame
  gsva.vs.exp.summary <- bind_rows(gsva.vs.exp.summary, gsva.exp.cor.df)
}
# Trim off the first row (leftover from the placeholder) and reorganize the data
gsva.vs.exp.summary <- gsva.vs.exp.summary[-1,]
gsva.vs.exp.summary <- gsva.vs.exp.summary %>% select(c(5,4,1:3)) %>% arrange(Signature, FDR)
# Export the complete summary table
write.csv(gsva.vs.exp.summary,"sig2target_gsva.vs.exp.csv",row.names = FALSE)


## Import the gene dependency data, do a quality check and calculate required summary objects
  ## NOT IN FUNCTION: import steps, which are part of the function `sigCheck` or done manually by the user
  ## Example data objects will be included in the package
dep.dat <- read.csv("CRISPR_gene_effect.csv", header = TRUE, row.names = 1)
colnames(dep.dat) <- word(colnames(dep.dat), 1, sep=fixed('..'))
# Check the structure of the expression data frame
if(any(sapply(dep.dat, is.numeric) == FALSE)){
  stop("Third parameter, 'dep.dat', must be a data frame with ONLY numeric values.")
}
#Calculate the mean dependency of each gene
avg.dep <- sapply(dep.dat, mean)


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
write.csv(gsva.vs.dep.summary,"sig2target_gsva.vs.dep.csv",row.names = FALSE)


## Compile the key summary data (GSVA scores, expression correlations, dependency correlations) into a list and return
summary.list <- list(gsva.scores, gsva.vs.exp.summary, gsva.vs.dep.summary)
names(summary.list) <- c("gsva.scores","gsva.vs.exp","gsva.vs.dep")
return(summary.list)
  # Signatures (or individual gene expression profiles) can be run against a target gene to compare their predictive accuracy using `sigAccuracy`



