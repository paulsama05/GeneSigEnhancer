## Inputs:
  # sample scores (e.g. single GSVA score)
  # dependency data set
  # target gene name
  # `sig.cutoff` for dependent/independent categorization
  # `window` size for scanning



# Start with `scored` data across samples (e.g. GSVA scores or expression values)
  # E.g., from the target2sig output when a `gene.sig` is provided
scores.in <- t2s.out$gsva.scores["top.overlap",]
  # OR alternately, can use:
    #scores.in <- sig.list$exp.dat$MX1
    #names(scores.in) <- rownames(sig.list$exp.dat)
if(!is.numeric(scores.in) | !is.vector(scores.in)){
  stop("`scores.in` must be a single vector of numerical values.")
}
if(is.null(names(scores.in))){
  stop("`scores.in` vector must contain `names`, e.g. sample names or cell line identifiers. ")
}

# Designate the dependency data set
dep.dat <- sig.list$dep.dat
if(any(sapply(dep.dat, is.numeric) == FALSE)){
  stop("`dep.dat` must be a data frame with ONLY numeric values.")
}

# Designate the target gene
gene <- "ADAR"
if(!(gene %in% colnames(dep.dat))){
  stop(paste0("Provided `gene`, \"",gene,"\", not found in dependency data set, `dep.dat`."))
}

#Set cutoff for dependency score significance (default -0.5)
sig.cutoff <- -0.5
if(!is.numeric(window)){
  stop("`sig.cutoff` must be a numerical value. Default is -0.5")
}

# Subset the target gene's dependency data
gene.dep <- dep.dat[[gene]]
names(gene.dep) <- rownames(dep.dat)

# Trim and align the GSVA data and the dependency data by sample (i.e. cell line ID)
scores.in.trim <- scores.in[names(scores.in) %in% names(gene.dep)]
scores.in.ordered <- scores.in.trim[order(names(scores.in.trim))]
gene.dep.trim <- gene.dep[names(gene.dep) %in% names (scores.in)]
gene.dep.ordered <- gene.dep.trim[order(names(gene.dep.trim))]

# Set the scanning window size (default 100)
window <- 100
if(!(length(window) == 1) | !is.numeric(window) | window <1 | window > length(scores.in.trim)){
  stop("`window` must be a single positive number smaller than the number of samples common to expression and dependency data sets.\nEntries rounded down to integers.")
}
window <- as.integer(window)

# Run linear model and extract the slope and intercept information
lm.result <- lm(gene.dep.ordered~scores.in.ordered)
lm.intercept <- lm.result$coefficients[1]
lm.slope <- lm.result$coefficients[2]

# Use the input data (e.g. sample-wise GSVA scores) the slope and intercept to calculate "predicted" dependency values for each cell line
predictions <- data.frame(scores.in.ordered, lm.slope*scores.in.ordered + lm.intercept, gene.dep.ordered)
names(predictions) <- c("Sample.Score", "Predicted.Dependency", "Actual.Dependency")
predictions <- arrange(predictions, Sample.Score)

# Assign predicted dependencies a value of dependent (<= sig.cutoff) or independent (> sig.cutoff)
predictions$Predicted.Binary <- "Independent"
predictions$Predicted.Binary[predictions$Predicted.Dependency<= sig.cutoff] <- "Dependent"

# Add actual ADAR dependencies as dependent/independent
predictions$Actual.Binary <- "Independent"
predictions$Actual.Binary[predictions$Actual.Dependency<= sig.cutoff] <- "Dependent"

# Generate a summary table for model accuracy (True Positives, False Positives...)
prediction.table <- predictions %>% group_by(Predicted.Binary,Actual.Binary) %>% summarize(Percentage = 100*n()/nrow(predictions))
write.csv(prediction.table,"sigAcc_prediction.table.csv", row.names = FALSE)

# Assign the accuracy of each individual sample
predictions$Is.Accurate <- FALSE
predictions$Is.Accurate[predictions$Predicted.Binary == predictions$Actual.Binary] <- TRUE
write.csv(predictions,"sigAcc_sample.predictions.csv")

## Scan the data in windows/"chunks" to determine where predictive accuracy is highest
  #Determine how many unique windows are possible with a given window size
scan.length <- nrow(predictions)-(window-1)
  # Create a placeholder vector to store chunk-wise score data
window.score <- rep(NA, scan.length)
  # Create a placeholder vector to store chunk-wise accuracy data
window.accuracy <- rep(NA, scan.length)
  # Walk along the accuracy column of `predictions` on by one and and return the average score and percent accuracy for each new window
for(i in 1:scan.length){
  window.score[i] <- mean(predictions$Sample.Score[i:(i+window-1)])
  window.accuracy[i] <- 100*sum(predictions$Is.Accurate[i:(i+window-1)])/window
}
  # Merge into a data frame with the corresponding "start" input score for each window
window.df <- data.frame(window.score, window.accuracy)
# Calculate the input score that would return the sig.cutoff value (i.e. the "pivot" between dependent and independent values)
pivot.pt <- (sig.cutoff - lm.intercept)/lm.slope

# Make a scatter plot of the input data
library(ggplot2)
gg.scatter <- ggplot(predictions, aes(Sample.Score, Actual.Dependency)) + geom_point(aes(colour = Actual.Binary)) + geom_vline(xintercept = pivot.pt) + geom_smooth(method = lm)

# Use the window df and the pivot point to plot the data
gg.windows <- ggplot(window.df, aes(window.score, window.accuracy)) + geom_point() + geom_vline(xintercept = pivot.pt) + ylim(c(0,100))

# Combine key outputs into a list and return
summary.list <- list(predictions, prediction.table, gg.scatter, gg.windows)
names(summary.list) <- c("sample.predictions","prediction.table","gg.scatter", "gg.windows")
return(summary.list)
