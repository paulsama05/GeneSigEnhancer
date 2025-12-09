#' Calculates how well sample scores predict target gene dependency
#'
#' User provides sample scores, target gene, and gene dependency data. Function compares predicted and actual dependencies for each sample and outputs both summary data and graphs. Data frame outputs are also written to files starting with "sigAcc".
#' @param scores.in Vector of scores (e.g. GSVA, expression data. etc.) with corresponding sample names (must be formatted like sample names in dependency data set).
#' @param dep.dat Data frame of gene dependency data organized with genes in columns, samples (e.g. cell lines) in rows.
#' @param gene Name of target gene. Must exactly match that in dependency data set.
#' @param sig.cutoff  Numerical value that sets the cutoff between "independent" and "dependent" values. Default is -0.5 (DepMap standard).
#' @param window Sets the window/chunk size used to scan sample scores and accuracy values that are used in final graph. Larger numbers return "smoother" graph, but less information.
#' @return A list 4 objects: 1) Data frame detailing scores and dependency outcomes for each sample, 2) summary data frame of prediction accuracy (True/false positive/negative), 3) GGplot scatter plot of scores vs dependency data, and 4) GGplot scatter plot of scores (by window) vs predictive accuracy.
#' @description
#' Runs linear model between sample scores and dependency values for target gene and uses slope + intercept data to predict dependency values for each sample. Compares predicted and actual dependency values to determine accuracy.
#'
#' @export
sigAccuracy <- function(scores.in, dep.dat, gene, sig.cutoff = -0.5, window = 100){
library(dplyr)
library(ggplot2)
if(!is.numeric(scores.in) | !is.vector(scores.in)){
  stop("`scores.in` must be a single vector of numerical values.")
}
if(is.null(names(scores.in))){
  stop("`scores.in` vector must contain `names`, e.g. sample names or cell line identifiers. ")
}
dep.dat <- sig.list$dep.dat
if(any(sapply(dep.dat, is.numeric) == FALSE)){
  stop("`dep.dat` must be a data frame with ONLY numeric values.")
}
gene <- "ADAR"
if(!(gene %in% colnames(dep.dat))){
  stop(paste0("Provided `gene`, \"",gene,"\", not found in dependency data set, `dep.dat`."))
}
sig.cutoff <- -0.5
if(!is.numeric(window)){
  stop("`sig.cutoff` must be a numerical value. Default is -0.5")
}
gene.dep <- dep.dat[[gene]]
names(gene.dep) <- rownames(dep.dat)
scores.in.trim <- scores.in[names(scores.in) %in% names(gene.dep)]
scores.in.ordered <- scores.in.trim[order(names(scores.in.trim))]
gene.dep.trim <- gene.dep[names(gene.dep) %in% names (scores.in)]
gene.dep.ordered <- gene.dep.trim[order(names(gene.dep.trim))]
window <- 100
if(!(length(window) == 1) | !is.numeric(window) | window <1 | window > length(scores.in.trim)){
  stop("`window` must be a single positive number smaller than the number of samples common to expression and dependency data sets.\nEntries rounded down to integers.")
}
window <- as.integer(window)
lm.result <- lm(gene.dep.ordered~scores.in.ordered)
lm.intercept <- lm.result$coefficients[1]
lm.slope <- lm.result$coefficients[2]
predictions <- data.frame(scores.in.ordered, lm.slope*scores.in.ordered + lm.intercept, gene.dep.ordered)
names(predictions) <- c("Sample.Score", "Predicted.Dependency", "Actual.Dependency")
predictions <- arrange(predictions, Sample.Score)
predictions$Predicted.Binary <- "Independent"
predictions$Predicted.Binary[predictions$Predicted.Dependency<= sig.cutoff] <- "Dependent"
predictions$Actual.Binary <- "Independent"
predictions$Actual.Binary[predictions$Actual.Dependency<= sig.cutoff] <- "Dependent"
prediction.table <- predictions %>% group_by(Predicted.Binary,Actual.Binary) %>% summarize(Percentage = 100*n()/nrow(predictions))
write.csv(prediction.table,"sigAcc_prediction.table.csv", row.names = FALSE)
predictions$Is.Accurate <- FALSE
predictions$Is.Accurate[predictions$Predicted.Binary == predictions$Actual.Binary] <- TRUE
write.csv(predictions,"sigAcc_sample.predictions.csv")
scan.length <- nrow(predictions)-(window-1)
window.score <- rep(NA, scan.length)
window.accuracy <- rep(NA, scan.length)
for(i in 1:scan.length){
  window.score[i] <- mean(predictions$Sample.Score[i:(i+window-1)])
  window.accuracy[i] <- 100*sum(predictions$Is.Accurate[i:(i+window-1)])/window
}
window.df <- data.frame(window.score, window.accuracy)
pivot.pt <- (sig.cutoff - lm.intercept)/lm.slope
library(ggplot2)
gg.scatter <- ggplot(predictions, aes(Sample.Score, Actual.Dependency)) + geom_point(aes(colour = Actual.Binary)) + geom_vline(xintercept = pivot.pt) + geom_smooth(method = lm)
gg.windows <- ggplot(window.df, aes(window.score, window.accuracy)) + geom_point() + geom_vline(xintercept = pivot.pt) + ylim(c(0,100))
summary.list <- list(predictions, prediction.table, gg.scatter, gg.windows)
names(summary.list) <- c("sample.predictions","prediction.table","gg.scatter", "gg.windows")
return(summary.list)
}
