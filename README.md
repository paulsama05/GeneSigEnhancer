
# GeneSigEnhancer

<!-- badges: start -->
<!-- badges: end -->

This package uses user-defined gene signatures, gene expression data (e.g. CCLE), and gene dependency data (e.g. DepMap) to explore the relationship between gene signatures and target genes.

## Installation

Installation is done in R using the following (or similar installation functions):

``` r
devtools::install_github("paulsama05/GeneSigEnhancer")
```


## Contents


This package contains 4 functions and 1 example data frame object, which can be viewed using:
``` r
library(GeneSigEnhancer)
ls("package:GeneSigEnhancer")
```
The functions can be more specifically viewed with:
``` r
library(GeneSigEnhancer)
lsf.str("package:GeneSigEnhancer")
```
For more details, please see the help documentation for each object.
 
The data frame object is derived from a raw .csv file, which can be retrieved using:
``` r
system.file("extdata","gene_signature_set.csv",package = "GeneSigEnhancer")
```
There are four other retrievable files: 1 raw .csv (gene_signature.csv) and 4 annotated .R files that describe how the functions work (sig2target_annotated.R, sigAccuracy_annotated.R, sigImport_annotated.R, target2sig_annotated.R). These can all be accessed using the `system.file(...)` annotation above.


## Requirements

The "GSVA" package must be installed from Bioconductor. The `sig2target` and `target2sig` functions will check for this and attempt to install it, but it is preferred to install beforehand.

For the `sigImport` function, all .csv files must be in the current working directory.

Gene signatures must have a header (signature name) and all genes must use a nomenclature consistent with that found in expression and dependency data sets. Example data is provided as an object and raw .csv file (see above), but user can provide their own.

Neither expression nor dependency data sets are provided in this package; it is recommended to download expression and dependency data as .csv files from DepMap, which will ensure consistency between gene and sample names. These files must be formatted with genes as column names and samples (e.g. cell line IDs) as rownames.


## Examples

The best way to understand the role of each function is to use them sequentially.

Although the user can do it manually, it is recommended to use `sigImport` to load the initial gene signature, gene expression and gene dependency files. For example:

``` r
library(GeneSigEnhancer)
sig.list <- sigImport()
```

For the purposes of this example, you should use the "gene_signature_set.csv" file in the extdat folder of this package (must be copied to working directory). After entering the complete target file names (without quotation marks), a list is assigned to the object `sig.list`. Take a peek in the list using:

``` r
sig.list$gene.sigs    # Data frame of gene signatures. If you used "my" data this should have 5 signatures labeled "A"-"E"
sig.list$exp.dat[1:10,1:10] # Data frame of expression data
sig.list$dep.dat[1:10,1:10] # Data frame of dependency data
```

These objects can then be used in `sig2target`:

``` r
s2t.out <- sig2target(sig.list$gene.sigs, sig.list$exp.dat, sig.list$dep.dat)
```

This function runs GSVA between each gene signature (gene.sigs) and the complete expression data set (exp.dat) then runs Pearson correlations between the resulting GSVA scores and the dependency data set (dep.dat).

Along the way .csv files are written to the working directory (all start with "sig2target") and a list is returned to the object `s2t.out`. Take a peek in the list using:

``` r
s2t.out$gsva.scores   # matrix of GSVA scores corresponding to each input signature
s2t.out$gsva.vs.exp   # data frame of Pearson correlation data between the GSVA scores and the expression data frame. Subset to only genes found in the input signatures so you can see which genes contributed most to the GSVA scores.
s2t.out$gsva.vs.dep   # data frame of Pearson correlation data between the GSVA scores and the dependency data frame. Subset to only genes that significantly correlate. Top hits (lowest FDRs) are good target candidates.
```

The output from `s2t.out$gsva.vs.dep` should give you a good idea of a target gene to go after. This gene can then be used to create new gene signatures with `target2sig`. For example, if you use the signature data frame provided in this package "ADAR" should be a top hit:

``` r
t2s.out <- target2sig(sig.list$dep.dat, gene = "ADAR", sig.list$exp.dat, cutoff = 100, gene.sig = sig.list$gene.sigs$A)
```

Note that I assigned `gene.sig = sig.list$gene.sigs$A`. This pulls a single signature list ("A") from sig.list$gene.sigs; in practice, you will have to know the name of your lists. Furthermore, the `gene.sig` field can be left empty

Also note that although I assigned a value for `cutoff`, this already defaults to 100.

Ultimately, this function runs a Pearson correlation between the dependency data of the target gene ("ADAR" here) and the expression data set, identifies the most highly correlated genes, and then keeps the top 100 (`cutoff`) as a new signature. Because we define `gene.sig`, the function creates 2 more signatures by further subsetting the full list of target-correlated genes AND the reduced list of 100 top genes.

Along the way .csv files are written to the working directory (all start with "target2sig") and a list is returned to the object `t2s.out`. Take a peek in the list using:

``` r
t2s.out$new.sigs   # list of new gene signatures. Includes 1) 100 top correlated genes, 2) the `gene.sig` full input signature, 3) the genes that overlap between `gene.sig` and ALL target-correlated genes, and 4) the genes that overlap between `gene.sig` and the top 100 correlated genes
t2s.out$gsva.scores   # matrix of GSVA scores corresponding to each input signature
t2s.out$gsva.vs.dep   # data frame of Pearson correlation data between the GSVA scores and the dependency data frame. Subset to only genes that significantly correlate. Expect to see input candidate gene near top of each list.
```

Scan the results of `t2s.out$gsva.vs.dep` to see that "ADAR" is near the top of the list. Also, see which list gives you the best Pearson correlation. In this example, it should be the "top.overlap" signature, which we can now take a deeper look at in `sigAccuracy`:

``` r
target.summary <- sigAccuracy(t2s.out$gsva.scores["top.overlap",], sig.list$dep.dat, "ADAR")
```

In this final function, we run a linear model between the "top.overlap" GSVA score data and the ADAR dependency data. By using the slope + intercept outputs we can predict the ADAR dependency of each sample/cell line and then compare to the actual value.

As above, .csv files are written to the working directory (all start with "sigAcc") and a list is returned to the object `target.summary`. Take a peek in the list using:

``` r
target.summary$sample.predictions   # A complete data frame documenting whether the input scores accurately predict the target gene dependency of each individual sample
target.summary$prediction.table   # A summary table of the above data documenting the percentage of True/False Positives/Negatives.
target.summary$gg.scatter   # A ggplot scatter plot comparing the input scores (GSVA here) and the dependency data. 
target.summary$gg.windows   # A ggplot scatter plot showing the predictive accuracy of each "window" of input scores
```

Note that higher input (GSVA) scores should correlate to higher dependence (lower DepMap score). In the final graph, we would want to see high predictive accuracy to the right of the vertical line (which indicates the border between "independent" and "dependent" samples). Unfortunately, the GSVA scores in this range are poor predictors of dependency!
