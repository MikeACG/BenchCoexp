# set working directory to the folder where this file is
# setwd()

# make sure these libraries are installed
lapply(c("Rcpp", "RcppParallel", "bigmemory", "data.table", "pracma"), require, character.only = TRUE)

# load source code
invisible(sapply(list.files("Rsrc/", full.names = TRUE), source))
invisible(sapply(list.files("CPPsrc/", full.names = TRUE), Rcpp::sourceCpp))

# set number of cores for parallel processes
RcppParallel::setThreadOptions(numThreads = floor(RcppParallel::defaultNumThreads() / 2))

# path to preprocessed (already normalized, corrected, etc.) expression data R object to benchmark
expPath <- "exampleExp.RData" # this is loaded inside the algorithm and not here for memory efficiency when working with many genes

# Load and parse geneset terms data
print(load("exampleTerms.RData")) # a subset of Gene Ontology Biological Process for this example
terms <- termSimil(terms) # this will convert the genes by terms matrix we just loaded to a gene-to-gene term similarity vector which is what is used for the benchmarking process
# NOTE: genes in the expression data and the terms data must be the same and appear in the same order. I've made sure this is true for the example data

# Source the association function to benchmark
Rcpp::sourceCpp("statG.cpp") # a G-statistic custom function as an example
f <- "statG" # indicate the name of the function as a string

# Source the preprocessing function
source("kmeansDisc.r") # will discretize the gene expression vectors before the G-statistic calculation
pref <- "kmeansDisc" # same as with the association function

# Indicate if wether or not the preprocessing function has a return value
prefret <- TRUE # if working with a large number of genes, setting this to FALSE in combination with a void cpp preprocessing function will optimize memory

# number of terms in common required for a gene pair to be considere d cofunctional
k <- 3 # increasing this means we require more evidence to consider a pair of genes as a "ground truth" positive

# Run the experiment
res <- benchCoexp(expPath, terms, f, k, pref, prefret)
res$timings # check how long each step of the algorithm took and total run time

# summarize result with AUC for example
tpr <- apply(res$cms, 3, function(m) m[1, 1] / (m[1, 1] + m[2, 1]))
fpr <- apply(res$cms, 3, function(m) m[1, 2] / (m[1, 2] + m[2, 2]))
auc <- pracma::trapz(rev(fpr), rev(tpr))
plot(fpr, tpr, type = "l", col = 3,
    xlab = "False Positive Rate", ylab = "True Positive Rate", main = "ROC Curve", sub = paste("AUC = ", signif(auc, 2))
)
abline(a = 0, b = 1, lty = 2, col = 8)
legend("topleft", c("Used association measure", "Expected at random"), lty = c(1, 2), col = c(3, 8))
