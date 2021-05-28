# BenchCoexp: automated benchmarking of association measures for coexpression analysis.

This is an R and C++ (through the Rcpp and RcppParallel environment) bundle of functions in order to evaluate the performance of any user provided association measure (e.g. Pearson Correlation Coefficient) at the task of predicting cofunctionality between gene pairs in a coexpression analysis context. More specifically, any {measure of association-expression dataset-geneset terms dataset} combination may be evaluated. Both the expression dataset and the terms dataset (e.g. A group of Gene Ontology terms) are user provided.

## Algorithm
The algorithm (`Rsrc/benchCoexp.r`) takes as inputs a gene expression matrix (genes by samples), a gene-to-gene term similarity vector (there is a helper function that takes care of producing this for the user given a genes by terms matrix with a 1 if a gene is part of a term and 0 otherwise), a function that computes the association measure, a preprocessing function for the association measure (if any) and the number of minimum required terms in common for two genes to be considered as cofunctioning.

The algorithm steps are as follows:
1. Apply the preprocessing function (if any) to the expression matrix. For example, we could use this if we need to discretize the gene vectors prior to a Mutual Information-based association measure calculation.
2. Calculate the co-association (co-expression) matrix of all possible combinations of two genes in the data using the association function provided.
3. For each gene `i`, rank the associations of `i` with every other possible gene `j`.
4. Compute the mutual rank (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2762411/) matrix. This will allow for different measures of associations to be comparable independently of their scale.
5. Conduct the benchmarking process (implemented with parallel computing). The algorithm used is the same as used in for example: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1700-9 "evaluation of co-expression networks", methods description.

The final return value is a 3-dimensional array with `s` slices. Each slice is a confusion matrix for a particular operation point (Mutual Rank value). At each operation point, this value is used as a threshold to predict if wether or not a pair of genes is cofunctional. Once all pairs of genes are predicted, the prediction is checked against the "Ground Truth" given by the terms data. This array of confusion matrices can be used in downstream analyses to determine the performance of the {association measure-expression dataset-terms dataset} combination in a similar fashion to what is done when evaluating a classifier in a Machine Learning binary classification task.

The file `example.r` walks through the usage of the program in detail. The final result (which in this case evaluates the measure of association using Area Under the Recieving Operator Curve) should look something like this:

![](demo.png?raw=true)
