#[edgeR]: http://www.bioconductor.org/packages/release/bioc/html/edgeR.html

#install
#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

#BiocManager::install("edgeR")
BiocManager::install("limma")

library("edgeR")

data_raw <- read.csv("gene_count_matrix.csv", header = TRUE)
