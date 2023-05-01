#script taken from https://rnnh.github.io/bioinfo-notebook/docs/DE_analysis_edgeR_script.html

#if (!requireNamespace("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")

#BiocManager::install("edgeR")

#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

#BiocManager::install("limma")

library(limma)
library(edgeR)

#read gene counts file
counts.df <- read.csv("gene_count_matrix.csv")

#make gene IDs as the row name and remove the first column (which is gene IDs)
rownames(counts.df) <- counts.df$gene_id
counts.df$gene_id <- NULL

#load experiment metadata
design.df <- read.table("Pheno_Data.txt", sep = "\t", header=TRUE)

counts_control_old.df  <- counts.df[,c("SRR6853319", "SRR6853326", "SRR6853335")]
counts_low_old.df  <- counts.df[,c("SRR6853322", "SRR6853324", "SRR6853325")]
counts_high_old.df  <- counts.df[,c("SRR6853320", "SRR6853321", "SRR6853323")]

counts_control_young.df  <- counts.df[,c("SRR6853328", "SRR6853331", "SRR6853334")]
counts_low_young.df  <- counts.df[,c("SRR6853327", "SRR6853330", "SRR6853333")]
counts_high_young.df  <- counts.df[,c("SRR6853329", "SRR6853332", "SRR6853336")]

#Testing standard deviation
#Define the function RSD.test
RSD.test <- function(dataframe){
  # This function tests whether the relative standard deviation (RSD) is less
  # than or equal to one for each row in a data frame.
  # It adds the result to a new variable in the data frame called "RSD.test".
  # For a given row, if data.frame$RSD.test is TRUE, that row has an RSD less
  # than or equal to one, i.e. RSD <= 1.
  # If data.frame$RSD.test is FALSE, that row has an RSD outside of this range.
  RSD_tests = dataframe[,1]
  for (row_index in 1:nrow(dataframe)){
    row = as.numeric(dataframe[row_index,])
    RSD = sd(row) / mean(row)
    RSD_tests[row_index] = RSD <= 1 || is.na(RSD)
  }
  dataframe$RSD.test <- as.factor(RSD_tests)
  levels(dataframe$RSD.test) <- c(FALSE, TRUE)
  return(dataframe)
}


#Run RSD.test on all
counts_control_old.df <- RSD.test(counts_control_old.df)
counts_low_old.df  <- RSD.test(counts_low_old.df)
counts_high_old.df  <- RSD.test(counts_high_old.df)

counts_control_young.df  <- RSD.test(counts_control_young.df)
counts_low_young.df  <- RSD.test(counts_low_young.df)
counts_high_young.df <- RSD.test(counts_high_young.df)

# Creating list of genes which failed RSD test
RSD_failed_genes <- rownames(counts_control_old.df[
  which(counts_control_old.df$RSD.test == FALSE),])

RSD_failed_genes <- append(RSD_failed_genes, rownames(counts_low_old.df[
  which(counts_low_old.df$RSD.test == FALSE),]))

RSD_failed_genes <- append(RSD_failed_genes, rownames(counts_high_old.df[
  which(counts_high_old.df$RSD.test == FALSE),]))

RSD_failed_genes <- append(RSD_failed_genes, rownames(counts_control_young.df[
  which(counts_control_young.df$RSD.test == FALSE),]))

RSD_failed_genes <- append(RSD_failed_genes, rownames(counts_low_young.df[
  which(counts_low_young.df$RSD.test == FALSE),]))

RSD_failed_genes <- append(RSD_failed_genes, rownames(counts_high_young.df[
  which(counts_high_young.df$RSD.test == FALSE),]))

RSD_failed_genes <- unique(RSD_failed_genes)
length(RSD_failed_genes)


# Filtering gene counts
filtered_counts.df <- counts.df[
  which(!rownames(counts.df) %in% RSD_failed_genes),]

# Checking that gene counts were correctly filtered
nrow(counts.df) - length(RSD_failed_genes) == nrow(filtered_counts.df)

# Creating a DGEList object using the filtered gene counts
counts.DGEList <- DGEList(counts = filtered_counts.df,
                          genes = rownames(filtered_counts.df))

# Confirming samples are in the same order in the gene counts and design table
summary(colnames(filtered_counts.df) == design.df$sample)

# Add grouping information to DGEList object
counts.DGEList$samples$group <- as.factor(design.df$treatment)
counts.DGEList
dim(counts.DGEList)

# Creating an object to filter genes with low expression
counts.keep <- filterByExpr(counts.DGEList)
summary(counts.keep)

# Filtering lowly expressed genes
counts.DGEList <- counts.DGEList[counts.keep, , keep.lib.sizes = FALSE]
dim(counts.DGEList)

# Confirming that the number of genes in counts.DGEList is the same as the
# number of TRUE values in counts.keep
length(counts.keep[counts.keep == TRUE]) == dim(counts.DGEList)[1]

# Printing the normalisation factors for the libraries
counts.DGEList$samples$norm.factors

# Calculating normalisation factors and applying them to counts.DGEList
counts.DGEList <- calcNormFactors(counts.DGEList)
counts.DGEList$samples$norm.factors

# Estimating common dispersion and tagwise dispersion
condition_ <- design.df$treatment
counts.DGEList <- estimateDisp(counts.DGEList, design = model.matrix(~condition_))
counts.DGEList

#View groups available for pairwise comparison
condition_

# Exact tests for differences between experimental conditions
std_Old_HighDosevsControl.DGEExact <- exactTest(counts.DGEList, pair = c("Control Old Mom",
                                                             "High Dose Old Mom"))

topTags(std_Old_HighDosevsControl.DGEExact)
#p.adj <- p.adjust(std_Old_HighDosevsControl.DGEExact, method="fdr")

#write.csv(std_Old_HighDosevsControl.DGEExact, file="Old_HighDosevsControl.csv") 


##Summary of DEGs based on FDR and MA plot

de1 <- decideTestsDGE(std_Old_HighDosevsControl.DGEExact, adjust.method="BH", p.value=0.05)
summary(de1)
de1tags12 <- rownames(de1)[as.logical(de1)] 
plotSmear(std_Old_HighDosevsControl.DGEExact, de.tags=de1tags12)
abline(h = c(-1, 1), col = "blue")

std_Old_LowDosevsControl.DGEExact <- exactTest(counts.DGEList, pair = c("Control Old Mom",
                                                                         "Low Dose Old Mom"))
#write.csv(std_Old_LowDosevsControl.DGEExact, file="Old_LowDosevsControl.csv") 

std_Old_HighDosevsLowDose.DGEExact <- exactTest(counts.DGEList, pair = c("Low Dose Old Mom",
                                                                        "High Dose Old Mom"))
#write.csv(std_Old_HighDosevsLowDose.DGEExact, file="Old_HighDosevsLowDose.csv") 



std_Young_HighDosevsControl.DGEExact <- exactTest(counts.DGEList, pair = c("Control Young Mom",
                                                                         "High Dose Young Mom"))
#write.csv(std_Young_HighDosevsControl.DGEExact, file="Young_HighDosevsControl.csv") 

de2 <- decideTestsDGE(std_Young_HighDosevsControl.DGEExact, adjust.method="BH", p.value=0.05)
summary(de2)
de2tags12 <- rownames(de2)[as.logical(de2)] 
plotSmear(std_Young_HighDosevsControl.DGEExact, de.tags=de2tags12)
abline(h = c(-2, 2), col = "blue")

std_Young_LowDosevsControl.DGEExact <- exactTest(counts.DGEList, pair = c("Control Young Mom",
                                                                        "Low Dose Young Mom"))
#write.csv(std_Young_LowDosevsControl.DGEExact, file="Young_LowDosevsControl.csv") 

std_Young_HighDosevsLowDose.DGEExact <- exactTest(counts.DGEList, pair = c("Low Dose Young Mom",
                                                                         "High Dose Young Mom"))
#write.csv(std_Young_HighDosevsLowDose.DGEExact, file="Young_HighDosevsLowDose.csv") 

