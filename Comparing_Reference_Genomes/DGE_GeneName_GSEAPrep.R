####DIFFERENTIAL GENE EXPRESSION ANALYSIS PREPARATION: GENE NAME MERGING

#In this script, we will merge the output file from our DEG analysis, in this case using the platform Limma-Voom, with 
#   the provided table translating the annotation FUN-ID's with the actual gene names from NCBI Blast databases. 
#INPUT: DGE (Limma-Voom) .csv file + NCBI Blast FUN-ID's and corresponding Gene Names .csv file 
#OUTPUT: DGE .csv file(s) of results form DGE (Limma-Voom) with their corresponding gene names; 
#   .csv file(s) with ranked genes based on log2fold changes and p-values 


### Merge 'gene names' with DGE results by Gene Model (gene_id)
#NOTE: for this to work, the DGE results .csv files must have a "gene_id" column header above the FUN-ID's
#NOTE: make sure to start this script in the proper working directory that holds all specified input files needed. 

## Import Annotation file with results from Blast to databases
Anno <- read.csv("Daphnia_pulex.annotations_Name.csv", stringsAsFactors = FALSE, na.strings=c(""))
summary(Anno)
dim(Anno)
#For this file, there are 30,370 Fun_IDs and their corresponding gene names, when known. 

## Import the DGE results file; make sure the gene model name (FUN_ID's) is 'gene_id' to match annotation file
DGEresults0 <- read.csv("PA42_limma-voom_MP_vsDose0.csv", stringsAsFactors = FALSE)
summary(DGEresults0)
dim(DGEresults0)

DGEresults102 <- read.csv("PA42_limma-voom_MP_vsDose102.csv", stringsAsFactors = FALSE)
summary(DGEresults102)
dim(DGEresults102)

DGEresults408 <- read.csv("PA42_limma-voom_MP_vsDose408.csv", stringsAsFactors = FALSE)
summary(DGEresults408)
dim(DGEresults408)

## Merge annotation file that contains the Gene names with DGE results; merging occurs by gene_id (the first column). 
# The all.y=FALSE will result in only the rows that are in common between the two files
DGE_Anno_0 <- merge(Anno,DGEresults0,by="gene_id",all.y = FALSE)
dim(DGE_Anno_0)
summary(DGE_Anno_0)

write.csv(as.data.frame(DGE_Anno_0), file="DGE_results_Dos0_GeneName.csv", row.names=FALSE) 

DGE_Anno_102 <- merge(Anno,DGEresults102,by="gene_id",all.y = FALSE)
dim(DGE_Anno_102)
summary(DGE_Anno_102)

write.csv(as.data.frame(DGE_Anno_102), file="DGE_results_Dos102_GeneName.csv", row.names=FALSE) 

DGE_Anno_408 <- merge(Anno,DGEresults408,by="gene_id",all.y = FALSE)
dim(DGE_Anno_408)
summary(DGE_Anno_408)

write.csv(as.data.frame(DGE_Anno_408), file="DGE_results_Dos408_GeneName.csv", row.names=FALSE)


############################# Make ranked list for GSEA ####################
## Here we are calculating a rank for each gene, and adding that column
## The rank is based on the log2fold change (how up or down regulated they are) and the pvalue (significance) 
DGE_Anno_Rank_0 <-  within(DGE_Anno_0, rank <- sign(logFC) * -log10(P.Value))
DGE_Anno_Rank_0 

DGE_Anno_Rank_102 <-  within(DGE_Anno_102, rank <- sign(logFC) * -log10(P.Value))
DGE_Anno_Rank_102 

DGE_Anno_Rank_408 <-  within(DGE_Anno_408, rank <- sign(logFC) * -log10(P.Value))
DGE_Anno_Rank_408 

#subset the results so only Gene Name and rank
DGErank_0 = subset(DGE_Anno_Rank_0, select = c(Name,rank) )
DGErank_0

DGErank_102 = subset(DGE_Anno_Rank_102, select = c(Name,rank) )
DGErank_102

DGErank_408 = subset(DGE_Anno_Rank_408, select = c(Name,rank) )
DGErank_408

#subset the results so only rows with a Name and rank
DGErank_withName_0 <- na.omit(DGErank_0)
DGErank_withName_0
dim(DGErank_withName_0)

DGErank_withName_102 <- na.omit(DGErank_102)
DGErank_withName_102
dim(DGErank_withName_102)

DGErank_withName_408 <- na.omit(DGErank_408)
DGErank_withName_408
dim(DGErank_withName_408)

## for use in GSEA pre-ranked analysis, the file must be tab delimitied and end in .rnk
  #Otherwise Could Use: write.csv(as.data.frame(DGErank_withName), file="DGErankName.csv", row.names=FALSE)

write.csv(as.data.frame(DGErank_withName_0), file="DGErankName0.csv", row.names=FALSE) 

write.csv(as.data.frame(DGErank_withName_102), file="DGErankName102.csv", row.names=FALSE) 

write.csv(as.data.frame(DGErank_withName_408), file="DGErankName408.csv", row.names=FALSE) 

write.table(DGErank_withName_0, file="DGErankName_0.rnk", sep = "\t", row.names=FALSE)  

write.table(DGErank_withName_102, file="DGErankName_102.rnk", sep = "\t", row.names=FALSE) 

write.table(DGErank_withName_408, file="DGErankName_408.rnk", sep = "\t", row.names=FALSE)  

##### We will load packages in order to set up the workspace
# Function to check for a package host on CRAN, then install (if needed) and library the package
prep_cranpack <- function (x){
  if (!requireNamespace(x, quietly = TRUE)) {
    install.packages(x)
    library(x, character.only = TRUE, quietly = TRUE)
  } else {
    library(x, character.only = TRUE, quietly = TRUE)
  }}

# Function to check for a package host on bioconductor, then install (if needed) and library the package
prep_biocpack <- function (x){
  if (!requireNamespace(x, quietly = TRUE)) {
    BiocManager::install(x)
    library(x, character.only = TRUE, quietly = TRUE)
  } else {
    library(x, character.only = TRUE, quietly = TRUE)
  }}

# loading list of CRAN packages
cranpacks <- c("BiocManager", "tools", "devtools", "tidyverse", "RColorBrewer", "stringr", "msigdbr")
invisible(lapply(cranpacks, prep_cranpack))

#BP-merge the Schwartz Lab annotation file with the additional gene names generated by Biomart
geneAnnoBP <- read.csv("~/Desktop/BPanno.csv")
biomart <- read.csv("~/Desktop/Biomartgenenames.csv")

resultAnnoBP <- merge(geneAnnoBP , biomart , by = "Gene.stable.ID", all = T)
head(resultAnnoBP)

write.csv(resultAnnoBP, "~/Desktop/resultAnnoBP.csv", row.names=FALSE)

#Merge the Schwartz Lab annotation file with Se_treatment 0 rank results 
resultsRanked0 <- read.csv("~/Desktop/DGErankName0.csv")
mergeresultAnnoBP <- read.csv("~/Desktop/mergeresultAnnoBP.csv")

resultAnno0 <- merge(resultsRanked0,mergeresultAnnoBP, by = "Name", all = T)
head(resultAnno0)

write.csv(resultAnno0, "~/Desktop/resultAnno.csv", row.names=FALSE)

# Merge the Schwartz Lab annotation file with Se_treatment 102 rank results
resultsRanked102 <- read.csv("~/Desktop/DGErankName102.csv")
mergeresultAnnoBP <- read.csv("~/Desktop/mergeresultAnnoBP.csv")

resultAnno102 <- merge(resultsRanked102,mergeresultAnnoBP, by = "Name", all = T)
head(resultAnno102)

write.csv(resultAnno102, "~/Desktop/resultAnno102.csv", row.names=FALSE)

#Merge the Schwartz Lab annotation with Se_treatment 408 rank results
resultsRanked408 <- read.csv("~/Desktop/DGErankName408.csv")
mergeresultAnnoBP <- read.csv("~/Desktop/mergeresultAnnoBP.csv")

resultAnno408 <- merge(resultsRanked408,mergeresultAnnoBP, by = "Name", all = T)
head(resultAnno408)

write.csv(resultAnno408, "~/Desktop/resultAnno408.csv", row.names=FALSE)

# drop NA's, select desired columns (rank and gene names), filter out missing gene names, and remove duplicated rows
library(dplyr)
#For 0 Se_treatment
rankFile0 <- resultAnno0 %>% drop_NA() %>% dplyr::select(rank, Name) %>% filter(Name != "") %>% distinct()
write.csv(rankFile0, "~/Desktop/rankFile0.csv", row.names=FALSE)

#For 102 Se_treatment
rankFile102 <- resultAnno102 %>% drop_NA() %>% dplyr::select(rank, Name) %>% filter(Name != "") %>% distinct()
write.csv(rankFile102, "~/Desktop/rankFile102.csv", row.names=FALSE)

#For 408 Se_treatment
rankFile408 <- resultAnno408 %>% drop_NA() %>% dplyr::select(rank, Name) %>% filter(Name != "") %>% distinct()
write.csv(rankFile408, "~/Desktop/rankFile408.csv", row.names=FALSE)

#For the original annotation of the PA42 genome comparison
PA42anno <- read.csv("~/Desktop/Daphnia_pulex.annotations_Name.csv")

# Merge PA42 annotation file with Se_treatment 102 rank results
PA42resultAnno102 <- merge(resultsRanked102,PA42anno, by = "Name", all = T)
head(PA42resultAnno102)

write.csv(PA42resultAnno102, "~/Desktop/PA42resultAnno102.csv", row.names=FALSE)

#Merge PA42 annotation file with Se_treatment 0 rank results 
PA42resultAnno0 <- merge(resultsRanked0,PA42anno, by = "Name", all = T)
head(PA42resultAnno0)

write.csv(PA42resultAnno0, "~/Desktop/PA42resultAnno.csv", row.names=FALSE)

#Merge PA42 annotation with Se_treatment 408 rank results
PA42resultAnno408 <- merge(resultsRanked408,PA42anno, by = "Name", all = T)
head(PA42resultAnno408)

write.csv(PA42resultAnno408, "~/Desktop/PA42resultAnno.csv", row.names=FALSE)

# drop NA's, select desired columns (rank and gene names), filter out missing gene names, and remove duplicated rows
#For 0 Se_treatment
PA42rankFile0 <- PA42resultAnno0 %>% drop() %>% dplyr::select(rank, Name) %>% filter(Name != "") %>% distinct()

write.csv(PA42rankFile0, "~/Desktop/PA42rankFile0.csv", row.names=FALSE)

#For 102 Se_treatment
PA42rankFile102 <- PA42resultAnno102 %>% drop_NA() %>% dplyr::select(rank, Name) %>% filter(Name != "") %>% distinct()

write.csv(PA42rankFile102, "~/Desktop/PA42rankFile102.csv", row.names=FALSE)

#For 408 Se_treatment
PA42rankFile408 <- PA42resultAnno408 %>% drop() %>% dplyr::select(rank, Name) %>% filter(Name != "") %>% distinct()

write.csv(PA42rankFile408, "~/Desktop/PA42rankFile408.csv", row.names=FALSE)

####### #At this point everything is prepped to be loaded into the GSEA software (v4.3.2) for analysis.

#Not completed: 
####  We also need the normalized count data. Here we are going back to the dds object
nt <- normTransform(dds) # defaults to log2(x+1)
head(assay(nt))
# compare to original count data
head(countdata)

# make it a new dataframe
NormTransExp<-assay(nt)
summary(NormTransExp)
head(NormTransExp)

gene_id <-gsub("^[^-]+-", "", rownames(NormTransExp))
NormTransExpIDs  <-cbind(gene_id,NormTransExp)
head(NormTransExpIDs)

#merge with name
Name_NormTransExp <- merge(Anno,NormTransExpIDs,by="gene_id",all.y = FALSE)
dim(Name_NormTransExp )
summary(Name_NormTransExp)

write.csv(as.data.frame(Name_NormTransExp), file="NormTransExpressionData.csv", row.names=FALSE)