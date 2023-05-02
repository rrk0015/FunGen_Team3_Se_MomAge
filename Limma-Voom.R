##Limma-Voom Protocol 
##BIOL5850 Final Project
##Jordan Chow

________________________________________________________________________________
#Set up

#Install edgeR package (installs limma as a dependency)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR")

# install the current version of biomartr on your system
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


#Load edgeR package
    library(edgeR)


#Read in the counts table
counts <- read.csv("gene_count_matrix.csv")
head(counts)
dim(counts)
View(counts)


#add this part to new genome data 
counts$gene_id  <- gsub("gene-", "", counts$gene_id)
#write to new file
write.csv(counts, file="corrected.csv") #Delete number column in excel before re-reading the file
#load new file
counts <- read.csv("corrected.csv", row.names="gene_id")
head(counts)
dim(counts)

#Create DDGEList object
         d0 <- DGEList(counts)
________________________________________________________________________________
#Preprocessing

#Calculate normalization factors
    d0 <- calcNormFactors(d0)
    d0
    
#filter low-expressed genes
    cutoff <- 1
    drop <- which(apply(cpm(d0), 1, max) < cutoff)
    d <- d0[-drop,] 
    dim(d) # number of genes left
##Sample names    
    snames <- colnames(counts) # Sample names
    snames
    
### Input the phenotype data
    #Note: The PHENO_DATA file contains information on each sample, e.g., sex or population. The exact way to import this depends on the format of the file.
    ##Make sure the individual names match between the count data and the metadata
    coldata <-(read.csv("Phenodata.csv", header=TRUE, row.names=1))
    dim(coldata)
    head(coldata)
##Check all sample IDs in colData are also in CountData and match their orders
    all(rownames(coldata) %in% colnames(counts))
    counts<- counts[, rownames(coldata)]
    all(rownames(coldata) == colnames(counts))
##format data 
  Maternalpop <- substr(coldata$MaternalPop,1,nchar(coldata$MaternalPop)) 
  Maternalpop
  Dose <-substr(coldata$SeMetDose,1,nchar(coldata$SeMetDose))
  Dose
  group <- interaction(Maternalpop,Dose)
  group

  
##MDS Plot
  plotMDS(d, col = as.numeric(group))
  
  
__________________________________________________________________________________________________________
#Voom transformation and calculation of variance weights
#fit to interact with Age and Dose factors 
mm1 <- model.matrix(~group, data = d)
rownames(mm1)
#voom: Mean-varience trend for filtered genes
y1 <- voom(d, design =  mm1, plot = T)
fit1 <- lmFit(y1, mm1)
fit1 <- eBayes(fit1)
head(coef(fit1))

_________________________________________________________________________________________________________________
#Fitting to a 2 factor model in limma
colnames(coef(fit1))

tmp1 <- contrasts.fit(fit1, coef =2) #directly contrast 2nd coefficient (group young)

tmp1 <- eBayes(tmp1, trend = TRUE)
tab1 <- topTable(tmp1, n=Inf)
head(tab1, 20)
length(which(tab1$adj.P.Val < 0.05)) # number of DE genes
write.csv(as.data.frame(tab1), file="limma-voom_Dose0.csv") 

mask <- tab1$adj.P.Val < 0.05

head(mask)

deGenes1 <- (tab1[mask, ])
head(deGenes1)
write.csv(as.data.frame(deGenes1), file="LV_Dose0_DEGs.csv") 

tmp2 <- contrasts.fit(fit1, coef = 4) #directly contrast 4th coefficient (MaternalpopYoung:Dose102)
tmp2 <- eBayes(tmp2, trend = TRUE)
tab2 <- topTable(tmp2, sort.by = "p", n = Inf)
head(tab2, 20)
length(which(tab2$adj.P.Val < 0.05)) # number of DE genes
write.csv(as.data.frame(tab2), file="limma-voom_Dose102.csv") 
mask <- tab2$adj.P.Val < 0.05

head(mask)

deGenes2 <- (tab2[mask, ])
head(deGenes2)
deGenes2
write.csv(as.data.frame(deGenes2), file="LV_Dose102_DEGs.csv") 

tmp3 <- contrasts.fit(fit1, coef = 6) #directly contrast 6th coefficient (MaternalpopYoung:Dose408)
tmp3 <- eBayes(tmp3, trend = TRUE)
tab3 <- topTable(tmp3, sort.by = "P", n = Inf)
mask <- tab3$adj.P.Val < 0.05
head(mask)
deGenes3 <- (tab3[mask, ])
deGenes3 <- head(deGenes3,20)
head(deGenes3)
dim(deGenes3)
write.csv(as.data.frame(deGenes3), file="LV_Dose408_DEGs.csv") 
head(tab3)


length(which(tab3$adj.P.Val < 0.05)) # number of DE genes
write.csv(as.data.frame(tab3), file="limma-voom_Dose408.csv") 

