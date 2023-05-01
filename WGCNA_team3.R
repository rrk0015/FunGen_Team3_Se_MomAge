#Link to tutorial: https://alexslemonade.github.io/refinebio-examples/04-advanced-topics/network-analysis_rnaseq_01_wgcna.html#3_Using_a_different_refinebio_dataset_with_this_analysis
#############################################################################
#WGCNA-RNASeq
#Se-toxicity_MaternalAge_Daphnia


#Getting set up! Specifying file paths to working data. 
# Create the data folder if it doesn't exist
if (!dir.exists("data")) {
  dir.create("data")
}

# Define the file path to the plots directory
plots_dir <- "plots"

# Create the plots folder if it doesn't exist
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}

# Define the file path to the results directory
results_dir <- "results"

# Create the results folder if it doesn't exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

# Define the file path to the data directory
# Replace with the path of the folder the files will be in
data_dir <- file.path("data", "rawdata")

# Declare the file path to the gene expression matrix file
# inside directory saved as `data_dir`
# Replace with the path to your dataset file
data_file <- file.path(data_dir, "gene_count_matrix.csv")

# Declare the file path to the metadata file
# inside the directory saved as `data_dir`
# Replace with the path to your metadata file
metadata_file <- file.path(data_dir, "PHENO_DATA.txt")

# Check if the gene expression matrix file is at the path stored in `data_file`
file.exists(data_file)

# Check if the metadata file is at the file path stored in `metadata_file`
file.exists(metadata_file)

#installing packages and loading libraries
if (!("DESeq2" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("DESeq2", update = FALSE)
}

if (!("impute" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("impute")
}

if (!("WGCNA" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  install.packages("WGCNA")
}

if (!("GO.db" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("GO.db")
}

if (!("preprocessCore" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("preprocessCore")
}

if (!("ggforce" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  install.packages("ggforce")
}

if (!("ComplexHeatmap" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("ComplexHeatmap")
}

if (!("limma" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("limma")
}


# Attach the DESeq2 library
library(DESeq2)

# We will need this so we can use the pipe: %>%
library(magrittr)

# We'll need this for finding gene modules
library(WGCNA)

# We'll be making some plots
library(ggplot2)

#Run linear model on each module
library(limma)

# Read in metadata TSV file
metadata <- readr::read_tsv(metadata_file)

# Read in data TSV file
df <- readr::read_csv(data_file) %>%
  # Here we are going to store the gene IDs as row names so that we can have a numeric matrix to perform calculations on later
  tibble::column_to_rownames("gene_id")

# Make the data in the order of the metadata
df <- df %>%
  dplyr::select(metadata$sample)

# Check if this is in the same order
all.equal(colnames(df), metadata$sample)

# The next DESeq2 functions need the values to be converted to integers
df <- round(df) %>%
  # The next steps require a data frame and round() returns a matrix
  as.data.frame() %>%
  # Only keep rows that have total counts above the cutoff
  dplyr::filter(rowSums(.) >= 50)

metadata <- metadata %>%
  dplyr::mutate(
    treatment = dplyr::case_when(
      # Create our new variable based on treatment containing 
      stringr::str_detect(treatment, "0ug SeMet/L Dose From OldMoms") ~ "Control Old Mom",
      stringr::str_detect(treatment, "102ug SeMet/L Dose From OldMoms") ~ "LowDose Old Mom",
      stringr::str_detect(treatment, "408ug SeMet/L Dose From OldMoms") ~ "HighDose Old Mom",
      stringr::str_detect(treatment, "0ug SeMet/L Dose From YoungMoms") ~ "Control Young Mom",
      stringr::str_detect(treatment, "102ug SeMet/L Dose From YoungMoms") ~ "LowDose Young Mom",
      stringr::str_detect(treatment, "408ug SeMet/L Dose From YoungMoms") ~ "HighDose Young Mom",
    ),
    # It's easier for future items if this is already set up as a factor
    treatment = as.factor(treatment)
  )

levels(metadata$treatment)

# Create a `DESeqDataSet` object
dds <- DESeqDataSetFromMatrix(
  countData = df, # Our prepped data frame with counts
  colData = metadata, # Data frame with annotation for our samples
  design = ~1 # Here we are not specifying a model
)


# Normalize and transform the data in the `DESeqDataSet` object using the `vst()`
# function from the `DESEq2` R package
dds_norm <- vst(dds)

gsg = goodSamplesGenes(df, verbose = 3);
gsg$allOK
#TRUE
#All genes have passed the cuts, no outliers to be removed
#If false, those genes should be removed

# Retrieve the normalized data from the `DESeqDataSet`
normalized_counts <- assay(dds_norm) %>%
  t() # Transpose this data

sft <- pickSoftThreshold(normalized_counts,
                         dataIsExpr = TRUE,
                         corFnc = cor,
                         networkType = "signed"
)

sft_df <- data.frame(sft$fitIndices) %>%
  dplyr::mutate(model_fit = -sign(slope) * SFT.R.sq)

#plotting
ggplot(sft_df, aes(x = Power, y = model_fit, label = Power)) +
  # Plot the points
  geom_point() +
  # We'll put the Power labels slightly above the data points
  geom_text(nudge_y = 0.1) +
  # We will plot what WGCNA recommends as an R^2 cutoff
  geom_hline(yintercept = 0.80, col = "red") +
  # Just in case our values are low, we want to make sure we can still see the 0.80 level
  ylim(c(min(sft_df$model_fit), 1.05)) +
  # We can add more sensible labels for our axis
  xlab("Soft Threshold (power)") +
  ylab("Scale Free Topology Model Fit, signed R^2") +
  ggtitle("Scale independence") +
  # This adds some nicer aesthetics to our plot
  theme_classic()

#Adjust power here
bwnet <- blockwiseModules(normalized_counts,
                          maxBlockSize = 5000, # What size chunks (how many genes) the calculations should be run in
                          TOMType = "signed", # topological overlap matrix
                          power = 20, # soft threshold for network construction
                          #numericLabels = TRUE, # Let's use numbers instead of colors for module labels
                          randomSeed = 1234, # there's some randomness associated with this calculation
                          # so we should set a seed
)

#filepath
readr::write_rds(bwnet,
                 file = file.path("results", "SRP140558_wgcna_results.RDS")
)

module_eigengenes <- bwnet$MEs

# Print out a preview
head(module_eigengenes)

all.equal(metadata$sample, rownames(module_eigengenes))

# Create the design matrix from the `time_point` variable
des_mat <- model.matrix(~ metadata$treatment)

# lmFit() needs a transposed version of the matrix
fit <- limma::lmFit(t(module_eigengenes), design = des_mat)

# Apply empirical Bayes to smooth standard errors
fit <- limma::eBayes(fit)

# Apply multiple testing correction and obtain stats
stats_df <- limma::topTable(fit, number = ncol(module_eigengenes)) %>%
  tibble::rownames_to_column("module")

write.csv(stats_df, file = "statsdf.csv")

head(stats_df)

module_t_df <- module_eigengenes %>%
  tibble::rownames_to_column("sample") %>%
  # Here we are performing an inner join with a subset of metadata
  dplyr::inner_join(metadata %>%
                      dplyr::select(sample, treatment),
                    by = c("sample" = "sample")
  )

#plot expression for a module
ggplot(
  module_t_df,
  aes(
    x = treatment,
    y = MEturquoise,
    color = treatment
  )
) +
  # a boxplot with outlier points hidden (they will be in the sina plot)
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  # A sina plot to show all of the individual data points
  ggforce::geom_sina(maxwidth = 0.3) +
  theme_classic()

gene_module_key <- tibble::enframe(bwnet$colors, name = "gene", value = "module") %>%
  # Let's add the `ME` part so its more clear what these numbers are and it matches elsewhere
  dplyr::mutate(module = paste0("ME", module))

gene_module_key %>%
  dplyr::filter(module == "MEturquoise")

readr::write_tsv(gene_module_key,
                 file = file.path("results", "SRP140558_wgcna_gene_to_module.tsv")
)

#establish heatmap function
make_module_heatmap <- function(module_name = "MEturquoise",
                                expression_mat = normalized_counts,
                                metadata_df = metadata,
                                gene_module_key_df = gene_module_key,
                                module_eigengenes_df = module_eigengenes) {
  # Create a summary heatmap of a given module.
  #
  # Args:
  # module_name: a character indicating what module should be plotted, e.g. "ME19"
  # expression_mat: The full gene expression matrix. Default is `normalized_counts`.
  # metadata_df: a data frame with refinebio_accession_code and time_point
  #              as columns. Default is `metadata`.
  # gene_module_key: a data.frame indicating what genes are a part of what modules. Default is `gene_module_key`.
  # module_eigengenes: a sample x eigengene data.frame with samples as row names. Default is `module_eigengenes`.
  #
  # Returns:
  # A heatmap of expression matrix for a module's genes, with a barplot of the
  # eigengene expression for that module.
  
  # Set up the module eigengene with its refinebio_accession_code
  module_eigengene <- module_eigengenes_df %>%
    dplyr::select(all_of(module_name)) %>%
    tibble::rownames_to_column("sample")
  
  # Set up column annotation from metadata
  col_annot_df <- metadata_df %>%
    # Only select the treatment and sample ID columns
    dplyr::select(treatment, treatment, sample) %>%
    # Add on the eigengene expression by joining with sample IDs
    dplyr::inner_join(module_eigengene, by = "sample") %>%
    # Arrange by patient and time point
    dplyr::arrange(treatment, sample) %>%
    # Store sample
    tibble::column_to_rownames("sample")
  
  # Create the ComplexHeatmap column annotation object
  col_annot <- ComplexHeatmap::HeatmapAnnotation(
    # Supply treatment labels
    treatment = col_annot_df$treatment,
    # Add annotation barplot
    module_eigengene = ComplexHeatmap::anno_barplot(dplyr::select(col_annot_df, module_name)),
    # Pick colors for each experimental group in time_point
    col = list(treatment = c("Control Old Mom" = "#f1a340", "HighDose Old Mom" = "#998ec3", "LowDose Old Mom" = "#56B4E9", "Control Young Mom" = "#2E9FDF", "HighDose Young Mom" = "#00AFBB", "LowDose Young Mom" = "#E7B800"))
  )
  
  # Get a vector of the Ensembl gene IDs that correspond to this module
  module_genes <- gene_module_key_df %>%
    dplyr::filter(module == module_name) %>%
    dplyr::pull(gene)
  
  # Set up the gene expression data frame
  mod_mat <- expression_mat %>%
    t() %>%
    as.data.frame() %>%
    # Only keep genes from this module
    dplyr::filter(rownames(.) %in% module_genes) %>%
    # Order the samples to match col_annot_df
    dplyr::select(rownames(col_annot_df)) %>%
    # Data needs to be a matrix
    as.matrix()
  
  # Normalize the gene expression values
  mod_mat <- mod_mat %>%
    # Scale can work on matrices, but it does it by column so we will need to
    # transpose first
    t() %>%
    scale() %>%
    # And now we need to transpose back
    t()
  
  # Create a color function based on standardized scale
  color_func <- circlize::colorRamp2(
    c(-2, 0, 2),
    c("#67a9cf", "#f7f7f7", "#ef8a62")
  )
  
  # Plot on a heatmap
  heatmap <- ComplexHeatmap::Heatmap(mod_mat,
                                     name = module_name,
                                     # Supply color function
                                     col = color_func,
                                     # Supply column annotation
                                     bottom_annotation = col_annot,
                                     # We don't want to cluster samples
                                     cluster_columns = FALSE,
                                     # We don't need to show sample or gene labels
                                     show_row_names = FALSE,
                                     show_column_names = FALSE
  )
  
  # Return heatmap
  return(heatmap)
}

#plot heatmap
mod_t_heatmap <- make_module_heatmap(module_name = "MEturquoise")

# Print out the plot=
mod_t_heatmap

#save png
png(file.path("results", "SRP140558_module_19_heatmap.png"))
mod_t_heatmap
dev.off()

# Print session info
sessioninfo::session_info()
