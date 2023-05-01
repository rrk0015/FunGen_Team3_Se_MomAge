################# Performing Stats Analysis on HISEQ2 Results ###############

# Load the data from file
data = read.csv("hiseq_comparisons_reads_oldvnew.csv")

# Assign column variables
Reads = data$Reads
Reference.Genome = data$Reference.Genome

# Run the actual T-test
t.test(Reads ~ Reference.Genome, var.equal = TRUE)


