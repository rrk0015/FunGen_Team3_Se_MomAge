#!/bin/sh

source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load star
#module load samtools

# Define file paths
fastq_dir="/scratch/mzr0111/CleanData/"
genome_dir="/scratch/mzr0111/DaphniaRefGenome_6/"
genome_fasta="GCF_021134715.1_ASM2113471v1_genomic.fasta"
gtf_file="GCF_021134715.1_ASM2113471v1_genomic.gtf"

# Create output directory for genome indexes
index_dir="${genome_dir}/genome_index"
#mkdir -p "${index_dir}"

# Generate genome indexes using STAR
#STAR --runThreadN 8 \
#     --runMode genomeGenerate \
#     --genomeDir "${index_dir}" \
#     --genomeFastaFiles "${genome_dir}/${genome_fasta}" \
#     --sjdbGTFfile "${genome_dir}/${gtf_file}"

# Perform read mapping using STAR
sample_name="SRR6853330"
fastq_file_1="${fastq_dir}/${sample_name}_1_paired.fastq"
fastq_file_2="${fastq_dir}/${sample_name}_2_paired.fastq"
output_dir="${fastq_dir}/star_output"
#mkdir -p "${output_dir}"

STAR --runThreadN 20 \
     --genomeDir "${index_dir}" \
     --readFilesIn "${fastq_file_1}" "${fastq_file_2}" \
     --outFileNamePrefix "${output_dir}/${sample_name}."

# Convert SAM to BAM using samtools
#output_sam="${output_dir}/${sample_name}.sam"
#output_bam="${output_dir}/${sample_name}.bam"
#samtools view -bS "${output_sam}" > "${output_bam}"

# Sort BAM file
#samtools sort "${output_bam}" -o "${output_bam%.bam}.sorted.bam"

# Index BAM file
#samtools index "${output_bam%.bam}.sorted.bam"

# Generate flagstat statistics
#samtools flagstat "${output_bam%.bam}.sorted.bam" > "${output_dir}/${sample_name}.flagstat.txt"

