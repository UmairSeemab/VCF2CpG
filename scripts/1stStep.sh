#!/bin/bash
# Convert all .bed.gz methylation files in the current folder to .vcf.gz

echo "Starting BED to VCF conversion..."

sudo apt install tabix
sudo apt install bgzip
python3 bed_to_vcf_auto.py

echo "Prepare and inspect the data"

sudo apt install bcftools


bcftools view GSE186458_blocks.s205.vcf.gz | head
bcftools stats GSE186458_blocks.s205.vcf.gz > vcf_summary.txt
bcftools view GSE186458_blocks.s207.hg38.vcf.gz | head
bcftools stats GSE186458_blocks.s207.hg38.vcf.gz > hg38vcf_summary.txt

echo "Extract methylation data from each VCF"

conda install -c conda-forge -c bioconda cyvcf2
python3 extract_methylation_vcf.py

echo "Merge s205 and s207 methylation data"
python3 merge_methylation_csv.py

echo "open R studio and set the path according to your data and run the R file's each command" 








