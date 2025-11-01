VCF2CpG
Simple and efficient Python tool to extract CpG methylation and coverage data from VCF files and merge them into a single matrix for downstream analysis.

Features
- Reads VCF files containing CpG methylation information.
- Extracts METH and COV fields using cyvcf2.
- Exports methylation data to clean CSV files.
- Merges multiple sample files into a combined CpG matrix.
- Ready for use in methylation analysis pipelines (e.g., methylKit, Minfi, custom R/Python scripts).

Requirements
- Python 3.8 or higher
- Packages:
  pip install cyvcf2 pandas

Usage
Step 1. Extract methylation from VCFs
  python extract_methylation_vcf.py
This creates individual .csv files for each sample.

Step 2. Merge extracted files
  python merge_methylation_csv.py
The merged file GSE186458_merged_methylation.csv will be saved in the output folder.

Step 3. Automate both steps (optional)
  bash 1stStep.sh

Output Example
Chromosome,Position,Methylation_s205,Coverage_s205,Methylation_s207,Coverage_s207
chr1,10012,0.83,15,0.78,18
chr1,10027,0.56,12,0.60,10

License
MIT License
