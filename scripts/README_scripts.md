VCF2CpG Scripts Folder

This folder contains all helper scripts used for the methylation extraction, merging, and downstream analysis workflow.

Files:
- 1stStep.sh : Runs the entire process automatically.
- extract_methylation_vcf.py : Extracts methylation and coverage from VCF files.
- merge_methylation_csv.py : Merges multiple sample CSVs into one CpG matrix.
- MethylationPipeline.R : Performs downstream R-based analysis such as correlation, visualization, or statistical comparison.

Usage:
1. Download all scripts in this folder.
2. Place them in the same directory as your data files (VCF or CSV).
3. Rename input file names inside scripts to match your data.
4. Run in order:
   bash 1stStep.sh
   Rscript MethylationPipeline.R

Note:
- Keep data files and scripts in the same folder.
- Modify file names only, not column headers or code logic.
- Output files are saved automatically in the working directory.
