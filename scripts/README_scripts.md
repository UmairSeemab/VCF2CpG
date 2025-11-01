VCF2CpG Scripts Folder

This folder contains helper scripts written in Bash, Python, and R to process DNA methylation data.

Folder structure
scripts/
│
├── 1stStep.sh
├── extract_methylation_vcf.py
├── merge_methylation_csv.py
├── your_analysis.R

Instructions
1. Download all scripts in this folder.
2. Place them in the same directory as your data files (for example, VCF or CSV files).
3. Rename the input file names inside the scripts to match your own data files.
   - Example in extract_methylation_vcf.py:
     extract_vcf("YourSample1.vcf.gz", "sample1.csv")
     extract_vcf("YourSample2.vcf.gz", "sample2.csv")
4. Run the scripts step by step or through the provided Bash file.

Typical workflow
1. Extract CpG data from VCF files
   python extract_methylation_vcf.py
2. Merge extracted CSV files
   python merge_methylation_csv.py
3. Run everything automatically
   bash 1stStep.sh
4. Proceed to downstream R analysis
   Rscript your_analysis.R

Notes
- All scripts must be in the same folder as the input data files.
- Modify file names only, not column headers or script logic.
- The merged output file will be created in the same directory or an output subfolder, depending on script settings.
