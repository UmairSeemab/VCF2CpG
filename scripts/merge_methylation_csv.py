import pandas as pd

# Load individual sample files
s205 = pd.read_csv("s205.csv")
s207 = pd.read_csv("s207.csv")

# Rename methylation columns for clarity
s205 = s205.rename(columns={"Methylation": "Methylation_s205", "Coverage": "Coverage_s205"})
s207 = s207.rename(columns={"Methylation": "Methylation_s207", "Coverage": "Coverage_s207"})

# Merge on Chromosome and Position
merged = pd.merge(s205, s207, on=["Chromosome", "Position"], how="inner")

# Drop missing or invalid methylation entries
merged = merged.dropna(subset=["Methylation_s205", "Methylation_s207"])

# Save the merged matrix
merged.to_csv("GSE186458_merged_methylation.csv", index=False)
print(f"Merged file saved as merged_methylation.csv with {len(merged)} CpG sites.")
