from cyvcf2 import VCF
import pandas as pd

def extract_vcf(vcf_file, output_csv):
    data = []
    vcf = VCF(vcf_file)
    for record in vcf:
        meth = record.INFO.get('METH')
        cov = record.INFO.get('COV')
        if meth is not None:
            data.append([record.CHROM, record.POS, meth, cov])
    df = pd.DataFrame(data, columns=['Chromosome', 'Position', 'Methylation', 'Coverage'])
    df.to_csv(output_csv, index=False)
    print(f"Saved {len(df)} CpG entries to {output_csv}")

extract_vcf("GSE186458_blocks.s205.vcf.gz", "s205.csv")
extract_vcf("GSE186458_blocks.s207.hg38.vcf.gz", "s207.csv")
