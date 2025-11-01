

setwd("/home/umair/Documents/ZiaBhaiTask/Data")


#####################  differential_methylation  ################

# Install methylKit

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

#BiocManager::install("methylKit")

library(methylKit)

# ---------- Load merged data ----------
data <- read.csv("GSE186458_merged_methylation.csv")

# ---------- Clean data (remove zero/NA coverage) ----------
data <- subset(data, Coverage_s205 > 0 & Coverage_s207 > 0)
data <- na.omit(data)

# ---------- Compute counts ----------
data$NumCs_s205 <- data$Methylation_s205 * data$Coverage_s205
data$NumTs_s205 <- data$Coverage_s205 - data$NumCs_s205
data$NumCs_s207 <- data$Methylation_s207 * data$Coverage_s207
data$NumTs_s207 <- data$Coverage_s207 - data$NumCs_s207

# ---------- Create methylRaw objects ----------
df_s205 <- data.frame(chr = data$Chromosome,
                      start = data$Position,
                      end = data$Position,
                      strand = "+",
                      coverage = data$Coverage_s205,
                      numCs = data$NumCs_s205,
                      numTs = data$NumTs_s205)

df_s207 <- data.frame(chr = data$Chromosome,
                      start = data$Position,
                      end = data$Position,
                      strand = "+",
                      coverage = data$Coverage_s207,
                      numCs = data$NumCs_s207,
                      numTs = data$NumTs_s207)

meth_s205 <- new("methylRaw", df_s205,
                 sample.id = "s205",
                 assembly = "hg38",
                 context = "CpG",
                 resolution = "base")

meth_s207 <- new("methylRaw", df_s207,
                 sample.id = "s207",
                 assembly = "hg38",
                 context = "CpG",
                 resolution = "base")

# ---------- Combine & filter ----------
meth_list <- new("methylRawList", list(meth_s205, meth_s207), treatment = c(0, 1))
meth_filt <- filterByCoverage(meth_list, lo.count = 10, hi.perc = 99.9)

# ---------- Unite common CpG sites ----------
meth_united <- unite(meth_filt, destrand = FALSE)
df_united <- getData(meth_united)

# ---------- Gentle numeric cleanup ----------
df_united$coverage1 <- round(df_united$coverage1)
df_united$coverage2 <- round(df_united$coverage2)
df_united$numCs1    <- pmax(0, round(df_united$numCs1))
df_united$numTs1    <- pmax(0, round(df_united$numTs1))
df_united$numCs2    <- pmax(0, round(df_united$numCs2))
df_united$numTs2    <- pmax(0, round(df_united$numTs2))

# Remove rows with missing or zero coverage in any sample
keep <- complete.cases(df_united) &
  df_united$coverage1 > 0 & df_united$coverage2 > 0
df_united <- df_united[keep, ]

cat("Remaining CpGs after filtering:", nrow(df_united), "\n")

if (nrow(df_united) < 100) {
  stop("Too few CpGs remain after filtering. Check input data scaling or filtering thresholds.")
}

# ---------- Rebuild methylBase ----------
meth_united <- new("methylBase",
                   df_united,
                   sample.ids = meth_united@sample.ids,
                   assembly   = meth_united@assembly,
                   context    = meth_united@context,
                   treatment  = meth_united@treatment,
                   destranded = FALSE,
                   resolution = "base")

cat("Running differential methylation test on", nrow(df_united), "CpGs...\n")

# ---------- Differential methylation ----------
diff_meth <- calculateDiffMeth(meth_united, test = "F", overdispersion = "none")

# ---------- Get significant CpGs ----------
sig_cpg <- getMethylDiff(diff_meth, difference = 25, qvalue = 0.01)

# ---------- Save output ----------
write.table(sig_cpg, "GSE186458_significant_CpGs.txt", sep="\t", quote=FALSE, row.names=FALSE)

cat("Differential methylation complete.\n")
cat("Significant CpGs:", nrow(sig_cpg), "\n")


# Fix scaling
sig_cpg <- read.table("GSE186458_significant_CpGs.txt", header=TRUE, sep="\t")

# Detect approximate scale factor (e.g. 100, 1000, 10000)
scale_factor <- max(abs(sig_cpg$meth.diff), na.rm=TRUE) / 100
sig_cpg$meth.diff <- sig_cpg$meth.diff / scale_factor

# Save corrected values
write.table(sig_cpg, "GSE186458_significant_CpGs_scaled.txt", sep="\t", quote=FALSE, row.names=FALSE)

summary(sig_cpg$meth.diff)
hist(sig_cpg$meth.diff, breaks=100, main="Methylation difference distribution", xlab="Methylation (%)")

pdf("GSE186458_methylation_diff_histogram.pdf", width=8, height=6)
hist(sig_cpg$meth.diff, breaks=100,
     main="Distribution of Differential Methylation",
     xlab="Methylation Difference (%)",
     col="lightblue", border="white")
dev.off()



########################  annotate_CpGs  ########################

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("annotatr", ask = FALSE, update = FALSE)
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
BiocManager::install("org.Hs.eg.db")




library(annotatr)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

# ---------- Load CpG data ----------
sig_cpg <- read.table("GSE186458_significant_CpGs_scaled.txt", header=TRUE, sep="\t")

# ---------- Convert to GRanges ----------
gr_cpg <- GRanges(seqnames = sig_cpg$chr,
                  ranges = IRanges(start = sig_cpg$start, end = sig_cpg$end),
                  strand = sig_cpg$strand,
                  meth.diff = sig_cpg$meth.diff,
                  qvalue = sig_cpg$qvalue)

# ---------- Build richer annotations ----------
builtin_annotations()
annots <- build_annotations(genome = "hg38",
                            annotations = c("hg38_basicgenes",
                                            "hg38_genes_promoters",
                                            "hg38_cpgs"))


# ---------- Annotate CpGs ----------
annotated_cpgs <- annotate_regions(
  regions = gr_cpg,
  annotations = annots,
  ignore.strand = TRUE,
  quiet = FALSE
)

annotated_df <- as.data.frame(annotated_cpgs)


# ---------- Convert to data frame ----------
annotated_df <- as.data.frame(annotated_cpgs)

# Inspect new columns
cat("Annotation columns available:\n")
print(colnames(annotated_df)[1:15])

annotated_df$gene_symbol <- annotated_df$annot.symbol
write.table(annotated_df,
            "GSE186458_annotated_CpGs_with_genes.txt",
            sep="\t", quote=FALSE, row.names=FALSE)



########################## Pathway and GO enrichment ######################################
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

# ---------- Load annotated CpG data ----------
annotated_df <- read.table("GSE186458_annotated_CpGs_with_genes.txt", header=TRUE, sep="\t")

# ---------- Extract unique gene symbols ----------
genes <- unique(na.omit(annotated_df$gene_symbol))

cat("Total unique genes for enrichment:", length(genes), "\n")

# ---------- Map symbols to Entrez IDs ----------
entrez_ids <- bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
entrez_gene_list <- entrez_ids$ENTREZID

# ---------- GO enrichment (Biological Process) ----------
ego <- enrichGO(gene          = entrez_gene_list,
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENTREZID",
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

# ---------- KEGG pathway enrichment ----------
ekegg <- enrichKEGG(gene         = entrez_gene_list,
                    organism     = "hsa",
                    pvalueCutoff = 0.05)

# ---------- Save results ----------
write.table(as.data.frame(ego), "GSE186458_GO_enrichment.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(as.data.frame(ekegg), "GSE186458_KEGG_enrichment.txt", sep="\t", quote=FALSE, row.names=FALSE)

# ---------- Visualization ----------
# ---------- GO Plot ----------
if (nrow(ego) > 0) {
  pdf("GSE186458_GO_barplot.pdf", width=8, height=6)
  print(dotplot(ego, showCategory=15, title="Top GO Biological Processes"))
  dev.off()
}
# ---------- KEGG Plot ----------
if (nrow(ekegg) > 0) {
  pdf("GSE186458_KEGG_barplot.pdf", width=8, height=6)
  print(dotplot(ekegg, showCategory=15, title="Top KEGG Pathways"))
  dev.off()
}

cat("Enrichment analysis complete.\n")
cat("Results saved as GSE186458_GO_enrichment.txt and GSE186458_KEGG_enrichment.txt\n")

################################ Combined enrichment summary figure ####################
library(ggplot2)
library(dplyr)
library(forcats)

# ---------- Load enrichment results ----------
go <- read.table("GSE186458_GO_enrichment.txt", header=TRUE, sep="\t", quote="")
kegg <- read.table("GSE186458_KEGG_enrichment.txt", header=TRUE, sep="\t", quote="")

# ---------- Prepare GO data ----------
go_plot <- go %>%
  mutate(Type = "GO",
         Description = as.character(Description),
         negLogP = -log10(p.adjust)) %>%
  arrange(p.adjust) %>%
  slice_head(n=10)

# ---------- Prepare KEGG data ----------
kegg_plot <- kegg %>%
  mutate(Type = "KEGG",
         Description = as.character(Description),
         negLogP = -log10(p.adjust)) %>%
  arrange(p.adjust) %>%
  slice_head(n=10)

# ---------- Combine both ----------
combined <- bind_rows(go_plot, kegg_plot)

# ---------- Order terms by significance ----------
combined <- combined %>%
  mutate(Description = fct_reorder(Description, negLogP))

# ---------- Plot ----------
p <- ggplot(combined, aes(x=Description, y=negLogP, fill=Type)) +
  geom_bar(stat="identity", position="dodge") +
  coord_flip() +
  labs(title="Top Enriched GO and KEGG Terms",
       x="Biological Process / Pathway",
       y="-log10(Adjusted p-value)") +
  theme_minimal(base_size=13) +
  scale_fill_manual(values=c("GO"="steelblue", "KEGG"="tomato")) +
  theme(legend.position="top",
        axis.text.y=element_text(size=10))

# ---------- Save ----------
ggsave("GSE186458_combined_GO_KEGG.pdf", plot=p, width=9, height=7)
ggsave("GSE186458_combined_GO_KEGG.png", plot=p, width=9, height=7, dpi=300)

cat("Combined enrichment figure saved as GSE186458_combined_GO_KEGG.pdf\n")

################   CpG-Level Visualizations  #################

################  Chromosomal Distribution  ###############

df$chr <- factor(df$chr, levels=paste0("chr",1:22))
barplot(table(df$chr), las=2, col="lightblue",
        main="Distribution of Significant CpGs per Chromosome",
        ylab="Number of CpGs")

################  Density Plot of Methylation Differences  ###############

plot(density(df$meth.diff), main="Density of Methylation Differences",
     xlab="Methylation Difference (%)", col="darkgreen", lwd=2)

################  Gene-Centric Heatmap  ###############

library(pheatmap)
ann <- read.table("GSE186458_annotated_CpGs_with_genes.txt", header=TRUE, sep="\t")
gene_meth <- aggregate(ann$meth.diff, by=list(ann$gene_symbol), mean, na.rm=TRUE)
top <- head(gene_meth[order(-abs(gene_meth$x)),], 50)
mat <- matrix(top$x, nrow=1)
colnames(mat) <- top$Group.1
pheatmap(mat, cluster_cols=TRUE, cluster_rows=FALSE,
         color=colorRampPalette(c("blue","white","red"))(50),
         main="Top Differentially Methylated Genes")

################  Genomic Context Distribution  ###############

library(ggplot2)
ggplot(ann, aes(x=annot.type, fill=annot.type)) +
  geom_bar() +
  coord_flip() +
  theme_minimal() +
  labs(title="CpGs by Genomic Region", x="Region Type", y="Count") +
  theme(legend.position="none")
ggsave("GSE186458_region_distribution.png", width=8, height=6, dpi=300)


################  Enrichment Network Plot  ###############

library(enrichplot)
cnetplot(ego, showCategory=10, foldChange=NULL, circular=TRUE, colorEdge=TRUE)



################  enrichment map  ###############

library(enrichplot)
library(clusterProfiler)

# Compute term similarity for GO enrichment result
ego_sim <- pairwise_termsim(ego)

# Now plot the enrichment map
emapplot(ego_sim, showCategory = 20, cex_label_category = 0.7)
pdf("GSE186458_GO_enrichment_map.pdf", width=9, height=7)
emapplot(ego_sim, showCategory=20, cex_label_category=0.7)
dev.off()

################  Word Cloud of GO Terms  ###############
install.packages(c("wordcloud", "tm", "slam"))
library(wordcloud)
words <- ego$Description
wordcloud(words, max.words=50, random.order=FALSE, colors=brewer.pal(8,"Dark2"))



################  Multi-Panel Publication Figure  ###############

library(patchwork)
p1 <- ggplot(df, aes(x=meth.diff, y=-log10(qvalue))) + geom_point(alpha=0.4)
p2 <- ggplot(ann, aes(x=annot.type, fill=annot.type)) + geom_bar()
p3 <- dotplot(ego, showCategory=10)
(p1 / p2) | p3
ggsave("GSE186458_summary_figure.pdf", width=12, height=8)

################  Correlation plot between samples  ###############

library(ggplot2)
merged <- read.csv("GSE186458_merged_methylation.csv")
ggplot(merged, aes(x=Methylation_s205, y=Methylation_s207)) +
  geom_hex() +
  geom_abline(color="red") +
  theme_minimal() +
  xlab("Sample s205") + ylab("Sample s207") +
  ggtitle("CpG-wise correlation between samples")


treeplot(ego, showCategory=20)


merged <- read.csv("GSE186458_merged_methylation.csv")
cor(merged$Methylation_s205, merged$Methylation_s207, use="complete.obs")

# Number of significant CpGs
nrow(sig_cpg)

# Number of unique genes affected
length(unique(annotated_df$gene_symbol))

# Examine top 10 hyper/hypo-methylated genes
head(annotated_df[order(-annotated_df$meth.diff), c("gene_symbol","meth.diff")], 10)
head(annotated_df[order(annotated_df$meth.diff), c("gene_symbol","meth.diff")], 10)

