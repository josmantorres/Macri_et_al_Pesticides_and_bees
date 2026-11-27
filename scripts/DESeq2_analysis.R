
# ============================================================
# RNA-seq differential expression analysis
# DESeq2 and apeglm
#
# Experimental comparison: After vs Before
# Experimental design: ~ Colony + Treatment
#
# Differential expression significance:
# Wald test with Benjamini-Hochberg correction
# DEG threshold: adjusted p-value < 0.05
#
# Log2 fold-change estimation:
# apeglm shrinkage without adaptive prior estimation
# (apeAdapt = FALSE)
# ============================================================


# ============================================================
# 0. LOAD PACKAGES AND CONFIGURE ANALYSIS
# ============================================================

library(DESeq2)

# Input and output directories
input_dir <- "input"
output_dir <- "output"

dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

# Statistical significance threshold
alpha <- 0.05

# ============================================================
# 1. IMPORT RAW COUNTS
# ============================================================

counts_data <- read.delim(
  file.path(input_dir, "Raw_counts.txt"),
  sep = "\t",
  header = TRUE,
  check.names = FALSE
)

# Extract gene identifiers and count matrix

gene_ids <- counts_data[, 1]

counts_matrix <- data.matrix(
  counts_data[, 2:ncol(counts_data)]
)

rownames(counts_matrix) <- gene_ids

# Check count matrix integrity

stopifnot(!anyDuplicated(rownames(counts_matrix)))
stopifnot(!anyDuplicated(colnames(counts_matrix)))
stopifnot(all(is.finite(counts_matrix)))
stopifnot(all(counts_matrix >= 0))
stopifnot(all(counts_matrix == round(counts_matrix)))

storage.mode(counts_matrix) <- "integer"


# ============================================================
# 2. IMPORT EXPERIMENTAL METADATA
# ============================================================

experiment <- read.delim(
  file.path(input_dir, "experiment.txt"),
  sep = "\t",
  header = TRUE
)

# Check metadata structure

stopifnot(
  all(c("Sample", "Treatment", "Colony") %in% colnames(experiment))
)

stopifnot(!anyDuplicated(experiment$Sample))

# Check correspondence between counts and metadata

stopifnot(
  setequal(colnames(counts_matrix), experiment$Sample)
)

# Reorder metadata to match the count matrix

experiment <- experiment[
  match(colnames(counts_matrix), experiment$Sample),
  ,
  drop = FALSE
]

stopifnot(
  identical(colnames(counts_matrix), experiment$Sample)
)

# Define experimental factors

experiment$Treatment <- factor(
  experiment$Treatment,
  levels = c("Before", "After")
)

experiment$Colony <- factor(experiment$Colony)

# Check sample distribution across colonies and treatments

print(table(experiment$Colony, experiment$Treatment))

# Confirm one sample per colony and treatment combination

stopifnot(
  all(table(experiment$Colony, experiment$Treatment) == 1)
)


# ============================================================
# 3. CONSTRUCT DESEQ2 DATASET AND FIT THE MODEL
# ============================================================

# Experimental design:
# Colony accounts for differences between colonies.
# Treatment estimates the After vs Before contrast.

dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData = experiment,
  design = ~ Colony + Treatment
)

# Fit the negative binomial model using the Wald test

dds <- DESeq(
  dds,
  test = "Wald"
)

# Inspect model coefficients

print(resultsNames(dds))

# Define the coefficient for the After vs Before comparison

coef_name <- "Treatment_After_vs_Before"

stopifnot(coef_name %in% resultsNames(dds))


# ============================================================
# 4. DIFFERENTIAL EXPRESSION ANALYSIS
# ============================================================

# Wald test:
# Null hypothesis: log2FC = 0
# No additional log2FC threshold

res_wald <- results(
  dds,
  name = coef_name,
  alpha = alpha,
  independentFiltering = TRUE,
  cooksCutoff = TRUE
)

# Display differential expression summary

summary(res_wald, alpha = alpha)


# ============================================================
# 5. LOG2 FOLD-CHANGE SHRINKAGE USING APEGLM
# ============================================================

# Estimate shrunken log2 fold changes using apeglm.
# Adaptive prior estimation is disabled to avoid the
# excessive shrinkage observed with apeAdapt = TRUE.

res_apeglm <- lfcShrink(
  dds,
  coef = coef_name,
  res = res_wald,
  type = "apeglm",
  apeAdapt = FALSE
)


# ============================================================
# 6. INTEGRATE WALD AND APEGLM RESULTS
# ============================================================

# Create a unified results table.
# Wald statistics are retained for differential expression
# testing, while apeglm estimates are used for effect sizes.

res_final <- data.frame(
  baseMean = res_wald$baseMean,

  # Original Wald estimates
  log2FC_Wald = res_wald$log2FoldChange,
  lfcSE_Wald = res_wald$lfcSE,
  stat_Wald = res_wald$stat,
  pvalue_Wald = res_wald$pvalue,
  padj_Wald = res_wald$padj,

  # Shrunken estimates from apeglm
  log2FC_apeglm = res_apeglm$log2FoldChange,
  lfcSD_apeglm = res_apeglm$lfcSE,

  row.names = rownames(res_wald)
)

# Count significant differentially expressed genes

n_DEGs <- sum(
  !is.na(res_final$padj_Wald) &
    res_final$padj_Wald < alpha
)

cat("Number of differentially expressed genes:", n_DEGs, "\n")


# ============================================================
# 7. EXPORT COMPLETE DIFFERENTIAL EXPRESSION RESULTS
# ============================================================

write.table(
  res_final,
  file = file.path(output_dir, "Deseq2_statistics.txt"),
  sep = "\t",
  dec = ".",
  row.names = TRUE,
  col.names = NA,
  quote = FALSE
)


# ============================================================
# 8. IDENTIFICATION AND EXPORT OF DEGs
# ============================================================

# DEG selection criterion:
# Benjamini-Hochberg adjusted Wald p-value < 0.05
# No additional log2FC threshold.

degs <- res_final[
  !is.na(res_final$padj_Wald) &
    res_final$padj_Wald < alpha,
  ,
  drop = FALSE
]

# Classify regulation using shrunken log2FC estimates

degs$Regulation <- ifelse(
  degs$log2FC_apeglm > 0,
  "Up_in_After",
  "Up_in_Before"
)

# Order DEGs by adjusted p-value

degs <- degs[
  order(degs$padj_Wald),
  ,
  drop = FALSE
]

# Display the number of genes in each regulation category

print(table(degs$Regulation))

# Export DEG results

write.table(
  degs,
  file = file.path(output_dir, "DEGs_After_vs_Before.txt"),
  sep = "\t",
  dec = ".",
  row.names = TRUE,
  col.names = NA,
  quote = FALSE
)


# ============================================================
# 9. MA PLOT
# ============================================================

tiff(
  file.path(output_dir, "MAplot_DESeq2.tif"),
  width = 15,
  height = 15,
  units = "cm",
  res = 300,
  pointsize = 10,
  compression = "lzw"
)

plotMA(
  res_apeglm,
  alpha = alpha,
  ylim = c(-10, 10),
  main = "After vs Before",
  cex = 0.6
)

abline(
  h = 0,
  col = "blue",
  lty = 2
)

dev.off()


# ============================================================
# 10. RLOG TRANSFORMATION
# ============================================================

# Apply regularized log transformation for exploratory analyses.
# blind = TRUE estimates the transformation independently
# of the experimental design.

rlog_data <- rlogTransformation(
  dds,
  blind = TRUE,
  fitType = "parametric"
)

# Export transformed expression values

write.table(
  data.frame(assay(rlog_data)),
  file = file.path(output_dir, "RLOG_Transformed_values.txt"),
  sep = "\t",
  dec = ".",
  row.names = TRUE,
  col.names = NA,
  quote = FALSE
)
