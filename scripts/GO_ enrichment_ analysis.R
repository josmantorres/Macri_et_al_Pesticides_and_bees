
# ============================================================
# Gene Ontology (GO) enrichment analysis using ermineR
# ============================================================

# Required packages
library(ermineR)
library(rJava)
library(httr)

# Increase timeout for internet requests
options(timeout = max(10000, getOption("timeout")))

# Define input and output directories
input_dir <- file.path(base_dir, "input")
output_dir <- file.path(base_dir, "output")

# Input files
go_annotation_file <- file.path(input_dir, "GO_annotations.txt")
raw_pvalues_file <- file.path(input_dir, "DESeq2_raw_pvalues.txt")

# Check that input files exist
if (!file.exists(go_annotation_file)) {
  stop("GO annotation file not found: ", go_annotation_file)
}

if (!file.exists(raw_pvalues_file)) {
  stop("DESeq2 raw p-values file not found: ", raw_pvalues_file)
}

# ------------------------------------------------------------
# 1. Import GO annotations
# ------------------------------------------------------------

# Expected columns: GeneID and GO
# Each row represents one gene-GO term association.
goAnnots <- read.delim(
  go_annotation_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!all(c("GeneID", "GO") %in% colnames(goAnnots))) {
  stop("GO_annotations.txt must contain columns named 'GeneID' and 'GO'.")
}

# Remove missing or empty gene identifiers and GO terms
goAnnots <- goAnnots[
  !is.na(goAnnots$GeneID) &
    !is.na(goAnnots$GO) &
    goAnnots$GeneID != "" &
    goAnnots$GO != "",
]

# Create a list of GO terms associated with each gene
GOList <- split(
  as.character(goAnnots$GO),
  as.character(goAnnots$GeneID)
)

# Format annotations for ermineR
annotationList <- makeAnnotation(
  GOList,
  return = TRUE
)

# ------------------------------------------------------------
# 2. Import DESeq2 raw p-values
# ------------------------------------------------------------

# Expected columns: GeneID and pvalue
# pvalue corresponds to the nominal, unadjusted DESeq2 p-value.
rawPvalues <- read.delim(
  raw_pvalues_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!all(c("GeneID", "pvalue") %in% colnames(rawPvalues))) {
  stop("DESeq2_raw_pvalues.txt must contain columns named 'GeneID' and 'pvalue'.")
}

# Remove missing values
rawPvalues <- rawPvalues[
  !is.na(rawPvalues$GeneID) &
    !is.na(rawPvalues$pvalue),
]

# Check p-value format
if (!is.numeric(rawPvalues$pvalue) ||
    any(rawPvalues$pvalue < 0 | rawPvalues$pvalue > 1)) {
  stop("The 'pvalue' column must contain numeric values between 0 and 1.")
}

# Ensure that gene identifiers are unique
if (anyDuplicated(rawPvalues$GeneID)) {
  stop("Duplicated gene identifiers found in DESeq2_raw_pvalues.txt.")
}

# Create the score table required by ermineR
scoreListRawP <- data.frame(
  pvalue = rawPvalues$pvalue,
  row.names = rawPvalues$GeneID
)

# ------------------------------------------------------------
# 3. Gene Score Resampling (GSR) analysis
# ------------------------------------------------------------

# Biological Process (BP)
ermineR(
  annotation = annotationList,
  expression = NULL,
  aspects = "B",
  scores = scoreListRawP,
  scoreColumn = 1,
  logTrans = TRUE,
  bigIsBetter = FALSE,
  test = "GSR",
  pAdjust = "FDR",
  iterations = 200000,
  geneReplicates = "mean",
  stats = "mean",
  return = FALSE,
  minClassSize = 10,
  maxClassSize = 200,
  output = file.path(output_dir, "GSR_pvalue_BP.txt")
)

# Molecular Function (MF)
ermineR(
  annotation = annotationList,
  expression = NULL,
  aspects = "M",
  scores = scoreListRawP,
  scoreColumn = 1,
  logTrans = TRUE,
  bigIsBetter = FALSE,
  test = "GSR",
  pAdjust = "FDR",
  iterations = 200000,
  geneReplicates = "mean",
  stats = "mean",
  return = FALSE,
  minClassSize = 10,
  maxClassSize = 200,
  output = file.path(output_dir, "GSR_pvalue_MF.txt")
)

# Cellular Component (CC)
ermineR(
  annotation = annotationList,
  expression = NULL,
  aspects = "C",
  scores = scoreListRawP,
  scoreColumn = 1,
  logTrans = TRUE,
  bigIsBetter = FALSE,
  test = "GSR",
  pAdjust = "FDR",
  iterations = 200000,
  geneReplicates = "mean",
  stats = "mean",
  return = FALSE,
  minClassSize = 10,
  maxClassSize = 200,
  output = file.path(output_dir, "GSR_pvalue_CC.txt")
)

