# Differential Gene Expression Analysis of Bee Brain RNA-seq Data

## Overview

This repository contains the R scripts, input and output files used to perform differential gene expression analysis of RNA-seq data from bee brains collected before and after fungicide and pesticide application.

The analysis was performed using **DESeq2**, accounting for the paired experimental design in which the same independent colonies were sampled at both time points. The main objective was to identify genes showing changes in expression between samples collected before and after the application.

## Experimental Design

The experiment included three independent bee colonies (C1, C3, and C5), each sampled before and after the combined fungicide and pesticide application. Each biological sample consisted of a pool of 10 bee brains. Thus, the colony was considered the experimental unit, with three biological replicates per sampling condition. The same colonies were sampled at both time points, although the individual bees used for each sampling were not necessarily the same.
The experimental design comprised three biological pairs (three colonies sampled at both time points), resulting in six RNA-seq libraries.

Because the same colonies were sampled before and after treatment, colony identity was included in the DESeq2 model to account for between-colony variability.

## Repository Structure

```
.
│   ├── Raw_counts.txt
│   ├── experiment.txt
│   ├── GO_annotations.txt
│   └── DESeq2_raw_pvalues.txt
│
├── scripts/
│   └── DESeq2_analysis.R
│   └── GO_ enrichment_ analysis
│
├── output/
│   ├── Deseq2_statistics.txt
│   ├── RLog_values.txt
│   ├── GSR_pvalue_BP.txt
│   ├── GSR_pvalue_MF.txt
│   └── GSR_pvalue_CC.txt
│
├── README.md
└── LICENSE
```

- **Raw_counts.txt:** Gene-level raw read count matrix. Rows represent genes and columns represent RNA-seq samples.
- **experiment.txt:** Sample metadata, including sample identifiers, treatment condition, and colony identity.
- **Deseq2_analysis.R:** R script used to perform differential gene expression analysis.
- **Deseq2_statistics.txt:** Table with the results of the differential expression analysis.
- **RLoG_Transformed_values:** Table with the rlog-transformed values.
- **README.md:** Repository documentation.

## Input Files

### 1. Raw_counts.txt

The count matrix must contain raw, non-normalized gene-level read counts. The first column should contain gene identifiers, while the remaining columns should correspond to individual RNA-seq samples.

Example:

| Gene_ID  | C1A | C3A | C5A | C1D | C3D | C5D |
| -------- | --- | --- | --- | --- | --- | --- |
| Gene_001 | 125 | 143 | 98  | 112 | 130 | 105 |
| Gene_002 | 45  | 32  | 51  | 67  | 40  | 58  |
| Gene_003 | 0   | 12  | 5   | 8   | 15  | 3   |

The values above are illustrative and do not represent the actual dataset.

### 2. experiment.txt

The metadata file must contain three columns:

- **Sample:** Sample identifier, matching the column names in `Raw_counts.txt`.
- **Treatment:** Sampling condition (`Before` or `After`).
- **Colony:** Colony identifier.

The experimental metadata file should follow this structure:

```
Sample    Treatment    Colony
C1A       Before       C1
C3A       Before       C3
C5A       Before       C5
C1D       After        C1
C3D       After        C3
C5D       After        C5
```

Each colony must have one sample for each treatment condition. Sample identifiers and their corresponding metadata must match exactly.

## Differential Expression Analysis

Differential gene expression analysis was performed using the **DESeq2** R package.

The statistical model was specified as:

```
design = ~ Colony + Treatment
```

This model accounts for colony-specific differences while estimating the effect associated with the sampling condition.

The analysis included the following steps:

1. Import and preparation of the raw count matrix and sample metadata.
2. Construction of a DESeq2 dataset using the paired experimental design.
3. Differential expression analysis using the Wald test.
4. Estimation of shrunken log2 fold changes using the apeglm method. Log2 fold changes were shrunk using apeglm with apeAdapt = FALSE, while statistical significance was assessed using the original Wald test p-values adjusted with the Benjamini–Hochberg procedure.
5. Identification and export of differentially expressed genes.

Genes with an adjusted p-value (`padj`) below 0.05 were considered differentially expressed.

Positive log2 fold-change values indicate higher expression after treatment, whereas negative values indicate higher expression before treatment.

## Requirements

The analysis was performed in R using the following packages:

- `DESeq2`
- `pheatmap`
- `RColorBrewer`
- `ggplot2`
- `ggrepel`

The required packages can be installed using:

```
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("DESeq2")

install.packages(c(
    "pheatmap",
    "RColorBrewer",
    "ggplot2",
    "ggrepel"
))
```

## Running the Analysis

1. Clone or download this repository.
2. Ensure that the input files are located in the input/ directory.
3. Open the DESeq2_analysis.R script located in the scripts/ directory.
4. Set the working directory to the repository root.
5. Run the script in R or RStudio.

The script automatically imports the count matrix and metadata, checks sample correspondence, fits the DESeq2 model, and exports the results.

## Output Files

The analysis generates the following main output files:

- **Deseq2_statistics.txt:** Complete differential expression results, including log2 fold changes, standard errors, test statistics, p-values, and adjusted p-values.
- **DEGs_After_vs_Before.txt:** List of genes identified as differentially expressed using the adjusted p-value threshold (`padj < 0.05`).
- **rlog_transformed_counts.txt:** Regularized log-transformed expression values.

## Important Considerations

The paired design accounts for differences among colonies by including colony identity as a blocking factor in the statistical model. However, because all colonies were sampled before and after treatment application, and no untreated colonies were included as a temporal control, changes in gene expression cannot be unequivocally attributed to the application itself. The results should therefore be interpreted as transcriptional changes associated with the post-application sampling period.

## Gene Ontology (GO) enrichment analysis

Gene Ontology (GO) enrichment analysis was performed using the **ermineR** R package and the Gene Score Resampling (GSR) method. This approach evaluates GO term enrichment using gene-level statistical scores rather than selecting genes based on a predefined significance threshold.

### GO annotation

Gene-to-GO term associations for *Apis mellifera* were obtained from the [Hymenoptera Genome Database (HGD)](https://hymenoptera.elsiklab.missouri.edu/hgd-go-annotation).

The annotation file was formatted as a two-column tab-delimited table containing gene identifiers (`GeneID`) and their associated GO terms (`GO`). Each gene–GO term association was represented in a separate row. Genes associated with multiple GO terms therefore appeared in multiple rows.

The annotation table was processed in R and converted into the format required by ermineR using the `makeAnnotation()` function.

### Gene scores

The input scores consisted of **raw (unadjusted) p-values generated by DESeq2**. Specifically, the nominal `pvalue` column from the DESeq2 differential expression results was used; Benjamini–Hochberg-adjusted p-values (`padj`) were not used as input.

No preliminary significance cutoff was applied to select genes for the GSR analysis. The analysis used the p-values available for the genes included in the input score file and their corresponding GO annotations.

The p-values were log-transformed by ermineR (`logTrans = TRUE`), so that smaller p-values corresponded to higher scores of statistical significance. Accordingly, `bigIsBetter = FALSE` was specified.

### GSR analysis parameters

The GSR analysis was performed independently for the three GO aspects:

- **Biological Process (BP)**
- **Molecular Function (MF)**
- **Cellular Component (CC)**

The following parameters were used:

| Parameter | Setting |
|---|---|
| Enrichment method | Gene Score Resampling (GSR) |
| Input scores | Raw DESeq2 p-values |
| Log transformation | TRUE |
| Higher scores are better | FALSE |
| Multiple-testing correction | False Discovery Rate (FDR) |
| Number of iterations | 200,000 |
| Multiple scores per gene | Mean |
| Gene set statistic | Mean |
| Minimum GO term size | 10 genes |
| Maximum GO term size | 200 genes |

FDR correction was applied by ermineR to the GO enrichment results to account for multiple testing. This correction is separate from the p-values used as input from DESeq2.

### Input and output files

The following files are required for the GO enrichment analysis:

- `GO_annotations.txt`: gene-to-GO term annotation table obtained from HGD.
- `DESeq2_raw_pvalues.txt`: gene identifiers and corresponding raw p-values from DESeq2.

The enrichment results are saved separately for each GO aspect:

- `GSR_pvalue_BP.txt`: Biological Process results.
- `GSR_pvalue_MF.txt`: Molecular Function results.
- `GSR_pvalue_CC.txt`: Cellular Component results.

The analysis requires the `ermineR`, `rJava`, and `httr` R packages, as well as a compatible Java installation.

## License

This repository is distributed under the MIT License. See the LICENSE file for details
