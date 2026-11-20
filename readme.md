# Differential Gene Expression Analysis of Bee Brain RNA-seq Data

## Overview

This repository contains the R scripts and input files used to perform differential gene expression analysis of RNA-seq data from bee brains collected before and after fungicide and pesticide application.

The analysis was performed using **DESeq2**, accounting for the paired experimental design in which the same independent colonies were sampled at both time points. The main objective was to identify genes showing changes in expression between samples collected before and after the application.

## Experimental Design

The experiment included samples collected from three independent bee colonies (C1, C3, and C5) at two time points:

- **Before:** samples collected prior to treatment application.
- **After:** samples collected following treatment application.

Each biological sample consisted of a pool of 10 bee brains. Therefore, the experimental unit was the colony, rather than the individual bee.

The experimental design comprised three biological pairs (three colonies sampled at both time points), resulting in six RNA-seq libraries.

Because the same colonies were sampled before and after treatment, colony identity was included in the DESeq2 model to account for between-colony variability.

## Repository Structure

```
.
├── Raw_counts.txt
├── experiment.txt
├── DESeq2_analysis.R
└── README.md
```

- **Raw_counts.txt:** Gene-level raw read count matrix. Rows represent genes and columns represent RNA-seq samples.
- **experiment.txt:** Sample metadata, including sample identifiers, treatment condition, and colony identity.
- **DESeq2_analysis.R:** R script used to perform differential gene expression analysis.
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
4. Estimation of moderated log2 fold changes using the `apeglm` shrinkage method.
5. Adjustment of p-values for multiple testing using the Benjamini–Hochberg procedure.
6. Identification and export of differentially expressed genes.

Genes with an adjusted p-value (`padj`) below 0.05 were considered differentially expressed.

Positive log2 fold-change values indicate higher expression in the After condition, whereas negative values indicate higher expression in the Before condition.

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
2. Place `Raw_counts.txt`, `experiment.txt`, and `DESeq2_analysis.R` in the same working directory.
3. Open the R script and update the input/output directory (`base_dir`) to match your local directory.
4. Run the script in R or RStudio.

The script automatically imports the count matrix and metadata, checks sample correspondence, fits the DESeq2 model, and exports the results.

## Output Files

The analysis generates the following main output files:

- **Deseq2_statistics.txt:** Complete differential expression results, including log2 fold changes, standard errors, test statistics, p-values, and adjusted p-values.
- **DEGs_After_vs_Before.txt:** List of genes identified as differentially expressed using the adjusted p-value threshold (`padj < 0.05`).
- **rlog_transformed_counts.txt:** Regularized log-transformed expression values.

Additional exploratory plots are generated by the script.

## Important Considerations

The paired design accounts for differences among colonies by including colony identity as a blocking factor in the statistical model. However, because all colonies were sampled before and after treatment application, and no untreated colonies were included as a temporal control, changes in gene expression cannot be unequivocally attributed to the application itself. The results should therefore be interpreted as transcriptional changes associated with the post-application sampling period.

## Citation

If you use this code or the associated analysis workflow, please cite this repository and the DESeq2 methodology:

Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15, 550. [https://doi.org/10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8)

Zhu, A., Ibrahim, J. G., & Love, M. I. (2019). Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences. *Bioinformatics*, 35(12), 2084–2092. [https://doi.org/10.1093/bioinformatics/bty895](https://doi.org/10.1093/bioinformatics/bty895)

## License

Please specify the license under which this repository is distributed.

**Una consideración estadística importante:** con solo tres pares de colonias, el modelo `~ Colony + Treatment` sigue siendo el apropiado para representar el diseño, pero la potencia estadística será limitada. Por eso conviene interpretar los resultados teniendo en cuenta el tamaño muestral y no solamente el umbral de `padj < 0.05`.
