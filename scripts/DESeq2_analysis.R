
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)

# =========================================================
# CONFIGURACION GENERAL
# =========================================================

base_dir <- "/home/josemlatorre/Documentos/Ivi_SClaudio/"

alpha <- 0.05
dec_out <- ","

# =========================================================
# 1. IMPORTAR COUNTS
# =========================================================

Counts <- read.delim(
  paste0(base_dir, "Raw_counts.txt"),
  sep = "\t",
  header = TRUE,
  check.names = FALSE
)

Counts_matrix <- data.matrix(Counts[, 2:ncol(Counts)])
rownames(Counts_matrix) <- Counts[, 1]

# =========================================================
# 2. IMPORTAR METADATA
# =========================================================

experiment <- read.delim(
  paste0(base_dir, "experiment.txt"),
  sep = "\t",
  header = TRUE
)

# Comprobar correspondencia entre counts y metadata

stopifnot(all(colnames(Counts_matrix) %in% experiment$Sample))

# Ordenar metadata según la matriz de counts

experiment <- experiment[
  match(colnames(Counts_matrix), experiment$Sample),
  ,
  drop = FALSE
]

stopifnot(identical(colnames(Counts_matrix), experiment$Sample))

# Definir factores experimentales

experiment$Treatment <- factor(
  experiment$Treatment,
  levels = c("Before", "After")
)

experiment$Colony <- factor(experiment$Colony)

# Comprobar número de muestras por colmena y condición

print(table(experiment$Colony, experiment$Treatment))

# Comprobar que cada colmena tenga una muestra Before y una After

stopifnot(
  all(table(experiment$Colony, experiment$Treatment) == 1)
)

# =========================================================
# 3. CONSTRUIR OBJETO DESEQ2 Y AJUSTAR EL MODELO
# =========================================================

# Diseño pareado:
# Colony controla las diferencias entre colmenas.
# Treatment estima el efecto Before vs After.

dds <- DESeqDataSetFromMatrix(
  countData = Counts_matrix,
  colData = experiment,
  design = ~ Colony + Treatment
)

# Ajustar modelo mediante Wald test

dds <- DESeq(dds, test = "Wald")

# Comprobar coeficientes del modelo

resultsNames(dds)

# Comparación After vs Before

coef_name <- "Treatment_After_vs_Before"

# =========================================================
# 4. ANALISIS DE EXPRESION DIFERENCIAL
# =========================================================

# Wald test convencional
# H0: log2FC = 0
# Sin umbral de log2FC

res <- results(
  dds,
  name = coef_name,
  alpha = alpha,
  independentFiltering = TRUE,
  cooksCutoff = TRUE
)

# =========================================================
# 5. SHRINKAGE DE LOG2FC CON APEGLM
# =========================================================

resLFC <- lfcShrink(
  dds,
  coef = coef_name,
  res = res,
  type = "apeglm"
)

# Incorporar estadísticos Wald a la tabla final

resLFC$pvalue_Wald <- res$pvalue
resLFC$padj_Wald <- res$padj

# Resumen general de resultados

summary(res, alpha = alpha)

# =========================================================
# 6. EXPORTAR RESULTADOS COMPLETOS
# =========================================================

res_final <- as.data.frame(resLFC)

write.table(
  res_final,
  file = paste0(base_dir, "Deseq2_statistics.txt"),
  sep = "\t",
  dec = dec_out,
  row.names = TRUE,
  col.names = NA,
  quote = FALSE
)

# =========================================================
# 7. IDENTIFICACION Y EXPORTACION DE DEGS
# =========================================================

# Criterio: adjusted p-value < 0.05
# Corrección Benjamini-Hochberg
# Sin umbral adicional de log2FC

degs <- res_final[
  !is.na(res_final$padj_Wald) &
    res_final$padj_Wald < alpha,
]

# Clasificar genes según dirección de expresión

degs$Regulation <- ifelse(
  degs$log2FoldChange > 0,
  "Up_in_After",
  "Up_in_Before"
)

# Ordenar por adjusted p-value

degs <- degs[order(degs$padj_Wald), ]

# Número de genes por dirección de regulación

print(table(degs$Regulation))

# Exportar DEGs

write.table(
  degs,
  file = paste0(base_dir, "DEGs_After_vs_Before.txt"),
  sep = "\t",
  dec = dec_out,
  row.names = TRUE,
  col.names = NA,
  quote = FALSE
)

# =========================================================
# 8. MA PLOT
# =========================================================

tiff(
  paste0(base_dir, "MAplot_DESeq2.tif"),
  width = 15,
  height = 15,
  units = "cm",
  res = 300,
  pointsize = 10,
  compression = "lzw"
)

plotMA(
  resLFC,
  alpha = alpha,
  ylim = c(-10, 10),
  main = "After vs Before",
  cex = 0.6
)

abline(h = 0, col = "blue", lty = 2)

dev.off()

# =========================================================
# 9. TRANSFORMACION RLOG
# =========================================================

# blind = TRUE: transformación independiente del diseño
# experimental, destinada al análisis exploratorio.

rlog_data <- rlogTransformation(
  dds,
  blind = TRUE,
  fitType = "parametric"
)

# Exportar valores transformados

write.table(
  data.frame(assay(rlog_data)),
  file = paste0(base_dir, "RLOG_Transformed_values.txt"),
  sep = "\t",
  dec = dec_out,
  row.names = TRUE,
  col.names = NA,
  quote = FALSE
)
