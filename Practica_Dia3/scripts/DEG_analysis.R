######
# Script : Analisis de expresion diferencial
# Author: Sofia Salazar y Evelia Coss
# Date: 03/03/2025
# Description: El siguiente script nos permite realiza el Analisis de expresion Diferencial
# a partir de los datos provenientes del alineamiento de STAR a R,
# Primero correr el script "import_txi_inR.R"
# Usage: Correr las lineas en un nodo de prueba en el cluster.
# Arguments:
#   - Input: Cargar la variable raw_counts.RData que contiene la matriz de cuentas y la metadata
#   - Output: DEG
#######

# qlogin 
# module load r/4.0.2
# R

# --- Load packages ----------
library(DESeq2)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(reshape2)

# --- Load data -----
# Cargar archivos
workdir <- 'RNA_lupus/'
load(file = paste0(workdir, 'counts/txi_filtered.RData'))
outdir = 'RNA_lupus/'

metadata <- read.csv(file = paste0(workdir, 'metadata/metadata.csv'), header = T)

# Solo control
control_df <- metadata[metadata$Group == "Ctrl",]
dim(control_df) # [1] 50 14


## ---- Construct DESEQDataSet Object (dds object) -----
# If you have performed transcript quantification (with Salmon, kallisto, RSEM, etc.) 
# you could import the data with tximport, which produces a list, and then you can 
# use DESeqDataSetFromTximport().

names(filtered_txi) # columnas en Txi
## [1] "abundance"           "counts"              "length"             
## [4] "countsFromAbundance"

dds <- DESeqDataSetFromTximport(txi = filtered_txi,
                                colData = metadata,
                                design = ~ Dose_category + Cell_type + SLEDAI_category)
dds

# using counts and average transcript lengths from tximport

#class: DESeqDataSet
#dim: 29744 165
#metadata(1): version
#assays(2): counts avgTxLength
#rownames(29744): A1BG A1BG-AS1 ... ZZEF1 ZZZ3
#rowData names(0):
#  colnames(165): QR011_0 QR011_1 ... QR108_3 QR108_4
#colData names(7): sample_ID Group ... Cell_type Internal

dim(dds)
# [1] 21816   165

## ---- Normalizacion de los datos -----
counts <- assays(dds)$counts

# > vsd
vsd <- varianceStabilizingTransformation(counts)
head(vsd)

# Save vsd
save(vsd, file = paste0(outdir, "vsd.RData"))

# > vsd all info
vsd_all <- vst(dds, blind=FALSE)

# Save vsd
save(vsd_all, file = paste0(outdir, "vsd_all.RData"))

# > Melting dataset for visualization
melted_norm_counts <-data.frame(melt(assay(vsd_all)))
# rename columns
colnames(melted_norm_counts) <- c("gene", "sample_ID", "normalized_counts")
# join metadata
melted_norm_counts <- left_join(melted_norm_counts, metadata, by = "sample_ID")

# > z-score
var_non_zero <- apply(counts, 1, var) !=0 #  filter out genes that have low variance across samples
filtered_counts <- counts[var_non_zero, ]

zscores <- t(scale(t(filtered_counts)))
dim(zscores) #[1] 18782   165

zscore_mat <- as.matrix(zscores)

# Save melted_norm_counts with z-scores
save(zscore_mat, file = paste0(outdir, "zscore_mat.RData"))  


# ------ Plot PCA -------
# Obtain matrix
mat <- as.matrix(vsd) # col = samples, rows  = genes

# PCA (principal components) analysis
pc <- prcomp(t(mat)) # col = genes , rows = samples

pca_df <- cbind(pc$x[,1:2], metadata[,c(2,7)]) %>% as.data.frame()
pca_df$sample <- rownames(pca_df)

# Acomodar levels
pca_df$Group <- as.factor(pca_df$Group)
levels(pca_df$Group)
# [1] "Ctrl" "SLE" 

levels(pca_df$Cell_type)
# [1] "moDC"     "moDCIMQ"  "monocyte" "tolDC"    "tolDCIMQ"

# Agregar varianza
variance_explained <- pc$sdev^2 / sum(pc$sdev^2) * 100

p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Cell_type)) +
  stat_ellipse() +  # Elipse por 'Cell_type'
  geom_point(aes(shape = Group), size = 2.5) +  # Mantener color por 'Group', pero no en la elipse
  scale_shape_manual(values=c(0, 16)) + 
  scale_color_manual(values = c('palevioletred','#256139','black','gold1','dodgerblue3')) +
  labs(x = sprintf('Principal Component 1 (%.2f%%) of variance', variance_explained[1]),
       y = sprintf('Principal Component 2 (%.2f%%) of variance',variance_explained[2]),
       color = "Cell type",
       shape = "Group") +
  theme_minimal() +
  theme(plot.background = element_rect(fill = "white"),
        text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.key.size = unit(1.5, "lines"))

p


##  -- Asignar la referencia y generar contrastes -----
# Select only control
dds_control <- dds[, control_df$sample_ID] # select columns
dim(dds_control)
# [1] 21816    50

copy_dds <- dds_control
## --- Differential expression analysis (Design 1) ----
# Design 1
# Add comparison
dds_control$Cell_type <-  relevel(dds_control$Cell_type, ref = 'monocyte')
# add design
design(dds_control) <- ~ Cell_type # only cell type because dose is the same
design(dds_control)

# Correr analisis
dds_control <- DESeq(dds_control)
# estimating size factors
# using 'avgTxLength' from assays(dds), correcting for library size
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing
# -- replacing outliers and refitting for 64 genes
# -- DESeq argument 'minReplicatesForReplace' = 7 
# -- original counts are preserved in counts(dds)
# estimating dispersions
# fitting model and testing

resultsNames(dds_control) # lists the coefficients
# [1] "Intercept"                      "Cell_type_moDC_vs_monocyte"    
# [3] "Cell_type_moDCIMQ_vs_monocyte"  "Cell_type_tolDC_vs_monocyte"   
# [5] "Cell_type_tolDCIMQ_vs_monocyte"

# save
save(metadata, dds_control, file = paste0(outdir, 'dds_objects/dds_control_design1.RData'))