# Subset txi counts
# Evelia Coss
# 3 de Marzo

library(tximport)
library(DESeq2)
library(dplyr)

# ----- Cargar archivos --------
load(file = 'Practica_Dia3/counts/txi_filtered.RData') # filtered_txi
# Cargar metadata
metadata <- read.csv("Practica_Dia3/metadata/metadata.csv", header = F)

# --- Extraer informacion (monocitos) -----

# Filtrar solo monocitos
monocytes_metadata <- metadata %>%
  filter(Cell_type == "monocyte")

set.seed(123)  # Fijar semilla para reproducibilidad
# Seleccionar 3 muestras aleatorias por grupo
selected_samples <- monocytes_metadata %>%
  group_by(Group) %>%
  slice_sample(n = 3) %>%
  pull(sample_ID)

# Crear subset en los datos de tximport
filtered_txi_subset <- list(
  counts = filtered_txi$counts[, selected_samples, drop=FALSE],
  abundance = filtered_txi$abundance[, selected_samples, drop=FALSE],
  length = filtered_txi$length[, selected_samples, drop=FALSE]
)

# Guardar informacion






