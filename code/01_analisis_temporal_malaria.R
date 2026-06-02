# ==============================================================================
# PROYECTO: Dinámica Transcriptómica Temporal en Malaria (GSE279789)
# AUTOR: Guillermo Ciesielski Calderón
# DESCRIPCIÓN: Pipeline completo de RNA-seq (edgeR) desde conteos crudos hasta
#              la extracción de firmas termporales y visualización UpSet.
# ==============================================================================

# 1. CARGA DE PAQUETES ---------------------------------------------------------

suppressPackageStartupMessages({
    library(edgeR)
    library(org.Mm.eg.db)
    library(AnnotationDbi)
    library(UpSetR)
})

# 2. IMPORTACIÓN DE DATOS Y GENERACIÓN DE METADATOS ----------------------------

cat("-> Cargando matriz de conteos crudos...\n")

# Revisar previamente formato de matriz de recuentos crudos para asignar
# correctamente colnames y rownames y garantizar que todas las columnas son conteos.
counts_matrix <- read.delim("../data/GSE279789_raw_counts.txt",
    header = TRUE,
    row.names = 1
)

nombres_muestras <- colnames(counts_matrix)
condiciones <- sub("_.*", "", nombres_muestras)

# Ordenamos cronologicamente para que el modelo tenga sentido biológico
metadata <- data.frame(
  Muestra = nombres_muestras,
  Grupo = factor(condiciones, levels = c("Naive", "D3", "D5", "D7", "D9", "D11", "D13"))
)

# 3. ANOTACIÓN GENÓMICA ESTRICTA (1:1) -----------------------------------------

cat("-> Anotando identificadores de Ensembl a Símbolos oficiales...\n")

genes_ids <- rownames(counts_matrix)

# Usamos mapIds con "first" para evitar la inflación de filas por duplicados
ann_clean <- data.frame(
  ENSEMBL = genes_ids,
  SYMBOL = mapIds(org.Mm.eg.db, keys = genes_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"),
  GENENAME = mapIds(org.Mm.eg.db, keys = genes_ids, column = "GENENAME", keytype = "ENSEMBL", multiVals = "first"),
  stringsAsFactors = FALSE
)

# 4. PREPROCESAMIENTO: FILTRADO Y NORMALIZACIÓN --------------------------------

cat("-> Construyendo objeto DGEList, filtrando y normalizando...\n")
y <- DGEList(counts = counts_matrix, group = metadata$Grupo)
y$genes <- ann_clean

# Eliminamos genes de baja expresión para aumentar el poder estadístico
# Recalculamos el tamaño de las librerías
keep <- filterByExpr(y)
y <- y[keep, keep.lib.sizes = FALSE]

# Normalización TMM (Trimmed Mean of M-values).
# Calcula factores de normalización para corregir el sesgo de composición de la
# librería (RNA composition bias) y obtener los tamaños efectivos de las librerías.
y <- calcNormFactors(y, method = "TMM")

# 5. MODELADO ESTADÍSTICOY CONTROL DE CALIDAD(GLM Quasi-Likelihood Robusto) ----

cat("-> Entrenando modelo estadístico y generando gráficos de QC...\n")

design <- model.matrix(~ 0 + Grupo, data = metadata)
colnames(design) <- levels(metadata$Grupo)
rownames(design) <- metadata$Muestra # Para comprobar asignacion correcta View()

# Parámetro robust = TRUE vital para modelos in vivo con alta variabilidad
y <- estimateDisp(y, design, robust = TRUE)

# --- QC 1: Gráfico BCV (Biological Coefficient of Variation) ---
# Evalúa la variabilidad biológica global frente a la abundancia del gen
png("../results/plots/QC_plot_BCV.png", width = 1000, height = 800, res = 120)
plotBCV(y, 
        xlab = "Abundancia media de expresión (log2 CPM)", 
        ylab = "Coeficiente de Variación Biológica (BCV)", 
        main = "Dispersión Biológica - Infección por Malaria")
dev.off()

# Ajuste del modelo QL (con protección contra outliers)
fit <- glmQLFit(y, design, robust = TRUE)

# --- QC 2: Gráfico de Dispersión QL (Quasi-Likelihood Fit) ---
# Verifica cómo el modelo encoge las varianzas empíricas hacia la tendencia global
png("../results/plots/QC_plot_QLDisp.png", width = 1000, height = 800, res = 120)
plotQLDisp(fit, 
           xlab = "Abundancia media de expresión (log2 CPM)", 
           ylab = "Dispersión QL (Raíz cuarta de la desviación)", 
           main = "Ajuste Robusto del Modelo QL")
dev.off()


# 6. EXTRACCIÓN TEMPORAL DE DEGs (Todos los días; Separando UP y DOWN) ---------

cat("-> Extrayendo firmas genéticas (UP y DOWN) significativas por día...\n")

# Función auxiliar para automatizar la extracción de cada contraste
# Acepta la dirección de la expresión (UP / DOWN)
extract_degs <- function(fit_model, contrast_name, design_matrix, direction = "up") {
  contraste <- makeContrasts(contrasts = contrast_name, levels = design_matrix)
  res <- glmQLFTest(fit_model, contrast = contraste)
  tabla <- topTags(res, n = Inf)$table
  
  # Filtro estricto: FDR < 0.05 y Log2FoldChange > 1 o < -1
  if (direction == "up") {
    significativos <- tabla[tabla$FDR < 0.05 & tabla$logFC > 1, ]
  } else if (direction == "down") {
    significativos <- tabla[tabla$FDR < 0.05 & tabla$logFC < -1, ]
  }
  return(na.omit(significativos$SYMBOL))
}

# Creamos la lista de Sobreexpresados (UP)
lista_up <- list(
  Dia_03 = extract_degs(fit, "D3 - Naive", design, "up"),
  Dia_05 = extract_degs(fit, "D5 - Naive", design, "up"),
  Dia_07 = extract_degs(fit, "D7 - Naive", design, "up"),
  Dia_09 = extract_degs(fit, "D9 - Naive", design, "up"),
  Dia_11 = extract_degs(fit, "D11 - Naive", design, "up"),
  Dia_13 = extract_degs(fit, "D13 - Naive", design, "up")
)

# Creamos la lista de Infraexpresados (DOWN)
lista_down <- list(
  Dia_03 = extract_degs(fit, "D3 - Naive", design, "down"),
  Dia_05 = extract_degs(fit, "D5 - Naive", design, "down"),
  Dia_07 = extract_degs(fit, "D7 - Naive", design, "down"),
  Dia_09 = extract_degs(fit, "D9 - Naive", design, "down"),
  Dia_11 = extract_degs(fit, "D11 - Naive", design, "down"),
  Dia_13 = extract_degs(fit, "D13 - Naive", design, "down")
)

# 7. VISUALIZACIÓN UPSET (Dos gráficos separados) ------------------------------

cat("-> Dibujando gráficos UpSet separados para UP y DOWN...\n")

# Gráfico 1: Genes Sobreexpresados (UP)
# Control de las intersecciones (Columnas): Mostramos el Top 30 de cruces con más genes
# nintersects = 30
png("../results/plots/upset_6_dias_UP.png", width = 1400, height = 800, res = 120)
upset(fromList(lista_up), nsets = 6, nintersects = 30, order.by = "freq", keep.order = TRUE,
      text.scale = 1.5, main.bar.color = "#d73027",
      mainbar.y.label = "Genes UP compartidos", sets.x.label = "Total Genes UP")
dev.off()

# Gráfico 2: Genes Sobreexpresados (UP)
# Control de las intersecciones (Columnas): Mostramos el Top 30 de cruces con más genes
# nintersects = 30
png("../results/plots/upset_6_dias_DOWN.png", width = 1400, height = 800, res = 120)
upset(fromList(lista_down), nsets = 6, nintersects = 30, order.by = "freq", keep.order = TRUE,
      text.scale = 1.5, main.bar.color = "#4575b4",
      mainbar.y.label = "Genes DOWN compartidos", sets.x.label = "Total Genes DOWN")
dev.off()

# 8. EXPORTACIÓN PARA ANÁLISIS DE ATRIBUTOS (Fase 2) ---------------------------

# Guardamos ambas listas limpias y el modelo ajustado
saveRDS(list(UP = lista_up, DOWN = lista_down), file = "../data/listas_degs_malaria.rds")
saveRDS(fit, file = "../data/modelo_estadistico_fit.rds")

# ==============================================================================