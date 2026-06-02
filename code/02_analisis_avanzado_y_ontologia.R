# ==============================================================================
# PROYECTO: Dinámica Transcriptómica Temporal en Malaria (GSE279789)
# AUTOR: Guillermo Ciesielski Calderón
# DESCRIPCIÓN: Análisis de Atributos y Ontología Génica.
# ==============================================================================

# 1. CARGA DE PAQUETES Y DATOS -------------------------------------------------

suppressPackageStartupMessages({
  library(UpSetR)
  library(edgeR)
})

cat("-> Cargando listas de genes y el modelo estadístico...\n")
listas <- readRDS("../data/listas_degs_malaria.rds")
fit <- readRDS("../data/modelo_estadistico_fit.rds")

# 2. PREPARACIÓN DE ATRIBUTOS (Contexto Biológico) -----------------------------
# Para que los diagramas de cajas (boxplots) funcionen, necesitamos extraer
# los valores numéricos (LogFC) de nuestro día más crítico: el día 7

cat("-> Extrayendo estadística del Día 7 para usarla como atributo...\n")

res_D7 <- glmQLFTest(fit, contrast = makeContrasts(D7 - Naive, levels = fit$design))
# Guardamos la tabla completa con todos los genes y sus LogFC
tabla_D7 <- topTags(res_D7, n = Inf)$table

# 3. FUNCIÓN CREADORA DE MATRIZ BINARIA + ATRIBUTO -----------------------------

construir_matriz_atributos <- function(lista_temporal, tabla_stats) {
  
  # PASO A: Extraer un vector con todos los genes únicos de la lista
  genes_unicos <- unique(unlist(lista_temporal))
  
  # PASO B: Crear esqueleto del data.frame
  df <- data.frame(SYMBOL = genes_unicos)
  
  # PASO C: Rellenar matriz binaria (1 = Está en ese día, 0 = No está)
  # El operador '%in%' comprueba si el gen está en la lista existente de ese día
  df$D03 <- ifelse(df$SYMBOL %in% lista_temporal$Dia_03, 1, 0)
  df$D05 <- ifelse(df$SYMBOL %in% lista_temporal$Dia_05, 1, 0)
  df$D07 <- ifelse(df$SYMBOL %in% lista_temporal$Dia_07, 1, 0)
  df$D09 <- ifelse(df$SYMBOL %in% lista_temporal$Dia_09, 1, 0)
  df$D11 <- ifelse(df$SYMBOL %in% lista_temporal$Dia_11, 1, 0)
  df$D13 <- ifelse(df$SYMBOL %in% lista_temporal$Dia_13, 1, 0)
  
  # PASO C: Pegar los atributos
  # match() busca la posición de nuestro gen en la tabla estadística del día 7
  posiciones <- match(df$SYMBOL, tabla_stats$SYMBOL)
  
  # Ahora que sabemos en qué fila está cada gen, nos traemos su LogFC exacto
  df$LogFC_D7 <- tabla_stats$logFC[posiciones]
  
  # Devolvemos matriz completa para UpSetR
  return(df)
}

# 4. APLICAMOS LA FUNCIÓN A NUESTRAS DOS REALIDADES (UP / DOWN) ----------------

cat("-> Construyendo matrices para UP y DOWN...\n")

# Fabricamos la matriz de inflamación (UP)
df_up <- construir_matriz_atributos(listas$UP, tabla_D7)

# Frabricamos la matriz de inflamación (DOWN)
df_down <- construir_matriz_atributos(listas$DOWN, tabla_D7)

# 5. GENERACIÓN DE GRÁFICOS ----------------------------------------------------

cat("-> Dibujando gráficos UpSet con diagramas de cajas...\n")

# -- GRÁFICO 1: LA CASCADA INFLAMATORIA (UP) --
png("../results/plots/upset_atributos_UP.png", width = 1600, height = 1000, res = 130)
upset(
  df_up, # Usamos nuestra matriz compleja en lugar de la función fromList()
  sets = c("D03", "D05", "D07", "D09", "D11", "D13"), 
  keep.order = TRUE, nintersects = 20, order.by = "freq", 
  main.bar.color = "#d73027", text.scale = 1.5,
  mainbar.y.label = "Genes UP Compartidos",
  
  # AÑADIMOS EL ATRIBUTO: Le decimos que coja la columna LogFC_D7 y dibuje un boxplot
  boxplot.summary = c("LogFC_D7")
)
dev.off()

# -- GRÁFICO 2: EL COLAPSO DEL ÓRGANO (DOWN) --
png("../results/plots/upset_atributos_DOWN.png", width = 1600, height = 1000, res = 130)
upset(
  df_down, 
  sets = c("D03", "D05", "D07", "D09", "D11", "D13"), 
  keep.order = TRUE, nintersects = 20, order.by = "freq", 
  main.bar.color = "#4575b4", text.scale = 1.5,
  mainbar.y.label = "Genes DOWN Compartidos",
  
  # En los genes reprimidos, el boxplot saldrá en valores negativos (LogFC < 0)
  boxplot.summary = c("LogFC_D7")
)
dev.off()

cat("Fase 1 completada - Gráficos avanzados exportados.\n")

# 6. ANÁLISIS DE ENRIQUECIMIENTO (Ontología Génica) ----------------------------

# Cargamos la librerías necesarias
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db) 
  library(ggplot2)      
})

# PASO 1: Aislamiento de los genes de interés
# La función Reduce() con intersect vusca los genes que están en TODAS nuestras listas UP.
# Aislamos exactamente los 307 genes que veíamos en el UpSetPlot
genes_core_up <- Reduce(intersect, listas$UP)

# PASO 2: Enrichment Analysis
ego_up <- enrichGO(
  gene = genes_core_up,
  keyType = "SYMBOL",     # Nuestros genes se llaman por su símbolo (ej. Spn)
  OrgDb = org.Mm.eg.db,   # Base de datos del ratón
  ont = "BP",             # Queremos saber los "Biological Processes" (Procesos Biológicos)
  pAdjustMethod = "BH",   # Filtro estadístico estricto
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

# PASO 3: Dibujar Dotplot
cat("-> Generando Dotplot de rutas inflamatorias...\n")

png("../results/plots/GO_Dotplot_Core_UP.png", width = 1200, height = 900, res = 130)

# Un dotplot dibuja las rutas más significativas. 
# El tamaño del punto es el número de genes, el color es el P-valor.
dotplot(ego_up, showCategory = 15, title = "Vías Biológicas Alteradas: Núcleo Inflamatorio") + 
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

dev.off()

# 7. ANÁLISIS DE LA RESPUESTA TARDÍA (Genes exclusivos a partir del Día 7) ------

cat("-> Aislando genes de la respuesta inflamatoria tardía...\n")

# Usamos la matriz df_up para filtrar la barra exacta del UpSet:
# Queremos genes apagados (0) en D03 y D05, pero encendidos (1) desde D07 hasta el final.
genes_tardios <- df_up$SYMBOL[df_up$D03 == 0 & df_up$D05 == 0 & 
                                df_up$D07 == 1 & df_up$D09 == 1 & 
                                df_up$D11 == 1 & df_up$D13 == 1]

cat(paste("Se han aislado", length(genes_tardios), "genes.\n"))

# Pasamos esta nueva lista por la ontología
ego_tardio <- enrichGO(
  gene = genes_tardios,
  keyType = "SYMBOL",
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

# Generamos el segundo Dotplot comparativo
png("../results/plots/GO_Dotplot_Respuesta_Tardia.png", width = 1200, height = 900, res = 130)

dotplot(ego_tardio, 
        showCategory = 15, 
        title = "Vías Alteradas: Tormenta Tardía (Día 7-13)",
        label_format = 45) + # Permite líneas más largas antes de cortarlas
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, color = "#d73027"),
    # AQUÍ ESTÁ EL ARREGLO: Reducimos el tamaño de fuente del eje Y y el interlineado
    axis.text.y = element_text(size = 11, lineheight = 0.7) 
  )

dev.off()

cat("Fase 2 completada y proyecto analítico cerrado\n")

# ==============================================================================