# ==============================================================================
# PROYECTO: Dinámica Transcriptómica Temporal en Malaria (GSE279789)
# AUTOR: Guillermo Ciesielski Calderon
# DESCRIPCIÓN: Instalación de Dependencias y Preparación del Entorno
# ==============================================================================

cat("-> Comprobando el entorno de paquetes...\n")

# 1. INSTALAR BIOCMANAGER (Gestor principal) -----------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# 2. DEFINIR LISTA DE HERRAMIENTAS NECESARIAS ----------------------------------
paquetes_necesarios <- c(
  "edgeR",           # Modelado estadístico de expresión diferencial
  "UpSetR",          # Visualización de intersecciones
  "clusterProfiler", # Análisis de ontología y rutas biológicas
  "org.Mm.eg.db",    # Diccionario del genoma del ratón
  "ggplot2"          # Motor gráfico subyacente
)

# 3. COMPROBAR CUALES FALTAN POR INSTALAR --------------------------------------
nuevos_paquetes <- paquetes_necesarios[!(paquetes_necesarios %in% installed.packages()[,"Package"])]

# 4. INSTALAR AUSENTES ---------------------------------------------------------
if(length(nuevos_paquetes) > 0) {
  cat(paste("-> Instalando", length(nuevos_paquetes), "paquete(s) faltante(s). Esto puede tardar unos minutos...\n"))
  BiocManager::install(nuevos_paquetes, update = FALSE, ask = FALSE)
  cat("-> Instalación completada.\n")
} else {
  cat("-> ¡Todo listo! Todos los paquetes ya están instalados en tu sistema.\n")
}

cat("-> Entorno preparado con éxito. Puedes proceder a ejecutar el Script 01.\n")