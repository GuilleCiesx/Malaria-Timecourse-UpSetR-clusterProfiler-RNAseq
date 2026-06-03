# ==============================================================================
# PROYECTO: Dinámica Transcriptómica Temporal en Malaria (GSE279789)
# AUTOR: Guillermo Ciesielski Calderón
# DESCRIPCIÓN: Validación Biológica de Genes Clave (Cruce con Literatura).
# ==============================================================================

cat("-> Cargando matrices de expresión termporal...\n")

# Cargamos el archivo .rds que guardamos anteriormente
listas <- readRDS("../data/listas_degs_malaria.rds")

# El paper original destaca:
# 1. Señalización de IFN-gamma (Ifng)
# 2. Receptores de células T (Cd8a, Cd8b1)
# 3. Marcadores de monocitos proinflamatorios (Ly6c1, Ly6c2)

genes_candidatos <- c("Ifng", "Cd8a", "Cd8b1", "Ly6c1", "Ly6c2")

cat("\n-> Buscando firmas genéticas del paper original...\n")

# Creamos una tabla limpia a medida
resultado <- data.frame(Gen = genes_candidatos)

# Comprobamos día a día si cada gen está en la lista de sobreexpresados (UP)
# as.integer convierte los TRUE/FALSE en 1/0
resultado$D03 <- as.integer(genes_candidatos %in% listas$UP$Dia_03)
resultado$D05 <- as.integer(genes_candidatos %in% listas$UP$Dia_05)
resultado$D07 <- as.integer(genes_candidatos %in% listas$UP$Dia_07)
resultado$D09 <- as.integer(genes_candidatos %in% listas$UP$Dia_09)
resultado$D11 <- as.integer(genes_candidatos %in% listas$UP$Dia_11)
resultado$D13 <- as.integer(genes_candidatos %in% listas$UP$Dia_13)

# Imprimimos el resultado en consola
print(resultado, row.names = FALSE)

cat("\n-> (Leyenda: 0 = No significativo, 1 = Sobreexpresado)\n")
