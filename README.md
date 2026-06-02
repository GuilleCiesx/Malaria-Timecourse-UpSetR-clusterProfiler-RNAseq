# Dinámica Transcriptómica Temporal en Malaria (GSE279789)

![R](https://img.shields.io/badge/R-276DC3?style=for-the-badge&logo=r&logoColor=white)
![edgeR](https://img.shields.io/badge/edgeR-BioConductor-green?style=for-the-badge)
![UpSetR](https://img.shields.io/badge/UpSetR-Data_Viz-orange?style=for-the-badge)
![clusterProfiler](https://img.shields.io/badge/clusterProfiler-Ontology-purple?style=for-the-badge)

## 📖 Descripción del Proyecto
Este proyecto reproduce y amplía el análisis transcriptómico del tejido pulmonar murino infectado por *Plasmodium chabaudi*, utilizando el dataset público **GSE279789**. 

El objetivo principal es rastrear la evolución temporal de la expresión génica desde la fase temprana de la infección (Día 3) hasta el pico letal de la enfermedad (Día 13), desentrañando la cascada biológica que conduce a la tormenta de citoquinas y al colapso del órgano.

---
## ⚙️ Metodología y Herramientas Analíticas
Para garantizar la robustez estadística y la interpretabilidad biológica, el pipeline se dividió en dos fases:

**Modelado Estadístico (`edgeR`):** Se ajustaron modelos lineales generalizados (GLM) basados en cuasi-verosimilitud (`glmQLFit`) con el parámetro `robust = TRUE` para controlar la alta variabilidad biológica intrínseca a los modelos *in vivo*. Se evaluó la calidad del ajuste mediante gráficos de dispersión biológica (BCV) y QL-Dispersion.

Dada la alta variabilidad intrínseca de los modelos de infección *in vivo*, el control de la dispersión es un paso crítico antes de evaluar la expresión diferencial. 

1. **Variabilidad Biológica (BCV):** El gráfico de dispersión (izquierda) muestra el Coeficiente de Variación Biológica. Como es esperable en tejido pulmonar completo sometido a una infección sistémica, la dispersión general es alta. 
2. **Ajuste de Cuasi-verosimilitud (QL):** Para evitar una tasa elevada de falsos positivos (genes que parecen alterados solo por ruido biológico), se aplicó el modelo estadístico robusto de `edgeR` (derecha). Este método "comprime" las estimaciones de dispersión individuales hacia una tendencia general (línea azul), penalizando los genes con variabilidad anómala.

<p align="center">
  <img src="results/plots/QC_plot_BCV.png" width="45%" />
  <img src="results/plots/QC_plot_QLDisp.png" width="45%" />
</p>
 
 **Visualización Compleja y Ontología (`UpSetR` y `clusterProfiler`):** Tradicionalmente, las intersecciones de conjuntos genéticos se han representado mediante diagramas de Venn o Euler. Sin embargo, estas representaciones resultan inadecuadas y difíciles de interpretar cuando se maneja un número elevado de conjuntos experimentales (Conway et al., 2017). Para superar esta limitación geométrica, se implementó el paquete `UpSetR`, empleando una visualización escalable basada en matrices. Esto permitió integrar el *Log2 Fold Change* del pico inflamatorio (Día 7) como diagramas de cajas superpuestos. Finalmente, se mapearon las firmas genéticas a la base de datos *Gene Ontology (GO)*.

---
## 📊 Resultados Clave: La evolución de la patología

### 1. La Respuesta Hiperinflamatoria Tardía (Sobreexpresión)
Al integrar los niveles de expresión (LogFC) del **Día 7** como *atributos* sobre las intersecciones temporales, observamos dos dinámicas clave:
* **El núcleo duro (307 genes):** Se activa en el Día 3 y se mantiene constitutivamente encendido, alcanzando los niveles más extremos de sobreexpresión durante el pico letal. Actúa como el motor destructivo de la enfermedad. 
* **La segunda oleada (599 genes):** Un grupo masivo que se activa a partir del Día 5, sumándose a la carga inflamatoria basal.

![UpSet Inflamación](results/plots/upset_atributos_UP.png)

### 2. El Colapso Pulmonar (Infraexpresión)
La represión génica (apagado de funciones basales del tejido) es un evento **estrictamente tardío**. Las intersecciones más masivas de genes infraexpresados pertenecen exclusivamente a los Días 9, 11 y 13, demostrando que el pulmón no falla hasta que la infiltración inmunitaria supera un umbral crítico.

![UpSet Represión](results/plots/upset_atributos_DOWN.png)

### 3. Ontología Génica (Dinámica Inmune)
El Análisis de Enriquecimiento (GO) reveló una clara transición en la fisiología del ratón:
1. **Fase de Expansión (Núcleo D3-D13):** Dominada por vías de proliferación celular y toxicidad mediada por leucocitos. El sistema inmune se multiplica para combatir el parásito.
2. **Fase de Tormenta (Exclusivos D7-D13):** El foco cambia hacia la *regulación de la respuesta inflamatoria* y la *inmunidad adaptativa*, evidenciando el fracaso del organismo para frenar el daño colateral.

<p align="center">
  <img src="results/plots/GO_Dotplot_Core_UP.png" width="45%" />
  <img src="results/plots/GO_Dotplot_Respuesta_Tardia.png" width="45%" />
</p>

---

## 💻 Reproducibilidad y Estructura del Repositorio
Todo el análisis ha sido programado en R. Para reproducir los resultados:
1. Clona este repositorio.
2. Ejecuta los scripts en la carpeta `code/` en orden:
   * `00_instalar_dependencias.R`   
   * `01_analisis_temporal_malaria.R`
   * `02_analisis_avanzado_y_ontologia.R`

---

## 📚 Referencias y Bibliografía
Este análisis se apoya en el desarrollo metodológico de las siguientes herramientas y publicaciones:
* **UpSetR:** Conway, J. R., Lex, A., & Gehlenborg, N. (2017). UpSetR: an R package for the visualization of intersecting sets and their properties. *Bioinformatics*, 33(18), 2938-2940.
* **clusterProfiler:** Yu, G., Wang, L. G., Han, Y., & He, Q. Y. (2012). clusterProfiler: an R package for comparing biological themes among gene clusters. *OMICS: A Journal of Integrative Biology*, 16(5), 284-287.
* **edgeR:** Robinson, M. D., McCarthy, D. J., & Smyth, G. K. (2010). edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. *Bioinformatics*, 26(1), 139-140.
* **Datos Originales (Malaria):** Chen, S. S., Yang, Q., Zhong, Y., Liu, D., Zhou, L., Wei, H. C., ... & Lin, J. W. (2026). Comprehensive immune profiling reveals IFN-γ signaling in T cells mediates parasite phagocytosis in a rodent malaria model. _Mbio_, 17(4), e03938-25.