# ==========================================================
#          Análisis GSEA con fgsea en R
# ==========================================================
# TCGA

# Cargar librerías necesarias
library(tidyverse) 
library(RColorBrewer)
library(ggplot2)
library(fgsea)
library(dplyr)
library(clusterProfiler)  # Para obtener rutas KEGG
library(ReactomePA)  # Para obtener rutas Reactome
library(msigdbr)  # Para usar conjuntos de genes de MSigDB

# Cargar datos de DEG ================================
# Reemplaza 'tu_archivo.csv' con el nombre de tu archivo de DEG
df <- readRDS("C:/Users/VALERIA/OneDrive - UCB/PDAC/Bases de Datos/TCGA/TCGA_R/res_df_TCGA_V3.rds") 

# Asegurar que las columnas necesarias están presentes
if (!all(c("hgnc_symbol", "log2FoldChange", "pvalue") %in% colnames(df))) {
  stop("⚠️ Error: Tu archivo debe contener las columnas 'gene', 'log2FoldChange' y 'pvalue'.")
}

#Eliminar Duplicados============================================================
# Contar los genes duplicados
duplicados <- df$hgnc_symbol[duplicated(df$hgnc_symbol)]

# Mostrar el número de genes duplicados únicos
cat("Número de genes duplicados únicos:", length(unique(duplicados)), "\n")

# Mostrar la lista de genes duplicados únicos
print(unique(duplicados))

# Filtrar genes que aparecen más de una vez
df_duplicados <- df %>%
  group_by(hgnc_symbol) %>%
  filter(n() > 1) %>%
  ungroup()

# Ver el dataframe reducido con solo los duplicados
print(df_duplicados)

# Agregar sufijos a los genes duplicados
df$hgnc_symbol <- ave(df$hgnc_symbol, df$hgnc_symbol, FUN = function(x) {
  if (length(x) > 1) {
    paste0(x, ".", seq_along(x))  # Agrega un sufijo numérico a los duplicados
  } else {
    x
  }
})

#HACER LO MISMO CON "GENE" PORQUE PUEDE PASAR QUE GENES ENSBML SEAN IGUALES PERO HGNC_SYMBOL NO

# Filtrar los genes no deseados por duplicación Ensembl y todas las varibles de DEGs
df <- df %>% filter(!hgnc_symbol %in% c("LINC00856", "ENTPD3-AS1", "RCAN3AS"))

#GSEA===========================================================================

# Crear ranking de genes para fgsea
rankings <- setNames(sign(df$log2FoldChange) * (-log10(df$pvalue)), df$hgnc_symbol)

sum(duplicated(names(rankings))) #Revisar si hay genes duplicados
table(names(rankings))[table(names(rankings)) > 1] #Si sí, cuáles son?

# Reemplazar valores infinitos para evitar errores en fgsea
max_ranking <- max(rankings[is.finite(rankings)], na.rm = TRUE)
min_ranking <- min(rankings[is.finite(rankings)], na.rm = TRUE)
rankings[!is.finite(rankings)] <- ifelse(rankings > 0, max_ranking * 10, min_ranking * 10)


#rankings <- rankings + rnorm(length(rankings), mean = 0, sd = 1e-5)
#rankings <- rankings[!duplicated(names(rankings))]

# Ordenar de mayor a menor
rankings <- sort(rankings, decreasing = TRUE)

# Obtener conjuntos de genes KEGG y Reactome
# Usando msigdbr para obtener términos KEGG y Reactome
kegg_sets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>%
  split(x = .$gene_symbol, f = .$gs_name)

reactome_sets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>%
  split(x = .$gene_symbol, f = .$gs_name)

# Ejecutar análisis GSEA con fgsea 
gsea_kegg <- fgsea(pathways = kegg_sets, stats = rankings, minSize = 10, maxSize = 500, nproc = 1)
gsea_reactome <- fgsea(pathways = reactome_sets, stats = rankings, minSize = 10, maxSize = 500, nproc = 1)

# Guardar resultados
write.csv(gsea_kegg, "GSEA_KEGG_results.csv", row.names = FALSE)
write.csv(gsea_reactome, "GSEA_Reactome_results.csv", row.names = FALSE)


# Visualización de resultados ==========================

# 10 términos más enriquecidos===========================
print(head(gsea_kegg[order(gsea_kegg$pval), ], 10))
print(head(gsea_reactome[order(gsea_reactome$pval), ], 10))

# Seleccionar las 10 vías más enriquecidas (según p-valor)
top_pathways_kegg <- gsea_kegg[order(gsea_kegg$pval), ][1:10, ]
top_pathways_reactome <- gsea_reactome[order(gsea_reactome$pval), ][1:10, ]

# Crear gráfico
ggplot(top_pathways_reactome, aes(x = reorder(pathway, -NES), y = NES, fill = pval)) +
  geom_bar(stat = "identity") +
  coord_flip() +  # Para que las barras sean horizontales
  scale_fill_gradient(low = "red", high = "blue") +
  labs(title = "Top 10 vías enriquecidas (GSEA KEGG)", x = "Vía", y = "NES") +
  theme_minimal()


# Graficar enriquecimiento de un pathway específico

#Añadimos una ruta donde nos guarde los plots
ruta_guardado <- "C:/Users/VALERIA/OneDrive - UCB/PDAC/Bases de Datos/TCGA/TCGA_Imagenes"

# Define el nombre del archivo
nombre_archivo <- paste0("Gsea_TCGA_", pathway_to_plot, ".png")



# Inflamación y microambiente tumoral - Grafica de enriquecimiento ===================================
# REACTOME_INTERFERON_SIGNALING ====
RIS_gsea_plot <- "REACTOME_INTERFERON_SIGNALING"  # Cambia esto por el pathway que te interese
if (RIS_gsea_plot %in% names(reactome_sets)) {
  plotEnrichment(reactome_sets[[RIS_gsea_plot]], rankings) +
    labs(title = paste("Reactome pathway:", RIS_gsea_plot)) +
    theme_classic()
} else {
  print("⚠️ El pathway seleccionado no está en los resultados.")
}

#Extraer los genes enriquecidos de un pathway
# Seleccionar el pathway de interés
RIS_gsea <- "REACTOME_INTERFERON_SIGNALING"  

# Extraer los genes core de este pathway
genes_enriched <- gsea_reactome %>%
  filter(pathway == RIS_gsea) %>%
  pull(leadingEdge) %>%
  unlist()

print(genes_enriched)

# Obtener el ranking de estos genes
# Filtrar los genes enriquecidos en el ranking original
RIS_ranked_genes <- rankings[names(rankings) %in% genes_enriched]

# Ordenar por ranking
RIS_ranked_genes <- RIS_ranked_genes[order(RIS_ranked_genes, decreasing = TRUE)]

print(RIS_ranked_genes)

# Crear un DataFrame con la información de los genes y el NES
# Extraer el enrichment score (NES) del pathway
nes_value <- gsea_reactome %>%
  filter(pathway == RIS_gsea) %>%
  pull(NES)

# Crear un DataFrame con la información
RIS_genes_info <- data.frame(
  Gene = names(RIS_ranked_genes),
  Ranking = RIS_ranked_genes,
  NES = nes_value
)

#Correlacionar con DEGs
# Cruzar los genes enriquecidos con los datos de DEGs (log2FoldChange)
RIS_genes_info <- RIS_genes_info %>%
  left_join(df %>% dplyr::select(hgnc_symbol, log2FoldChange), by = c("Gene" = "hgnc_symbol"))

#Verificar si los genes se unieron correctamente
genes_no_encontrados <- RIS_genes_info$Gene[!(RIS_genes_info$Gene %in% df$hgnc_symbol)]
print(genes_no_encontrados)

df <- df %>%
  mutate(hgnc_symbol = tolower(trimws(hgnc_symbol)))

RIS_genes_info <- RIS_genes_info %>%
  mutate(Gene = tolower(trimws(Gene)))

RIS_genes_info <- RIS_genes_info %>%
  left_join(df %>% dplyr::select(hgnc_symbol, log2FoldChange), by = c("Gene" = "hgnc_symbol"))

# Ver qué genes están realmente upregulated o downregulated
RIS_genes_info <- RIS_genes_info %>%
  mutate(Expression = ifelse(log2FoldChange > 0, "Upregulated", "Downregulated"))

# Mostrar resultados
print(RIS_genes_info)

# Mostrar en pantalla o guardar
print(genes_info)
write.csv(genes_info, paste0("GSEA_", RIS_gsea, "_genes_info.csv"), row.names = FALSE)

# KEGG_IL6_SIGNALING_PATHWAY ====
IL6_gsea_plot <- "KEGG_IL6_SIGNALING_PATHWAY"  # Cambia esto por el pathway que te interese
if (IL6_gsea_plot %in% names(kegg_sets)) {
  plotEnrichment(kegg_sets[[IL6_gsea_plot]], rankings) +
    labs(title = paste("KEGG pathway:", IL6_gsea_plot)) +
    theme_classic()
} else {
  print("⚠️ El pathway seleccionado no está en los resultados.")
}


# Vías clave en PDAC - Grafica de enriquecimiento ===================================
# REACTOME_RAS_SIGNALING_PATHWAY ====
pathway_to_plot <- "REACTOME_RAS_SIGNALING_PATHWAY"  # Cambia esto por el pathway que te interese
if (pathway_to_plot %in% names(reactome_sets)) {
  plotEnrichment(reactome_sets[[pathway_to_plot]], rankings) +
    labs(title = paste("Reactome pathway:", pathway_to_plot)) +
    theme_classic()
} else {
  print("⚠️ El pathway seleccionado no está en los resultados.")
}


# KEGG_TGF_BETA_SIGNALING_PATHWAY =====
TGFB_plot <- "KEGG_TGF_BETA_SIGNALING_PATHWAY"  # Cambia esto por el pathway que te interese
if (TGFB_plot  %in% names(kegg_sets)) {
  plotEnrichment(kegg_sets[[TGFB_plot ]], rankings) +
    labs(title = paste("KEGG pathway:", TGFB_plot )) +
    theme_classic()
} else {
  print("⚠️ El pathway seleccionado no está en los resultados.")
}

# Seleccionar el pathway de interés
TGFB_gsea <- "KEGG_TGF_BETA_SIGNALING_PATHWAY"  

# Extraer los genes core de este pathway
TGFB_genes_enriched <- gsea_kegg %>%
  filter(pathway == TGFB_gsea) %>%
  pull(leadingEdge) %>%
  unlist()

print(TGFB_genes_enriched)

# Obtener el ranking de estos genes
# Filtrar los genes enriquecidos en el ranking original
TGFB_ranked_genes <- rankings[names(rankings) %in% TGFB_genes_enriched]

# Ordenar por ranking
TGFB_ranked_genes <- TGFB_ranked_genes[order(TGFB_ranked_genes, decreasing = TRUE)]

print(TGFB_ranked_genes)

# Crear un DataFrame con la información de los genes y el NES
# Extraer el enrichment score (NES) del pathway
TGFB_nes_value <- gsea_kegg %>%
  filter(pathway == TGFB_gsea) %>%
  pull(NES)

# Crear un DataFrame con la información
TGFB_genes_info <- data.frame(
  Gene = names(TGFB_ranked_genes),
  Ranking = TGFB_ranked_genes,
  NES = TGFB_nes_value
)

#Correlacionar con DEGs
#Convertiremos todo a minúsculas y eliminaremos espacios ocultos entre los nombres de los genes
df <- df %>%
  mutate(hgnc_symbol = tolower(trimws(hgnc_symbol)))

TGFB_genes_info <- TGFB_genes_info %>%
  mutate(Gene = tolower(trimws(Gene)))

# Cruzar los genes enriquecidos con los datos de DEGs (log2FoldChange)
TGFB_genes_info <- TGFB_genes_info %>%
  left_join(df %>% dplyr::select(hgnc_symbol, log2FoldChange), by = c("Gene" = "hgnc_symbol"))

# Ver qué genes están realmente upregulated o downregulated
TGFB_genes_info <- TGFB_genes_info %>%
  mutate(Expression = ifelse(log2FoldChange > 0, "Upregulated", "Downregulated"))


# Respuesta inmunitaria general - Grafica de enriquecimiento ===================================
# REACTOME_INNATE_IMMUNE_SYSTEM ====
RIIN_plot <- "REACTOME_INNATE_IMMUNE_SYSTEM"  # Cambia esto por el pathway que te interese
if (RIIN_plot %in% names(reactome_sets)) {
  plotEnrichment(reactome_sets[[RIIN_plot]], rankings) +
    labs(title = paste("Reactome pathway:", RIIN_plot)) +
    theme_classic()
} else {
  print("⚠️ El pathway seleccionado no está en los resultados.")
}

#Extraer los genes enriquecidos de un pathway
# Seleccionar el pathway de interés
RIIN_gsea <- "REACTOME_INNATE_IMMUNE_SYSTEM"  

# Extraer los genes core de este pathway
RIIN_genes_enriched <- gsea_reactome %>%
  filter(pathway == RIIN_gsea) %>%
  pull(leadingEdge) %>%
  unlist()

print(RIIN_genes_enriched)

# Obtener el ranking de estos genes
# Filtrar los genes enriquecidos en el ranking original
RIIN_ranked_genes <- rankings[names(rankings) %in% RIIN_genes_enriched]

# Ordenar por ranking
RIIN_ranked_genes <- RIIN_ranked_genes[order(RIIN_ranked_genes, decreasing = TRUE)]

print(RIIN_ranked_genes)

# Crear un DataFrame con la información de los genes y el NES
# Extraer el enrichment score (NES) del pathway
RIIN_nes_value <- gsea_reactome %>%
  filter(pathway == RIIN_gsea) %>%
  pull(NES)

# Crear un DataFrame con la información
RIIN_genes_info <- data.frame(
  Gene = names(RIIN_ranked_genes),
  Ranking = RIIN_ranked_genes,
  NES = RIIN_nes_value
)

#Correlacionar con DEGs
# Convertir genes en minúsculas
df <- df %>%
  mutate(hgnc_symbol = tolower(trimws(hgnc_symbol)))


RIIN_genes_info <- RIIN_genes_info %>%
  mutate(Gene = tolower(trimws(Gene)))

# Cruzar los genes enriquecidos con los datos de DEGs (log2FoldChange)
RIIN_genes_info <- RIIN_genes_info %>%
  left_join(df %>% dplyr::select(hgnc_symbol, log2FoldChange), by = c("Gene" = "hgnc_symbol"))

#Verificar si los genes se unieron correctamente
genes_no_encontrados <- RIIN_genes_info$Gene[!(RIIN_genes_info$Gene %in% df$hgnc_symbol)]
print(genes_no_encontrados)

# Ver qué genes están realmente upregulated o downregulated
RIIN_genes_info <- RIIN_genes_info %>%
  mutate(Expression = ifelse(log2FoldChange > 0, "Upregulated", "Downregulated"))

# Mostrar resultados
print(RIIN_genes_info)

# Mostrar en pantalla o guardar
print(RIIN_genes_info)
write.csv(genes_info, paste0("GSEA_", RIIN_gsea, "_genes_info.csv"), row.names = FALSE)






#

pathway_to_plot <- "KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY"  # Cambia esto por el pathway que te interese
if (pathway_to_plot %in% names(kegg_sets)) {
  plotEnrichment(kegg_sets[[pathway_to_plot]], rankings) +
    labs(title = paste("KEGG pathway:", pathway_to_plot)) +
    theme_classic()
} else {
  print("⚠️ El pathway seleccionado no está en los resultados.")
}

#Vías inmunosupresoras
pathway_to_plot <- "KEGG_PDL1_PD1_CHECKPOINT_PATHWAY_IN_CANCER"  # Cambia esto por el pathway que te interese
if (pathway_to_plot %in% names(kegg_sets)) {
  plotEnrichment(kegg_sets[[pathway_to_plot]], rankings) +
    labs(title = paste("KEGG pathway:", pathway_to_plot)) +
    theme_classic()
} else {
  print("⚠️ El pathway seleccionado no está en los resultados.")
}

pathway_to_plot <- "REACTOME_TGF_BETA_RECEPTOR_SIGNALING"  # Cambia esto por el pathway que te interese
if (pathway_to_plot %in% names(reactome_sets)) {
  plotEnrichment(reactome_sets[[pathway_to_plot]], rankings) +
    labs(title = paste("Reactome pathway:", pathway_to_plot)) +
    theme_classic()
} else {
  print("⚠️ El pathway seleccionado no está en los resultados.")
}


# Infiltración de linfocitos T citotóxicos - Grafica de enriquecimiento ===================================
pathway_to_plot <- "KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY"  # Cambia esto por el pathway que te interese
if (pathway_to_plot %in% names(kegg_sets)) {
  plotEnrichment(kegg_sets[[pathway_to_plot]], rankings) +
    labs(title = paste("KEGG pathway:", pathway_to_plot)) +
    theme_classic()
} else {
  print("⚠️ El pathway seleccionado no está en los resultados.")
}

pathway_to_plot <- "REACTOME_INTERLEUKIN_SIGNALING"  # Cambia esto por el pathway que te interese
if (pathway_to_plot %in% names(reactome_sets)) {
  plotEnrichment(reactome_sets[[pathway_to_plot]], rankings) +
    labs(title = paste("Reactome pathway:", pathway_to_plot)) +
    theme_classic()
} else {
  print("⚠️ El pathway seleccionado no está en los resultados.")
}










