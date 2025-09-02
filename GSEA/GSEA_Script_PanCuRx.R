# ==========================================================
#          Análisis GSEA con fgsea en R
# ==========================================================
# PanCuRx

# Cargar librerías necesarias
library(tidyverse) 
library(RColorBrewer)
library(ggplot2)
library(fgsea)
library(clusterProfiler)  # Para obtener rutas KEGG
library(ReactomePA)  # Para obtener rutas Reactome
library(msigdbr)  # Para usar conjuntos de genes de MSigDB

# 📌 2. Cargar datos de DEG ================================
# Reemplaza 'tu_archivo.csv' con el nombre de tu archivo de DEG
df <- readRDS("C:/Users/VALERIA/OneDrive - UCB/PDAC/Bases de Datos/PanCuRx/PanCuRx_R/res_df_PanCuRx_V1.rds") 

# Asegurar que las columnas necesarias están presentes
if (!all(c("gene", "log2FoldChange", "pvalue") %in% colnames(df))) {
  stop("⚠️ Error: Tu archivo debe contener las columnas 'gene', 'log2FoldChange' y 'pvalue'.")
}

# Crear ranking de genes para fgsea
rankings <- setNames(sign(df$log2FoldChange) * (-log10(df$pvalue)), df$gene)

# Reemplazar valores infinitos para evitar errores en fgsea
max_ranking <- max(rankings[is.finite(rankings)], na.rm = TRUE)
min_ranking <- min(rankings[is.finite(rankings)], na.rm = TRUE)
rankings[!is.finite(rankings)] <- ifelse(rankings > 0, max_ranking * 10, min_ranking * 10)

# Ordenar de mayor a menor
rankings <- sort(rankings, decreasing = TRUE)

# 📌 3. Obtener conjuntos de genes KEGG y Reactome ==========
# Usando msigdbr para obtener términos KEGG y Reactome
kegg_sets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>%
  split(x = .$gene_symbol, f = .$gs_name)

reactome_sets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>%
  split(x = .$gene_symbol, f = .$gs_name)

# 📌 4. Ejecutar análisis GSEA con fgsea =====================
gsea_kegg <- fgsea(pathways = kegg_sets, stats = rankings, minSize = 10, maxSize = 500, nproc = 1)
gsea_reactome <- fgsea(pathways = reactome_sets, stats = rankings, minSize = 10, maxSize = 500, nproc = 1)

# 📌 5. Guardar resultados ===================================
write.csv(gsea_kegg, "GSEA_KEGG_results.csv", row.names = FALSE)
write.csv(gsea_reactome, "GSEA_Reactome_results.csv", row.names = FALSE)

# 📌 6. Visualización de resultados ==========================
# Mostrar los 10 términos más enriquecidos
print(head(gsea_kegg[order(gsea_kegg$pval), ], 10))
print(head(gsea_reactome[order(gsea_reactome$pval), ], 10))

# Seleccionar las 10 vías más enriquecidas (según p-valor)
top_pathways_kegg <- gsea_kegg[order(gsea_kegg$pval), ][1:10, ]
top_pathways_reactome <- gsea_reactome[order(gsea_reactome$pval), ][1:10, ]

# Crear gráfico
ggplot(top_pathways_kegg, aes(x = reorder(pathway, -NES), y = NES, fill = pval)) +
  geom_bar(stat = "identity") +
  coord_flip() +  # Para que las barras sean horizontales
  scale_fill_gradient(low = "red", high = "blue") +
  labs(title = "Top 10 vías enriquecidas (GSEA KEGG)", x = "Vía", y = "NES") +
  theme_minimal()



# Vías de reparación del ADN - Grafica de enriquecimiento ===================================
REACTOME_DNA_REPAIR <- "REACTOME_DNA_REPAIR"  # Cambia esto por el pathway que te interese
if (REACTOME_DNA_REPAIR %in% names(reactome_sets)) {
  plotEnrichment(reactome_sets[[REACTOME_DNA_REPAIR]], rankings) +
    labs(title = paste("Reactome pathway:", REACTOME_DNA_REPAIR)) +
    theme_classic()
} else {
  print("⚠️ El pathway seleccionado no está en los resultados.")
}


#
pathway_to_plot <- "KEGG_MISMATCH_REPAIR"  # Cambia esto por el pathway que te interese
if (pathway_to_plot %in% names(kegg_sets)) {
  plotEnrichment(kegg_sets[[pathway_to_plot]], rankings) +
    labs(title = paste("KEGG pathway:", pathway_to_plot)) +
    theme_classic()
} else {
  print("⚠️ El pathway seleccionado no está en los resultados.")
}

pathway_to_plot <- "KEGG_HOMOLOGOUS_RECOMBINATION"  # Cambia esto por el pathway que te interese
if (pathway_to_plot %in% names(kegg_sets)) {
  plotEnrichment(kegg_sets[[pathway_to_plot]], rankings) +
    labs(title = paste("KEGG pathway:", pathway_to_plot)) +
    theme_classic()
} else {
  print("⚠️ El pathway seleccionado no está en los resultados.")
}

# Vías de Detoxificación de fármacos y metabolismo - Grafica de enriquecimiento ===================================
pathway_to_plot <- "REACTOME_PHASE_I_FUNCTIONALIZATION_OF_COMPOUNDS"  # Cambia esto por el pathway que te interese
if (pathway_to_plot %in% names(reactome_sets)) {
  plotEnrichment(reactome_sets[[pathway_to_plot]], rankings) +
    labs(title = paste("Reactome pathway:", pathway_to_plot)) +
    theme_classic()
} else {
  print("⚠️ El pathway seleccionado no está en los resultados.")
}

pathway_to_plot <- "KEGG_DRUG_METABOLISM_CYTOCHROME_P450"  # Cambia esto por el pathway que te interese
if (pathway_to_plot %in% names(kegg_sets)) {
  plotEnrichment(kegg_sets[[pathway_to_plot]], rankings) +
    labs(title = paste("KEGG pathway:", pathway_to_plot)) +
    theme_classic()
} else {
  print("⚠️ El pathway seleccionado no está en los resultados.")
}

# Resistencia a apoptosis - Grafica de enriquecimiento ===================================
pathway_to_plot <- "KEGG_P53_SIGNALING_PATHWAY"  # Cambia esto por el pathway que te interese
if (pathway_to_plot %in% names(kegg_sets)) {
  plotEnrichment(kegg_sets[[pathway_to_plot]], rankings) +
    labs(title = paste("KEGG pathway:", pathway_to_plot)) +
    theme_classic()
} else {
  print("⚠️ El pathway seleccionado no está en los resultados.")
}

# Inflamación y microambiente tumoral - Grafica de enriquecimiento ===================================
REACTOME_INTERFERON_SIGNALING <- "REACTOME_INTERFERON_SIGNALING"  # Cambia esto por el pathway que te interese
if (REACTOME_INTERFERON_SIGNALING %in% names(reactome_sets)) {
  plotEnrichment(reactome_sets[[REACTOME_INTERFERON_SIGNALING]], rankings) +
    labs(title = paste("Reactome pathway:", REACTOME_INTERFERON_SIGNALING)) +
    theme_classic()
} else {
  print("⚠️ El pathway seleccionado no está en los resultados.")
}

# Conocer los genes enriquecidos
names(gsea_reactome$leadingEdge) <- gsea_reactome$pathway
# Asegurar que la vía está en los resultados de fgsea
if ("REACTOME_INTERFERON_SIGNALING" %in% gsea_reactome$pathway) {
  
  # Obtener genes enriquecidos
  enriched_genes <- gsea_reactome$leadingEdge[["REACTOME_INTERFERON_SIGNALING"]]
  
  # Verificar que se obtuvieron genes
  if (!is.null(enriched_genes) && length(enriched_genes) > 0) {
    
    # Crear dataframe con rank y enrichment score
    df_enriquecidos <- data.frame(
      gene = enriched_genes,
      rank = match(enriched_genes, names(rankings)),  # Ranking en la lista ordenada
      enrichment_score = rankings[enriched_genes]  # Obtener ES desde rankings
    )
    
    # Ordenar por ranking
    df_enriquecidos <- df_enriquecidos %>% arrange(rank)
    
    # Mostrar primeros genes enriquecidos
    print(head(df_enriquecidos))
  } else {
    print("⚠️ No hay genes enriquecidos en esta vía.")
  }
} else {
  print("⚠️ La vía seleccionada no está en los resultados de fgsea.")
}




#





pathway_to_plot <- "KEGG_IL6_SIGNALING_PATHWAY"  # Cambia esto por el pathway que te interese
if (pathway_to_plot %in% names(kegg_sets)) {
  plotEnrichment(kegg_sets[[pathway_to_plot]], rankings) +
    labs(title = paste("KEGG pathway:", pathway_to_plot)) +
    theme_classic()
} else {
  print("⚠️ El pathway seleccionado no está en los resultados.")
}

# Vías clave en PDAC - Grafica de enriquecimiento ===================================
pathway_to_plot <- "KEGG_PANCREATIC_CANCER"  # Cambia esto por el pathway que te interese
if (pathway_to_plot %in% names(kegg_sets)) {
  plotEnrichment(kegg_sets[[pathway_to_plot]], rankings) +
    labs(title = paste("KEGG pathway:", pathway_to_plot)) +
    theme_classic()
} else {
  print("⚠️ El pathway seleccionado no está en los resultados.")
}

pathway_to_plot <- "KEGG_TGF_BETA_SIGNALING_PATHWAY"  # Cambia esto por el pathway que te interese
if (pathway_to_plot %in% names(kegg_sets)) {
  plotEnrichment(kegg_sets[[pathway_to_plot]], rankings) +
    labs(title = paste("KEGG pathway:", pathway_to_plot)) +
    theme_classic()
} else {
  print("⚠️ El pathway seleccionado no está en los resultados.")
}


# Vías de metástasis y migración celular - Grafica de enriquecimiento ===================================
pathway_to_plot <- "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION"  # Cambia esto por el pathway que te interese
if (pathway_to_plot %in% names(reactome_sets)) {
  plotEnrichment(reactome_sets[[pathway_to_plot]], rankings) +
    labs(title = paste("Reactome pathway:", pathway_to_plot)) +
    theme_classic()
} else {
  print("⚠️ El pathway seleccionado no está en los resultados.")
}


pathway_to_plot <- "KEGG_CELL_ADHESION_MOLECULES_CAMS"  # Cambia esto por el pathway que te interese
if (pathway_to_plot %in% names(kegg_sets)) {
  plotEnrichment(kegg_sets[[pathway_to_plot]], rankings) +
    labs(title = paste("KEGG pathway:", pathway_to_plot)) +
    theme_classic()
} else {
  print("⚠️ El pathway seleccionado no está en los resultados.")
}


pathway_to_plot <- " KEGG_EPIDERMAL_GROWTH_FACTOR_RECEPTOR_SIGNALING_PATHWAY"  # Cambia esto por el pathway que te interese
if (pathway_to_plot %in% names(kegg_sets)) {
  plotEnrichment(kegg_sets[[pathway_to_plot]], rankings) +
    labs(title = paste("KEGG pathway:", pathway_to_plot)) +
    theme_classic()
} else {
  print("⚠️ El pathway seleccionado no está en los resultados.")
}

# Respuesta inmunitaria general - Grafica de enriquecimiento ===================================
pathway_to_plot <- "REACTOME_INNATE_IMMUNE_SYSTEM"  # Cambia esto por el pathway que te interese
if (pathway_to_plot %in% names(reactome_sets)) {
  plotEnrichment(reactome_sets[[pathway_to_plot]], rankings) +
    labs(title = paste("Reactome pathway:", pathway_to_plot)) +
    theme_classic()
} else {
  print("⚠️ El pathway seleccionado no está en los resultados.")
}

pathway_to_plot <- "KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY"  # Cambia esto por el pathway que te interese
if (pathway_to_plot %in% names(kegg_sets)) {
  plotEnrichment(kegg_sets[[pathway_to_plot]], rankings) +
    labs(title = paste("KEGG pathway:", pathway_to_plot)) +
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


# Conocer los genes enriquecidos =======================================
genes_enriquecidos <- gsea_kegg$leadingEdge[["KEGG_TGF_BETA_SIGNALING_PATHWAY"]]

# Convertir a dataframe
df_enriquecidos <- data.frame(
  gene = genes_enriquecidos,
  rank = match(genes_enriquecidos, names(ranked_genes)),  # Obtener su ranking
  enrichment_score = ranked_genes[genes_enriquecidos]     # Obtener su ES
)

# Ordenar por ranking
df_enriquecidos <- df_enriquecidos %>% arrange(rank)

# Ver los primeros genes enriquecidos
print(head(df_enriquecidos))









