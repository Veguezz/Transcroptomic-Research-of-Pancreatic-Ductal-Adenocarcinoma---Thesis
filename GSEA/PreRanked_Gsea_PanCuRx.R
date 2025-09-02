#Pre-Ranked para GSEA

# Cargar datos de DEG ==========================================================
# Reemplaza 'tu_archivo.csv' con el nombre de tu archivo de DEG
df <- readRDS("C:/Users/VALERIA/OneDrive - UCB/PDAC/Bases de Datos/PanCuRx/PanCuRx_R/res_df_PanCuRx_V4.rds") 

#Definir ruta de guardado
# Especificar la ruta donde quieres guardar el archivo .rnk
ruta_guardado <- "C:/Users/VALERIA/OneDrive - UCB/PDAC/Bases de Datos/PanCuRx/PanCuRx_R/PanCuRx_V4_Preranked.rnk"

#Eliminar Duplicados (Si aplica) ===============================================
# Contar los genes duplicados
duplicados <- df$gene[duplicated(df$gene)]

# Mostrar el número de genes duplicados únicos
cat("Número de genes duplicados únicos:", length(unique(duplicados)), "\n")

# Mostrar la lista de genes duplicados únicos
print(unique(duplicados))

# Filtrar genes que aparecen más de una vez
df_duplicados <- df %>%
  group_by(gene) %>%
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
df <- df %>% filter(!gene %in% c("LINC00856", "ENTPD3-AS1", "RCAN3AS"))

# Ranking para GSEA ========================================================
# Crear ranking de genes para fgsea
rankings <- setNames(sign(df$log2FoldChange) * (-log10(df$pvalue)), df$gene)

# Reemplazar valores infinitos para evitar errores en fgsea (Opcional)

#Revisa si tienes infinitos
sum(is.infinite(rankings)) 

# Si sí acomódalos
max_ranking <- max(rankings[is.finite(rankings)], na.rm = TRUE)
min_ranking <- min(rankings[is.finite(rankings)], na.rm = TRUE)
rankings[!is.finite(rankings)] <- ifelse(rankings > 0, max_ranking * 10, min_ranking * 10)

# Ordenar de mayor a menor
rankings <- sort(rankings, decreasing = TRUE)

# Guardar el archico .rnK para GSEA ============================================
# Crear un DataFrame en formato .rnk
rnk_file <- data.frame(
  Gene = names(rankings),  # Nombre del gen
  Score = rankings         # Ranking basado en log2FC y p-valor
)

# Guardar como archivo .rnk (sin encabezados, delimitado por tabulaciones)
write.table(rnk_file, ruta_guardado, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Verificación
readLines("ranking_GSEA_NoResponders.rnk", n = 10)










