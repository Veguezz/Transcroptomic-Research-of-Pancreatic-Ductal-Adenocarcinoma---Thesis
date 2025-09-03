#ORA

if (!requireNamespace("clusterProfiler", quietly = TRUE)) install.packages("clusterProfiler")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db")
if (!requireNamespace("ReactomePA", quietly = TRUE)) BiocManager::install("ReactomePA")


#Cargar librerías
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(dplyr)
library(enrichplot)

#Cargar datos
df <- readRDS("C:/Users/VALERIA/OneDrive - UCB/PDAC/Bases de Datos/PanCuRx/PanCuRx_R/res_df_PanCuRx_V1.rds") 

# Asegurar que los genes sean caracteres
df$gene <- as.character(df$gene)

#Crear una columna que contenga los nombres en ENTREZ ID para que pueda asociarse a las bases de datos Reactome y Kegg
# Convertir símbolos de genes a Entrez ID
df$entrez_id <- mapIds(org.Hs.eg.db, 
                       keys = df$gene,       # Columna con símbolos de genes
                       column = "ENTREZID",  # Queremos obtener Entrez ID
                       keytype = "SYMBOL",   # Nuestros genes están en formato HGNC
                       multiVals = "first")  # Si hay varios, toma el primero

#Verificar si hay genes no mapeados
sum(is.na(df$entrez_id))  # Cuenta cuántos NA hay en entrez_id

# Remover genes sin conversión exitosa
df <- df[!is.na(df$entrez_id), ]

# Filtrar genes significativos
degs_filtered <- df %>%
  filter(padj < 0.1 & abs(log2FoldChange) > 1)  

# Extraer los genes en formato de vector (como caracteres)
gene_list <- as.character(degs_filtered$entrez_id)

# Definir el universo de genes expresados (opcional pero recomendable)
gene_universe <- as.character(df$entrez_id)


# Análisis de enriquecimiento en KEGG y Reactome
#KEGG
kegg_results <- enrichKEGG(gene = gene_list, 
                           organism = "hsa", 
                           universe = gene_universe,
                           pAdjustMethod = "BH", 
                           pvalueCutoff = 0.05)

#Les colocamos nuevamente el símbolo hgnc
gene_symbols <- mapIds(org.Hs.eg.db, 
                       keys = unique(unlist(strsplit(kegg_results@result$geneID, "/"))),  
                       column = "SYMBOL", 
                       keytype = "ENTREZID", 
                       multiVals = "first")

# Reemplazar los números por nombres de genes en el objeto de resultados
kegg_results@result$geneID <- sapply(strsplit(kegg_results@result$geneID, "/"), 
                                     function(x) paste(na.omit(gene_symbols[x]), collapse = "/"))
#REACTOME
reactome_results <- enrichPathway(gene = gene_list, 
                                  organism = "human",
                                  universe = gene_universe,
                                  pAdjustMethod = "BH", 
                                  pvalueCutoff = 0.05)

# Convertir de Entrez ID a HGNC Symbol
gene_symbols <- mapIds(org.Hs.eg.db, 
                       keys = unique(unlist(strsplit(reactome_results@result$geneID, "/"))), 
                       column = "SYMBOL", 
                       keytype = "ENTREZID", 
                       multiVals = "first")

# Reemplazar los números por nombres de genes en el objeto de resultados
reactome_results@result$geneID <- sapply(strsplit(reactome_results@result$geneID, "/"), 
                                         function(x) paste(na.omit(gene_symbols[x]), collapse = "/"))

dim(reactome_results@result)  # Devuelve filas y columnas del dataframe

# GRAFICAS ######################################################################
# Graficar resultados
dotplot(kegg_results, showCategory = 15)  # KEGG pathways
dotplot(reactome_results, showCategory = 15)  # Reactome pathways

# De barras
# KEGG
barplot(kegg_results, showCategory = 10, title = "Top 10 KEGG Pathways")
# Reactome
barplot(reactome_results, showCategory = 10, title = "Top 10 Reactome Pathways")


#Diagrama de vías o Scatter Plot
# Scatter plot (cnetplot) para KEGG - Interacciones de genes en una vía específica
cnetplot(kegg_results, showCategory = 5, categorySize = "pvalue", foldChange = degs_filtered$log2FoldChange)

# Scatter plot para Reactome - Relaciones entre genes en una vía específica
cnetplot(reactome_results, showCategory = 5, categorySize = "pvalue", foldChange = degs_filtered$log2FoldChange)



