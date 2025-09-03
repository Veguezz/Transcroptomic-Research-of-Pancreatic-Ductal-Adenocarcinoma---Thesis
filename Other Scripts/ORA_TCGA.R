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
df <- readRDS("C:/Users/VALERIA/OneDrive - UCB/PDAC/Bases de Datos/TCGA/TCGA_R/res_df_TCGA_V3.rds") 



# Asegurar que los genes sean caracteres
df$hgnc_symbol <- as.character(df$hgnc_symbol)

#Crear una columna que contenga los nombres en ENTREZ ID para que pueda asociarse a las bases de datos Reactome y Kegg
# Convertir símbolos de genes a Entrez ID
df$entrez_id <- mapIds(org.Hs.eg.db, 
                       keys = df$hgnc_symbol,       # Columna con símbolos de genes
                       column = "ENTREZID",  # Queremos obtener Entrez ID
                       keytype = "SYMBOL",   # Nuestros genes están en formato HGNC
                       multiVals = "first")  # Si hay varios, toma el primero

#Verificar si hay genes no mapeados
sum(is.na(df$entrez_id))  # Cuenta cuántos NA hay en entrez_id

# Remover genes sin conversión exitosa
df <- df[!is.na(df$entrez_id), ]

