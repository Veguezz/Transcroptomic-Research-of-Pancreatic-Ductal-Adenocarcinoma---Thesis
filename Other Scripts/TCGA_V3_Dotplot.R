# Instalación de paquetes
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Paquetes de Bioconductor
BiocManager::install(c("DESeq2", "biomaRt", "pheatmap", "ComplexHeatmap", "EnhancedVolcano"))

# Paquetes de CRAN
install.packages(c("tidyverse", "readr"))

# Cargar librerías
library(tidyverse)
library(readr)
library(DESeq2)
library(biomaRt)
library(pheatmap)
library(ComplexHeatmap)
library(EnhancedVolcano)

# Definir el directorio de trabajo
setwd("C:/Users/VALERIA/OneDrive - UCB/PDAC/RNAseq")
# Cargar el archivo CSV en un DataFrame
conteos <- read.csv("Conteos_preparados.csv")

# Verificar la carga del archivo
head(conteos)

########## Agregamos los símbolos a los genes.
conteos$GENE_ID <- sub("\\.\\d+$", "", conteos$GENE_ID)  # Remover sufijos numéricos
conteos$GENE_ID <- gsub("_.*", "", conteos$GENE_ID)  # Remover sufijos no estándar

# Usa la columna GENE_ID ya limpia para crear
ensg_ids_clean <- conteos$GENE_ID

# Verificar si todavía hay identificadores con sufijos largos
genes_with_suffixes <- ensg_ids_clean[grepl("_", ensg_ids_clean)]
if (length(genes_with_suffixes) > 0) {
  cat("Hay genes con sufijos restantes:\n")
  print(genes_with_suffixes)
} else {
  cat("No hay genes con sufijos restantes.\n")
}

# Contar el número de duplicados en la columna GENE_ID
num_duplicados <- sum(duplicated(conteos$GENE_ID))

# Imprimir el número de duplicados
cat("Número de identificadores GENE_ID duplicados:", num_duplicados, "\n")


# Contar el número de genes únicos en el DataFrame
num_genes_unicos <- length(unique(conteos$GENE_ID))

# Función para conectar a Ensembl probando diferentes espejos
connect_to_ensembl <- function(mirrors) {
  for (mirror in mirrors) {
    message(paste("Intentando conectar al espejo:", mirror))
    ensembl <- tryCatch(
      useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = mirror),
      error = function(e) {
        message(paste("Error al conectar con el espejo:", mirror, "- Probando otro..."))
        return(NULL)
      }
    )
    if (!is.null(ensembl)) {
      message(paste("Conectado exitosamente al espejo:", mirror))
      return(ensembl)
    }
  }
  stop("No se pudo conectar a ningún servidor de Ensembl. Intenta nuevamente más tarde.")
}

# Espejos a probar
mirrors <- c("www", "useast", "asia")

# Conectar a Ensembl
ensembl <- connect_to_ensembl(mirrors)

# Realizar el mapeo de genes
gene_ids <- conteos$GENE_ID  # Usar la columna GENE_ID del DataFrame "conteos"

gene_mapping <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = gene_ids,  # Los identificadores ya deben estar limpios
  mart = ensembl
)

# Imprimir cuántos genes fueron mapeados
cat("Número de genes mapeados:", nrow(gene_mapping), "\n")

# Identificar cuántos no fueron mapeados
genes_no_mapeados <- setdiff(gene_ids, gene_mapping$ensembl_gene_id)
cat("Número de genes no mapeados:", length(genes_no_mapeados), "\n")

# Mostrar los genes no mapeados
cat("Lista de genes no mapeados:\n")
print(genes_no_mapeados)

# Unir el resultado de mapeo con el DataFrame original para agregar la columna de nombres de genes Hugo
conteos <- merge(conteos, gene_mapping, by.x = "GENE_ID", by.y = "ensembl_gene_id", all.x = TRUE)

# Mostrar un resumen del DataFrame actualizado
head(conteos)

# Verificar si hay duplicados en la columna de nombres de genes mapeados (hgnc_symbol)
genes_mapeados_duplicados <- conteos$hgnc_symbol[duplicated(conteos$hgnc_symbol) & !is.na(conteos$hgnc_symbol)]

# Imprimir los genes mapeados que están duplicados
if (length(genes_mapeados_duplicados) > 0) {
  cat("Genes mapeados duplicados:\n")
  print(unique(genes_mapeados_duplicados))
} else {
  cat("No hay genes mapeados duplicados.\n")
}

# Eliminar las filas donde el nombre del gen mapeado (hgnc_symbol) es NA
conteos <- conteos[!is.na(conteos$hgnc_symbol), ]

# Verificar el número de genes restantes después de eliminar los no mapeados
cat("Número de genes después de eliminar los no mapeados:", nrow(conteos), "\n")


# Crear un nuevo DataFrame con "hgnc_symbol" como la primera columna y sin la columna "GENE_ID"
conteos2 <- conteos[, !names(conteos) %in% "GENE_ID"]  # Eliminar la columna "GENE_ID"
conteos2 <- conteos2[, c("hgnc_symbol", setdiff(names(conteos2), "hgnc_symbol"))]  # Reordenar para que "hgnc_symbol" sea la primera columna

# Verificar que el nuevo DataFrame esté en el formato correcto
head(conteos2)

# Cargar la librería DESeq2
library(DESeq2)

# Crear una matriz de conteos sin la columna de nombres de genes (hgnc_symbol)
conteos_matriz <- as.matrix(conteos2[, -which(names(conteos2) == "hgnc_symbol")])
rownames(conteos_matriz) <- conteos2$hgnc_symbol  # Usar los nombres de genes como nombres de filas


# Crear un objeto de diseño para DESeq2 indicando los grupos (Sensible o Resistente)
# Ajusta este vector según tus datos (por ejemplo, las primeras 8 columnas como "Sensible" y el resto como "Resistente")
grupo <- c(rep("Respondedor", 8), rep("No Respondedor", ncol(conteos_matriz) - 8))

is.factor(grupo) #Verificamos en que formato se encuentra

grupo <- as.factor(grupo) #Lo volvemos un factor para definir el grupo de referencia

grupo <- factor(grupo, levels = c("Respondedor", "No Respondedor"))

levels(grupo) #El que aparezca primero es el grupo de referencia. 

# Crear un objeto DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = conteos_matriz, colData = data.frame(grupo), design = ~ grupo)

# Realizar la normalización y análisis de DESeq2
dds <- DESeq(dds)

# Obtener los resultados normalizados
conteos_normalizados <- counts(dds, normalized = TRUE)

# Ver un resumen de los datos normalizados
head(conteos_normalizados)

# Convertir la matriz de conteos normalizados a un DataFrame
conteos_normalizados_df <- as.data.frame(conteos_normalizados)

# Ver un resumen del DataFrame normalizado
head(conteos_normalizados_df)

# Asegurarnos de que tenemos el objeto DESeq2 y los resultados
res <- results(dds)

# Convertir los resultados a un DataFrame si no lo has hecho aún
res_df <- as.data.frame(res)

# Asegurarte de que tienes una columna con los nombres de los genes
res_df$gene <- rownames(res_df)

# Filtrar los genes que tienen nombres identificados y eliminar aquellos que no fueron mapeados correctamente
res_df_filtrado <- res_df %>%
  filter(!grepl("^\\.", gene) & gene != "")  # Elimina genes que comienzan con un punto o que están vacíos

# Verificar el nuevo DataFrame filtrado
head(res_df_filtrado)

# Instalar ggplot2 si no lo tienes ya instalado
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

# Cargar la librería ggplot2
library(ggplot2)

# Paso 1: Especificar el gen de interés
gen_de_interes <- "CYSLTR1"  # Cambia este nombre por el gen que deseas graficar

# Paso 2: Extraer los valores de expresión normalizados del gen de interés
# Asegúrate de que el gen de interés esté en el data frame de conteos normalizados
if (!gen_de_interes %in% rownames(conteos_normalizados_df)) {
  stop(paste("El gen", gen_de_interes, "no se encuentra en los conteos normalizados. Verifica el nombre del gen."))
}

# Extraer los conteos normalizados del gen de interés
conteos_gen_interes <- conteos_normalizados_df[gen_de_interes, ]

# Paso 3: Crear un data frame que contenga los valores de expresión y los grupos
# Convertimos los valores de conteos a un data frame con la información del grupo
df_grafico <- data.frame(
  Expresion = as.numeric(conteos_gen_interes),  # Convertir los conteos a numérico si es necesario
  Grupo = grupo
)

# Verificar los primeros elementos del DataFrame
head(df_grafico)

# Paso 4: Crear el gráfico de dispersión (scatter plot) con ggplot2, eliminando el título y ajustando la proximidad entre grupos
ggplot(df_grafico, aes(x = Grupo, y = Expresion)) +
  geom_jitter(width = 0.15, size = 3, color = "dodgerblue", alpha = 0.7) +  # Ajuste de tamaño, color y transparencia de puntos
  theme_minimal() +
  labs(y = paste("Conteos"), x = NULL) +  # Eliminar título del eje X
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotar etiquetas del eje X
    axis.title.y = element_text(size = 14),  # Ajustar tamaño del texto del eje Y
    axis.text = element_text(size = 12)  # Ajustar tamaño de las etiquetas
  ) +
  ggtitle(paste(gen_de_interes)) +  # Agregar un título descriptivo
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))  # Centrar y estilizar el título


