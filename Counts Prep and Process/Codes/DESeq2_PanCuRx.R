#Instalar y cargar los paquetes necesarios

##BiocManager permite a los usuarios instalar y gestionar paquetes del proyecto Bioconductor. Bioconductor se centra en el análisis estadístico y la comprensión de datos genómicos de alto rendimiento.
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

##DESeq2 se utiliza para estimar la varianza-dependencia media en los datos de recuento a partir de ensayos de secuenciación de alto rendimiento y prueba para la expresión diferencial basada en un modelo que utiliza la distribución binomial negativa.
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  BiocManager::install("DESeq2")
}

##readxl se utiliza para poder leer los datos existentes en un archivo excel
if (!requireNamespace("readxl", quietly = TRUE)) {
  install.packages("readxl")
}

# Instalar y cargar biomaRt si no lo tienes ya
if (!requireNamespace("biomaRt", quietly = TRUE)) {
  install.packages("biomaRt")
}

if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
  install.packages("AnnotationDbi")
}
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db")
}

if (!requireNamespace("annotables", quietly = TRUE)) {
  remotes::install_github("stephenturner/annotables")
}

if (!requireNamespace("HGNChelper", quietly = TRUE)) {
  install.packages("HGNChelper")
}

#librerías
library(DESeq2)
library(readxl)
library(biomaRt)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(annotables)
library(HGNChelper)


#Lectura de datos
file_path <- "C:/Users/VALERIA/OneDrive - UCB/PDAC/Bases de Datos/PanCuRx/PanCuRx_Excel/PanCuRx_V4.xlsx"

# Verificar si el archivo existe
if (!file.exists(file_path)) {
  stop("El archivo no existe en la ruta especificada.")
} else {
  message("El archivo fue encontrado en la ruta especificada.")
}

#Verificar los nombres de las hojas
sheets <- excel_sheets(file_path)
print(sheets)

raw_data <- read_excel(file_path, sheet = "Counts")

# Verifica que tu columna de genes esté presente
print(colnames(raw_data))

#Inspecciona las primeras filas de tus datos para verificar su estrucutra
head(raw_data)

#ELIMINAMOS SUFIJOS DE LOS GENES#################################################
raw_data$Gene <- sub("\\.\\d+$", "", raw_data$Gene)  # Remover sufijos numéricos
raw_data$Gene <- gsub("_.*", "", raw_data$Gene)  # Remover sufijos no estándar

# Usa la columna GENE_ID ya limpia para crear
clean_gene <- raw_data$Gene

# Verificar si todavía hay identificadores con sufijos largos
genes_with_suffixes <- clean_gene[grepl("_", clean_gene)]
if (length(genes_with_suffixes) > 0) {
  cat("Hay genes con sufijos restantes:\n")
  print(genes_with_suffixes)
} else {
  cat("No hay genes con sufijos restantes.\n")
}

#Verifica si hay Valores duplicados en la columna de Gene
sum(duplicated(raw_data$Gene))

# Agrega sufijos a los genes duplicados (la primera ocurrencia se mantiene igual)
raw_data$Gene <- make.unique(as.character(raw_data$Gene))

#Elimina los genes duplicados
#raw_data <- raw_data[!duplicated(raw_data$Gene), ]

#Convertir a data.frame para modificar rownames
raw_data <- as.data.frame(raw_data)

#Asignar los rownames
rownames(raw_data) <- raw_data$Gene

#Elimina la columna de los genes
raw_data$Gene <- NULL

#Convierte en una matriz
counts_matrix <- as.matrix(raw_data)

#Verificar si hay valores NA en count_matrix y reemplazarlos con 0 (opcional)
if (any(is.na(counts_matrix))) {
  cat("⚠️ Se encontraron valores NA en count_matrix. Reemplazando con 0...\n")
  counts_matrix[is.na(count_matrix)] <- 0
} else {
  cat("✅ No se encontraron valores NA en count_matrix.\n")
}

# Verificar las primeras filas de la matriz de counts para asegurar su corrección
head(counts_matrix)
head(rownames(counts_matrix)) 

# Obtener el número de muestras
num_samples <- ncol(counts_matrix)

#Crear un dataframe de metadata para DESeq2 (ejemplo: dos condiciones)
#Ajustar la información de la muestra según sea necesario
condition <- c(rep("No Respondedores", 28), rep("Respondedores", 20))
sample_info <- data.frame(
  row.names = colnames(counts_matrix),
  condition = condition  # Ajustar según tus datos
)

#Chequeamos que el dataframe esté generado de manera correcta
print(sample_info)

#Empezamos con nuestro análisis de expresión diferencial con DESeq2
# Creamos un DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = counts_matrix, colData = sample_info, design = ~condition)
dds$condition <- relevel(dds$condition, ref = "Respondedores")

# Corremos el análisis de DESeq2 
dds <- DESeq(dds)

# Obtener los resultados de DESeq2
res <- results(dds)

# Convertir a data frame y agregar los nombres de genes
res_df <- as.data.frame(res)

res_df$gene <- rownames(res_df)  # Esto agrega la fila de nombres como una columna

# Ordenar los resultados por p-value
res_df <- res_df[order(res_df$pvalue), ]

#Verifica si hay filas con valores NA
sum(is.na(res_df$log2FoldChange) | is.na(res_df$pvalue)) 
res_df <- na.omit(res_df)  # Elimina filas con NA en cualquier columna

# Verificar los resultados
head(res_df)

#Guarda tus resultados
saveRDS(res_df, file = "C:/Users/VALERIA/OneDrive - UCB/PDAC/Bases de Datos/PanCuRx/PanCuRx_R/res_df_PanCuRx_V4.rds")

#Guardamos para CIBERSORT
norm_counts <- counts(dds, normalized = TRUE)

saveRDS(norm_counts, file = "C:/Users/VALERIA/OneDrive - UCB/PDAC/Cibersort/norm_dds_counts_PanCuRx_V4.rds")

