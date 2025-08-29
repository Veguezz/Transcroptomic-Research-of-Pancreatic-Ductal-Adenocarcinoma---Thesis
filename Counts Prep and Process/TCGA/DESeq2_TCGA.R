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
library(EnhancedVolcano)
library(readxl)
library(biomaRt)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(annotables)
library(HGNChelper)

### PREPARACIÓN DE DATOS ##############

#Lectura de datos
file_path <- "C:/Users/VALERIA/OneDrive - UCB/PDAC/Bases de Datos/TCGA/TCGA_Excel/TCGA_V3.xlsx"

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



#ELIMINAMOS SUFIJOS DE LOS GENES
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

#Elimina los genes duplicados
raw_data <- raw_data[!duplicated(raw_data$Gene), ]

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

### DESeq2 ######

#Crear un dataframe de metadata para DESeq2 (ejemplo: dos condiciones)
#Ajustar la información de la muestra según sea necesario
condition <- c(rep("No Respondedores", 15), rep("Respondedores", 8))
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

### AGREGAR SÍMBOLO HGNC ########
# Definir la función para conectar a Ensembl usando diferentes espejos
connect_to_ensembl <- function(mirrors) {
  for (mirror in mirrors) {
    message(paste("Intentando conectar al espejo:", mirror))
    ensembl <- tryCatch(
      useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = mirror),
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


# Conectar a Ensembl
mirrors <- c("www", "useast", "asia")
ensembl <- connect_to_ensembl(mirrors)

# Extraer los IDs Ensembl de las rownames
ensembl_ids <- rownames(res_df)

# Obtener información de los símbolos HGNC
gene_info <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = ensembl
)

# Convertir rownames en una columna para hacer merge()
res_df$ensembl_gene_id <- rownames(res_df)

# Unir la información obtenida con el dataframe de DESeq2
res_df <- merge(res_df, gene_info, by = "ensembl_gene_id", all.x = TRUE)

# Reemplazar los rownames con los símbolos HGNC, manteniendo los genes sin símbolo
new_rownames <- ifelse(is.na(res_df$hgnc_symbol) | res_df$hgnc_symbol == "", 
                       res_df$ensembl_gene_id, 
                       res_df$hgnc_symbol)

# Asegurar que los rownames sean únicos
rownames(res_df) <- make.unique(new_rownames)

# Eliminar la columna ensembl_gene_id (ya está en rownames)
res_df$gene <- NULL

# Verificar los resultados
head(res_df)

#Eliminar genes NA o sin mapear
#Verificar si hay genes sin mapear
unmapped_genes <- res_df$hgnc_symbol[is.na(res_df$hgnc_symbol)]
if (length(unmapped_genes) > 0) {
  cat("Genes no mapeados encontrados:\n")
  print(unmapped_genes)
} else {
  cat("Todos los genes fueron mapeados correctamente.\n")
}

#Elimina los genes no mapeados o con NA
res_df_mapped <- res_df %>% filter(!is.na(hgnc_symbol) & hgnc_symbol != "")

#Guarda tus resultados
saveRDS(res_df_mapped, file = "C:/Users/VALERIA/OneDrive - UCB/PDAC/Bases de Datos/TCGA/TCGA_R/res_df_TCGA_V3.2.rds")
write.csv(res_df_mapped,"C:/Users/VALERIA/OneDrive - UCB/PDAC/Bases de Datos/TCGA/TCGA_R/res_df_TCGA_V3.cvs", sep = "\t", row.names = FALSE, quote = FALSE)


### Heatmap ##############

# Defino mis genes de interés

genes_of_interest <- c("CD3D", "GRAP", "ITGAL", "LTB", "LYVE1", "NCF1",
                       "PLA2G2A", "PLA2G2D", "TNFRSF9", "GBP6", "IL20RB",
                       "MME", "WNT11")

# Cálculo de transformación de expresión (rlog)

rld <- rlogTransformation(dds, blind = FALSE) # a partir del dds
rld_mat <- assay(rld)  # filas: genes, columnas: muestras

### Símbolos HGNC #####
connect_to_ensembl <- function(mirrors) {
  for (mirror in mirrors) {
    message(paste("Intentando conectar al espejo:", mirror))
    ensembl <- tryCatch(
      useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = mirror),
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


# Conectar a Ensembl
mirrors <- c("www", "useast", "asia")
ensembl <- connect_to_ensembl(mirrors)

# Extraer los IDs Ensembl de las rownames
ensembl_ids <- rownames(rld_mat)

# Obtener información de los símbolos HGNC
gene_info <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = ensembl
)

# Convertir rownames en una columna para hacer merge()
rld_mat$ensembl_gene_id <- rownames(rld_mat)

# Unir la información obtenida con el dataframe de DESeq2
rld_mat <- merge(rld_mat, gene_info, by = "ensembl_gene_id", all.x = TRUE)

# Reemplazar los rownames con los símbolos HGNC, manteniendo los genes sin símbolo
new_rownames <- ifelse(is.na(rld_mat$hgnc_symbol) | rld_mat$hgnc_symbol == "", 
                       rld_mat$ensembl_gene_id, 
                       rld_mat$hgnc_symbol)

# Asegurar que los rownames sean únicos
rownames(rld_mat) <- make.unique(new_rownames)

# Eliminar la columna ensembl_gene_id (ya está en rownames)
rld_mat$gene <- NULL

# Verificar los resultados
head(rld_mat)

#Eliminar genes NA o sin mapear
#Verificar si hay genes sin mapear
unmapped_genes <- rld_mat$hgnc_symbol[is.na(rld_mat$hgnc_symbol)]
if (length(unmapped_genes) > 0) {
  cat("Genes no mapeados encontrados:\n")
  print(unmapped_genes)
} else {
  cat("Todos los genes fueron mapeados correctamente.\n")
}

#Elimina los genes no mapeados o con NA
res_df_mapped <- rld_mat %>% filter(!is.na(hgnc_symbol) & hgnc_symbol != "")





















