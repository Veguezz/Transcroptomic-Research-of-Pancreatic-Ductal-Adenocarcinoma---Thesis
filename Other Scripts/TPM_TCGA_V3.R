library(biomaRt)
library(readxl)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(annotables)
library(HGNChelper)
library(dplyr)

# MIS DATOS ####
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

counts <- read_excel(file_path, sheet = "Counts")

counts_df <- data.frame(counts)

# Convertimos ensembl a Hgnc symbol ####
colnames(counts_df)[1] <- "ensembl_gene_id"

# ELIMINAMOS SUFIJOS DE LOS GENES
counts_df$ensembl_gene_id <- sub("\\.\\d+$", "", counts_df$ensembl_gene_id)  # Remover sufijos numéricos
counts_df$ensembl_gene_id <- gsub("_.*", "", counts_df$ensembl_gene_id)  # Remover sufijos no estándar

#Verifica si hay Valores duplicados en la columna de Gene
sum(duplicated(counts_df$ensembl_gene_id))

# Estudiamos estos genes
duplicated_genes <- counts_df %>% 
  filter(ensembl_gene_id %in% ensembl_gene_id[duplicated(ensembl_gene_id)])

#Elimina los genes duplicados
counts_df <- counts_df[!duplicated(counts_df$ensembl_gene_id), ]

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
ensembl_ids <- counts$Gene

# Obtener información de los símbolos HGNC
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters    = "ensembl_gene_id",
  values     = counts_df$ensembl_gene_id,
  mart       = ensembl
)

# Unir recuentos + símbolos y filtrar no mapeados
counts_hgnc <- counts_df %>%
  left_join(gene_map, by = "ensembl_gene_id") %>%
  filter(!is.na(hgnc_symbol) & hgnc_symbol != "") 

#Verifica si hay Valores duplicados en la columna de hgnc
sum(duplicated(counts_hgnc$hgnc_symbol))

# Estudiamos estos genes
duplicated_genes <- counts_hgnc %>% 
  filter(hgnc_symbol %in% hgnc_symbol[duplicated(hgnc_symbol)])

# Sumamos los valores de genes duplicados
counts_hgnc <- counts_hgnc %>%
  group_by(hgnc_symbol) %>%
  summarise(across(-ensembl_gene_id, sum), .groups = "drop")

# colocamos hgnc_symbol como row names
counts_hgnc <- as.data.frame(counts_hgnc)
rownames(counts_hgnc) <- counts_hgnc$hgnc_symbol
counts_hgnc$hgnc_symbol <- NULL

# Data frame 'counts' = tu matriz de counts (genes x muestras)
gene_symbols <- rownames(counts_hgnc)

# OBTENEMOS LONGITUDES ####

# Conéctate al servidor Ensembl (genes humanos)
mart <- useEnsembl("genes", dataset = "hsapiens_gene_ensembl")

# Descarga la longitud colapsada de exones (gene_length)
len_tbl <- getBM(
  attributes = c("hgnc_symbol", "start_position", "end_position"),
  filters    = "hgnc_symbol",
  values     = gene_symbols,          # tus símbolos
  mart       = mart
)

len_tbl$genomic_length <- len_tbl$end_position - len_tbl$start_position + 1

# vector con longitudes y nombres = símbolo
gene_len <- len_tbl$genomic_length
names(gene_len) <- len_tbl$hgnc_symbol

# Alinea al orden de los genes de tu matriz
gene_len <- gene_len[match(rownames(counts_hgnc), names(gene_len))]

# Elimina genes sin longitud o con length = 0
keep <- !is.na(gene_len) & gene_len > 0
counts_ok <- counts_hgnc[keep, ]
len_ok    <- gene_len[keep]

########  (b) Reads-Per-Kilobase  ########
rpk <- counts_ok / (len_ok / 1e3)     # divide por kb

########  (c) TPM por muestra ########
tpm <- t( t(rpk) / colSums(rpk) ) * 1e6

########  (d) Verifica que cada columna sume ~1e6 ########
stopifnot(all(abs(colSums(tpm) - 1e6) < 1e-3))

# Añade una columna Gene al inicio (CIBERSORTx la necesita)
tpm_out <- cbind(Gene = rownames(tpm), as.data.frame(tpm))

# Guardamos
write.table(
  tpm_out,
  file      = "C:/Users/VALERIA/OneDrive - UCB/PDAC/Bases de Datos/TCGA/TCGA_R/TPM_TCGA_V3.txt",
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)
