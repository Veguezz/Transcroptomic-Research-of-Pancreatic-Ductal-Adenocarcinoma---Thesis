library(biomaRt)
library(readxl)

# Conéctate al servidor Ensembl (genes humanos)
mart <- useEnsembl("genes", dataset = "hsapiens_gene_ensembl")

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

counts <- read_excel(file_path, sheet = "Counts")

counts <- data.frame(counts)

rownames(counts) <- counts$Gene

# Data frame 'counts' = tu matriz de counts (genes x muestras)
gene_symbols <- rownames(counts)

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
gene_len <- gene_len[match(rownames(counts), names(gene_len))]

# Elimina genes sin longitud o con length = 0
keep <- !is.na(gene_len) & gene_len > 0
counts_ok <- counts[keep, ]
counts_ok$Gene <- NULL
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
  file      = "C:/Users/VALERIA/OneDrive - UCB/PDAC/Bases de Datos/PanCuRx/PanCuRx_R/TPM_PanCuRx_V4.txt",
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)
