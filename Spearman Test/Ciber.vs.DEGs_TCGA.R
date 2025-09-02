library(Hmisc)
library(dplyr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(grid)


# PREAPAR DATOS ####

# Matríz de Expresión Diferencial
ddsT <- readRDS("C:/Users/VALERIA/OneDrive - UCB/PDAC/Cibersort/norm_dds_counts_TCGA_V3.rds")

ddsT <- as.data.frame(ddsT)

# Nombre de muestras originales
ruta <- "C:/Users/VALERIA/OneDrive - UCB/PDAC/Bases de Datos/TCGA/TCGA_Excel/TCGA_V3.xlsx"

sheets <- excel_sheets(ruta)
print(sheets)

NombresT <- read_excel(ruta, sheet = "Pacientes")

Nombres <- NombresT$Nombres

# Remplazo los nombbres en las filas de CiberP
# Verificamos que el número de muestras coincida
stopifnot(length(Nombres) == ncol(ddsT))

# Reasignamos los nombres de las filas (muestras)
colnames(ddsT) <- Nombres

# Matríz de Infiltración Inmunitaria
CiberT <- readRDS("C:/Users/VALERIA/OneDrive - UCB/PDAC/Cibersort/FirmaMod_TPM_TCGA_V3.rds")


genes_interest <- c(
  "CD3D", "GRAP", "ITGAL", "LTB", "LYVE1", "NCF1", 
  "PLA2G2A", "PLA2G2D", "TNFRSF9", "GBP6", 
  "IL20RB", "MME", "WNT11")

# Subset de expresión sólo para esos genes
DEGsT <- ddsT[genes_interest, , drop = FALSE]

# Alineamos muestras
common_samps <- intersect(colnames(DEGsT), rownames(CiberT))
# Expr: filas = genes, columnas = muestras
expr_sub <- DEGsT[, common_samps, drop = FALSE]  

# CIBERSORT: filas = muestras, columnas = tipos celulares
cib_sub  <- CiberT[common_samps, , drop = FALSE]

# SPEARMAN ####

# Convertimos en Matriz
if(is.data.frame(expr_sub)) {
  expr_sub <- as.matrix(expr_sub)
}

if(is.data.frame(cib_sub)) {
  cib_sub <- as.matrix(cib_sub)
}
#Corremos Spearman

res <- rcorr(
  x = t(expr_sub),  # transpose a M observaciones × 13 genes
  y = cib_sub,      # M observaciones × 7 células
  type = "spearman"
)


# VISUALIZACIÓN ####

# Extraer las matrices de coeficientes (ρ) y p-values
r_mat <- res$r        # matriz completa de correlaciones
p_mat <- res$P        # matriz completa de p-values

# Definir filas y columnas de interés
genes <- rownames(expr_sub)     # los 13 genes
cells <- colnames(cib_sub)      # las 7 poblaciones celulares

# Subset a genes × células
r_gc <- r_mat[genes, cells]     # 13 × 7 de ρ
p_gc <- p_mat[genes, cells]     # 13 × 7 de p-values

# Calcula –log10(p-value) para colorear
logp_gc <- -log10(p_gc)

# Define dos “paletas” de color
#   – Para ρ: azul (–1) → blanco (0) → rojo (+1)
col_fun_rho <- colorRamp2(c(-1, 0, 1), c("navyblue", "white", "red3"))
#   – Para –log10(p): desde p = 0.05 (–log10≈1.3) hasta el máximo que tengas
max_logp <- ceiling(max(logp_gc, na.rm = TRUE))
col_fun_p   <- colorRamp2(c(0, max_logp/2, max_logp),
                          c("white", "powderblue", "steelblue3"))

# Dibuja el heatmap dividido
ht <- Heatmap(
  matrix = r_gc,
  name   = "rho",                 # nombre para la leyenda de rho
  col    = col_fun_rho,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  cell_fun = function(j, i, x, y, w, h, fill) {
    # j = índice columna, i = índice fila
    # x,y = centro de la celda; w,h = ancho/alto de la celda
    
    # triángulo superior (ρ)
    grid.polygon(
      x = c(x - w/2, x + w/2, x - w/2),
      y = c(y + h/2, y + h/2, y - h/2),
      gp = gpar(fill = col_fun_rho(r_gc[i, j]), col = NA)
    )
    
    # triángulo inferior (–log10 p-value)
    grid.polygon(
      x = c(x - w/2, x + w/2, x + w/2),
      y = c(y - h/2, y + h/2, y - h/2),
      gp = gpar(fill = col_fun_p(logp_gc[i, j]), col = NA)
    )
    
    # texto de ρ en la esquina superior izquierda
    grid.text(
      label = sprintf("%.2f", r_gc[i, j]),
      x     = x - w * 0.15,
      y     = y + h * 0.15,
      gp    = gpar(fontsize = 8)
    )
    
    # texto de –log10(p) en la esquina inferior derecha
    grid.text(
      label = sprintf("%.1f", logp_gc[i, j]),
      x     = x + w * 0.15,
      y     = y - h * 0.15,
      gp    = gpar(fontsize = 8)
    )
  }
)

# Añade también una segunda leyenda para el –log10(p)
lgd_p <- Legend(
  at    = c(0, max_logp/2, max_logp),
  col_fun = col_fun_p,
  title = "-log10(p-value)"
)

# Dibuja todo junto
draw(ht, annotation_legend_list = list(lgd_p), heatmap_legend_side = "right")
