library(readxl)
library(tidyr)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(ggpubr)
library(limma)
library(ggrepel)
library(tibble)
library(corrplot)

# DATOS ####
CiberT <- readRDS("C:/Users/VALERIA/OneDrive - UCB/PDAC/Cibersort/Res_Ciber_TPM_TCGA_V3.rds")

CiberT_df <- as.data.frame(CiberT)

CiberT_dfF <- CiberT_df %>%
  transmute(
    Linfocitos_B = `B cells naive` + `B cells memory` + `Plasma cells`,
    Linfocitos_T   = `T cells CD8` + `T cells CD4 memory resting` + `T cells CD4 memory activated`,
    Treg      = `T cells regulatory (Tregs)`,
    Células_NK  = `NK cells resting` + `NK cells activated`,
    Macrófagos_M2 = `Macrophages M2`,
    Macrófagos_M0M1 = `Macrophages M0` + `Macrophages M1`,
    Células_Dendtríticas  = `Dendritic cells resting` + `Dendritic cells activated`
  )

# AGREGAMOS CATEGORÍAS ####

# Nombre de muestras originales

ruta2 <- "C:/Users/VALERIA/OneDrive - UCB/PDAC/Bases de Datos/TCGA/TCGA_Excel/TCGA_V3.xlsx"

sheets2 <- excel_sheets(ruta2)
print(sheets2)

NombresT <- read_excel(ruta2, sheet = "Pacientes")

NombresT <- NombresT$Nombres

# Remplazo los nombbres en las filas de CiberP
# Verificamos que el número de muestras coincida
stopifnot(length(NombresT) == nrow(CiberT_dfF))

# Reasignamos los nombres de las filas (muestras)
rownames(CiberT_dfF) <- NombresT

#Crear dataframe con clasificación Respondedor y No Respondedor
CiberT_Catego <- as.data.frame(CiberT_dfF)

# Vector con las categorías
categorias <- c(
  rep("No Respondedor", 15),
  rep("Respondedor", 8)
)

# Paso 2: Verifica que el número coincide
length(categorias) == nrow(CiberT_Catego)  # Esto debe dar TRUE

# 2) Asignar directamente el vector como nueva columna
CiberT_Catego$Categoria <- categorias

# 3) Pasar a formato largo, manteniendo Categoria
long <- CiberT_Catego %>%
  rownames_to_column("Sample") %>%               # si necesitas nombre de muestra
  pivot_longer(
    cols      = -c(Sample, Categoria),
    names_to  = "CellType",
    values_to = "Fraction"
  )


Cibersort_Catego_Spearman <- CiberT_Catego
Cibersort_Catego_Spearman$Categoria <- NULL
saveRDS(Cibersort_Catego_Spearman, file = "C:/Users/VALERIA/OneDrive - UCB/PDAC/Cibersort/FirmaMod_TPM_TCGA_V3.rds")

# DIAGRAMA DE BARRAS ####
# Convertir a formato largo
# 1) Construye tu data.frame asegurándote de conservar los nombres
dfT <- as.data.frame(CiberT_Catego
                     , check.names = FALSE)
dfT$Sample <- rownames(dfT)

# Graficamos
ggplot(long, aes(x = Sample, y = Fraction, fill = CellType)) +
  geom_col() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Muestra",
    y = "Fracción estimada",
    fill = "Tipo Celular",
    title = "Composición inmunitaria por muestra (TCGA)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )


# BOXPLOTS ####

# Boxplots
ggplot(long, aes(x = Categoria, y = Fraction, fill = Categoria)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 0.8, alpha = 0.4) +
  facet_wrap(~ CellType, scales = "free_y") +
  labs(
    title = "Infiltración inmunitaria por grupo",
    x = "Grupo de pacientes",
    y = "Fracción estimada",
    fill = "Grupo"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 0),
    strip.text = element_text(size = 9)
  )



# WILCOX.TEST ####
# Comparación entre grupos por tipo celular (Wilcoxon)
test_resultados <- long %>%
  group_by(CellType) %>%
  summarise(
    p_value = wilcox.test(Fraction ~ Categoria)$p.value,
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))  # Ajuste por FDR

# Ver resultados ordenados por significancia
test_resultados %>% arrange(p_adj)

# Gráfico de barras del log10(p-valor ajustado)

test_resultados %>%
  mutate(significativo = ifelse(p_adj < 0.05, "Sí", "No")) %>%
  ggplot(aes(x = reorder(CellType, -log10(p_adj)), y = -log10(p_adj), fill = significativo)) +
  geom_col() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  coord_flip() +
  scale_fill_manual(values = c("Sí" = "tomato", "No" = "grey70")) +
  labs(
    title = "Significancia de diferencias por célula (Wilcoxon)",
    x = "Tipo Celular",
    y = expression(-log[10](p[adj])),
    fill = "Significativo"
  ) +
  theme_minimal()


# GRÁFICO DE VIOLIN ####

# Primero calcula la altura máxima de cada CellType para que ahí vayan los P-values
pval_heights <- long %>%
  group_by(CellType) %>%
  summarise(max_frac = max(Fraction, na.rm = TRUE)) %>%
  mutate(label_y = max_frac + 0.05)  # Agregamos un pequeño espacio extra arriba


# Graficamos
ggplot(long, aes(x = CellType, y = Fraction, fill = Categoria)) +
  geom_violin(trim = FALSE, scale = "width", alpha = 0.6) +
  geom_boxplot(width = 0.1,
               position = position_dodge(0.9),
               outlier.shape = NA,
               alpha = 0.8) +
  stat_compare_means(aes(group = Categoria),
                     method = "wilcox.test",
                     label = "p.format",
                     size = 3,
                     hide.ns = TRUE,
                     label.y = pval_heights$label_y) +
  scale_fill_manual(
    values = c(
      "No Respondedor" = "#CD0000",   # rojo
      "Respondedor"  = "#1874CD"    # verde
    )
  ) +
  labs(title = "Comparación de infiltración inmune por tipo celular (TCGA)",
       x     = "Tipo celular",
       y     = "Fracción estimada",
       fill  = "Grupo") +
  theme_bw() +
  theme(
    axis.text.x    = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.title     = element_text(hjust = 0.5, face = "bold"),
    legend.title   = element_text(face = "bold"),
    legend.position= "top"
  )


# HEATMAP & DENDOGRAMA ####

# Sacar solo las columnas de fracciones (las primeras 22 columnas)
matriz_heatmap <- as.matrix(CiberT_Catego[, 1:6])

# Asegúrate que las filas son las muestras
rownames(matriz_heatmap) <- rownames(CiberT_Catego)

# Creamos un data frame para las anotaciones de las muestras
anotacion_muestras <- data.frame(
  Grupo = CiberT_Catego$Categoria
)
rownames(anotacion_muestras) <- rownames(CiberT_Catego)

#Define los colores de los grupos
ann_colors <- list(
  Grupo = c("Respondedor" = "skyblue", "No Respondedor" = "salmon")
)

pheatmap(
  t(matriz_heatmap),             # Transponemos para que células sean filas
  scale = "row",
  name = "Escala Z",
  clustering_distance_cols = "euclidean",
  clustering_distance_rows = "euclidean",
  clustering_method = "complete",
  annotation_col = anotacion_muestras,
  annotation_colors = ann_colors,
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize_row = 7,
  fontsize_col = 5,
  main = "Infiltración inmune en muestras (CIBERSORT)",
  color = colorRampPalette(c("green", "lightgrey", "red"))(100),
  legend_title = "Escala Z"
)


# WILCOXON & LIMMA ####

# 1. Fracciones de CIBERSORT
fractions <- CiberT_Catego[, 1:7]  # Solo columnas de tipos celulares (fracciones)

# 2. Vector de grupos
group <- factor(CiberT_Catego$Categoria)  # De tu columna "Categoria"

# 3. Comprobación rápida

# --- Análisis Wilcoxon ---

wilcox_results <- apply(fractions, 2, function(cell){
  wilcox.test(cell ~ group)$p.value
})

wilcox_df <- data.frame(
  CellType = names(wilcox_results),
  pvalue   = wilcox_results,
  padj     = p.adjust(wilcox_results, method = "BH"),
  MedianDiff = apply(fractions, 2, function(x)
    median(x[group == "Respondedor"]) - median(x[group == "No Respondedor"]))
) %>% arrange(padj)

# --- Visualización rápida de Wilcoxon ---
ggplot(wilcox_df, aes(MedianDiff, -log10(padj), label = CellType)) +
  geom_point() +
  ggrepel::geom_text_repel(max.overlaps = 20) +
  theme_bw() +
  labs(
    title = "Comparación de Infiltración Inmune - Wilcoxon (TCGA)",
    x = "Diferencia de medianas (Respondedor - No Respondedor)",
    y = "-log10(p ajustado)"
  )

# --- Análisis Limma-voom ---

# 1. Transformación de fracciones
y <- t(asin(sqrt(as.matrix(fractions))))  # t() = transponer

# 2. Diseño del modelo
design <- model.matrix(~ group)

# 3. Aplicar voom
v <- voom(y, design, plot = TRUE)

# 4. Ajuste lineal
fit <- lmFit(v, design)
fit <- eBayes(fit)

# 5. Resultados limma
limma_df <- topTable(fit, coef = "groupRespondedor", number = Inf)

# --- Visualización rápida de limma ---
limma_df <- limma_df %>% rownames_to_column(var = "CellType")

ggplot(limma_df, aes(logFC, -log10(adj.P.Val), label = CellType)) +
  geom_point() +
  ggrepel::geom_text_repel(max.overlaps = 20) +
  theme_bw() +
  labs(
    title = "Comparación de Infiltración Inmune - Limma-voom (TCGA)",
    x = "Log Fold Change (Respondedor vs No Respondedor)",
    y = "-log10(FDR)"
  )



# HEATMAP & DENDOGRAMA ####


# Calcular matriz de correlación
matriz_cor <- cor(dfT[, 1:7], method = "spearman")

# Dibujar heatmap de correlación
corrplot(matriz_cor, method = "color", type = "upper", 
         col = colorRampPalette(c("blue", "white", "red"))(200),
         addCoef.col = "black", # valores de correlación
         tl.col = "black", tl.srt = 45, # texto de las etiquetas
         diag = TRUE)
