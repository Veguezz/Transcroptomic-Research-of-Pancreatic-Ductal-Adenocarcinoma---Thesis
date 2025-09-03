library(EPIC)

# EPIC ####

MatrizP <- readRDS("C:/Users/VALERIA/OneDrive - UCB/PDAC/Cibersort/norm_dds_counts_PanCuRx_V4.rds")

# Asegúrate de que esté en formato matrix
MatrizP <- as.matrix(MatrizP)

# Corre Epic

results_epic <- EPIC(bulk = MatrizP)

results_epic_df <- as.data.frame(results_epic)

# Adaptamos EPIC

results_epic_df <- results_epic_df %>%
  transmute(
    Linfocitos_B = `cellFractions.Bcells`,
    Linfocitos_T   = `cellFractions.CD8_Tcells` + `cellFractions.CD4_Tcells`, 
    Células_NK  = `cellFractions.NKcells`,
    Macrófagos = `cellFractions.Macrophages`,
    Otras_Células  = `cellFractions.otherCells`,
    Código_de_convergencia_numérica = `fit.gof.convergeCode`,
    RMSE = `fit.gof.RMSE_weighted`,
    RMSE_Genes = `fit.gof.Root_mean_squared_geneExpr_weighted`,
    Correlación_Spearman  = `fit.gof.spearmanR`,
    ValorP_Spearman = `fit.gof.spearmanP`,
    Correlación_Pearson = `fit.gof.pearsonR`,
    ValorP_Pearson = `fit.gof.pearsonP`,
    Pendiente_regresión = `fit.gof.regline_a_x`,
    Intercepto_regresión = `fit.gof.regline_b`,
    Suma_Total_Fracciones = `fit.gof.sum_mRNAProportions`,
    
  )

# CATEGORÍAS ####

# Nombre de muestras originales
ruta <- "C:/Users/VALERIA/OneDrive - UCB/PDAC/Bases de Datos/PanCuRx/PanCuRx_Excel/PanCuRx_V4.xlsx"

sheets <- excel_sheets(ruta)
print(sheets)

NombresP <- read_excel(ruta, sheet = "Pacientes")

Nombres <- NombresP$Nombres

# Remplazo los nombbres en las filas de CiberP
# Verificamos que el número de muestras coincida
stopifnot(length(Nombres) == nrow(results_epic_df))

# Reasignamos los nombres de las filas (muestras)
rownames(results_epic_df) <- Nombres

#Crear dataframe con clasificación Respondedor y No Respondedor
EPIC_Catego <- as.data.frame(results_epic_df)

# Vector con las categorías
categorias <- c(
  rep("No Respondedor", 28),
  rep("Respondedor", 20)
)

# Paso 2: Verifica que el número coincide
length(categorias) == nrow(EPIC_Catego)  # Esto debe dar TRUE

# Paso 3: Agrega la columna al data.frame
EPIC_Catego$Categoria <- categorias


# DIAGRAMA DE BARRAS ####
# Convertir a formato largo
# 1) Construye tu data.frame asegurándote de conservar los nombres
df <- as.data.frame(EPIC_Catego
                    , check.names = FALSE)
df$Sample <- rownames(df)

# 2) Pivotar únicamente las primeras 22 columnas (tipos celulares)
long <- df %>%
  pivot_longer(
    cols = 1:5,                 # aquí seleccionas las columnas 1 a 22
    names_to = "CellType",
    values_to = "Fraction"
  )

# Graficamos
ggplot(long, aes(x = Sample, y = Fraction, fill = CellType)) +
  geom_col() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Muestra",
    y = "Fracción estimada",
    fill = "Tipo Celular",
    title = "Composición inmunitaria por muestra (CIBERSORT)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )


# BOXPLOTS ####
# Pivotar las primeras 22 columnas (las fracciones celulares)
long <- Cibersort_Catego %>%
  pivot_longer(
    cols = 1:5,  # Ajusta si tus columnas de tipos celulares están en otro rango
    names_to = "CellType",
    values_to = "Fraction"
  )

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

# Convertir tu Cibersort con categoría a formato largo
long <- df %>%
  mutate(Sample = rownames(df)) %>%
  pivot_longer(
    cols = 1:5,          # columnas de fracción celular
    names_to = "CellType",
    values_to = "Fraction"
  )

# Primero calcula la altura máxima de cada CellType para que ahí vayan los P-values
pval_heights <- long %>%
  group_by(CellType) %>%
  summarise(max_frac = max(Fraction, na.rm = TRUE)) %>%
  mutate(label_y = max_frac + 0.05)  # Agregamos un pequeño espacio extra arriba


# Graficamos
ggplot(long, aes(x = CellType, y = Fraction, fill = Categoria)) +
  geom_violin(trim = FALSE, scale = "width", alpha = 0.6) +
  geom_boxplot(width = 0.1, position = position_dodge(0.9), outlier.shape = NA, alpha = 0.8) +
  stat_compare_means(
    aes(group = Categoria),
    method = "wilcox.test",
    label = "p.format",
    size = 3,
    hide.ns = TRUE,
    label.y = pval_heights$label_y  # <- Ahora usa alturas dinámicas
  ) +
  labs(
    title = "Comparación de infiltración inmune por tipo celular",
    x = "Tipo celular",
    y = "Fracción estimada",
    fill = "Grupo"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.position = "top"
  )

# HEATMAP & DENDOGRAMA ####

# Sacar solo las columnas de fracciones (las primeras 22 columnas)
matriz_heatmap <- as.matrix(EPIC_Catego[, 1:5])

# Asegúrate que las filas son las muestras
rownames(matriz_heatmap) <- rownames(EPIC_Catego)

# Creamos un data frame para las anotaciones de las muestras
anotacion_muestras <- data.frame(
  Grupo = EPIC_Catego$Categoria
)
rownames(anotacion_muestras) <- rownames(EPIC_Catego)

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
  main = "Infiltración inmune en muestras (EPIC)",
  color = colorRampPalette(c("green", "lightgrey", "red"))(100),
  legend_title = "Escala Z"
)


# WILCOXON & LIMMA ####

# 1. Fracciones de CIBERSORT
fractions <- EPIC_Catego[, 1:5]  # Solo columnas de tipos celulares (fracciones)

# 2. Vector de grupos
group <- factor(EPIC_Catego$Categoria)  # De tu columna "Categoria"

# 2½. Renormalizar para que cada fila sume 1
fractions <- fractions / rowSums(fractions)

# 3. Comprobación rápida
# (opcional en tu caso porque tus fracciones NO suman a 1 exactamente, pero igual checamos)
stopifnot(all(abs(rowSums(fractions) - 1) < 1e-6))

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
    title = "Comparación de Infiltración Inmune (Wilcoxon)",
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
    title = "Comparación de Infiltración Inmune (Limma-voom)",
    x = "Log Fold Change (Respondedor vs No Respondedor)",
    y = "-log10(FDR)"
  )



# HEATMAP DE CORRELACIÓN ####


# Calcular matriz de correlación
matriz_cor <- cor(EPIC_Catego[, 1:5], method = "spearman")

# Dibujar heatmap de correlación
corrplot(matriz_cor, method = "color", type = "upper", 
         col = colorRampPalette(c("blue", "white", "red"))(200),
         addCoef.col = "black", # valores de correlación
         tl.col = "black", tl.srt = 45, # texto de las etiquetas
         diag = TRUE)






