library(readxl)
library(dplyr)
library(ggplot2)
library(writexl)

# CARGAR DATOS ####

# Cargar las hojas del archivo
vias_tcga <- read_excel("C:/Users/VALERIA/OneDrive - UCB/PDAC/C Pearson/Vías Disminuidas.xlsx", sheet = "VíasT")  # "VíasT" → Vías de TCGA
vias_pancurx <- read_excel("C:/Users/VALERIA/OneDrive - UCB/PDAC/C Pearson/Vías Disminuidas.xlsx", sheet = "VíasP")  # "VíasP" → Vías de PanCuRx

# Leer archivo de Leading Edge genes
genes_tcga <- read_excel("C:/Users/VALERIA/OneDrive - UCB/PDAC/C Pearson/Vías Disminuidas.xlsx", sheet = "LEGsT")   # Gene, Pathway (esta hoja la creaste tú)
genes_pancurx <- read_excel("C:/Users/VALERIA/OneDrive - UCB/PDAC/C Pearson/Vías Disminuidas.xlsx", sheet = "LEGsP") # Gene, Pathway

# Leer los datos de DEGs para TCGA y PanCuRx
degs_tcga <- readRDS("C:/Users/VALERIA/OneDrive - UCB/PDAC/Bases de Datos/TCGA/TCGA_R/res_df_TCGA_V3.rds")    # columnas: Gene, log2FoldChange
degs_pancurx <- readRDS("C:/Users/VALERIA/OneDrive - UCB/PDAC/Bases de Datos/PanCuRx/PanCuRx_R/res_df_PanCuRx_V1.rds") 




# MEDIA DE Log2FC POR VÍA ####

# TCGA
# Renombrar hgnc_symbol a Gene
degs_tcga <- degs_tcga %>%
  rename(Gene = hgnc_symbol)

# Unir genes Leading Edge con log2FoldChange
genes_tcga_joined <- genes_tcga %>%
  inner_join(degs_tcga, by = "Gene")

# Calcular media de log2FoldChange por vía
mean_log2fc_tcga <- genes_tcga_joined %>%
  group_by(Pathway) %>%
  summarise(mean_log2FC = mean(log2FoldChange, na.rm = TRUE))

# Unir con tabla de NES
final_tcga <- vias_tcga %>%
  inner_join(mean_log2fc_tcga, by = c("Vías" = "Pathway"))  # aquí 'Vías' es tu columna de nombre de la vía

# PANCURX
# Renombrar hgnc_symbol a Gene
degs_pancurx <- degs_pancurx %>%
  rename(Gene = gene)

# Unir genes Leading Edge con log2FoldChange
genes_pancurx_joined <- genes_pancurx %>%
  inner_join(degs_pancurx, by = "Gene")

# Calcular media de log2FoldChange por vía
mean_log2fc_pancurx <- genes_pancurx_joined %>%
  group_by(Pathway) %>%
  summarise(mean_log2FC = mean(log2FoldChange, na.rm = TRUE))

# Unir con tabla de NES
final_pancurx <- vias_pancurx %>%
  inner_join(mean_log2fc_pancurx, by = c("Vías" = "Pathway"))


#CORRELACIÓN DE PEARSON ####

# TCGA
correlation_tcga <- cor.test(final_tcga$mean_log2FC, final_tcga$NES, method = "pearson")
print(correlation_tcga)

# PANCURX
correlation_pancurx <- cor.test(final_pancurx$mean_log2FC, final_pancurx$NES, method = "pearson")
print(correlation_pancurx)


# VISUALIZACIÓN ####
# Gráfico TCGA
ggplot(final_tcga, aes(x = mean_log2FC, y = NES)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  theme_minimal() +
  labs(title = "Correlación NES vs Media log2FoldChange (TCGA)",
       x = "Media log2FoldChange",
       y = "Normalized Enrichment Score (NES)")

# Gráfico PanCuRx
ggplot(final_pancurx, aes(x = mean_log2FC, y = NES)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  theme_minimal() +
  labs(title = "Correlación NES vs Media log2FoldChange (PanCuRx)",
       x = "Media log2FoldChange",
       y = "Normalized Enrichment Score (NES)")

# y


