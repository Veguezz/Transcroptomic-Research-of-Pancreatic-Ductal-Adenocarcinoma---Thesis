# Cargar librerías
library(tidyverse)
library(readr)
library(pheatmap)
library(ComplexHeatmap)

#Cargamos datos
res_dfT <- readRDS("C:/Users/VALERIA/OneDrive - UCB/PDAC/Bases de Datos/TCGA/TCGA_R/res_df_TCGA_V3.2.rds") 

# Vector con los 13 genes de interés
genes_significativos <- c("CD3D", "GRAP", "ITGAL", "LTB", "LYVE1", "NCF1", 
                          "PLA2G2A", "PLA2G2D", "TNFRSF9", "GBP6", "IL20RB", 
                          "MME", "WNT11")

# Filtrar la matriz de datos normalizados para los genes significativos
res_dfT_heatmap <- res_dfT[rownames(res_dfT) %in% genes_significativos, ]

# Normalizar los datos con escala Z si es necesario
conteos_heatmap_scaled <- t(scale(t(res_dfT_heatmap)))

# Crear un vector que identifique los grupos de las muestras
grupo_colores <- c(rep("Respondedores", 8), rep("No Respondedores", 15))  # Ajusta estos números según tus datos

# Asignar colores a los grupos
colores_grupo <- setNames(c("lightgreen", "khaki"), c("No Respondedores", "Respondedores"))


# Crear una anotación para los grupos con una cuadrícula
anotacion_colores <- HeatmapAnnotation(
  Grupo = grupo_colores,
  col = list(Grupo = colores_grupo),
  show_legend = TRUE,  # Mostrar la leyenda de la anotación
  annotation_legend_param = list(title = NULL),  # Ocultar el título de la leyenda
  show_annotation_name = FALSE,  # Ocultar el nombre de la fila "Grupo"
  gp = gpar(col = "black", lwd = 0.5)  # Agregar la cuadrícula a la anotación
)

# Crear un esquema de colores más suaves para el heatmap
color_palette <- colorRampPalette(c("navyblue", "white", "red3"))(100)

# Crear el heatmap con ComplexHeatmap
Heatmap50<- Heatmap(
  conteos_heatmap_scaled,
  name = "Escala Z",  # Título de la escala de colores
  col = color_palette,  # Aplicar la paleta de colores suaves
  top_annotation = anotacion_colores,
  show_row_names = TRUE,
  show_column_names = TRUE,
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  row_names_side = "right",  # Mostrar los nombres de los genes en la derecha
  row_names_gp = gpar(fontsize = 8),
  column_names_side = "bottom",
  column_names_gp = gpar(fontsize = 10),
  rect_gp = gpar(col = "black", lwd = 0.5),  # Agregar cuadriculado para divisiones de filas y columnas
  heatmap_legend_param = list(title = "Escala Z")  # Agregar título a la leyenda de la escala de colores
)

dev.new()
print(Heatmap50)
