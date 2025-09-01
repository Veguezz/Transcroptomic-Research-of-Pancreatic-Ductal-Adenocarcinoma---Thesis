# Cargar paquetes
library(readxl)
library(dplyr)

# Leer el archivo de Excel (ajustar la ruta del archivo)
archivo <- "C:/Users/VALERIA/OneDrive - UCB/PDAC/Bases de Datos/PanCuRx/PanCuRx_Excel/PanCuRx_V1_Variables.xlsx"  # Cambia esto por la ubicación real del archivo
datos <- read_excel(archivo, sheet = "Edad")

# Calcular la media y desviación estándar por grupo
resumen <- datos %>%
  group_by(Grupo) %>%
  summarise(
    Media_Edad = mean(Edad, na.rm = TRUE),
    Desviacion_Estandar = sd(Edad, na.rm = TRUE)
  )

# Mostrar los resultados
print(resumen)
