# Cargar paquetes
library(readxl)
library(dplyr)

# Leer el archivo de Excel (ajustar la ruta del archivo)
archivo <- "C:/Users/VALERIA/OneDrive - UCB/PDAC/Bases de Datos/TCGA/TCGA_Excel/TCGA_V3_Variables.xlsx"  # Cambia esto por la ubicación real del archivo
datos <- read_excel(archivo, sheet = "DispEdad")

# Calcular la media
media_edad <- mean(datos$Edad, na.rm = TRUE)

# Calcular la desviación estándar
desviacion_estandar_edad <- sd(datos$Edad, na.rm = TRUE)

# Mostrar los resultados
cat("Media de las edades:", media_edad, "\n")
cat("Desviación estándar de las edades:", desviacion_estandar_edad, "\n")

# Mostrar los resultados
print(resumen)
