#La prueba de U de Mann-Whitney es una prueba para comparar
#dos grupos en una variable continua es recomedables para datos con pocas muestras.
#No hacen suposiciones sobre la distribución de los datos.
#Comparan la mediana y la distribución general de los grupos.

library(readxl)

# Lectura de datos
file_path <- "C:/Users/VALERIA/OneDrive - UCB/PDAC/Bases de Datos/PanCuRx/PanCuRx_Excel/PanCuRx_V4_Variables.xlsx"

# Carga la hoja "Edad" y guarda los datos en una variable
Edad <- read_excel(file_path, sheet = "Edad")

# Verifica los datos cargados
print(colnames(Edad))

# Asegúrate de que las columnas estén correctamente nombradas, por ejemplo:
# Columna 'Edad' para la edad y 'Grupo' para el grupo de respondedores y no respondedores

# Asegurar que Grupo sea un factor
Edad$Grupo <- as.factor(Edad$Grupo)

# Crear vectores separados para cada grupo
respondedores <- Edad$Edad[Edad$Grupo == "Respondedor"]
no_respondedores <- Edad$Edad[Edad$Grupo == "No Respondedor"]

# Realiza la prueba de U de Mann-Whitney
resultado_mann_whitney <- wilcox.test(respondedores, no_respondedores)

# Imprime los resultados de la prueba
print(resultado_mann_whitney)
