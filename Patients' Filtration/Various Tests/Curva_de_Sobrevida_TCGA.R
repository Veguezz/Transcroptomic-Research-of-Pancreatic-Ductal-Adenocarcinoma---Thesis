#Instalar paquetes
install.packages("survival")   # Análisis de supervivencia
install.packages("survminer")  # Visualización de curvas

#Cargar librerías

library(readxl)
library(survival)
library(survminer)
library(extrafont)

#Lectura de datos
file_path <- "C:/Users/VALERIA/OneDrive - UCB/PDAC/Bases de Datos/TCGA/TCGA_Excel/TCGA_V3_Variables.xlsx"

sheets <- excel_sheets(file_path)

print(sheets)

# Carga la hoja específica y guarda los datos en una variable
datos <- read_excel(file_path, sheet = "Curva de Sobrevida")

print(datos)

# Clasificar en grupos de edad
datos$group <- ifelse(datos$age >= 40 & datos$age <= 60, "40-60 años",
                      ifelse(datos$age >= 70 & datos$age <= 80, "70-80 años", NA))

# Filtrar solo los rangos de edad de interés
datos <- subset(datos, group %in% c("40-60 años", "70-80 años"))

# Crear objeto de supervivencia
surv_object <- Surv(time = datos$time, event = datos$status)

# Ajustar el modelo Kaplan-Meier por grupo de edad
fit_group <- survfit(surv_object ~ group, data = datos)

# Graficar las curvas
library(survminer)
ggsurvplot(
  fit_group,
  data = datos,
  conf.int = TRUE,               # Mostrar intervalos de confianza
  pval = TRUE,                   # Valor p para comparar grupos
  legend.title = "Rango de edad", # Título de la leyenda
  legend.labs = c("40-60 años", "70-80 años"), # Etiquetas de la leyenda
  xlab = "Tiempo (días)",        # Cambiar "Time" por "Tiempo (días)"
  ylab = "Probabilidad de sobrevida", # Cambiar "Survival Probability" por "Probabilidad de sobrevida"
  title = "Curva de Sobrevida (pacientes de 40-60 años y 70-80 años, TCGA)", # Definir título correctamente
  risk.table = TRUE,             # Mostrar tabla de riesgo
  risk.table.col = "strata",      # Colorear la tabla de riesgo por grupo
  risk.table.title = "Número en riesgo", # Mantener el título de la tabla de riesgo
  risk.table.fontsize = 4,        # Ajustar tamaño de la fuente en la tabla de riesgo
  ggtheme = theme_minimal() + 
    theme(
      text = element_text(family = "Times New Roman", size = 12), # Aplicar Times New Roman
      plot.title = element_text(hjust = 0.5, face = "bold", family = "Times New Roman"), # Título centrado en Times New Roman
      axis.title = element_text(size = 12, family = "Times New Roman"), 
      axis.text = element_text(size = 12, family = "Times New Roman"), 
      legend.text = element_text(size = 12, family = "Times New Roman"), 
      legend.title = element_text(size = 12, face = "bold", family = "Times New Roman") 
    )
)


