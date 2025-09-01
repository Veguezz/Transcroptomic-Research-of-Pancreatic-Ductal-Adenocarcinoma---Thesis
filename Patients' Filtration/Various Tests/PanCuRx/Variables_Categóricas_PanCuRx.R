#Para variables categóricas (como sexo, radioterapia, etc.), 
#usa la prueba exacta de Fisher.
#Esta es la mejor opción para  categorías con pocos pacientes (por ejemplo, solo 5 respondedores)
#Es más precisa que la prueba de chi-cuadrado cuando hay frecuencias bajas en las celdas de la tabla.

library(readxl)

#Lectura de datos
file_path <- "C:/Users/VALERIA/OneDrive - UCB/PDAC/Bases de Datos/PanCuRx/PanCuRx_Excel/PanCuRx_V4_Variables.xlsx"

sheets <- excel_sheets(file_path)

print(sheets)

#############################SEXO
# Carga la hoja específica y guarda los datos en una variable
Sexo <- read_excel(file_path, sheet = "Sexo")

# Verifica los datos cargados
print(Sexo)

# Convierte la tabla en una matriz de contingencia si es necesario (por ejemplo, si vas a hacer la prueba de Fisher)
tabla_contingencia_sexo <- as.matrix(Sexo[, -1]) 

tabla_contingencia_sexo <- apply(tabla_contingencia_sexo, c(1, 2), as.numeric)

# Verifica que todos los valores sean numéricos y no haya valores faltantes
if (any(is.na(tabla_contingencia_sexo))) {
  cat("Advertencia: Se encontraron valores NA en la matriz.\n")
  # Puedes reemplazar los NA por 0 si es adecuado para tu análisis (opcional)
  tabla_contingencia_sexo[is.na(tabla_contingencia_sexo)] <- 0
}
print(tabla_contingencia_sexo)

# Realiza la prueba exacta de Fisher
resultado_fisher_sexo <- fisher.test(tabla_contingencia_sexo)

# Imprime los resultados de la prueba
print(resultado_fisher_sexo)


#######################SUBDIVISIÓN ANATÓMICA
# Carga la hoja específica y guarda los datos en una variable
SubDiv <- read_excel(file_path, sheet = "SubDiv")

# Verifica los datos cargados
print(SubDiv)

# Convierte la tabla en una matriz de contingencia si es necesario (por ejemplo, si vas a hacer la prueba de Fisher)
tabla_contingencia_SubDiv <- as.matrix(SubDiv[, -1]) 

tabla_contingencia_SubDiv <- apply(tabla_contingencia_SubDiv, c(1, 2), as.numeric)

# Verifica que todos los valores sean numéricos y no haya valores faltantes
if (any(is.na(tabla_contingencia_SubDiv))) {
  cat("Advertencia: Se encontraron valores NA en la matriz.\n")
  # Puedes reemplazar los NA por 0 si es adecuado para tu análisis (opcional)
  tabla_contingencia_SubDiv[is.na(tabla_contingencia_SubDiv)] <- 0
}
print(tabla_contingencia_SubDiv)

# Realiza la prueba exacta de Fisher
resultado_fisher_SubDiv <- fisher.test(tabla_contingencia_SubDiv)

# Imprime los resultados de la prueba
print(resultado_fisher_SubDiv)


##########################ESTADO DEL PACIENTE
# Carga la hoja específica y guarda los datos en una variable
Estadio <- read_excel(file_path, sheet = "Estadio") 

# Verifica los datos cargados
print(Estadio)

# Convierte la tabla en una matriz de contingencia si es necesario (por ejemplo, si vas a hacer la prueba de Fisher)
TC_Estadio <- as.matrix(Estadio[, -1]) 

print(TC_Estadio)

TC_Estadio <- apply(TC_Estadio, c(1, 2), as.numeric)

# Verifica que todos los valores sean numéricos y no haya valores faltantes
if (any(is.na(TC_Estadio))) {
  cat("Advertencia: Se encontraron valores NA en la matriz.\n")
  # Puedes reemplazar los NA por 0 si es adecuado para tu análisis (opcional)
  TC_Estadio[is.na(TC_Estadio)] <- 0
}
print(TC_Estadio)

# Realiza la prueba exacta de Fisher
Res_Estadio <- fisher.test(TC_Estadio)

# Imprime los resultados de la prueba
print(Res_Estadio)


##########################RADIOTERAPIA
# Carga la hoja específica y guarda los datos en una variable
Radio <- read_excel(file_path, sheet = "Radio") 

# Verifica los datos cargados
print(Radio)

# Convierte la tabla en una matriz de contingencia si es necesario (por ejemplo, si vas a hacer la prueba de Fisher)
TC_Radio <- as.matrix(Radio[, -1]) 

print(TC_Radio)

TC_Radio <- apply(TC_Radio, c(1, 2), as.numeric)

# Verifica que todos los valores sean numéricos y no haya valores faltantes
if (any(is.na(TC_Radio))) {
  cat("Advertencia: Se encontraron valores NA en la matriz.\n")
  # Puedes reemplazar los NA por 0 si es adecuado para tu análisis (opcional)
  TC_Radio[is.na(TC_Radio)] <- 0
}
print(TC_Radio)

# Realiza la prueba exacta de Fisher
Res_Radio <- fisher.test(TC_Radio)

# Imprime los resultados de la prueba
print(Res_Radio)


#RELATIVO CON CANCER#########################
# Carga la hoja específica y guarda los datos en una variable
Rel_Cancer <- read_excel(file_path, sheet = "Rel_Cancer") 

# Verifica los datos cargados
print(Rel_Cancer)

# Convierte la tabla en una matriz de contingencia si es necesario (por ejemplo, si vas a hacer la prueba de Fisher)
TC_Rel_Cancer <- as.matrix(Rel_Cancer[, -1]) 

print(TC_Rel_Cancer)

TC_Rel_Cancer <- apply(Rel_Cancer, c(1, 2), as.numeric)

# Verifica que todos los valores sean numéricos y no haya valores faltantes
if (any(is.na(TC_Rel_Cancer))) {
  cat("Advertencia: Se encontraron valores NA en la matriz.\n")
  # Puedes reemplazar los NA por 0 si es adecuado para tu análisis (opcional)
  TC_Rel_Cancer[is.na(TC_Rel_Cancer)] <- 0
}
print(TC_Rel_Cancer)

# Realiza la prueba exacta de Fisher
Res_Rel_Cancer <- fisher.test(TC_Rel_Cancer)

# Imprime los resultados de la prueba
print(Res_Rel_Cancer)



#ESTADO VITAL#########################
# Carga la hoja específica y guarda los datos en una variable
EstadoV <- read_excel(file_path, sheet = "EstadoV") 

# Verifica los datos cargados
print(EstadoV)

# Convierte la tabla en una matriz de contingencia si es necesario (por ejemplo, si vas a hacer la prueba de Fisher)
TC_EstadoV <- as.matrix(EstadoV[, -1]) 

print(TC_EstadoV)

TC_EstadoV <- apply(EstadoV, c(1, 2), as.numeric)

# Verifica que todos los valores sean numéricos y no haya valores faltantes
if (any(is.na(TC_EstadoV))) {
  cat("Advertencia: Se encontraron valores NA en la matriz.\n")
  # Puedes reemplazar los NA por 0 si es adecuado para tu análisis (opcional)
  TC_EstadoV[is.na(TC_EstadoV)] <- 0
}
print(TC_EstadoV)

# Realiza la prueba exacta de Fisher
Res_EstadoV <- fisher.test(TC_EstadoV)

# Imprime los resultados de la prueba
print(Res_EstadoV)



#RESPONSE#########################
# Carga la hoja específica y guarda los datos en una variable
Response <- read_excel(file_path, sheet = "Response") 

# Verifica los datos cargados
print(Response)

# Convierte la tabla en una matriz de contingencia si es necesario (por ejemplo, si vas a hacer la prueba de Fisher)
TC_Response <- as.matrix(Response[, -1]) 

print(TC_Response)

TC_Response <- apply(Response, c(1, 2), as.numeric)

# Verifica que todos los valores sean numéricos y no haya valores faltantes
if (any(is.na(TC_Response))) {
  cat("Advertencia: Se encontraron valores NA en la matriz.\n")
  # Puedes reemplazar los NA por 0 si es adecuado para tu análisis (opcional)
  TC_Response[is.na(TC_Response)] <- 0
}
print(TC_Response)

# Realiza la prueba exacta de Fisher
Res_Response <- fisher.test(TC_Response)

# Imprime los resultados de la prueba
print(Res_Response)

