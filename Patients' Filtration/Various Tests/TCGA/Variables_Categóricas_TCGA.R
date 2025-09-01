#Para variables categóricas (como sexo, radioterapia, etc.), 
#usa la prueba exacta de Fisher.
#Esta es la mejor opción para  categorías con pocos pacientes (por ejemplo, solo 5 respondedores)
#Es más precisa que la prueba de chi-cuadrado cuando hay frecuencias bajas en las celdas de la tabla.

library(readxl)

#Lectura de datos
file_path <- "C:/Users/VALERIA/OneDrive - UCB/PDAC/Bases de Datos/TCGA/TCGA_Excel/TCGA_V3_Variables.xlsx"

sheets <- excel_sheets(file_path)

print(sheets)

#TUMOR RESIDUAL############################
# Carga la hoja específica y guarda los datos en una variable
ResidualT <- read_excel(file_path, sheet = "ResidualT")

# Verifica los datos cargados
print(ResidualT)

# Convierte la tabla en una matriz de contingencia si es necesario (por ejemplo, si vas a hacer la prueba de Fisher)
TC_ResidualT <- as.matrix(ResidualT[, -1]) 

TC_ResidualT <- apply(TC_ResidualT, c(1, 2), as.numeric)

# Verifica que todos los valores sean numéricos y no haya valores faltantes
if (any(is.na(TC_ResidualT))) {
  cat("Advertencia: Se encontraron valores NA en la matriz.\n")
  # Puedes reemplazar los NA por 0 si es adecuado para tu análisis (opcional)
  TC_ResidualT[is.na(TC_ResidualT)] <- 0
}
print(TC_ResidualT)

# Realiza la prueba exacta de Fisher
res_ResidualT <- fisher.test(TC_ResidualT)

# Imprime los resultados de la prueba
print(res_ResidualT)










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
SubDiv <- read_excel(file_path, sheet = "Subdivisión")

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

##########################RADIOTERAPIA
# Carga la hoja específica y guarda los datos en una variable
RadT <- read_excel(file_path, sheet = "Radioterapia")

# Verifica los datos cargados
print(RadT)

# Convierte la tabla en una matriz de contingencia si es necesario (por ejemplo, si vas a hacer la prueba de Fisher)
TC_RadT <- as.matrix(RadT[, -1]) 

TC_RadT <- apply(TC_RadT, c(1, 2), as.numeric)

# Verifica que todos los valores sean numéricos y no haya valores faltantes
if (any(is.na(TC_RadT))) {
  cat("Advertencia: Se encontraron valores NA en la matriz.\n")
  # Puedes reemplazar los NA por 0 si es adecuado para tu análisis (opcional)
  TC_RadT[is.na(TC_RadT)] <- 0
}
print(TC_RadT)

# Realiza la prueba exacta de Fisher
Res_TC_RadT <- fisher.test(TC_RadT)

# Imprime los resultados de la prueba
print(Res_TC_RadT)

#######################FÁRMACO
# Carga la hoja específica y guarda los datos en una variable
Drug <- read_excel(file_path, sheet = "Fármaco")

# Verifica los datos cargados
print(Drug)

# Convierte la tabla en una matriz de contingencia si es necesario (por ejemplo, si vas a hacer la prueba de Fisher)
tabla_contingencia_Drug <- as.matrix(Drug[, -1]) 

tabla_contingencia_Drug <- apply(tabla_contingencia_Drug, c(1, 2), as.numeric)

# Verifica que todos los valores sean numéricos y no haya valores faltantes
if (any(is.na(tabla_contingencia_Drug))) {
  cat("Advertencia: Se encontraron valores NA en la matriz.\n")
  # Puedes reemplazar los NA por 0 si es adecuado para tu análisis (opcional)
  tabla_contingencia_Drug[is.na(tabla_contingencia_Drug)] <- 0
}
print(tabla_contingencia_Drug)

# Realiza la prueba exacta de Fisher
resultado_fisher_Drug <- fisher.test(tabla_contingencia_Drug)

# Imprime los resultados de la prueba
print(resultado_fisher_Drug)

##########################ESTADIO
# Carga la hoja específica y guarda los datos en una variable
TStage <- read_excel(file_path, sheet = "Estadio")

# Verifica los datos cargados
print(TStage)

# Convierte la tabla en una matriz de contingencia si es necesario (por ejemplo, si vas a hacer la prueba de Fisher)
tabla_contingencia_TStage <- as.matrix(TStage[, -1]) 

tabla_contingencia_TStage <- apply(tabla_contingencia_TStage, c(1, 2), as.numeric)

# Verifica que todos los valores sean numéricos y no haya valores faltantes
if (any(is.na(tabla_contingencia_TStage))) {
  cat("Advertencia: Se encontraron valores NA en la matriz.\n")
  # Puedes reemplazar los NA por 0 si es adecuado para tu análisis (opcional)
  tabla_contingencia_TStage[is.na(tabla_contingencia_TStage)] <- 0
}
print(tabla_contingencia_TStage)

# Realiza la prueba exacta de Fisher
resultado_fisher_TStage <- fisher.test(tabla_contingencia_TStage)

# Imprime los resultados de la prueba
print(resultado_fisher_TStage)

##########################ESTADO DEL PACIENTE
# Carga la hoja específica y guarda los datos en una variable
Estado <- read_excel(file_path, sheet = "Estado") 

# Verifica los datos cargados
print(Estado)

# Convierte la tabla en una matriz de contingencia si es necesario (por ejemplo, si vas a hacer la prueba de Fisher)
TC_Estado <- as.matrix(Estado[, -1]) 

print(TC_Estado)

TC_Estado <- apply(TC_Estado, c(1, 2), as.numeric)

# Verifica que todos los valores sean numéricos y no haya valores faltantes
if (any(is.na(TC_Estado))) {
  cat("Advertencia: Se encontraron valores NA en la matriz.\n")
  # Puedes reemplazar los NA por 0 si es adecuado para tu análisis (opcional)
  TC_Estado[is.na(TC_Estado)] <- 0
}
print(TC_Estado)

# Realiza la prueba exacta de Fisher
Res_Estado <- fisher.test(TC_Estado)

# Imprime los resultados de la prueba
print(Res_Estado)


#CIRUGÍA#########################
# Carga la hoja específica y guarda los datos en una variable
Cirugia <- read_excel(file_path, sheet = "Cirugía") 

# Verifica los datos cargados
print(Cirugia)

# Convierte la tabla en una matriz de contingencia si es necesario (por ejemplo, si vas a hacer la prueba de Fisher)
TC_Cirugia <- as.matrix(Cirugia[, -1]) 

print(TC_Cirugia)

TC_Cirugia <- apply(TC_Cirugia, c(1, 2), as.numeric)

# Verifica que todos los valores sean numéricos y no haya valores faltantes
if (any(is.na(TC_Cirugia))) {
  cat("Advertencia: Se encontraron valores NA en la matriz.\n")
  # Puedes reemplazar los NA por 0 si es adecuado para tu análisis (opcional)
  TC_Cirugia[is.na(TC_Cirugia)] <- 0
}
print(TC_Cirugia)

# Realiza la prueba exacta de Fisher
Res_Cirugia <- fisher.test(TC_Cirugia)

# Imprime los resultados de la prueba
print(Res_Cirugia)


#CIRUGÍA#########################
# Carga la hoja específica y guarda los datos en una variable
Pancreatitis <- read_excel(file_path, sheet = "Pancreatitis") 

# Verifica los datos cargados
print(Pancreatitis)

# Convierte la tabla en una matriz de contingencia si es necesario (por ejemplo, si vas a hacer la prueba de Fisher)
TC_Pancreatitis <- as.matrix(Pancreatitis[, -1]) 

print(TC_Pancreatitis)

TC_Pancreatitis <- apply(Pancreatitis, c(1, 2), as.numeric)

# Verifica que todos los valores sean numéricos y no haya valores faltantes
if (any(is.na(TC_Pancreatitis)) {
  cat("Advertencia: Se encontraron valores NA en la matriz.\n")
  # Puedes reemplazar los NA por 0 si es adecuado para tu análisis (opcional)
  TC_Pancreatitis[is.na(TC_Pancreatitis)] <- 0
}
print(TC_Pancreatitis)

# Realiza la prueba exacta de Fisher
Res_Pancreatitis <- fisher.test(TC_Pancreatitis)

# Imprime los resultados de la prueba
print(Res_Pancreatitis)
