# Instalar los paquetes si no los tienes
# install.packages("readxl")
# install.packages("ggplot2")
# install.packages("tidyr")

# Cargar los paquetes
library(readxl)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggrepel)


# Cargar los datos desde el archivo de Excel
file_path <- "C:/Users/VALERIA/OneDrive - UCB/PDAC/Bases de Datos/TCGA/Edades TCGAV3.xlsx"
data <- read_excel(file_path)

# Convertir los datos a formato largo para ggplot
data_long <- pivot_longer(data, cols = everything(), names_to = "Dataset", values_to = "Age")

# Eliminar filas con valores NA
data_long <- na.omit(data_long)


# Identificar el valor máximo y mínimo para cada grupo
max_tcga <- max(data$TCGA, na.rm = TRUE)
min_tcga <- min(data$TCGA, na.rm = TRUE)

cat("El valor máximo en TCGA es:", max_tcga, "\n")
cat("El valor mínimo en TCGA es:", min_tcga, "\n")

# Seleccionar los 3 valores más altos y los 3 más bajos por grupo
extreme_values <- data_long %>%
  group_by(Dataset) %>%
  arrange(Dataset, Age) %>%
  slice(c(1:3, (n()-2):n())) %>%
  ungroup()

# Mostrar los resultados
print(extreme_values)

# Calcular la mediana de cada grupo y guardarla como vector
median_values <- data_long %>%
  group_by(Dataset) %>%
  summarize(Mediana = median(Age, na.rm = TRUE))


# Convertir la columna de medianas a un vector con nombres de cada grupo
mediana <- setNames(median_values$Mediana, median_values$Dataset)

# Mostrar el vector resultante
print(mediana)

# Calcular la media de cada grupo y guardarla como vector
mean_values <- data_long %>%
  group_by(Dataset) %>%
  summarize(Media = mean(Age, na.rm = TRUE))

# Convertir la columna de medias a un vector con nombres de cada grupo
media <- setNames(mean_values$Media, mean_values$Dataset)

# Mostrar el vector resultante
print(media)

# Calcular la desviación estándar de cada grupo y guardarla como vector
sd_values <- data_long %>%
  group_by(Dataset) %>%
  summarize(DesviacionEstandar = sd(Age, na.rm = TRUE))

# Convertir la columna de desviaciones estándar a un vector con nombres de cada grupo
D_estadar <- setNames(sd_values$DesviacionEstandar, sd_values$Dataset)
print(D_estadar)


# Graficar el dispersión con la línea de la mediana para cada grupo
ggplot(data_long, aes(x = Dataset, y = Age)) +
  geom_jitter(width = 0.05, size = 1.5) +  # Agrega dispersión a los puntos
  labs(title = "Distribución de Edades entre TCGA", x = "Dataset", y = "Edad") +
  ylim(0, 100) +  # Limitar el rango del eje y de 0 a 100
  theme_minimal() +
  geom_hline(data = median_values, aes(yintercept = Mediana, color = Dataset), 
             linetype = "dashed", size = 1) +  # Agregar la línea de mediana para cada grupo
  scale_color_manual(values = c("TCGA" = "blue", "PanCuRx" = "red")) +  # Colores para las medianas
  theme(legend.position = "none")  # Eliminar la leyenda si no es necesaria

# Seleccionar los 3 valores más altos y los 3 más bajos por grupo para etiquetarlos
extreme_values <- data_long %>%
  group_by(Dataset) %>%
  arrange(Age) %>%
  slice(c(1:3, (n()-2):n())) %>%
  ungroup()

# Dividir en valores bajos y altos para ambos grupos
extreme_values_low <- extreme_values %>%
  group_by(Dataset) %>%
  slice(1:3) %>% 
  ungroup()

extreme_values_high <- extreme_values %>%
  group_by(Dataset) %>%
  slice((n()-2):n()) %>% 
  ungroup()

# Graficar el dispersión con la línea de la mediana para cada grupo y etiquetas en los extremos
ggplot(data_long, aes(x = Dataset, y = Age)) +
  geom_jitter(width = 0.05, size = 1.5) +  # Agrega dispersión a los puntos
  labs(title = "Distribución de Edades entre TCGA y PanCuRx", x = "Dataset", y = "Edad") +
  ylim(0, 100) +  # Limitar el rango del eje y de 0 a 100
  theme_minimal() +
  geom_hline(data = median_values, aes(yintercept = Mediana, color = Dataset), 
             linetype = "dashed", size = 1) +  # Agregar la línea de mediana para cada grupo
  # Etiquetas para valores bajos, apuntando hacia abajo
  geom_text_repel(data = extreme_values_low, aes(label = round(Age, 1)), 
                  size = 3, color = "black", nudge_y = -10, 
                  direction = "x", segment.color = "black", segment.size = 0.5) +
  # Etiquetas para valores altos, apuntando hacia arriba
  geom_text_repel(data = extreme_values_high, aes(label = round(Age, 1)), 
                  size = 3, color = "black", nudge_y = 10, 
                  direction = "x", segment.color = "black", segment.size = 0.5) +
  scale_color_manual(values = c("TCGA" = "blue", "PanCuRx" = "red")) +  # Colores para las medianas
  theme(legend.position = "none")  # Eliminar la leyenda si no es necesaria



#############
# Agregar las medias y desviaciones estándar a data_long
data_long <- data_long %>%
  left_join(mean_values, by = "Dataset") %>%
  left_join(sd_values, by = "Dataset") %>%
  mutate(InRange = ifelse(Age >= (Media - DesviacionEstandar) & Age <= (Media + DesviacionEstandar), "Dentro", "Fuera"))

# Filtrar los valores que están dentro de la desviación estándar
ages_in_range <- data_long %>%
  filter(InRange == "Dentro") %>%
  select(Dataset, Age)

# Mostrar las edades dentro de la desviación estándar
print(ages_in_range)


# Graficar la dispersión con colores según la desviación estándar
ggplot(data_long, aes(x = Dataset, y = Age, color = InRange)) +
  geom_jitter(width = 0.05, size = 1.5) +  # Agregar dispersión a los puntos
  labs(title = "Distribución de Edades", x = "Dataset", y = "Edad") +
  ylim(0, 100) +  # Limitar el rango del eje y de 0 a 100
  theme_minimal() +
  geom_hline(data = median_values, aes(yintercept = Mediana), 
             linetype = "dashed", size = 1, color = "black") +  # Línea de mediana
  scale_color_manual(values = c("Dentro" = "green", "Fuera" = "gray")) +  # Colores
  theme(legend.position = "right")  # Mostrar leyenda


