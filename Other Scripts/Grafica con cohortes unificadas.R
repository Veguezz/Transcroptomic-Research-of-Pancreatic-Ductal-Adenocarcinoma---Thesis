library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)

# CARGAMOS DATOS ####
CiberP <- readRDS("C:/Users/VALERIA/OneDrive - UCB/PDAC/Cibersort/Res_Ciber_TPM_PanCuRx.rds")
CiberP_df <- as.data.frame(CiberP)
CiberP_dfF <- CiberP_df %>%
  transmute(
    Linfocitos_B = `B cells naive` + `B cells memory` + `Plasma cells`,
    Linfocitos_T   = `T cells CD8`+`T cells CD4 memory resting` + `T cells CD4 memory activated`, 
    Treg      = `T cells regulatory (Tregs)`,
    Células_NK  = `NK cells resting` + `NK cells activated`,
    Macrófagos_M2 = `Macrophages M2`,
    Macrófagos_M0M1 = `Macrophages M0` + `Macrophages M1`,
    Células_Dendtríticas  = `Dendritic cells resting` + `Dendritic cells activated`
  )

CiberT <- readRDS("C:/Users/VALERIA/OneDrive - UCB/PDAC/Cibersort/Res_Ciber_TPM_TCGA_V3.rds")
CiberT_df <- as.data.frame(CiberT)
CiberT_dfF <- CiberT_df %>%
  transmute(
    Linfocitos_B = `B cells naive` + `B cells memory` + `Plasma cells`,
    Linfocitos_T   = `T cells CD8` + `T cells CD4 memory resting` + `T cells CD4 memory activated`,
    Treg      = `T cells regulatory (Tregs)`,
    Células_NK  = `NK cells resting` + `NK cells activated`,
    Macrófagos_M2 = `Macrophages M2`,
    Macrófagos_M0M1 = `Macrophages M0` + `Macrophages M1`,
    Células_Dendtríticas  = `Dendritic cells resting` + `Dendritic cells activated`
  )

# MODIFICACIÓN DE MATRICES PARA SU UNIÓN ####
# Agregamos columna de cohorte
CiberT_dfF$Cohorte <- "TCGA"
CiberP_dfF$Cohorte <- "PanCuRx"

# Convertimos rownames en una columna
CiberT_dfF <- CiberT_dfF %>%
  rownames_to_column(var = "SampleID")
CiberP_dfF <- CiberP_dfF %>%
  rownames_to_column(var = "SampleID")

# Agrego la categoría de los grupos

#TCGA
# Vector con las categorías
categoriasT <- c(
  rep("No Respondedor", 15),
  rep("Respondedor", 8)
)

# Paso 2: Verifica que el número coincide
length(categoriasT) == nrow(CiberT_dfF)  # Esto debe dar TRUE

# 2) Asignar directamente el vector como nueva columna
CiberT_dfF$Grupo <- categoriasT

# PANCURX
# Vector con las categorías
categoriasP <- c(
  rep("No Respondedor", 28),
  rep("Respondedor", 20)
)

# Paso 2: Verifica que el número coincide
length(categoriasP) == nrow(CiberP_dfF)  # Esto debe dar TRUE

# 2) Asignar directamente el vector como nueva columna
CiberP_dfF$Grupo <- categoriasP

# UNIÓN ####
# Transformar de wide a long (Esto reorganiza tu tabla para que cada fila sea una combinación única de muestra, célula y fracción estimada)
tcga_long <- CiberT_dfF %>%
  pivot_longer(cols = -c(SampleID, Grupo, Cohorte),
               names_to = "TipoCelular",
               values_to = "Fraccion")

pancurx_long <- CiberP_dfF %>%
  pivot_longer(cols = -c(SampleID, Grupo, Cohorte),
               names_to = "TipoCelular",
               values_to = "Fraccion")

datos_combinados <- bind_rows(tcga_long, pancurx_long)

# GRÁFICA ####

# Crear variable de eje x (tipo celular + grupo de respuesta)
datos_combinados <- datos_combinados %>%
  mutate(EjeX = paste(TipoCelular, Grupo, sep = "_"),
         GrupoCohorte = paste(Grupo, Cohorte, sep = "_"))

# Orden lógico en eje x
datos_combinados$EjeX <- factor(datos_combinados$EjeX,
                                levels = c(
                                  as.vector(sapply(unique(datos_combinados$TipoCelular), function(tc) {
                                    c(paste(tc, "No Respondedor", sep = "_"), paste(tc, "Respondedor", sep = "_"))
                                  }))
                                )
)

# Orden en colores
datos_combinados$GrupoCohorte <- factor(datos_combinados$GrupoCohorte,
                                        levels = c(
                                          "No Respondedor_PanCuRx",
                                          "No Respondedor_TCGA",
                                          "Respondedor_PanCuRx",
                                          "Respondedor_TCGA"
                                        )
)

# Crear gráfico
ggplot(datos_combinados, aes(x = EjeX, y = Fraccion, fill = GrupoCohorte)) +
  geom_violin(alpha = 0.4, position = position_identity(), color = NA, trim = TRUE) +
  
  # Añadir punto de media
  stat_summary(fun = mean, geom = "point", shape = 21, size = 2,
               color = "black", fill = "white", position = position_identity()) +
  
  # Añadir barra de desviación estándar
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
               geom = "errorbar", width = 0.2, color = "black",
               position = position_identity()) +
  
  # Paleta personalizada (tonos azul y rojo)
  scale_fill_manual(values = c(
    "Respondedor_TCGA" = "#4A90E2",      # Azul claro
    "Respondedor_PanCuRx" = "#003f7d",   # Azul oscuro
    "No Respondedor_TCGA" = "#E94F4F",   # Rojo claro
    "No Respondedor_PanCuRx" = "#7d0000" # Rojo oscuro
  )) +
  
  # Limpiar etiquetas del eje X
  scale_x_discrete(labels = function(x) gsub("_(No Respondedor|Respondedor)", "", x)) +
  
  labs(
    title = "Infiltración Inmunitaria según Respuesta y Cohorte",
    x = "Tipo celular",
    y = "Fracción estimada (CIBERSORTx)",
    fill = "Grupo y Cohorte"
  ) +
  
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  )
