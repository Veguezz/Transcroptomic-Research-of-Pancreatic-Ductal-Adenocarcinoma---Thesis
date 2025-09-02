library(fmsb)
library(scales)

valores <- data.frame(
  Drogabilidad = c(3, 3, 0),
  Toxicidad    = c(3, 3, 3),
  M   = c(0, 3, 0),
  A = c(3, 3, 3),
  Competencia = c(2, 2, 3),
  P = c(3, 3, 3),
  Novedad     = c(2, 0, 1),
  row.names   = c("CD3D","ITGAL","LYVE1")
)

orden_deseado <- c("Toxicidad", "Novedad", "P", "Competencia", "A", "M", "Drogabilidad")
valores <- valores[, orden_deseado]

# Define rango (ajusta los máximos si necesitas otro rango)
max_val <- rep(3, ncol(valores))
min_val <- rep(0, ncol(valores))
df_radar <- rbind(max_val, min_val, valores)

# 2) Ajusta el margen derecho
#    (mar = c(abajo, izquierda, arriba, derecha))
par(mar = c(1, 1, 1, 8))


# Parámetros de estilo (puedes ajustarlos)
radarchart(
  df_radar,
  axistype    = 1,               # etiquetas en cada eje radial
  seg         = 3,               # subdivisiones → círculos en 0,1,2,3
  pcol        = c("red","blue","green"),
  pfcol       = alpha(c("red","blue","green"), 0.3),
  plwd        = 2,
  cglcol      = "grey80",        # color de líneas concéntricas
  cglty       = 1,
  cglwd       = 0.8,             
  axislabcol  = "grey30",
  caxislabels = as.character(0:3),  # **aquí eliminas el “(%)”**
  vlcex       = 0.8               # tamaño de texto en vértices
)

# 4) Dibuja la leyenda fuera del área de trazado
par(xpd = NA)  
legend(
  "topleft",               # punto de referencia
  inset   = c(-0.25, 0),    # desplaza 25% hacia la derecha fuera del plot
  legend  = rownames(valores),
  col     = c("red","blue","green"),
  pch     = 20,
  pt.cex  = 1.5,
  bty     = "n",
  cex     = 0.9
)

