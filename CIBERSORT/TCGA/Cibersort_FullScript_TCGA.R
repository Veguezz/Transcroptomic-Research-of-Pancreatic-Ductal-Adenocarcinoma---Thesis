# LIBRERÍAS ####
library(e1071)
library(parallel)
library(preprocessCore)

# LECTURA DE GATOS ####

#Lectura Matriz LM22
LM22 <- read.table("C:/Users/VALERIA/OneDrive - UCB/PDAC/Cibersort/LM22.txt", header=TRUE, sep="\t", row.names=1, check.names=FALSE)

dim(LM22) #Verificamos dimensiones

head(colnames(LM22))

# Cargar la matriz de expresión normalizada
matriz_expres_PV4 <- readRDS("C:/Users/VALERIA/OneDrive - UCB/PDAC/Cibersort/norm_dds_counts_TCGA_V3.rds")

dim(matriz_expres_PV4) #Verificamos dimesiones

head(matriz_expres_PV4)

#Asegurar la Consistencia de genes entre ambas

# Extraer genes en común entre LM22 y la matriz de expresión
genes_comunes <- intersect(rownames(LM22), rownames(matriz_expres_PV4))

# Filtrar ambas matrices para que contengan solo esos genes
LM22 <- LM22[genes_comunes, ]
matriz_expres_PV4 <- matriz_expres_PV4[genes_comunes, ]

# CORE ALGORITHM ####

# Función CoreAlg (se mantiene igual que en el script original)
# Función CoreAlg adaptada
CoreAlg <- function(X, y, absolute, abs_method) {
  # Probar 3 valores de nu para el modelo SVM
  svn_itor <- 3
  
  # Función interna para entrenar el modelo con un valor específico de nu
  res <- function(i) {
    if (i == 1) { nus <- 0.25 }
    if (i == 2) { nus <- 0.5 }
    if (i == 3) { nus <- 0.75 }
    model <- svm(X, y, type = "nu-regression", kernel = "linear", nu = nus, scale = FALSE)
    model
  }
  
  # Entrenar los 3 modelos (en paralelo si el sistema lo permite)
  if (Sys.info()['sysname'] == 'Windows') {
    out <- mclapply(1:svn_itor, res, mc.cores = 1)
  } else {
    out <- mclapply(1:svn_itor, res, mc.cores = svn_itor)
  }
  
  # Inicializar vectores para almacenar las métricas de cada modelo
  nusvm <- rep(0, svn_itor)
  corrv <- rep(0, svn_itor)
  
  # CIBERSORT ####
  
  t <- 1
  while (t <= svn_itor) {
    weights <- t(out[[t]]$coefs) %*% out[[t]]$SV
    weights[weights < 0] <- 0
    w <- weights / sum(weights)
    
    u <- sweep(X, MARGIN = 2, w, "*")
    k <- apply(u, 1, sum)
    
    nusvm[t] <- sqrt(mean((k - y)^2))
    corrv[t] <- cor(k, y)
    t <- t + 1
  }
  
  # Seleccionar el modelo con menor RMSE
  rmses <- nusvm
  mn <- which.min(rmses)
  best_model <- out[[mn]]
  
  # Calcular y normalizar los coeficientes del mejor modelo
  q <- t(best_model$coefs) %*% best_model$SV
  q[q < 0] <- 0
  if (!absolute || abs_method == 'sig.score')
    w_final <- q / sum(q)  # Espacio relativo (fracciones)
  if (absolute && abs_method == 'no.sumto1')
    w_final <- q          # Espacio absoluto (scores)
  
  mix_rmse <- rmses[mn]
  mix_r <- corrv[mn]
  
  newList <- list("w" = w_final, "mix_rmse" = mix_rmse, "mix_r" = mix_r)
  return(newList)
}

# Función doPerm: realiza permutaciones para estimar una distribución nula
doPerm <- function(perm, X, Y, absolute, abs_method) {
  itor <- 1
  Ylist <- as.list(data.matrix(Y))
  dist <- matrix()
  
  while (itor <= perm) {
    # Generar una mezcla aleatoria a partir de Y
    yr <- as.numeric(Ylist[sample(length(Ylist), dim(X)[1])])
    yr <- (yr - mean(yr)) / sd(yr)
    
    # Ejecutar CoreAlg para la mezcla aleatoria
    result <- CoreAlg(X, yr, absolute, abs_method)
    mix_r <- result$mix_r
    
    if (itor == 1) { 
      dist <- mix_r 
    } else {
      dist <- rbind(dist, mix_r)
    }
    
    itor <- itor + 1
  }
  newList <- list("dist" = dist)
  return(newList)
}

# Función principal CIBERSORT #----
CIBERSORT <- function(sig_matrix, mixture_file, perm = 0, QN = TRUE, absolute = FALSE, abs_method = 'sig.score') {
  if (absolute && abs_method != 'no.sumto1' && abs_method != 'sig.score')
    stop("abs_method must be set to either 'sig.score' or 'no.sumto1'")
  
  # Utilizar los datos ya cargados (evitando read.table para acelerar el proceso)
  X <- sig_matrix
  Y <- mixture_file
  
  # Corregir posibles duplicados en los nombres de genes en Y
  dups <- dim(Y)[1] - length(unique(rownames(Y)))
  if (dups > 0) {
    warning(paste(dups, " duplicated gene symbol(s) found in mixture file!", sep = ""))
    rownames(Y) <- make.names(rownames(Y), unique = TRUE)
  }
  
  X <- data.matrix(X)
  Y <- data.matrix(Y)
  
  # Ordenar filas por nombre de gen
  X <- X[order(rownames(X)), ]
  Y <- Y[order(rownames(Y)), ]
  
  P <- perm
  
  # Si los valores son bajos (p.ej., log-scale), revertir transformándolos exponencialmente
  if (max(na.omit(as.numeric(Y))) < 50) { Y <- 2^Y }
  
  # Si se solicita, normalización cuantílica
  if (QN) {
    tmpc <- colnames(Y)
    tmpr <- rownames(Y)
    Y <- normalize.quantiles(Y)
    colnames(Y) <- tmpc
    rownames(Y) <- tmpr
  }
  
  # Guardar la matriz original y calcular la mediana
  Yorig <- Y
  Ymedian <- max(median(Yorig), 1)
  
  # Intersección de genes entre X y Y
  Xgns <- rownames(X)
  Ygns <- rownames(Y)
  YintX <- Ygns %in% Xgns
  Y <- Y[YintX, ]
  XintY <- Xgns %in% rownames(Y)
  X <- X[XintY, ]
  
  # Estandarizar la matriz de firmas
  X <- (X - mean(X)) / sd(as.vector(X))
  
  # Cálculo de la distribución nula por permutaciones (si perm > 0)
  if (P > 0) { nulldist <- sort(doPerm(P, X, Y, absolute, abs_method)$dist) }
  
  header <- c('Mixture', colnames(X), "P-value", "Correlation", "RMSE")
  if (absolute) header <- c(header, paste('Absolute score (', abs_method, ')', sep = ""))
  
  output <- matrix()
  itor <- 1
  mixtures <- dim(Y)[2]
  pval <- 9999
  
  while (itor <= mixtures) {
    y <- Y[, itor]
    # Estandarizar la muestra (mixture)
    y <- (y - mean(y)) / sd(y)
    
    # Ejecutar CoreAlg sobre la muestra
    result <- CoreAlg(X, y, absolute, abs_method)
    w <- result$w
    mix_r <- result$mix_r
    mix_rmse <- result$mix_rmse
    
    if (absolute && abs_method == 'sig.score') {
      w <- w * median(Y[, itor]) / Ymedian
    }
    
    if (P > 0) { pval <- 1 - (which.min(abs(nulldist - mix_r)) / length(nulldist)) }
    
    out <- c(colnames(Y)[itor], w, pval, mix_r, mix_rmse)
    if (absolute) out <- c(out, sum(w))
    
    if (itor == 1) {
      output <- out
    } else {
      output <- rbind(output, out)
    }
    itor <- itor + 1
  }
  
  # Preparar la salida final (tabla numérica)
  obj <- rbind(header, output)
  obj <- obj[, -1]
  obj <- obj[-1, ]
  obj <- matrix(as.numeric(unlist(obj)), nrow = nrow(obj))
  rownames(obj) <- colnames(Y)
  if (!absolute) {
    colnames(obj) <- c(colnames(X), "P-value", "Correlation", "RMSE")
  } else {
    colnames(obj) <- c(colnames(X), "P-value", "Correlation", "RMSE", paste('Absolute score (', abs_method, ')', sep = ""))
  }
  return(obj)
}

# ANALISIS POR CHUNKS ----


# Si deseas procesar todas las muestras a la vez usando la función CIBERSORT:
Y_final <- matriz_expres_PV4  # O bien, si la función CIBERSORT ya hace el preprocesamiento, usa esa variable.

# Definir el número total de muestras
total_muestras <- ncol(Y_final)

# Establecer el tamaño de chunk a 4
chunk_size <- 4

# Dividir los índices de las muestras en chunks de tamaño 4
indices_chunks <- split(1:total_muestras, ceiling(seq_along(1:total_muestras) / chunk_size))

# Lista para almacenar los resultados parciales (opcional)
resultados_chunks <- list()

# Procesar cada chunk de 4 muestras
for (chunk in indices_chunks) {
  cat("Procesando el chunk con muestras:", chunk, "\n")
  
  # Ejecutar la función CIBERSORT para el subset de muestras en este chunk.
  # Se asume que CIBERSORT acepta matrices de mezcla con varias muestras (columnas).
  partial_result <- CIBERSORT(LM22, Y_final[, chunk], perm = 50, QN = TRUE, absolute = FALSE, abs_method = 'sig.score')
  
  # Guardar el resultado parcial en una lista (clave basada en el rango de índices)
  chunk_name <- paste0("chunk_", min(chunk), "_to_", max(chunk))
  resultados_chunks[[chunk_name]] <- partial_result
  
  # Guardar a disco el resultado parcial
  saveRDS(partial_result, file = paste0("C:/Users/VALERIA/OneDrive - UCB/PDAC/Cibersort/Deconvolución_Celular_PanCuRx_V4", chunk_name, ".rds"))
  
  # Opcional: Pausa para dar "respiro" a la computadora
  Sys.sleep(300)  # Pausa de 10 segundos
}

# GUARDAR RDS ----

# Verifica la estructura de cada elemento (por ejemplo, todas deben ser matrices con el mismo número de columnas).
str(resultados_chunks)

# Combina todos los resultados parciales en una sola matriz utilizando do.call y rbind:
all_results <- do.call(rbind, resultados_chunks)

# Verifica el resultado combinado:
print(dim(all_results))
head(all_results)

# Guarda el objeto combinado en un archivo RDS para poder cargarlo en sesiones futuras sin tener que rehacer el análisis
saveRDS(all_results, file = "C:/Users/VALERIA/OneDrive - UCB/PDAC/Cibersort/Perfil_Deconvolución_Celular_TCGA_V3.rds")
