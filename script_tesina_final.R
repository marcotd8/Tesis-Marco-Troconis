# =====================================================================
# TESINA: Analisis de la sostenibilidad financiera de esquemas de retiro
# privados mediante pruebas de estres
# Script completo (version final)
# =====================================================================
# Antes de correr:
#   1) Ajustar ruta_tabla y carpeta_salida (mas abajo).
#   2) La tabla de mortalidad es un CSV con separador ";" y decimal ",",
#      con dos columnas (Edad ; q(x)) para las edades 65 a 99.
# Se recomienda correr por bloques (Cmd/Ctrl + Enter) verificando la salida.
# =====================================================================

ruta_tabla     <- "Tabla Mortalidad R.csv"
carpeta_salida <- "."

# Formato de numeros grandes (miles con punto)
fmt <- function(x) format(round(x), big.mark = ".", decimal.mark = ",", scientific = FALSE)

# =====================================================================
# BLOQUE 1: PARAMETROS GENERALES
# =====================================================================

edad_retiro     <- 65
anos_proyeccion <- 35
n_asegurados    <- 500
tasa_tecnica    <- 0.04
renta           <- 2608125.49     # B = 30.000.000 / a(65)

# Parametros de Vasicek
kappa <- 0.3
theta <- 0.04
sigma <- 0.01
r0    <- 0.04

# Simulacion (no agregar nada entre set.seed y el bloque 2)
n_simulaciones <- 10000
set.seed(123)

# =====================================================================
# BLOQUE 2: SIMULACION DE VASICEK (MONTE CARLO)
# =====================================================================

tasas <- matrix(0, nrow = n_simulaciones, ncol = anos_proyeccion)

for (sim in 1:n_simulaciones) {
  r <- r0
  for (t in 1:anos_proyeccion) {
    r <- r + kappa * (theta - r) + sigma * rnorm(1)
    tasas[sim, t] <- r
  }
}

cat("Simulacion de Vasicek completada\n")
cat("Tasa promedio ano 1:", round(mean(tasas[, 1]) * 100, 4), "%\n")

# =====================================================================
# BLOQUE 3: TABLA DE MORTALIDAD
# =====================================================================

tabla_mort <- read.csv(ruta_tabla, sep = ";", dec = ",")
names(tabla_mort) <- c("Edad", "q.x.")

# Se agrega la edad limite (100) con qx = 1
tabla_mort <- rbind(tabla_mort, data.frame(Edad = 100, q.x. = 1))

edades <- tabla_mort$Edad
qx     <- tabla_mort$q.x.
px     <- 1 - qx

cat("Tabla de mortalidad cargada. Edades:", min(edades), "a", max(edades), "\n")

# =====================================================================
# BLOQUE 4: SOBREVIVIENTES (Ecuacion 2)
# =====================================================================

l <- numeric(length(edades))
l[1] <- n_asegurados

for (i in 2:length(edades)) {
  l[i] <- l[i - 1] * px[i - 1]
}

cat("l_65:", round(l[1], 2), "| l_100:", round(l[length(l)], 4), "\n")

# =====================================================================
# BLOQUE 5: VALOR ACTUAL ACTUARIAL a(x) ANTICIPADO (Ecuacion 5)
# =====================================================================

ax <- numeric(length(edades))
ax[length(edades)] <- 1

for (i in (length(edades) - 1):1) {
  ax[i] <- 1 + px[i] * (1 / (1 + tasa_tecnica)) * ax[i + 1]
}

cat("a_65:", format(ax[1], digits = 10), "\n")

# =====================================================================
# BLOQUE 6: ACTIVO INICIAL Y PASIVO ACTUARIAL (Ecuaciones 8 y 9)
# =====================================================================

activo_inicial <- renta * ax[1] * l[1]

pasivo <- numeric(length(edades))
for (i in 1:length(edades)) {
  pasivo[i] <- l[i] * renta * ax[i]
}

cat("Activo inicial:", fmt(activo_inicial), "\n")
cat("Pasivo ano 0  :", fmt(pasivo[1]), "\n")

# =====================================================================
# BLOQUE 7: PROYECCION DEL ACTIVO EN LAS 10.000 SIMULACIONES (Ecuacion 12)
# activos[sim, k] = activo al inicio del ano k, antes del pago de la renta
# =====================================================================

activos <- matrix(0, nrow = n_simulaciones, ncol = length(edades))

for (sim in 1:n_simulaciones) {
  activo <- activo_inicial
  for (t in 1:length(edades)) {
    activos[sim, t] <- activo                 # antes del pago
    activo <- activo - renta * l[t]           # se paga la renta
    if (t < length(edades)) {
      activo <- activo * (1 + tasas[sim, t])  # se capitaliza
    }
  }
}

cat("Proyeccion de activos completada\n")

# =====================================================================
# BLOQUE 8: DISTRIBUCION DEL RESULTADO FINAL (Cuadro 8)
# =====================================================================

resultado_final <- activos[, length(edades)]

cat("\n--- DISTRIBUCION DEL ACTIVO FINAL ---\n")
cat("Percentil 1% :", fmt(as.numeric(quantile(resultado_final, 0.01))), "\n")
cat("Percentil 5% :", fmt(as.numeric(quantile(resultado_final, 0.05))), "\n")
cat("Mediana      :", fmt(as.numeric(quantile(resultado_final, 0.50))), "\n")
cat("Media        :", fmt(mean(resultado_final)), "\n")
cat("Percentil 95%:", fmt(as.numeric(quantile(resultado_final, 0.95))), "\n")
cat("Probabilidad de deficit:", round(mean(resultado_final < 0) * 100, 2), "%\n")

# Simulacion mas cercana al percentil 5% (escenario de estres de tasa)
percentil_5 <- quantile(resultado_final, 0.05)
indice_p5   <- which.min(abs(resultado_final - percentil_5))

cat("\nTasa promedio trayectoria P5, anos 1-35:",
    round(mean(tasas[indice_p5, ]) * 100, 2), "%\n")
cat("Tasa promedio trayectoria P5, anos 1-10:",
    round(mean(tasas[indice_p5, 1:10]) * 100, 2), "%\n")

# =====================================================================
# BLOQUE 9: ESCENARIO BASE Y ESTRES DE TASA
# =====================================================================

tasas_promedio <- colMeans(tasas)
anos <- 1:length(edades)

# Escenario base: tasas promedio de Vasicek
activo_base_vasicek <- numeric(length(edades))
activo_base_vasicek[1] <- activo_inicial
for (t in 1:(length(edades) - 1)) {
  activo_base_vasicek[t + 1] <- (activo_base_vasicek[t] - renta * l[t]) *
                                (1 + tasas_promedio[t])
}

# Activo con tasa tecnica fija (verificacion: termina en el ultimo pago)
activo_base <- numeric(length(edades))
activo_base[1] <- activo_inicial
for (t in 1:(length(edades) - 1)) {
  activo_base[t + 1] <- (activo_base[t] - renta * l[t]) * (1 + tasa_tecnica)
}

# Estres de tasa: trayectoria del percentil 5%
activo_estres    <- activos[indice_p5, ]
ano_quiebra_tasa <- anos[which(activo_estres < 0)[1]]

# =====================================================================
# BLOQUE 10: ESTRES DE LONGEVIDAD (Ecuacion 16)
# =====================================================================

lambda <- 0.015

qx_estres_long <- numeric(length(edades))
for (i in 1:length(edades)) {
  t <- i - 1
  qx_estres_long[i] <- tabla_mort$q.x.[i] * (1 - lambda)^t
}
px_estres_long <- 1 - qx_estres_long

l_estres_long <- numeric(length(edades))
l_estres_long[1] <- n_asegurados
for (i in 2:length(edades)) {
  l_estres_long[i] <- l_estres_long[i - 1] * px_estres_long[i - 1]
}

activo_estres_long <- numeric(length(edades))
activo_estres_long[1] <- activo_inicial
for (t in 1:(length(edades) - 1)) {
  activo_estres_long[t + 1] <- (activo_estres_long[t] - renta * l_estres_long[t]) *
                               (1 + tasa_tecnica)
}
ano_quiebra_long <- anos[which(activo_estres_long < 0)[1]]

# =====================================================================
# BLOQUE 11: SENSIBILIDAD A LAMBDA Y TASA DE EQUILIBRIO (Ecuacion 17)
# =====================================================================

lambdas <- c(0.01, 0.015, 0.02)
colores <- c("green", "red", "purple")
activos_lambda <- list()
l_fin_lambda   <- numeric(length(lambdas))
tasa_eq_lambda <- numeric(length(lambdas))

for (j in 1:length(lambdas)) {
  qx_lambda <- numeric(length(edades))
  for (i in 1:length(edades)) {
    t <- i - 1
    qx_lambda[i] <- tabla_mort$q.x.[i] * (1 - lambdas[j])^t
  }
  px_lambda <- 1 - qx_lambda

  l_lambda <- numeric(length(edades))
  l_lambda[1] <- n_asegurados
  for (i in 2:length(edades)) {
    l_lambda[i] <- l_lambda[i - 1] * px_lambda[i - 1]
  }

  activo_lambda <- numeric(length(edades))
  activo_lambda[1] <- activo_inicial
  for (t in 1:(length(edades) - 1)) {
    activo_lambda[t + 1] <- (activo_lambda[t] - renta * l_lambda[t]) *
                            (1 + tasa_tecnica)
  }
  activos_lambda[[j]] <- activo_lambda
  l_fin_lambda[j]     <- l_lambda[length(edades)]

  # Tasa de equilibrio i*: activo final = 0, incluyendo el ultimo pago
  activo_final_lambda <- function(tasa) {
    activo <- activo_inicial
    for (t in 1:(length(edades) - 1)) {
      activo <- (activo - renta * l_lambda[t]) * (1 + tasa)
    }
    activo <- activo - renta * l_lambda[length(edades)]
    return(activo)
  }
  tasa_eq_lambda[j] <- uniroot(activo_final_lambda,
                               interval = c(0.01, 0.20),
                               tol = 1e-12)$root
}

# =====================================================================
# BLOQUE 12: ESCENARIO COMBINADO
# (trayectoria de tasas del P5 obtenido con mortalidad base + mortalidad
#  con mejora lambda = 1,5%)
# =====================================================================

activo_combinado <- numeric(length(edades))
activo_combinado[1] <- activo_inicial
for (t in 1:(length(edades) - 1)) {
  activo_combinado[t + 1] <- (activo_combinado[t] - renta * l_estres_long[t]) *
                             (1 + tasas[indice_p5, t])
}
ano_quiebra_combinado <- anos[which(activo_combinado < 0)[1]]

# =====================================================================
# BLOQUE 13: INDICADORES DE RESULTADO (seccion 3.8)
#   Ano de quiebra k*  : primer ano cuyo inicio tiene activo < 0
#   Edad al inicio     : 64 + k*
#   Deficit final D    : max(0, -A_T), A_T = activo al inicio del ultimo
#                        ano (t = 35), antes del pago de esa renta
#   Interaccion I      : D_comb - (D_tasa + D_long)
# =====================================================================

T_final <- length(edades)

D_tasa  <- max(0, -activo_estres[T_final])
D_long  <- max(0, -activo_estres_long[T_final])
D_comb  <- max(0, -activo_combinado[T_final])
I_inter <- D_comb - (D_tasa + D_long)

cat("\n--- INDICADORES (seccion 3.8) ---\n")
cat("Ano de quiebra: tasa", ano_quiebra_tasa,
    "| longevidad", ano_quiebra_long,
    "| combinado", ano_quiebra_combinado, "\n")
cat("Edad al inicio de ese ano: tasa", 64 + ano_quiebra_tasa,
    "| longevidad", 64 + ano_quiebra_long,
    "| combinado", 64 + ano_quiebra_combinado, "\n")

cat("\nDeficit final D (pesos):\n")
cat("  Estres de tasa   :", fmt(D_tasa), "|", round(100 * D_tasa / activo_inicial, 1), "% de A0\n")
cat("  Estres longevidad:", fmt(D_long), "|", round(100 * D_long / activo_inicial, 1), "% de A0\n")
cat("  Estres combinado :", fmt(D_comb), "|", round(100 * D_comb / activo_inicial, 1), "% de A0\n")

cat("\nEn miles de millones (3 decimales):",
    round(D_tasa / 1e9, 3), "(tasa) +", round(D_long / 1e9, 3), "(long) | combinado:",
    round(D_comb / 1e9, 3), "\n")
cat("Interaccion I:", fmt(I_inter), "| % de la suma individual:",
    round(100 * I_inter / (D_tasa + D_long), 1), "% | % de A0:",
    round(100 * I_inter / activo_inicial, 1), "%\n")
cat("Cocientes: comb/tasa =", round(D_comb / D_tasa, 2),
    "| comb/long =", round(D_comb / D_long, 2),
    "| tasa/long =", round(D_tasa / D_long, 2), "\n")
cat("Participacion en D_comb: tasa", round(100 * D_tasa / D_comb, 1),
    "% | longevidad", round(100 * D_long / D_comb, 1), "%\n")

cat("\n--- SENSIBILIDAD A LAMBDA (Cuadro 9) ---\n")
for (j in 1:length(lambdas)) {
  ano_q <- anos[which(activos_lambda[[j]] < 0)[1]]
  D_j   <- max(0, -activos_lambda[[j]][T_final])
  cat("lambda =", lambdas[j] * 100, "% | ano quiebra:", ano_q,
      "| edad:", 64 + ano_q,
      "| deficit final:", fmt(D_j), "(", round(100 * D_j / activo_inicial, 1), "% de A0)\n")
  cat("   sobrevivientes a los 100:", round(l_fin_lambda[j], 2),
      "| extra vs base:", round(l_fin_lambda[j] - l[length(l)], 2),
      "| i*:", round(100 * tasa_eq_lambda[j], 4), "%",
      "| spread:", round(10000 * (tasa_eq_lambda[j] - tasa_tecnica), 1), "pb\n")
}
cat("Sobrevivientes a los 100 anos, escenario base:", round(l[length(l)], 2), "\n")

# =====================================================================
# BLOQUE 14: DISTRIBUCION DEL ANO DE QUIEBRA (Grafico 9)
# =====================================================================

ano_quiebra_sim <- numeric(n_simulaciones)
for (sim in 1:n_simulaciones) {
  quiebra <- which(activos[sim, ] < 0)[1]
  if (is.na(quiebra)) {
    ano_quiebra_sim[sim] <- NA
  } else {
    ano_quiebra_sim[sim] <- quiebra
  }
}

cat("\nSimulaciones con quiebra:", round(mean(!is.na(ano_quiebra_sim)) * 100, 2), "%\n")
cat("Simulaciones sin quiebra:", round(mean(is.na(ano_quiebra_sim)) * 100, 2), "%\n")

# Solo anos dentro del horizonte (1 a 35)
anos_con_quiebra <- ano_quiebra_sim[!is.na(ano_quiebra_sim) & ano_quiebra_sim <= 35]

# =====================================================================
# BLOQUE 15: GRAFICOS (se guardan en carpeta_salida)
# Los titulos usan codigos \u para tildes y enies.
# =====================================================================

setwd(carpeta_salida)

# GRAFICO 1
png("grafico1_trayectorias_vasicek.png", width = 800, height = 600)
plot(1:anos_proyeccion, tasas[1, ], type = "l", col = rgb(0, 0, 1, 0.1),
     ylim = c(min(tasas), max(tasas)),
     xlab = "A\u00f1o", ylab = "Tasa de inter\u00e9s",
     main = "Gr\u00e1fico 1: Trayectorias simuladas - Modelo de Vasicek")
for (sim in 2:100) {
  lines(1:anos_proyeccion, tasas[sim, ], col = rgb(0, 0, 1, 0.1))
}
lines(1:anos_proyeccion, colMeans(tasas), col = "red", lwd = 2)
abline(h = theta, col = "black", lty = 2)
legend("topright", legend = c("Trayectorias", "Promedio", "Theta"),
       col = c("blue", "red", "black"), lty = c(1, 1, 2))
dev.off()

# GRAFICO 2
png("grafico2_distribucion_activo_final.png", width = 800, height = 600)
hist(resultado_final / 1e9, breaks = 100,
     main = "Gr\u00e1fico 2: Distribuci\u00f3n del Activo Final",
     xlab = "Activo final (miles de millones)",
     col = "lightblue")
abline(v = percentil_5 / 1e9, col = "red", lwd = 2)
legend("topright", legend = "Percentil 5%", col = "red", lwd = 2)
dev.off()

# GRAFICO 3
png("grafico3_trayectoria_estres_tasa.png", width = 800, height = 600)
plot(1:anos_proyeccion, tasas[indice_p5, ], type = "l", col = "red", lwd = 2,
     xlab = "A\u00f1o", ylab = "Tasa",
     main = "Gr\u00e1fico 3: Trayectoria de Tasa - Escenario de Estr\u00e9s (P5%)")
abline(h = theta, col = "black", lty = 2)
lines(1:anos_proyeccion, colMeans(tasas), col = "blue", lwd = 2)
legend("topright", legend = c("Estr\u00e9s (P5%)", "Promedio", "Theta"),
       col = c("red", "blue", "black"), lty = c(1, 1, 2))
dev.off()

# GRAFICO 4
png("grafico4_activo_base_vs_estres_tasa.png", width = 800, height = 600)
plot(anos[1:35], activo_base_vasicek[1:35] / 1e9, type = "l", col = "blue", lwd = 2,
     xlab = "A\u00f1o", ylab = "Activo (miles de millones)",
     main = "Gr\u00e1fico 4: Evoluci\u00f3n del Activo: Base vs Estr\u00e9s de Tasa")
lines(anos[1:35], activo_estres[1:35] / 1e9, col = "red", lwd = 2)
abline(h = 0, col = "black", lty = 2)
legend("topright", legend = c("Base (promedio Vasicek)", "Estr\u00e9s tasa (P5%)"),
       col = c("blue", "red"), lwd = 2)
dev.off()

# GRAFICO 5
png("grafico5_sobrevivientes_base_vs_estres_long.png", width = 800, height = 600)
plot(anos[1:35], l[1:35], type = "l", col = "blue", lwd = 2,
     xlab = "A\u00f1o", ylab = "Sobrevivientes",
     main = "Gr\u00e1fico 5: Sobrevivientes: Base vs Estr\u00e9s de Longevidad")
lines(anos[1:35], l_estres_long[1:35], col = "red", lwd = 2)
legend("topright", legend = c("Base", "Estr\u00e9s longevidad"),
       col = c("blue", "red"), lwd = 2)
dev.off()

# GRAFICO 6
png("grafico6_activo_base_vs_estres_long.png", width = 800, height = 600)
plot(anos[1:35], activo_base_vasicek[1:35] / 1e9, type = "l", col = "blue", lwd = 2,
     xlab = "A\u00f1o", ylab = "Activo (miles de millones)",
     main = "Gr\u00e1fico 6: Evoluci\u00f3n del Activo: Base vs Estr\u00e9s de Longevidad")
lines(anos[1:35], activo_estres_long[1:35] / 1e9, col = "red", lwd = 2)
abline(h = 0, col = "black", lty = 2)
legend("topright", legend = c("Base (promedio Vasicek)", "Estr\u00e9s longevidad"),
       col = c("blue", "red"), lwd = 2)
dev.off()

# GRAFICO 7
png("grafico7_sensibilidad_lambda.png", width = 800, height = 600)
plot(anos[1:35], activo_base_vasicek[1:35] / 1e9, type = "l", col = "blue", lwd = 2,
     xlab = "A\u00f1o", ylab = "Activo (miles de millones)",
     main = "Gr\u00e1fico 7: An\u00e1lisis de Sensibilidad - Factor de Mejora de Mortalidad")
for (j in 1:length(lambdas)) {
  lines(anos[1:35], activos_lambda[[j]][1:35] / 1e9, col = colores[j], lwd = 2)
}
abline(h = 0, col = "black", lty = 2)
legend("topright",
       legend = c("Base (promedio Vasicek)", "\u03bb=1%", "\u03bb=1.5%", "\u03bb=2%"),
       col = c("blue", colores), lwd = 2)
dev.off()

# GRAFICO 8
png("grafico8_escenario_combinado.png", width = 800, height = 600)
plot(anos[1:35], activo_base_vasicek[1:35] / 1e9, type = "l", col = "blue", lwd = 2,
     ylim = c(min(activo_combinado[1:35] / 1e9), max(activo_base_vasicek[1:35] / 1e9)),
     xlab = "A\u00f1o", ylab = "Activo (miles de millones)",
     main = "Gr\u00e1fico 8: Evoluci\u00f3n del Activo: Todos los Escenarios")
lines(anos[1:35], activo_estres[1:35] / 1e9, col = "orange", lwd = 2)
lines(anos[1:35], activo_estres_long[1:35] / 1e9, col = "red", lwd = 2)
lines(anos[1:35], activo_combinado[1:35] / 1e9, col = "purple", lwd = 2)
abline(h = 0, col = "black", lty = 2)
legend("topright",
       legend = c("Base (promedio Vasicek)", "Estr\u00e9s tasa (P5%)",
                  "Estr\u00e9s longevidad (\u03bb=1.5%)", "Estr\u00e9s combinado"),
       col = c("blue", "orange", "red", "purple"), lwd = 2)
dev.off()

# GRAFICO 9
png("grafico9_distribucion_ano_quiebra.png", width = 800, height = 600)
hist(anos_con_quiebra, breaks = 35,
     main = "Gr\u00e1fico 9: Distribuci\u00f3n del a\u00f1o de quiebra",
     xlab = "A\u00f1o de quiebra",
     ylab = "Frecuencia",
     col = "lightcoral",
     xlim = c(1, 35))
abline(v = ano_quiebra_tasa, col = "red", lwd = 2, lty = 2)
legend("topright",
       legend = paste0("Estr\u00e9s tasa P5%: a\u00f1o ", ano_quiebra_tasa),
       col = "red", lwd = 2, lty = 2)
dev.off()

cat("\nGraficos guardados en:", getwd(), "\n")

# =====================================================================
# BLOQUE 16: EXPORTAR TABLA RESUMEN
# =====================================================================

tabla_resultados <- data.frame(
  Escenario = c("Base", "Estres de tasa (P5%)",
                "Estres de longevidad (lambda=1.5%)", "Estres combinado"),
  Ano_quiebra = c(NA, ano_quiebra_tasa, ano_quiebra_long, ano_quiebra_combinado),
  Edad_inicio_ano_quiebra = c(NA, 64 + ano_quiebra_tasa, 64 + ano_quiebra_long,
                              64 + ano_quiebra_combinado),
  Deficit_final = c(0, round(D_tasa), round(D_long), round(D_comb)),
  Deficit_pct_A0 = c(0, round(100 * D_tasa / activo_inicial, 1),
                     round(100 * D_long / activo_inicial, 1),
                     round(100 * D_comb / activo_inicial, 1))
)

write.csv2(tabla_resultados, "tabla_resumen_completa.csv", row.names = FALSE)

cat("\n=== SCRIPT COMPLETADO ===\n")
