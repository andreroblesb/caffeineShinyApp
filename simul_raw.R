library(shiny)
library(deSolve)
library(ggplot2)
library(bslib)

# Modelo Estocastico
# -----------------------------------------------------------
t <- 1:365
N <- 100
S <- rep(0, 365)
I <- rep(0, 365) 
R <- rep(0, 365) 
n <- rep(0, 365)

pI = 1
pR <- 0.01 # Saber si estaba bien antes e enfermarse o no etc.
pC <- 0.5 # Aqui esta establecido que AL MENOS la mitas de los INFECTADOS no los SUCEPTIBLES usan cubrebocas.

# Que pasa si los hospitales estan llenos?
# El primer dia no debe de haber los tres. 
# Debemos de cambiar las probabilidades conforme avanza el tiempo.
# Probabilidad condicional?
# pD <- 0.001 # Enfermedades graves antes de morir
# pV <- 0.02 # Probabilidad de que este vacunado
# pF <- 0.0001 # Dias festivos
# pHB <- 0.00001 # Probabilidad de cumpleanos con fiesta, sin subestimar los eventos.
# pvivienda <- 0# Promedio de personas que viven en cada casa 3.6 en mexico a comparacion de estados unidos 2.5. 
# Mundial?
# psecotrinformal <- 0#  Personas que trabajan en el mercado etc. Trabajadores que dependen del dia a dia. 
# Susptibilidad por edades. 
# Fuma o no.
# Clima.
# Densidad poblacional
# Region.
# Como cambias la probabilidad de ser infectado respecto al tiempo?
# Como cambias la probabilidad de recuperarte estando vacunado o no?

I[1] <- 1             
R[1] <- 0              
S[1] <- N - I[1] - R[1] 
# -----------------------------------------------------------


# Modelo Numerico
# -----------------------------------------------------------
mu <- 0.05    
beta <- 0.5   
gamma <- 0.1    
Nn <- 100           
S0 <- 100           
I0 <- 1            
R0 <- 0 
time <- seq(0, 100, by = 1)

sir_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- mu * N - (beta * I * S / N) - mu * S
    dI <- (beta * I * S / N) - gamma * I - mu * I
    dR <- gamma * I - mu * R
    return(list(c(dS, dI, dR)))
  })
}

rkf45 <- function(t0, X0, h, n, parameters) {
  t <- numeric(n + 1)
  S <- numeric(n + 1)
  I <- numeric(n + 1)
  R <- numeric(n + 1)
  
  t[1] <- t0
  S[1] <- X0[1]
  I[1] <- X0[2]
  R[1] <- X0[3]
  tol <- 1e-5
  i <- 1
  
  while (i <= n) {
    state <- c(S = S[i], I = I[i], R = R[i])
    
    k1 <- h * unlist(sir_model(t[i], state, parameters))
    k2 <- h * unlist(sir_model(t[i] + h / 4, state + k1 / 4, parameters))
    k3 <- h * unlist(sir_model(t[i] + 3 * h / 8, state + 3 * k1 / 32 + 9 * k2 / 32, parameters))
    k4 <- h * unlist(sir_model(t[i] + 12 * h / 13, state + 1932 * k1 / 2197 - 7200 * k2 / 2197 + 7296 * k3 / 2197, parameters))
    k5 <- h * unlist(sir_model(t[i] + h, state + 439 * k1 / 216 - 8 * k2 + 3680 * k3 / 513 - 845 * k4 / 4104, parameters))
    k6 <- h * unlist(sir_model(t[i] + h / 2, state - 8 * k1 / 27 + 2 * k2 - 3544 * k3 / 2565 + 1859 * k4 / 4104 - 11 * k5 / 40, parameters))
    
    y5 <- state + (16 * k1 / 135 + 6656 * k3 / 12825 + 28561 * k4 / 56430 - 9 * k5 / 50 + 2 * k6 / 55)
    y4 <- state + (25 * k1 / 216 + 1408 * k3 / 2565 + 2197 * k4 / 4104 - k5 / 5)
    error <- max(abs(y5 - y4))
    
    if (error < tol) {
      t[i + 1] <- t[i] + h
      S[i + 1] <- y5[1]
      I[i + 1] <- y5[2]
      R[i + 1] <- y5[3]
      i <- i + 1
      if (error < tol / 2) h <- h * 2  
    } else {
      h <- h / 2  
    }
    if (i == n) break  
  }
  
  return(data.frame(t = t, S = S, I = I, R = R))
}

# Definir los parámetros del modelo
parameters <- c(mu = mu, beta = beta, gamma = gamma)

# Estado inicial
t0 <- 0
X0 <- c(S = S0, I = I0, R = R0)
h <- 1
n <- length(time) - 1

# Resolver el sistema de ecuaciones diferenciales usando RKF45
sir_output <- rkf45(t0, X0, h, n, parameters)

# Graficar los resultados
library(ggplot2)
ggplot(sir_output, aes(x = t)) +
  geom_line(aes(y = S, color = "Susceptibles")) +
  geom_line(aes(y = I, color = "Infectados")) +
  geom_line(aes(y = R, color = "Recuperados")) +
  labs(title = "Modelo SIR Determinista con RKF45", x = "Tiempo", y = "Población") +
  scale_color_manual(values = c("Susceptibles" = "blue", "Infectados" = "red", "Recuperados" = "green")) +
  theme_minimal()
# -----------------------------------------------------------

# Loop Modelo estocastico
# -----------------------------------------------------------
for (day in 2:365) {
  
  face_mask <- rbinom(1, S[day - 1], pC)
  
  pI = runif(1) # Nueva probabilidad de infeccion para cada dia que pasa. 
  
  if (S[day - 1] > 0) {
    effective_pI <- pI * (1 - face_mask / S[day - 1])
  } else {
    effective_pI <- 0
  }

  new_infected <- rbinom(1, S[day - 1], effective_pI * I[day - 1] / N)
  new_recovered <- rbinom(1, I[day - 1], pR)
  
  S[day] <- S[day - 1] - new_infected
  I[day] <- I[day - 1] + new_infected - new_recovered
  R[day] <- R[day - 1] + new_recovered
  n[day] <- S[day] + I[day] + R[day]
}
# -----------------------------------------------------------


resultados <- data.frame(Dia = t, Susceptibles = S, Infectados = I, Recuperados = R, Poblacion = n)
print(resultados)
# Sacar una grafica para cada uno ademas de la combinada.
# Sacar una tabla para observar el comportamiento.

plot(t, S, type = "l", col = "blue", ylim = c(0, N), xlab = "Dia", ylab = "Poblacion", main = "Modelo SIR")
lines(t, I, col = "red")
lines(t, R, col = "green")
lines(t, n, col = "purple")