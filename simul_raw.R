t <- 1:200 # Cambiar el tiempo a 365 dias empezando desde 1 de enero.
N <- 100 # Mas personas? 1000? 2000? 10000?
S <- rep(0, 200)  # Actualizar los dias para cada uno de SIR.
I <- rep(0, 200) 
R <- rep(0, 200) 

pI <- 0.3 # Variable aleatoria usando runif(i) para cada dia que pasa, cambiando la proobabilidad ed infeccion.
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


for (day in 2:200) {
  face_mask <- rbinom(1, S[day - 1], pC)
  
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
}

data.frame(Dia = t, Susceptibles = S, Infectados = I, Recuperados = R)

# Sacar una grafica para cada uno ademas de la combinada.
# Sacar una tabla para observar el comportamiento.

plot(t, S, type = "l", col = "blue", ylim = c(0, N), xlab = "Dia", ylab = "Poblacion", main = "Modelo SIR")
lines(t, I, col = "red")
lines(t, R, col = "green")