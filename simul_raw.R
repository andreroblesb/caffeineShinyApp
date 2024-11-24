t <- 1:200
N <- 100 
S <- rep(0, 200)  
I <- rep(0, 200) 
R <- rep(0, 200) 

pI <- 0.3
pR <- 0.01 
pC <- 0.5

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

plot(t, S, type = "l", col = "blue", ylim = c(0, N), xlab = "Dia", ylab = "Poblacion", main = "Modelo SIR")
lines(t, I, col = "red")
lines(t, R, col = "green")