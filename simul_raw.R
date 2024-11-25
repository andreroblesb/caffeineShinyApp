library(shiny)
library(deSolve)
library(ggplot2)
library(bslib)

# Interfaz Shiny para el Modelo SIR
ui <- fluidPage(
  theme = bs_theme(version = 4, bootswatch = "flatly"),
  titlePanel("Modelos SIR: Estocástico y Numérico"),
  sidebarLayout(
    sidebarPanel(
      h4("Parámetros del Modelo Numérico"),
      numericInput("mu", "Tasa de natalidad/mortalidad (mu)", value = 0.05, min = 0, max = 1, step = 0.01),
      numericInput("beta", "Tasa de transmisión (beta)", value = 0.5, min = 0, max = 2, step = 0.01),
      numericInput("gamma", "Tasa de recuperación (gamma)", value = 0.1, min = 0, max = 1, step = 0.01),
      numericInput("N", "Tamaño de la población", value = 3000, min = 1, max = 10000, step = 1),
      actionButton("run_model", "Ejecutar Modelo")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Modelo Numérico", plotOutput("sirPlotNum"), tableOutput("sirTableNum")),
        tabPanel("Modelo Estocástico", plotOutput("sirPlotEst"), tableOutput("sirTableEst"))
      )
    )
  )
)

server <- function(input, output) {
  observeEvent(input$run_model, {
    # Modelo Numérico
    mu <- input$mu
    beta <- input$beta
    gamma <- input$gamma
    N <- input$N
    
    S0 <- 2700
    I0 <- 300
    R0 <- 0
    
    time <- seq(0, 100, by = 1)
    
    sir_model <- function(t, state, parameters) {
      with(as.list(c(state, parameters)), {
        dS <- mu * N - (beta * I * S / N) - mu * S
        dI <- (beta * I * S / N) - gamma * I - mu * I
        dR <- gamma * I - mu * R
        list(c(dS, dI, dR))
      })
    }
    
    state <- c(S = S0, I = I0, R = R0)
    parameters <- c(mu = mu, beta = beta, gamma = gamma, N = N)
    
    sir_output <- tryCatch({
      ode(y = state, times = time, func = sir_model, parms = parameters)
    }, error = function(e) {
      NULL
    })
    
    if (!is.null(sir_output)) {
      sir_output <- as.data.frame(sir_output)
      
      output$sirPlotNum <- renderPlot({
        plot(sir_output$time, sir_output$S, type = "l", col = "blue", ylim = c(0, N), xlab = "Tiempo", ylab = "Población", lwd = 2, main = "Modelo Numérico SIR")
        lines(sir_output$time, sir_output$I, col = "red", lwd = 2)
        lines(sir_output$time, sir_output$R, col = "green", lwd = 2)
        legend("right", legend = c("Susceptible", "Infectado", "Recuperado"), col = c("blue", "red", "green"), lty = 1, lwd = 2)
      })
      
      output$sirTableNum <- renderTable({
        sir_output
      })
    } else {
      output$sirPlotNum <- renderPlot({
        plot(1, type = "n", xlab = "", ylab = "", xlim = c(0, 1), ylim = c(0, 1))
        text(0.5, 0.5, "Error en la simulación del modelo numérico", cex = 1.5)
      })
      
      output$sirTableNum <- renderTable({
        data.frame(Mensaje = "Error en la simulación del modelo numérico")
      })
    }
    
    # Modelo Estocástico
    t <- 1:365
    S <- rep(0, 365)
    I <- rep(0, 365)
    R <- rep(0, 365)
    n <- rep(0, 365)
    
    pI <- 1
    pR <- 0.01
    pC <- 0.5
    
    N <- 100
    I[1] <- 1
    R[1] <- 0
    S[1] <- N - I[1] - R[1]
    
    for (day in 2:365) {
      face_mask <- rbinom(1, S[day - 1], pC)
      pI <- runif(1)
      
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
    
    resultados <- data.frame(Dia = t, Susceptibles = S, Infectados = I, Recuperados = R, Población = n)
    
    output$sirPlotEst <- renderPlot({
      plot(t, S, type = "l", col = "blue", ylim = c(0, N), xlab = "Día", ylab = "Población", main = "Modelo Estocástico SIR")
      lines(t, I, col = "red")
      lines(t, R, col = "green")
      lines(t, n, col = "purple")
      legend("right", legend = c("Susceptible", "Infectado", "Recuperado", "Población Total"), col = c("blue", "red", "green", "purple"), lty = 1)
    })
    
    output$sirTableEst <- renderTable({
      resultados
    })
  })
}

shinyApp(ui = ui, server = server)
