library(shiny)
library(deSolve)
library(ggplot2)
library(bslib)

# INTERFAZ
# Define the UI
# Define la UI
ui <- navbarPage(
  theme = bs_theme(version = 4, bootswatch = "flatly"),
  title = "Simulación de Modelos SIR",
  
  # Primera pestaña: Inicio
  tabPanel("Inicio",
           fluidPage(
             h3("Bienvenido a la Simulación de Modelos SIR"),
             p("Esta aplicación permite simular modelos SIR (Susceptibles, Infectados, Recuperados) 
                usando enfoques deterministas y estocásticos."),
             p("Para empezar:"),
             tags$ul(
               tags$li("Seleccione los parámetros del modelo en la sección de Modelaje."),
               tags$li("Elija entre un modelo determinista o estocástico."),
               tags$li("Presione el botón 'Simular' para visualizar los resultados.")
             ),
             p("Navegue por las pestañas para explorar los modelos o descargar resultados.")
           )
  ),
  
  # Segunda pestaña: Modelaje
  tabPanel("Modelaje",
           sidebarLayout(
             sidebarPanel(
               h4("Parámetros del Modelo"),
               numericInput("mu", "Tasa de natalidad/mortalidad (mu)", value = 0.05, min = 0, max = 1, step = 0.01),
               numericInput("beta", "Tasa de transmisión (beta)", value = 0.5, min = 0, max = 2, step = 0.01),
               numericInput("gamma", "Tasa de recuperación (gamma)", value = 0.1, min = 0, max = 1, step = 0.01),
               numericInput("N", "Tamaño de la población", value = 100, min = 1, max = 10000, step = 1),
               radioButtons("model_type", "Selecciona el Modelo:",
                            choices = list("Determinista" = "determinista", 
                                           "Estocástico" = "estocastico", "Ambos" = "Ambos"),
                            selected = "determinista"),
               actionButton("simular", "Simular")
             ),
             
             mainPanel(
               tabsetPanel(
                 id = "model_tabs",  # ID para controlar las pestañas internas
                 tabPanel("Modelo SIR", plotOutput("sirPlot")),
                 tabPanel("Descarga", 
                          h4("Descarga de Resultados"),
                          downloadButton("downloadData", "Descargar Resultados")
                 )
               )
             )
           )
  )
)



server <- function(input, output) {
  observeEvent(input$simular, {
    mu <- input$mu
    beta <- input$beta
    gamma <- input$gamma
    N <- input$N
    
    if (input$model_type == "determinista" || input$model_type == "Ambos") {
      # Deterministic Model using RKF45
      S0 <- N - 1
      I0 <- 1
      R0 <- 0
      max_time <- 365
      
      sir_model <- function(t, state, parameters) {
        with(as.list(c(state, parameters)), {
          dS <- mu * N - (beta * I * S / N) - mu * S
          dI <- (beta * I * S / N) - gamma * I - mu * I
          dR <- gamma * I - mu * R
          return(list(c(dS, dI, dR)))
        })
      }
      
      parameters <- c(mu = mu, beta = beta, gamma = gamma)
      state <- c(S = S0, I = I0, R = R0)
      time <- seq(0, max_time, by = 1)
      
      sir_deterministic <- ode(y = state, times = time, func = sir_model, parms = parameters)
      sir_deterministic <- as.data.frame(sir_deterministic)
    }
    
    if (input$model_type == "estocastico" || input$model_type == "Ambos") {
      # Stochastic Model using Random Variables
      t <- 1:365
      S <- rep(0, 365)
      I <- rep(0, 365)
      R <- rep(0, 365)
      pR <- 0.01
      pC <- 0.5
      
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
      }
      
      sir_stochastic <- data.frame(Dia = t, S = S, I = I, R = R)
    }
    
    output$sirPlot <- renderPlot({
      if (input$model_type == "determinista") {
        # Plot Deterministic Model
        plot(sir_deterministic$time, sir_deterministic$S, type = "l", col = "blue", ylim = c(0, N), xlab = "Tiempo", ylab = "Población", lwd = 2, main = "Modelo Determinista SIR")
        lines(sir_deterministic$time, sir_deterministic$I, col = "red", lwd = 2)
        lines(sir_deterministic$time, sir_deterministic$R, col = "green", lwd = 2)
        legend("right", legend = c("Susceptible", "Infectado", "Recuperado"), col = c("blue", "red", "green"), lty = 1, lwd = 2)
      } else if (input$model_type == "estocastico") {
        # Plot Stochastic Model
        plot(sir_stochastic$Dia, sir_stochastic$S, type = "l", col = "blue", ylim = c(0, N), xlab = "Día", ylab = "Población", main = "Modelo Estocástico SIR")
        lines(sir_stochastic$Dia, sir_stochastic$I, col = "red")
        lines(sir_stochastic$Dia, sir_stochastic$R, col = "green")
        legend("right", legend = c("Susceptible", "Infectado", "Recuperado"), col = c("blue", "red", "green"), lty = 1)
      } else if (input$model_type == "Ambos") {
        # Plot Both Models Side-by-Side
        library(gridExtra)
        p1 <- ggplot(sir_deterministic, aes(x = time)) +
          geom_line(aes(y = S, color = "Susceptibles")) +
          geom_line(aes(y = I, color = "Infectados")) +
          geom_line(aes(y = R, color = "Recuperados")) +
          labs(title = "Modelo Determinista", x = "Tiempo", y = "Población") +
          scale_color_manual(values = c("Susceptibles" = "blue", "Infectados" = "red", "Recuperados" = "green")) +
          theme_minimal()
        
        p2 <- ggplot(sir_stochastic, aes(x = Dia)) +
          geom_line(aes(y = S, color = "Susceptibles")) +
          geom_line(aes(y = I, color = "Infectados")) +
          geom_line(aes(y = R, color = "Recuperados")) +
          labs(title = "Modelo Estocástico", x = "Día", y = "Población") +
          scale_color_manual(values = c("Susceptibles" = "blue", "Infectados" = "red", "Recuperados" = "green")) +
          theme_minimal()
        
        grid.arrange(p1, p2, ncol = 2)
      }
    })
  })
}


# Run the app
shinyApp(ui = ui, server = server)
