library(shiny)
library(deSolve)
library(ggplot2)
library(bslib)
library(gridExtra)

# INTERFAZ
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
               tags$li("Seleccione el tipo de modelo que desea simular."),
               tags$li("Ingrese los parámetros específicos para el modelo seleccionado."),
               tags$li("Presione el botón 'Simular' para visualizar los resultados.")
             ),
             p("Navegue por las pestañas para explorar los modelos o descargar resultados.")
           )
  ),
  
  # Segunda pestaña: Modelaje
  tabPanel("Modelaje",
           sidebarLayout(
             sidebarPanel(
               h4("Selecciona el Modelo"),
               radioButtons("model_type", "Selecciona el Modelo:",
                            choices = list("Determinista" = "determinista", 
                                           "Estocástico" = "estocastico", 
                                           "Ambos" = "ambos"),
                            selected = "determinista"),
               
               # Parámetros visibles solo para los modelos individuales
               conditionalPanel(
                 condition = "input.model_type == 'determinista'",
                 h4("Parámetros para el Modelo Determinista"),
                 numericInput("mu_d", "Tasa de natalidad/mortalidad (mu)", value = 0.05, min = 0, max = 1, step = 0.01),
                 numericInput("beta_d", "Tasa de transmisión (beta)", value = 0.5, min = 0, max = 2, step = 0.01),
                 numericInput("gamma_d", "Tasa de recuperación (gamma)", value = 0.1, min = 0, max = 1, step = 0.01),
                 numericInput("N_d", "Tamaño de la población", value = 100, min = 1, max = 10000, step = 1)
               ),
               
               conditionalPanel(
                 condition = "input.model_type == 'estocastico'",
                 h4("Parámetros para el Modelo Estocástico"),
                 numericInput("mu_s", "Tasa de natalidad/mortalidad (mu)", value = 0.05, min = 0, max = 1, step = 0.01),
                 numericInput("pR_s", "Probabilidad de recuperación (pR)", value = 0.01, min = 0, max = 1, step = 0.01),
                 numericInput("pC_s", "Probabilidad de usar cubrebocas (pC)", value = 0.5, min = 0, max = 1, step = 0.01),
                 numericInput("N_s", "Tamaño de la población", value = 100, min = 1, max = 10000, step = 1)
               ),
               
               actionButton("simular", "Simular")
             ),
             
             mainPanel(
               tabsetPanel(
                 id = "model_tabs",
                 tabPanel("Modelo SIR", plotOutput("sirPlot")),
                 tabPanel("Descarga", 
                          fluidRow(
                            column(6,
                                   h4("Tabla Modelo Determinista"),
                                   conditionalPanel(
                                     condition = "input.model_type == 'determinista' || input.model_type == 'ambos'",
                                     tableOutput("resultsTableDet"),
                                     downloadButton("downloadDetTable", "Descargar Determinista")
                                   )
                            ),
                            column(6,
                                   h4("Tabla Modelo Estocástico"),
                                   conditionalPanel(
                                     condition = "input.model_type == 'estocastico' || input.model_type == 'ambos'",
                                     tableOutput("resultsTableSto"),
                                     downloadButton("downloadStoTable", "Descargar Estocástico")
                                   )
                            )
                          )
                 )
               )
             )
           )
  )
)

# SERVER
server <- function(input, output) {
  results_deterministic <- reactiveVal(NULL)
  results_stochastic <- reactiveVal(NULL)
  
  observeEvent(input$simular, {
    if (input$model_type == "determinista" || input$model_type == "ambos") {
      mu <- input$mu_d
      beta <- input$beta_d
      gamma <- input$gamma_d
      N <- input$N_d
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
      results_deterministic(as.data.frame(sir_deterministic))
    }
    
    if (input$model_type == "estocastico" || input$model_type == "ambos") {
      mu <- input$mu_s
      pR <- input$pR_s
      pC <- input$pC_s
      N <- input$N_s
      t <- 1:365
      S <- rep(0, 365)
      I <- rep(0, 365)
      R <- rep(0, 365)
      
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
      
      results_stochastic(data.frame(Dia = t, S = S, I = I, R = R))
    }
  })
  
  # Gráficas
  output$sirPlot <- renderPlot({
    if (input$model_type == "determinista" && !is.null(results_deterministic())) {
      df <- results_deterministic()
      plot(df$time, df$S, type = "l", col = "blue", ylim = c(0, max(df$S, df$I, df$R)), xlab = "Tiempo", ylab = "Población", main = "Modelo Determinista SIR")
      lines(df$time, df$I, col = "red")
      lines(df$time, df$R, col = "green")
      legend("topright", legend = c("Susceptibles", "Infectados", "Recuperados"), col = c("blue", "red", "green"), lty = 1)
    } else if (input$model_type == "estocastico" && !is.null(results_stochastic())) {
      df <- results_stochastic()
      plot(df$Dia, df$S, type = "l", col = "blue", ylim = c(0, max(df$S, df$I, df$R)), xlab = "Día", ylab = "Población", main = "Modelo Estocástico SIR")
      lines(df$Dia, df$I, col = "red")
      lines(df$Dia, df$R, col = "green")
      legend("topright", legend = c("Susceptibles", "Infectados", "Recuperados"), col = c("blue", "red", "green"), lty = 1)
    } else if (input$model_type == "ambos" && !is.null(results_deterministic()) && !is.null(results_stochastic())) {
      df_det <- results_deterministic()
      df_sto <- results_stochastic()
      
      p1 <- ggplot(df_det, aes(x = time)) +
        geom_line(aes(y = S, color = "Susceptibles")) +
        geom_line(aes(y = I, color = "Infectados")) +
        geom_line(aes(y = R, color = "Recuperados")) +
        labs(title = "Modelo Determinista", x = "Tiempo", y = "Población") +
        theme_minimal()
      
      p2 <- ggplot(df_sto, aes(x = Dia)) +
        geom_line(aes(y = S, color = "Susceptibles")) +
        geom_line(aes(y = I, color = "Infectados")) +
        geom_line(aes(y = R, color = "Recuperados")) +
        labs(title = "Modelo Estocástico", x = "Día", y = "Población") +
        theme_minimal()
      
      grid.arrange(p1, p2, ncol = 2)
    }
  })
  
  # Tablas y descargas
  output$resultsTableDet <- renderTable({
    results_deterministic()
  })
  
  output$resultsTableSto <- renderTable({
    results_stochastic()
  })
  
  output$downloadDetTable <- downloadHandler(
    filename = "Resultados_Determinista.csv",
    content = function(file) {
      write.csv(results_deterministic(), file, row.names = FALSE)
    }
  )
  
  output$downloadStoTable <- downloadHandler(
    filename = "Resultados_Estocastico.csv",
    content = function(file) {
      write.csv(results_stochastic(), file, row.names = FALSE)
    }
  )
}

# Ejecuta la app
shinyApp(ui = ui, server = server)
