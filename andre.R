library(shiny)
library(deSolve)
library(ggplot2)
library(bslib)
library(dplyr)

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
             p("Navegue por las pestañas para explorar los modelos o realizar análisis comparativos.")
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
               
               # Parámetros para los modelos
               numericInput("mu_d", "Tasa de natalidad/mortalidad (mu)", value = 0.05, min = 0, max = 1, step = 0.01),
               numericInput("beta_d", "Tasa de transmisión (beta)", value = 0.5, min = 0, max = 2, step = 0.01),
               numericInput("gamma_d", "Tasa de recuperación (gamma)", value = 0.1, min = 0, max = 1, step = 0.01),
               numericInput("N_d", "Tamaño de la población", value = 100, min = 1, max = 10000, step = 1),
               
               actionButton("simular", "Simular")
             ),
             mainPanel(
               tabsetPanel(
                 id = "model_tabs",
                 tabPanel("Modelo SIR", plotOutput("sirPlot")),
                 tabPanel("Descarga", 
                          fluidRow(
                            column(6,
                                   h4("Modelo Determinista"),
                                   conditionalPanel(
                                     condition = "input.model_type == 'determinista' || input.model_type == 'ambos'",
                                     downloadButton("downloadDetTable", "Descargar Determinista"),
                                     tableOutput("resultsTableDet")
                                   )
                            ),
                            column(6,
                                   h4("Modelo Estocástico"),
                                   conditionalPanel(
                                     condition = "input.model_type == 'estocastico' || input.model_type == 'ambos'",
                                     downloadButton("downloadStoTable", "Descargar Estocástico"),
                                     tableOutput("resultsTableSto")
                                   )
                            )
                          )
                 ),
                 tabPanel("Resultados", 
                          h4("Mean Squared Error (MSE)"),
                          verbatimTextOutput("mseOutput"),
                          fileInput("real_data_file", "Sube el archivo CSV con datos reales:",
                                    accept = c(".csv")),
                          plotOutput("errorComparisonPlot") # Nueva gráfica de comparación de errores
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
  mse_value <- reactiveVal(NULL)
  #
  #
  #
  # Reactive to read the uploaded CSV
  real_data <- reactive({
    req(input$real_data_file)
    read.csv(input$real_data_file$datapath) %>%
      mutate(fecha = as.Date(fecha, format = "%Y-%m-%d"))
  })
  
  observeEvent(input$simular, {
    mu <- input$mu_d
    beta <- input$beta_d
    gamma <- input$gamma_d
    N <- input$N_d
    S0 <- N - 1
    I0 <- 1
    R0 <- 0
    max_time <- 365
    
    # Modelo Determinista
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
    
    sir_deterministic <- as.data.frame(ode(y = state, times = time, func = sir_model, parms = parameters))
    results_deterministic(sir_deterministic)
    
    if (input$model_type == "ambos" || input$model_type == "estocastico") {
      S <- rep(0, max_time + 1)
      I <- rep(0, max_time + 1)
      R <- rep(0, max_time + 1)
      
      S[1] <- S0
      I[1] <- I0
      R[1] <- R0
      
      for (t in 2:(max_time + 1)) {
        new_infected <- rbinom(1, S[t - 1], beta * I[t - 1] / N)
        new_recovered <- rbinom(1, I[t - 1], gamma)
        
        S[t] <- S[t - 1] - new_infected
        I[t] <- I[t - 1] + new_infected - new_recovered
        R[t] <- R[t - 1] + new_recovered
      }
      
      results_stochastic(data.frame(time = time, S = S, I = I, R = R))
      
      if (input$model_type == "ambos") {
        mse <- mean((sir_deterministic$I - I)^2)
        mse_value(mse)
      }
    }
  })
  
  # Gráfico principal
  output$sirPlot <- renderPlot({
    if (input$model_type == "ambos" && !is.null(results_deterministic()) && !is.null(results_stochastic())) {
      df_det <- results_deterministic() %>% mutate(model = "Determinista")
      df_sto <- results_stochastic() %>% mutate(model = "Estocástico")
      
      combined_data <- bind_rows(
        df_det %>% select(time, S, I, R, model),
        df_sto %>% select(time, S, I, R, model)
      )
      
      ggplot(combined_data, aes(x = time)) +
        geom_line(aes(y = S, color = interaction(model, "Susceptibles"))) +
        geom_line(aes(y = I, color = interaction(model, "Infectados"))) +
        geom_line(aes(y = R, color = interaction(model, "Recuperados"))) +
        labs(title = "Comparación de Modelos Determinista y Estocástico", 
             x = "Tiempo", y = "Población") +
        scale_color_manual(values = c(
          "Determinista.Susceptibles" = "blue",
          "Determinista.Infectados" = "red",
          "Determinista.Recuperados" = "green",
          "Estocástico.Susceptibles" = "darkblue",
          "Estocástico.Infectados" = "darkred",
          "Estocástico.Recuperados" = "darkgreen"
        )) +
        theme_minimal()
    } else if (input$model_type == "determinista" && !is.null(results_deterministic())) {
      df <- results_deterministic()
      ggplot(df, aes(x = time)) +
        geom_line(aes(y = S, color = "Susceptibles")) +
        geom_line(aes(y = I, color = "Infectados")) +
        geom_line(aes(y = R, color = "Recuperados")) +
        labs(title = "Modelo Determinista SIR", x = "Tiempo", y = "Población") +
        scale_color_manual(values = c("Susceptibles" = "blue", "Infectados" = "red", "Recuperados" = "green")) +
        theme_minimal()
    } else if (input$model_type == "estocastico" && !is.null(results_stochastic())) {
      df <- results_stochastic()
      ggplot(df, aes(x = time)) +
        geom_line(aes(y = S, color = "Susceptibles")) +
        geom_line(aes(y = I, color = "Infectados")) +
        geom_line(aes(y = R, color = "Recuperados")) +
        labs(title = "Modelo Estocástico SIR", x = "Día", y = "Población") +
        scale_color_manual(values = c("Susceptibles" = "blue", "Infectados" = "red", "Recuperados" = "green")) +
        theme_minimal()
    }
  })
  
  # Nueva gráfica de comparación de errores
  output$errorComparisonPlot <- renderPlot({
    if (input$model_type == "ambos" &&
        !is.null(results_deterministic()) &&
        !is.null(results_stochastic()) &&
        !is.null(real_data())) {
      
      real <- real_data()
      df_det <- results_deterministic() %>% mutate(fecha = seq.Date(from = min(real$fecha), by = "day", length.out = n()))
      df_sto <- results_stochastic() %>% mutate(fecha = seq.Date(from = min(real$fecha), by = "day", length.out = n()))
      
      comparison <- real %>%
        left_join(df_det, by = "fecha") %>%
        left_join(df_sto, by = "fecha", suffix = c("_det", "_sto")) %>%
        mutate(
          error_det = abs(casos - I_det),
          error_sto = abs(casos - I_sto)
        ) %>%
        pivot_longer(cols = c(error_det, error_sto), names_to = "modelo", values_to = "error")
      
      ggplot(comparison, aes(x = fecha, y = error, color = modelo)) +
        geom_line() +
        labs(title = "Comparación de Errores Absolutos", x = "Fecha", y = "Error Absoluto") +
        scale_color_manual(values = c("error_det" = "red", "error_sto" = "blue")) +
        theme_minimal()
    }
  })
  
  
  # Mostrar el MSE
  output$mseOutput <- renderText({
    if (!is.null(mse_value())) {
      paste("El error cuadrado (MSE) entre el modelo determinista y estocástico es:", round(mse_value(), 4))
    } else {
      "Seleccione 'Ambos' para calcular el MSE."
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


shinyApp(ui = ui, server = server)
