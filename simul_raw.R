library(shiny)
library(deSolve)
library(ggplot2)
library(bslib)
library(gridExtra)
library(dplyr)
library(tidyr)

# INTERFAZ
ui <- navbarPage(
  theme = bs_theme(version = 4, bootswatch = "flatly"),
  title = "Simulación de Modelos SIR",
  
  # Primera pestaña: Inicio
  tabPanel("Inicio",
           fluidPage(
             # Encabezado principal con estilo
             div(style = "text-align: center; margin-bottom: 10px; padding: 10px; background-color: #f8f9fa; border-radius: 10px;",
                 h2("📊 Simulación de Modelos SIR", style = "color: #2c3e50; font-weight: bold;")
             ),
             
             # Contenido principal
             div(style = "margin: 0 auto; max-width: 900px;",
                 # Primera sección: Descripción de la aplicación
                 div(style = "margin-bottom: 10px;",
                     h3("¿Qué es esta aplicación?", style = "color: #16a085; font-weight: bold;"),
                     p("Esta herramienta interactiva permite modelar y simular el comportamiento de epidemias utilizando el modelo SIR.", 
                       style = "font-size: 16px; line-height: 1.6; color: #2c3e50;"),
                     p("Incluye opciones para explorar modelos deterministas y estocásticos, facilitando el entendimiento de cómo las variables afectan la propagación de una enfermedad.", 
                       style = "font-size: 16px; line-height: 1.6; color: #2c3e50;")
                 ),
                 
                 # Segunda sección: Cómo empezar
                 div(style = "margin-bottom: 10px;",
                     h3("¿Cuáles fueron los parametros?", style = "color: #e74c3c; font-weight: bold;"),
                     tags$div(
                       p("Para el modelo determinista, utilizaremos las siguientes ecuaciones diferenciales:")
                     ),
                     withMathJax(
                       tags$div(
                         "$$\\frac{dS}{dt} = \\mu N - \\frac{\\beta I S}{N} - \\mu S \\quad \\text{(natalidad, infección, muerte)}$$",
                         "$$\\frac{dI}{dt} = \\frac{\\beta I S}{N} - \\gamma I - \\mu I \\quad \\text{(infección, recuperación, muerte)}$$",
                         "$$\\frac{dR}{dt} = \\gamma I - \\mu R \\quad \\text{(recuperación, muerte)}$$",
                         p("La información para los parametros mu y gamma se obtuvieron a partir de una simulación que también utilizaba el modelo SIR, en un estudio realizado por estudiantes de la Universidad de Sevilla.", style = "font-size: 16px; color: #34495e; line-height: 1.8;")
                       )
                     ),
                     tags$div(
                       p("Por otro lado, para el modelo estocastico, se utilizaron variables aleatorias para generar la simulacion. Para que la simulacion sea mas acertada decidimos involucrar otras variables que modifiquen el comportamiento de la simulacion conforme avanzan los dias."),
                       tags$ul(
                         tags$li("Tasa de infeccion generada por una variable aleatoria."),
                         tags$li("Probabilidad de que los infectados usen cubrebocas."),
                         tags$li("Probabilidad de recuperacion."),
                         tags$li("Probabilidad de que la poblacion este vacunada.")
                       ),
                       p("Las condiciones iniciales para cada uno de los modelos ceran una poblacion que se puede seleccionar, y un individuo que este infectado con el viruz. Ambas simulaciones avanzaran al paso de un dia, donde el tiempo maximo sera de un año o 365 días.")
                     )
                 ),
                 
                 # Tercera sección: Beneficios
                 div(style = "margin-bottom: 10px; padding: 20px; background-color: #f9f9f9; border-radius: 10px; box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);",
                     h3("¿Qué puedes hacer aquí?", style = "color: #3498db; font-weight: bold;"),
                     tags$ul(
                       tags$li("Generar gráficos interactivos para analizar la propagación de enfermedades."),
                       tags$li("Comparar modelos deterministas y estocásticos."),
                       tags$li("Descargar resultados en formato CSV para análisis adicionales."),
                       tags$li("Comparar el resultado de los errores de ambos modelos.")
                     )
                 )
             ),
             
             # Pie de página con un estilo moderno
             div(style = "text-align: center; margin-top: 10px; padding: 10px; background-color: #e9ecef; border-top: 2px solid #dcdcdc;",
                 p("Desarrollado con ❤️ por PURO PADRE | 2024", style = "font-size: 14px; color: #7f8c8d;")
             )
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
  
  # Reactive to read the uploaded CSV
  real_data <- reactive({
    req(input$real_data_file)
    read.csv(input$real_data_file$datapath) %>%
      mutate(fecha = as.Date(fecha, format = "%Y-%m-%d"))
  })
  
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
      pR <- 0.01
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
        face_mask <- rbinom(1, I[day - 1], pC)
        pI <- runif(1)
        
        if (I[day - 1] > 0) {
          effective_pI <- pI * (1 - face_mask / I[day - 1])
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
    
    if (input$model_type == "ambos") {
      df_det <- results_deterministic()
      df_sto <- results_stochastic()
      mse <- mean((df_det$I - df_sto$I)^2)
      mse_value(mse)
    }
  })
  
  # Gráfico principal
  output$sirPlot <- renderPlot({
    if (input$model_type == "determinista" && !is.null(results_deterministic())) {
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
      ggplot(df, aes(x = Dia)) +
        geom_line(aes(y = S, color = "Susceptibles")) +
        geom_line(aes(y = I, color = "Infectados")) +
        geom_line(aes(y = R, color = "Recuperados")) +
        labs(title = "Modelo Estocástico SIR", x = "Día", y = "Población") +
        scale_color_manual(values = c("Susceptibles" = "blue", "Infectados" = "red", "Recuperados" = "green")) +
        theme_minimal()
    } else if (input$model_type == "ambos" && !is.null(results_deterministic()) && !is.null(results_stochastic())) {
      df_det <- results_deterministic()
      df_det$model <- "Determinista"
      df_sto <- results_stochastic()
      df_sto$model <- "Estocástico"
      
      # Asegurarse de que los nombres de columnas sean consistentes
      colnames(df_sto)[colnames(df_sto) == "Dia"] <- "time"
      
      combined_data <- rbind(
        df_det %>% mutate(model = "Determinista"),
        df_sto %>% mutate(model = "Estocástico")
      )
      
      ggplot(combined_data, aes(x = time)) +
        geom_line(aes(y = S, color = interaction(model, "Susceptibles")), size = 1) +
        geom_line(aes(y = I, color = interaction(model, "Infectados")), size = 1) +
        geom_line(aes(y = R, color = interaction(model, "Recuperados")), size = 1) +
        labs(title = "Comparación de Modelos Determinista y Estocástico", 
             x = "Tiempo", 
             y = "Población") +
        scale_color_manual(
          values = c(
            "Determinista.Susceptibles" = "blue",
            "Determinista.Infectados" = "red",
            "Determinista.Recuperados" = "green",
            "Estocástico.Susceptibles" = "darkblue",
            "Estocástico.Infectados" = "darkred",
            "Estocástico.Recuperados" = "darkgreen"
          ),
          labels = c(
            "Determinista.Susceptibles" = "Determinista - Susceptibles",
            "Determinista.Infectados" = "Determinista - Infectados",
            "Determinista.Recuperados" = "Determinista - Recuperados",
            "Estocástico.Susceptibles" = "Estocástico - Susceptibles",
            "Estocástico.Infectados" = "Estocástico - Infectados",
            "Estocástico.Recuperados" = "Estocástico - Recuperados"
          )
        ) +
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