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
  title = "SARS-COV2 visto desde SIR",
  
  # Primera pestaña: Inicio
  tabPanel("Inicio",
           fluidPage(
             div(style = "text-align: center; margin-bottom: 10px; padding: 10px; background-color: #f8f9fa; border-radius: 10px;",
                 h2("📊 Simulación de Modelos SIR - Caso de estudio: SARS-COV2", style = "color: #2c3e50; font-weight: bold;")
             ),
             # Contenido de la pestaña Inicio...
             # Contenido principal
             div(style = "margin: 0 auto; max-width: 900px;",
                 # Primera sección: Descripción de la aplicación
                 div(style = "margin-bottom: 10px;",
                     h3("¿Qué es esta aplicación?", style = "color: #16a085; font-weight: bold;"),
                     p("Esta herramienta interactiva permite modelar y simular el comportamiento de epidemias utilizando el modelo SIR.", 
                       style = "font-size: 16px; line-height: 1.6; color: #2c3e50;"),
                     p("Incluye opciones para explorar modelos deterministas y estocásticos, facilitando el entendimiento de cómo las variables afectan la propagación de una enfermedad.", 
                       style = "font-size: 16px; line-height: 1.6; color: #2c3e50;"),
                     p("En este caso, la enfermedad más relevante para la época es sin duda el COVID'19. La enfermedad por coronavirus (sARS-COV2) es una enfermedad infecciosa causada por el virus SARS-CoV-2. 
                       La mayoría de las personas infectadas por el virus experimentarán una enfermedad respiratoria de leve a moderada. Sin embargo, algunas enfermarán gravemente y requerirán atención médica (WHO). Por esto, es menestér generar aplicaciones enfocadas en el análisis profundo de las causas y como estas producen efectos en poblaciones humanas.",
                       style = "font-size: 16px; line-height: 1.6; font-style: italic; color: #7f8c8d;")
                 ),
                 
                 # Segunda sección: Cómo empezar
                 div(style = "margin-bottom: 10px;",
                     h3("¿Cuáles fueron los parámetros?", style = "color: #e74c3c; font-weight: bold;"),
                     tags$div(
                       p("Para el modelo determinista, utilizaremos las siguientes ecuaciones diferenciales:")
                     ),
                     withMathJax(
                       tags$div(
                         "$$\\frac{dS}{dt} = \\mu N - \\frac{\\beta I S}{N} - \\mu S \\quad \\text{(natalidad, infección, muerte)}$$",
                         "$$\\frac{dI}{dt} = \\frac{\\beta I S}{N} - \\gamma I - \\mu I \\quad \\text{(infección, recuperación, muerte)}$$",
                         "$$\\frac{dR}{dt} = \\gamma I - \\mu R \\quad \\text{(recuperación, muerte)}$$",
                         p("La información para los parámetros mu y gamma se obtuvieron a partir de una simulación que también utilizaba el modelo SIR, en un estudio realizado por estudiantes de la Universidad de Sevilla.", style = "font-size: 16px; color: #34495e; line-height: 1.8;")
                       )
                     ),
                     tags$div(
                       p("Por otro lado, para el modelo estocástico, se utilizaron variables aleatorias para generar la simulación. Para que la simulación sea más acertada decidimos involucrar otras variables que modifiquen el comportamiento de la simulación conforme avanzan los días."),
                       tags$ul(
                         tags$li("Tasa de infección generada por una variable aleatoria."),
                         tags$li("Probabilidad de que los infectados usen cubrebocas."),
                         tags$li("Probabilidad de recuperación."),
                         tags$li("Probabilidad de que la población esté vacunada.")
                       ),
                       p("Las condiciones iniciales para cada uno de los modelos serán una población que se puede seleccionar, y un individuo que esté infectado con el virus. Ambas simulaciones avanzarán al paso de un día, donde el tiempo máximo será de un año o 365 días.")
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
                 p("Desarrollado con ❤️ por PURO PADRE | 2024 | A01706832 | A01707339 | A01701234", style = "font-size: 14px; color: #7f8c8d;")
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
                          verbatimTextOutput("mseOutput"),
                          fileInput("real_data_file", "Sube el archivo CSV con datos reales:",
                                    accept = c(".csv")),
                          plotOutput("errorComparisonPlot")
                 ),
                 
                 # Pestaña del Asistente Virtual
                 tabPanel("Asistente Virtual",
                          fluidPage(
                            div(style = "text-align: center; margin-bottom: 10px; padding: 10px; background-color: #f8f9fa; border-radius: 10px;",
                                h2("Asistente Virtual", style = "color: #2c3e50; font-weight: bold;")
                            ),
                            div(style = "margin: 0 auto; max-width: 600px;",
                                textAreaInput("user_input", label = "Etoy aquí para complementar tu análisis, haz preguntas acerca de los valores de la ventana \"Resultados\".", placeholder = "¿Qué modelo se adecúa mejor los primeros 100 días?"),
                                actionButton("ask_button", "Preguntar"),
                                br(),
                                h4("Respuesta del Asistente:"),
                                verbatimTextOutput("assistant_response")
                            )
                          )
                 )
               )
             )
           )
  )
)


# SERVER
server <- function(input, output, session) {
  # Placeholder for the full response and the current typed response
  full_response <- reactiveVal("") # Full response from the API
  typed_response <- reactiveVal("") # Gradually typed response
  
  # Timer for typing effect
  typing_timer <- reactiveVal(NULL)
  
  # Handle the assistant response logic
  observeEvent(input$ask_button, {
    req(input$user_input)  # Ensure there is a question
    
    # Reset typed response and timer
    typed_response("")
    typing_timer(NULL)
    
    # Collect data for the context
    mse <- if (!is.null(mse_value())) paste("El MSE entre los modelos es:", round(mse_value(), 4)) else "No hay un MSE disponible."
    det_results <- results_deterministic()
    sto_results <- results_stochastic()
    real_data_df <- real_data()
    
    # Convert real data to text
    real_data_text <- if (!is.null(real_data_df)) {
      paste(capture.output(print(head(real_data_df, 5))), collapse = "\n")
    } else {
      "No se proporcionaron datos reales."
    }
    
    # Build the context
    context <- list(
      list(role = "system", content = "Eres un asistente virtual para un análisis epidemiológico basado en modelos SIR."),
      list(role = "user", content = paste(
        "Pregunta:", input$user_input, "\n",
        "Información:\n",
        mse, "\n",
        "Resultados del modelo determinista (primeros 5 registros):", 
        if (!is.null(det_results)) paste(capture.output(print(head(det_results, 5))), collapse = "\n") else "No disponible.", "\n",
        "Resultados del modelo estocástico (primeros 5 registros):", 
        if (!is.null(sto_results)) paste(capture.output(print(head(sto_results, 5))), collapse = "\n") else "No disponible.", "\n",
        "Datos reales del CSV (primeros 5 registros):", real_data_text
      ))
    )
    
    # Make the API call
    response <- tryCatch({
      httr::POST(
        url = "https://api.openai.com/v1/chat/completions",
        add_headers(Authorization = paste("Bearer", "YOURAPIKEY")), # Replace with your API key
        content_type_json(),
        body = toJSON(list(
          model = "gpt-4",
          messages = context
        ), auto_unbox = TRUE)
      )
    }, error = function(e) {
      full_response(paste("Error al conectar con el API de ChatGPT:", e$message))
      return(NULL)
    })
    
    # Handle response or errors
    if (is.null(response) || http_type(response) != "application/json") {
      full_response("Error al conectar con el API de ChatGPT. Revisa tu conexión o la configuración del API Key.")
    } else {
      result <- fromJSON(content(response, as = "text", encoding = "UTF-8"))
      if (!is.null(result$choices) && length(result$choices) > 0) {
        full_response(result$choices[[1]]$message$content)
      } else {
        full_response("No se recibió una respuesta válida del asistente.")
      }
    }
    
    # Start typing animation
    response_chars <- strsplit(full_response(), "")[[1]] # Split full response into characters
    counter <- reactiveVal(1) # Initialize counter
    
    typing_timer(invalidateLater(50, session)) # Start the timer (50 ms per character)
    
    observeEvent(typing_timer(), {
      if (counter() <= length(response_chars)) {
        # Append the next character to the typed response
        typed_response(paste0(typed_response(), response_chars[counter()]))
        counter(counter() + 1) # Move to the next character
        typing_timer(invalidateLater(50, session)) # Continue timer
      } else {
        # Stop the timer when the response is fully typed
        typing_timer(NULL)
      }
    })
  })
  
  # Render the animated response
  output$assistant_response <- renderText({
    typed_response()
  })
  
  # Reactive values for models
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
  # Nueva gráfica de comparación de errores con cálculo del área bajo la curva (AUC)
  output$errorComparisonPlot <- renderPlot({
    if (input$model_type == "ambos" &&
        !is.null(results_deterministic()) &&
        !is.null(results_stochastic()) &&
        !is.null(real_data())) {
      
      real <- real_data()
      df_det <- results_deterministic() %>% mutate(fecha = seq.Date(from = min(real$fecha), by = "day", length.out = n()))
      df_sto <- results_stochastic() %>% mutate(fecha = seq.Date(from = min(real$fecha), by = "day", length.out = n()))
      
      # Comparación de errores absolutos
      comparison <- real %>%
        left_join(df_det, by = "fecha") %>%
        left_join(df_sto, by = "fecha", suffix = c("_det", "_sto")) %>%
        mutate(
          error_det = abs(casos - I_det),
          error_sto = abs(casos - I_sto)
        )
      
      # Calcular AUC para cada modelo
      auc_det <- sum(diff(comparison$fecha) * (head(comparison$error_det, -1) + tail(comparison$error_det, -1)) / 2, na.rm = TRUE)
      auc_sto <- sum(diff(comparison$fecha) * (head(comparison$error_sto, -1) + tail(comparison$error_sto, -1)) / 2, na.rm = TRUE)
      
      # Pivotar datos para graficar
      comparison_long <- comparison %>%
        pivot_longer(cols = c(error_det, error_sto), names_to = "modelo", values_to = "error")
      
      # Graficar comparación de errores
      ggplot(comparison_long, aes(x = fecha, y = error, color = modelo)) +
        geom_line() +
        labs(
          title = paste(
            "Comparación de Errores Absolutos\n",
            sprintf("AUC Determinista: %.2f", auc_det),
            sprintf(" | AUC Estocástico: %.2f", auc_sto)
          ),
          x = "Fecha",
          y = "Error Absoluto"
        ) +
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