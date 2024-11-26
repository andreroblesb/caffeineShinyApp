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
               
               conditionalPanel(
                 condition = "input.model_type == 'determinista'",
                 h4("Parámetros para el Modelo Determinista"),
                 numericInput("mu_d", "Tasa de natalidad/mortalidad (mu)", value = 0.001, min = 0, max = 1, step = 0.01),
                 numericInput("beta_d", "Tasa de transmisión (beta)", value = 0.2, min = 0, max = 2, step = 0.01),
                 numericInput("gamma_d", "Tasa de recuperación (gamma)", value = 0.01, min = 0, max = 1, step = 0.01),
                 sliderInput("N_d", "Tamaño de la población", value = 100, min = 100, max = 10000000, step = 1),
                 numericInput("tiempo", "Tiempo de simulacion", value = 365, min = 100, max = 1000, step = 1)
               ),
               
               conditionalPanel(
                 condition = "input.model_type == 'estocastico'",
                 h4("Parámetros para el Modelo Estocástico"),
                 numericInput("Inf_I", "Numero de Infectados", value = 1, min = 1, max = 10000, step = 1),
                 numericInput("pC_s", "Probabilidad de usar cubrebocas (pC)", value = 0, min = 0, max = 1, step = 0.01),
                 sliderInput("N_s", "Tamaño de la población", value = 100, min = 100, max = 10000000, step = 1),
                 numericInput("time", "Tiempo de simulacion", value = 365, min = 100, max = 1000, step = 1),
                 numericInput("pV_s", "Probabilidad de Vacunacion", value = 0, min = 0, max = 1, step = 0.01)
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
                          plotOutput("errorComparisonPlot")
                 ),
                 
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
server <- function(input, output) {
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
  results_deterministic <- reactiveVal(NULL)
  results_stochastic <- reactiveVal(NULL)
  mse_value <- reactiveVal(NULL)
  
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
      max_time <- input$tiempo
      
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
      pV <- input$pV_s  # Probabilidad inicial de vacunación
      N <- input$N_s
      max_time <- input$time
      
      S <- rep(0, max_time + 1)
      I <- rep(0, max_time + 1)
      R <- rep(0, max_time + 1)
      
      I[1] <- input$Inf_I
      R[1] <- 0
      S[1] <- N - I[1] - R[1]
      
      for (t in 2:(max_time + 1)) {
        
        if (t > max_time * 1 / 3 && pV > 0) {
          pV <- min(pV + runif(1) / 10, 1)
        }
        if (t > max_time * 1 / 3) {
          vaccine_I <- ifelse(I[t - 1] > 0 && pV > 0, rbinom(1, I[t - 1], pV), 0)
          vaccine_S <- ifelse(S[t - 1] > 0 && pV > 0, rbinom(1, S[t - 1], pV), 0)
        } else {
          vaccine_I <- 0
          vaccine_S <- 0
        }
        
        face_mask <- ifelse(I[t - 1] > 0, rbinom(1, I[t - 1], pC), 0)
        
        pI <- runif(1)
        if (S[t - 1] > 0 && I[t - 1] > 0) {
          effective_pI <- pI * (1 - face_mask / max(I[t - 1], 1)) * 
            (1 - vaccine_S / max(S[t - 1], 1))
          effective_pI <- pmin(pmax(effective_pI, 0), 1)
        } else {
          effective_pI <- 0
        }
        
        new_infected <- ifelse(S[t - 1] > 0, rbinom(1, S[t - 1], effective_pI * I[t - 1] / N), 0)
        pR_effective <- ifelse(vaccine_I > 0, pR + 0.01, pR)
        new_recovered <- ifelse(I[t - 1] > 0, rbinom(1, I[t - 1], pR_effective), 0)
        
        S[t] <- max(S[t - 1] - new_infected, 0)
        I[t] <- max(I[t - 1] + new_infected - new_recovered, 0)
        R[t] <- max(R[t - 1] + new_recovered, 0)
      }
      
      results_stochastic(data.frame(time = 0:max_time, S = S, I = I, R = R))
    }
    
    
    if (input$model_type == "ambos") {
      
      df_det <- results_deterministic()
      df_sto <- results_stochastic()
      
      print("Data frame determinista:")
      print(head(df_det))
      print("Data frame estocástico:")
      print(head(df_sto))
      
      if (!is.null(df_det) && !is.null(df_sto) && nrow(df_det) == nrow(df_sto)) {
        # Calcular el error cuadrado medio (MSE)
        mse <- mean((df_det$I - df_sto$I)^2, na.rm = TRUE)
        mse_value(mse)  # Guardar el valor del MSE en la variable reactiva
      } else {
        print("Error: Los datos no son compatibles para calcular el MSE.")
        mse_value(NULL)  # Si los datos no coinciden, establecer MSE como NULL
      }
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
      ggplot(df, aes(x = time)) +
        geom_line(aes(y = S, color = "Susceptibles")) +
        geom_line(aes(y = I, color = "Infectados")) +
        geom_line(aes(y = R, color = "Recuperados")) +
        labs(title = "Modelo Estocástico SIR", x = "Día", y = "Población") +
        scale_color_manual(values = c("Susceptibles" = "blue", "Infectados" = "red", "Recuperados" = "green")) +
        theme_minimal()
    } else if (input$model_type == "ambos") {
      # Datos del modelo determinista y estocástico
      df_det <- results_deterministic()
      df_sto <- results_stochastic()
      
      if (!is.null(df_det) && !is.null(df_sto)) {
        # Preparar datos de los modelos
        df_det <- df_det %>% mutate(model = "Determinista")
        df_sto <- df_sto %>% mutate(model = "Estocástico")
        
        combined_data <- rbind(
          df_det %>% mutate(model = "Determinista"),
          df_sto %>% mutate(model = "Estocástico")
        )
        
        # Si no hay datos reales cargados, mostrar solo los modelos
        if (is.null(input$real_data_file)) {
          return(
            ggplot(combined_data, aes(x = time)) +
              geom_line(aes(y = S, color = interaction(model, "Susceptibles")), size = 1) +
              geom_line(aes(y = I, color = interaction(model, "Infectados")), size = 1) +
              geom_line(aes(y = R, color = interaction(model, "Recuperados")), size = 1) +
              scale_color_manual(
                values = c(
                  "Determinista.Susceptibles" = "blue",
                  "Determinista.Infectados" = "red",
                  "Determinista.Recuperados" = "green",
                  "Estocástico.Susceptibles" = "darkblue",
                  "Estocástico.Infectados" = "darkred",
                  "Estocástico.Recuperados" = "darkgreen"
                )
              ) +
              labs(
                title = "Comparación de Modelos Determinista y Estocástico", 
                x = "Tiempo", 
                y = "Población"
              ) +
              theme_minimal()
          )
        }
        
        # Incorporar datos reales cuando se suban
        real <- real_data()
        real <- real %>%
          mutate(model = "Datos Reales") %>%
          rename(I = casos) %>%
          mutate(time = as.numeric(difftime(fecha, min(fecha), units = "days")))
        
        combined_data <- bind_rows(
          combined_data,
          real %>% select(time, I, model)
        )
        
        # Graficar los tres conjuntos de datos
        return(
          ggplot(combined_data, aes(x = time, y = I, color = model, linetype = model)) +
            geom_line(size = 1) +
            scale_color_manual(values = c(
              "Determinista" = "red",
              "Estocástico" = "blue",
              "Datos Reales" = "purple"
            )) +
            scale_linetype_manual(values = c(
              "Determinista" = "solid",
              "Estocástico" = "solid",
              "Datos Reales" = "solid"
            )) +
            labs(
              title = "Comparación de Modelos y Datos Reales",
              x = "Tiempo (días)",
              y = "Número de Infectados",
              color = "Modelo",
              linetype = "Modelo"
            ) +
            theme_minimal()
        )
      }
    }
  })
  
  
  # Mostrar el MSE
  output$mseOutput <- renderText({
    if (!is.null(mse_value())) {
      paste("El error cuadrado medio (MSE) entre los modelos determinista y estocástico es:", round(mse_value(), 4))
    } else {
      "El MSE no pudo ser calculado. Asegúrese de que los modelos tengan datos comparables."
    }
  })
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