library(shiny)
library(deSolve)
library(ggplot2)
library(bslib)
library(dplyr)
library(gridExtra)  # For side-by-side plots

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
               h6("Configuración (Predisp. SARS-COV2)"),
               # Radio button to select model
               radioButtons("model_type", "Selecciona el Modelo:",
                            choices = list("Determinista" = "determinista", 
                                           "Estocástico" = "estocastico", 
                                           "Ambos" = "ambos"),
                            selected = "determinista"),
               
               # Conditional UI elements based on model selection
               conditionalPanel(
                 condition = "input.model_type == 'determinista' || input.model_type == 'ambos'",
                 h4("Parámetros para el Modelo Determinista"),
                 numericInput("mu_d", "Tasa de natalidad/mortalidad (mu)", value = 0.05, min = 0, max = 1, step = 0.01),
                 numericInput("beta_d", "Tasa de transmisión (beta)", value = 0.5, min = 0, max = 2, step = 0.01),
                 numericInput("gamma_d", "Tasa de recuperación (gamma)", value = 0.25, min = 0, max = 1, step = 0.01),
                 numericInput("N_d", "Tamaño de la población", value = 100, min = 1, max = 10000, step = 1)
               ),
               
               conditionalPanel(
                 condition = "input.model_type == 'estocastico' || input.model_type == 'ambos'",
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

# SERVER
server <- function(input, output) {
  observeEvent(input$simular, {
    if (input$model_type == "determinista" || input$model_type == "ambos") {
      # Deterministic Model
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
      sir_deterministic <- as.data.frame(sir_deterministic)
      sir_deterministic$model <- "Determinista"  # Add model label
    }
    
    if (input$model_type == "estocastico" || input$model_type == "ambos") {
      # Stochastic Model
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
      
      sir_stochastic <- data.frame(
        time = t,
        S = S,
        I = I,
        R = R,
        model = "Estocástico"  # Add model label
      )
    }
    
    output$sirPlot <- renderPlot({
      if (input$model_type == "determinista") {
        # Plot Deterministic Model
        ggplot(sir_deterministic, aes(x = time)) +
          geom_line(aes(y = S, color = "Susceptibles")) +
          geom_line(aes(y = I, color = "Infectados")) +
          geom_line(aes(y = R, color = "Recuperados")) +
          labs(title = "Modelo Determinista SIR", x = "Tiempo", y = "Población") +
          scale_color_manual(values = c("Susceptibles" = "blue", "Infectados" = "red", "Recuperados" = "green")) +
          theme_minimal()
      } else if (input$model_type == "estocastico") {
        # Plot Stochastic Model
        ggplot(sir_stochastic, aes(x = time)) +
          geom_line(aes(y = S, color = "Susceptibles")) +
          geom_line(aes(y = I, color = "Infectados")) +
          geom_line(aes(y = R, color = "Recuperados")) +
          labs(title = "Modelo Estocástico SIR", x = "Día", y = "Población") +
          scale_color_manual(values = c("Susceptibles" = "blue", "Infectados" = "red", "Recuperados" = "green")) +
          theme_minimal()
      } else if (input$model_type == "ambos") {
        # Combine both models into one data frame
        combined_data <- rbind(
          sir_deterministic %>% select(time, S, I, R, model),
          sir_stochastic %>% select(time, S, I, R, model)
        )
        
        # Plot both models together
        ggplot(combined_data, aes(x = time)) +
          geom_line(aes(y = S, color = interaction(model, "Susceptibles"))) +
          geom_line(aes(y = I, color = interaction(model, "Infectados"))) +
          geom_line(aes(y = R, color = interaction(model, "Recuperados"))) +
          labs(title = "Comparación de Modelos Determinista y Estocástico", x = "Tiempo", y = "Población") +
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
          theme_minimal()
      }
    })
  })
}

# Run the app
shinyApp(ui = ui, server = server)
