library(shiny)
library(deSolve)
library(ggplot2)

# holeeee

# Definir la ecuación diferencial para la eliminación de cafeína
f_cafeina <- function(t, X, a) {
  return(-0.15 * X + a)
}

# Método de Euler
metodo_euler <- function(t0, X0, h, n) {
  t <- numeric(n + 1)
  X <- numeric(n + 1)
  t[1] <- t0
  X[1] <- X0
  
  for (i in 1:n) {
    X[i + 1] <- X[i] + h * f_cafeina(t[i], X[i])
    t[i + 1] <- t[i] + h
  }
  
  return(data.frame(t = t, X = X))
}

# Método RK4
runge_kutta4 <- function(t0, X0, h, n) {
  t <- numeric(n + 1)
  X <- numeric(n + 1)
  t[1] <- t0
  X[1] <- X0
  
  for (i in 1:n) {
    k1 <- h * f_cafeina(t[i], X[i])
    k2 <- h * f_cafeina(t[i] + h / 2, X[i] + k1 / 2)
    k3 <- h * f_cafeina(t[i] + h / 2, X[i] + k2 / 2)
    k4 <- h * f_cafeina(t[i] + h, X[i] + k3)
    
    X[i + 1] <- X[i] + (k1 + 2 * k2 + 2 * k3 + k4) / 6
    t[i + 1] <- t[i] + h
  }
  
  return(data.frame(t = t, X = X))
}

# Ecuación diferencial para RK45 usando el paquete deSolve
ed_cafeina <- function(t, state, parms) {
  with(as.list(state), {
    dxdt <- -0.15 * X
    return(list(c(dxdt)))
  })
}

# UI (Interfaz de usuario)
ui <- fluidPage(
  theme = bslib::bs_theme(bootswatch = "journal"),
  titlePanel("☕ Eliminación de Cafeína: Comparación de Métodos Numéricos"),
  
  sidebarLayout(
    sidebarPanel(
      h3("Parámetros de Entrada"),
      numericInput("h", "Paso de tiempo (h):", 0.1, min = 0.01, step = 0.01),
      numericInput("n", "Número de pasos (n):", 30, min = 1),
      actionButton("solve_euler", "Resolver con Euler", icon = icon("chart-line"), class = "btn-primary"),
      actionButton("solve_rk4", "Resolver con RK4", icon = icon("chart-line"), class = "btn-danger"),
      actionButton("solve_rk45", "Resolver con RK45", icon = icon("chart-line"), class = "btn-success"),
      actionButton("solve_comparar", "Comparar Métodos", icon = icon("balance-scale"), class = "btn-info")
    ),
    
    mainPanel(
      plotOutput("plot_result", height = "600px"),
      h5("Este simulador permite analizar la eliminación de cafeína en el cuerpo humano usando distintos métodos numéricos. ¡Explora cómo evolucionan los niveles de cafeína con el tiempo!")
    )
  )
)

# Server (Lógica del servidor)
server <- function(input, output) {
  # Resolver con Euler
  observeEvent(input$solve_euler, {
    h <- input$h
    n <- input$n
    t0 <- 0
    X0 <- 90  # Dosis inicial de cafeína (mg)
    
    res <- metodo_euler(t0, X0, h, n)
    
    output$plot_result <- renderPlot({
      ggplot(res, aes(x = t, y = X)) +
        geom_line(color = "blue", size = 1.2) +
        geom_point(color = "blue", size = 2) +
        labs(title = "Método de Euler", x = "Tiempo (hr)", y = "Cafeína (mg)") +
        theme_minimal(base_size = 15)
    })
  })
  
  # Resolver con RK4
  observeEvent(input$solve_rk4, {
    h <- input$h
    n <- input$n
    t0 <- 0
    X0 <- 90
    
    res <- runge_kutta4(t0, X0, h, n)
    
    output$plot_result <- renderPlot({
      ggplot(res, aes(x = t, y = X)) +
        geom_line(color = "red", size = 1.2) +
        geom_point(color = "red", size = 2) +
        labs(title = "Método de RK4", x = "Tiempo (hr)", y = "Cafeína (mg)") +
        theme_minimal(base_size = 15)
    })
  })
  
  # Resolver con RK45
  observeEvent(input$solve_rk45, {
    h <- input$h
    n <- input$n
    t <- seq(0, n * h, by = h)
    init <- c(X = 90)  # Dosis inicial de cafeína (mg)
    
    res <- ode(y = init, times = t, func = ed_cafeina, parms = NULL, method = "ode45")
    
    output$plot_result <- renderPlot({
      res_df <- as.data.frame(res)
      ggplot(res_df, aes(x = time, y = X)) +
        geom_line(color = "green", size = 1.2) +
        labs(title = "Método de RK45", x = "Tiempo (hr)", y = "Cafeína (mg)") +
        theme_minimal(base_size = 15)
    })
  })
  
  # Comparar los tres métodos (Euler, RK4, RK45)
  observeEvent(input$solve_comparar, {
    h <- input$h
    n <- input$n
    t0 <- 0
    X0 <- 90  # Dosis inicial
    
    # Resultados de Euler
    res_euler <- metodo_euler(t0, X0, h, n)
    res_euler$Metodo <- "Euler"
    
    # Resultados de RK4
    res_rk4 <- runge_kutta4(t0, X0, h, n)
    res_rk4$Metodo <- "RK4"
    
    # Resultados de RK45
    t <- seq(0, n * h, by = h)
    res_rk45 <- as.data.frame(ode(y = c(X = 90), times = t, func = ed_cafeina, parms = NULL, method = "ode45"))
    res_rk45$Metodo <- "RK45"
    
    # Combinar resultados
    res_combined <- rbind(
      data.frame(t = res_euler$t, X = res_euler$X, Metodo = res_euler$Metodo),
      data.frame(t = res_rk4$t, X = res_rk4$X, Metodo = res_rk4$Metodo),
      data.frame(t = res_rk45$time, X = res_rk45$X, Metodo = res_rk45$Metodo)
    )
    
    output$plot_result <- renderPlot({
      ggplot(res_combined, aes(x = t, y = X, color = Metodo)) +
        geom_line(size = 1.2) +
        labs(title = "Comparación de Métodos Numéricos", x = "Tiempo (hr)", y = "Cafeína (mg)", color = "Método") +
        theme_minimal(base_size = 15) +
        scale_color_manual(values = c("Euler" = "blue", "RK4" = "red", "RK45" = "green"))
    })
  })
}

# Ejecutar la aplicación
shinyApp(ui = ui, server = server)