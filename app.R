library(shiny)
library(deSolve)
library(ggplot2)

# ec diferencial para la eliminación de cafeína
f_cafeina <- function(t, X, a) {
  return(-0.15 * X + a)
}

# Método de Euler
metodo_euler <- function(t0, X0, h, n, a) {
  t <- numeric(n + 1)
  X <- numeric(n + 1)
  t[1] <- t0
  X[1] <- X0
  
  for (i in 1:n) {
    X[i + 1] <- X[i] + h * f_cafeina(t[i], X[i], a)
    t[i + 1] <- t[i] + h
  }
  
  return(data.frame(t = t, X = X))
}

# Método RK4
runge_kutta4 <- function(t0, X0, h, n, a) {
  t <- numeric(n + 1)
  X <- numeric(n + 1)
  t[1] <- t0
  X[1] <- X0
  
  for (i in 1:n) {
    k1 <- h * f_cafeina(t[i], X[i], a)
    k2 <- h * f_cafeina(t[i] + h / 2, X[i] + k1 / 2, a)
    k3 <- h * f_cafeina(t[i] + h / 2, X[i] + k2 / 2, a)
    k4 <- h * f_cafeina(t[i] + h, X[i] + k3, a)
    
    X[i + 1] <- X[i] + (k1 + 2 * k2 + 2 * k3 + k4) / 6
    t[i + 1] <- t[i] + h
  }
  
  return(data.frame(t = t, X = X))
}

# El rk45
rk45 <- function(t0, X0, h, n, a) {
  t <- numeric(n + 1)
  y <- numeric(n + 1)
  t[1] <- t0
  y[1] <- X0
  
  for (i in 1:n) {
    k1 <- h * f_cafeina(t[i], y[i], a)
    k2 <- h * f_cafeina(t[i] + h / 4, y[i] + k1 / 4, a)
    k3 <- h * f_cafeina(t[i] + 3 * h / 8, y[i] + 3 * k1 / 32 + 9 * k2 / 32, a)
    k4 <- h * f_cafeina(t[i] + 12 * h / 13, y[i] + 1932 * k1 / 2197 - 7200 * k2 / 2197 + 7296 * k3 / 2197, a)
    k5 <- h * f_cafeina(t[i] + h, y[i] + 439 * k1 / 216 - 8 * k2 + 3680 * k3 / 513 - 845 * k4 / 4104, a)
    k6 <- h * f_cafeina(t[i] + h / 2, y[i] - 8 * k1 / 27 + 2 * k2 - 3544 * k3 / 2565 + 1859 * k4 / 4104 - 11 * k5 / 40, a)
    
    y[i + 1] <- y[i] + (16 * k1 / 135 + 6656 * k3 / 12825 + 28561 * k4 / 56430 - 9 * k5 / 50 + 2 * k6 / 55)
    t[i + 1] <- t[i] + h
  }
  
  return(data.frame(t = t,X=y))
}

# El rk45 fehlberg, ajuste de paso
rkf45 <- function(t0, X0, h, n, a) {
  t <- numeric(n + 1)
  y <- numeric(n + 1)
  t[1] <- t0
  y[1] <- X0
  tol <- 1e-5
  
  i <- 1
  while (i <= n) {
    k1 <- h * f_cafeina(t[i], y[i], a)
    k2 <- h * f_cafeina(t[i] + h / 4, y[i] + k1 / 4, a)
    k3 <- h * f_cafeina(t[i] + 3 * h / 8, y[i] + 3 * k1 / 32 + 9 * k2 / 32, a)
    k4 <- h * f_cafeina(t[i] + 12 * h / 13, y[i] + 1932 * k1 / 2197 - 7200 * k2 / 2197 + 7296 * k3 / 2197, a)
    k5 <- h * f_cafeina(t[i] + h, y[i] + 439 * k1 / 216 - 8 * k2 + 3680 * k3 / 513 - 845 * k4 / 4104, a)
    k6 <- h * f_cafeina(t[i] + h / 2, y[i] - 8 * k1 / 27 + 2 * k2 - 3544 * k3 / 2565 + 1859 * k4 / 4104 - 11 * k5 / 40, a)

    
    y5 = y[i] + (16 * k1 / 135 + 6656 * k3 / 12825 + 28561 * k4 / 56430 - 9 * k5 / 50 + 2 * k6 / 55)
    y4 = y[i] + (25 * k1 / 216 + 1408 * k3 / 2565 + 2197 * k4 / 4104 - k5 / 5)
    error = abs(y5 - y4)
    
    if (error < tol) {
      t[i + 1] = t[i] + h
      y[i + 1] = y5
      i <- i + 1
      if (error < tol / 2) h <- h * 2  
    } else {
      h <- h / 2  
    }
    if (i == n) break  
  }
  
  return(data.frame(t = t, X = y))
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
      actionButton("solve_rkf45", "Resolver con RKF45", icon = icon("chart-line"), class = "btn-warning"),
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
  a <- 50
  # Resolver con Euler
  observeEvent(input$solve_euler, {
    h <- input$h
    n <- input$n
    t0 <- 0
    X0 <- 90  # Dosis inicial de cafeína (mg)
    
    res <- metodo_euler(t0, X0, h, n, a)
    
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
    
    res <- runge_kutta4(t0, X0, h, n, a)
    
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
    t0 <- 0
    X0 <- 90
    # t <- seq(0, n * h, by = h)
    init <- c(X = 90)  # Dosis inicial de cafeína (mg)
    
    res <- rk45(t0, X0, h, n, a)
    
    output$plot_result <- renderPlot({
      # res_df <- as.data.frame(res)
      ggplot(res, aes(x = t, y = X)) +
        geom_line(color = "orange", size = 1.2) +
        labs(title = "Método de RK45", x = "Tiempo (hr)", y = "Cafeína (mg)") +
        theme_minimal(base_size = 15)
    })
  })
  
  # Resolver con RK45 fehlberg
  observeEvent(input$solve_rkf45, {
    h <- input$h
    n <- input$n
    t0 <- 0
    X0 <- 90
    
    res <- rkf45(t0, X0, h, n, a)
    
    output$plot_result <- renderPlot({
      ggplot(res, aes(x = t, y = X)) +
        geom_line(color = "green", size = 1.2) +
        labs(title = "Método de RKF45", x = "Tiempo (hr)", y = "Cafeína (mg)") +
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
    res_euler <- metodo_euler(t0, X0, h, n, a)
    res_euler$Metodo <- "Euler"
    
    # Resultados de RK4
    res_rk4 <- runge_kutta4(t0, X0, h, n, a)
    res_rk4$Metodo <- "RK4"
    
    # Resultados de RK45
    t <- seq(0, n * h, by = h)
    res_rk45 <- rk45(t0, X0, h, n, a)
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