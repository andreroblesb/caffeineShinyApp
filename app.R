library(shiny)
library(deSolve)
library(ggplot2)
library(bslib)

# Definir la ecuación diferencial para la eliminación de cafeína
f_cafeina <- function(t, X) {
  return(-0.15 * X)
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

# Método RK45
rk45 <- function(f_cafeina, t0, y0, h, n) {
  t <- numeric(n + 1)
  y <- numeric(n + 1)
  t[1] <- t0
  y[1] <- y0
  
  for (i in 1:n) {
    k1 <- h * f_cafeina(t[i], y[i])
    k2 <- h * f_cafeina(t[i] + h / 4, y[i] + k1 / 4)
    k3 <- h * f_cafeina(t[i] + 3 * h / 8, y[i] + 3 * k1 / 32 + 9 * k2 / 32)
    k4 <- h * f_cafeina(t[i] + 12 * h / 13, y[i] + 1932 * k1 / 2197 - 7200 * k2 / 2197 + 7296 * k3 / 2197)
    k5 <- h * f_cafeina(t[i] + h, y[i] + 439 * k1 / 216 - 8 * k2 + 3680 * k3 / 513 - 845 * k4 / 4104)
    k6 <- h * f_cafeina(t[i] + h / 2, y[i] - 8 * k1 / 27 + 2 * k2 - 3544 * k3 / 2565 + 1859 * k4 / 4104 - 11 * k5 / 40)
    
    y[i + 1] <- y[i] + (16 * k1 / 135 + 6656 * k3 / 12825 + 28561 * k4 / 56430 - 9 * k5 / 50 + 2 * k6 / 55)
    t[i + 1] <- t[i] + h
  }
  
  return(data.frame(Tiempo = t, Cafeina = y))
}

ui <- navbarPage(
  title = "Eliminación de Cafeína",
  theme = bslib::bs_theme(bootswatch = "journal"),

  tabPanel("Yay",
           p("Cosas pdres")
    
  ),
  
  tabPanel("Simulador",
           sidebarLayout(
             sidebarPanel(
               h3("Parámetros de Entrada"),
               sliderInput("h", "Paso de tiempo (h):", 0.1, min = 0.01, max = 1, step = 0.01),
               sliderInput("n", "Número de pasos (n):", 1, min = 1, max = 200),
               actionButton("solve_euler", "Resolver con Euler", icon = icon("chart-line"), class = "btn-primary", width = '215px', style = 'margin-bottom: 10px'),
               actionButton("solve_rk4", "Resolver con RK4", icon = icon("chart-line"), class = "btn-danger", width = '215px', style = 'margin-bottom: 10px'),
               actionButton("solve_rk45", "Resolver con RK45", icon = icon("chart-line"), class = "btn-success", width = '215px', style = 'margin-bottom: 10px'),
               actionButton("solve_comparar", "Comparar Métodos", icon = icon("balance-scale"), class = "btn-info", width = '215px', style = 'margin-bottom: 10px')
             ),
             
             mainPanel(
               plotOutput("plot_result_simulador", height = "600px"),
             )
           )
  ),
  
  tabPanel("Opciones avanzadas",
           sidebarLayout(
             sidebarPanel(
               h3("Parámetros de Entrada"),
               sliderInput("h", "Paso de tiempo (h):", 0.1, min = 0.01, max = 1, step = 0.01),
               sliderInput("n", "Número de pasos (n):", 0, min = 1, max = 200),
               sliderInput("x", "Cantidad de cafeína:", 90, min = 1, max = 200),
               actionButton("solve_euler2", "Resolver con Euler", icon = icon("chart-line"), class = "btn-primary", width = '215px', style = 'margin-bottom: 10px'),
               actionButton("solve_rk42", "Resolver con RK4", icon = icon("chart-line"), class = "btn-danger", width = '215px', style = 'margin-bottom: 10px'),
               actionButton("solve_rk452", "Resolver con RK45", icon = icon("chart-line"), class = "btn-success", width = '215px', style = 'margin-bottom: 10px'),
               actionButton("solve_comparar2", "Comparar Métodos", icon = icon("balance-scale"), class = "btn-info", width = '215px', style = 'margin-bottom: 10px')
             ),
             
             mainPanel(
               plotOutput("plot_result_avanzadas", height = "600px"),
             )
           )
  )
)


server <- function(input, output) {

  observeEvent(input$solve_euler, {
    h <- input$h
    n <- input$n
    t0 <- 0
    X0 <- 90
    
    res <- metodo_euler(t0, X0, h, n)
    
    output$plot_result_simulador <- renderPlot({
      ggplot(res, aes(x = t, y = X)) +
        geom_line(color = "blue", size = 1.2) +
        geom_point(color = "blue", size = 2) +
        labs(title = "Método de Euler", x = "Tiempo (hr)", y = "Cafeína (mg)") +
        theme_minimal(base_size = 15)
    })
  })
  
  observeEvent(input$solve_rk4, {
    h <- input$h
    n <- input$n
    t0 <- 0
    X0 <- 90
    
    res <- runge_kutta4(t0, X0, h, n)
    
    output$plot_result_simulador <- renderPlot({
      ggplot(res, aes(x = t, y = X)) +
        geom_line(color = "red", size = 1.2) +
        geom_point(color = "red", size = 2) +
        labs(title = "Método de RK4", x = "Tiempo (hr)", y = "Cafeína (mg)") +
        theme_minimal(base_size = 15)
    })
  })
  
  observeEvent(input$solve_rk45, {
    h <- input$h
    n <- input$n
    t0 <- 0
    y0 <- 90
    
    res <- rk45(f_cafeina, t0, y0, h, n)
    
    output$plot_result_simulador <- renderPlot({
      res_df <- as.data.frame(res)
      ggplot(res_df, aes(x = Tiempo, y = Cafeina)) +
        geom_line(color = "green", size = 1.2) +
        geom_point(color = "green", size = 2) +
        labs(title = "Método de RK45", x = "Tiempo (hr)", y = "Cafeína (mg)") +
        theme_minimal(base_size = 15)
    })
  })
  
  observeEvent(input$solve_comparar, {
    h <- input$h
    n <- input$n
    t0 <- 0
    X0 <- 90
    

    res_euler <- metodo_euler(t0, X0, h, n)
    res_euler$Metodo <- "Euler"
    

    res_rk4 <- runge_kutta4(t0, X0, h, n)
    res_rk4$Metodo <- "RK4"
    

    res_rk45 <- rk45(f_cafeina, t0, X0, h, n)
    res_rk45$Metodo <- "RK45"
    

    res_combined <- rbind(
      data.frame(t = res_euler$t, X = res_euler$X, Metodo = res_euler$Metodo),
      data.frame(t = res_rk4$t, X = res_rk4$X, Metodo = res_rk4$Metodo),
      data.frame(t = res_rk45$Tiempo, X = res_rk45$Cafeina, Metodo = res_rk45$Metodo)
    )
    
    output$plot_result_simulador <- renderPlot({
      ggplot(res_combined, aes(x = t, y = X, color = Metodo)) +
        geom_line(size = 1.2) +
        labs(title = "Comparación de Métodos Numéricos", x = "Tiempo (hr)", y = "Cafeína (mg)", color = "Método") +
        theme_minimal(base_size = 15) +
        scale_color_manual(values = c("Euler" = "blue", "RK4" = "red", "RK45" = "green"))
    })
  })
  

  observeEvent(input$solve_euler2, {
    h <- input$h
    n <- input$n
    t0 <- 0
    X0 <- 90
    
    res <- metodo_euler(t0, X0, h, n)
    
    output$plot_result_avanzadas <- renderPlot({
      ggplot(res, aes(x = t, y = X)) +
        geom_line(color = "blue", size = 1.2) +
        geom_point(color = "blue", size = 2) +
        labs(title = "Método de Euler", x = "Tiempo (hr)", y = "Cafeína (mg)") +
        theme_minimal(base_size = 15)
    })
  })
  
  observeEvent(input$solve_rk42, {
    h <- input$h
    n <- input$n
    t0 <- 0
    X0 <- 90
    
    res <- runge_kutta4(t0, X0, h, n)
    
    output$plot_result_avanzadas <- renderPlot({
      ggplot(res, aes(x = t, y = X)) +
        geom_line(color = "red", size = 1.2) +
        geom_point(color = "red", size = 2) +
        labs(title = "Método de RK4", x = "Tiempo (hr)", y = "Cafeína (mg)") +
        theme_minimal(base_size = 15)
    })
  })
  
  observeEvent(input$solve_rk452, {
    h <- input$h
    n <- input$n
    t0 <- 0
    y0 <- 90
    
    res <- rk45(f_cafeina, t0, y0, h, n)
    
    output$plot_result_avanzadas <- renderPlot({
      res_df <- as.data.frame(res)
      ggplot(res_df, aes(x = Tiempo, y = Cafeina)) +
        geom_line(color = "green", size = 1.2) +
        geom_point(color = "green", size = 2) +
        labs(title = "Método de RK45", x = "Tiempo (hr)", y = "Cafeína (mg)") +
        theme_minimal(base_size = 15)
    })
  })
  
  observeEvent(input$solve_comparar2, {
    h <- input$h
    n <- input$n
    t0 <- 0
    X0 <- 90
    
    res_euler <- metodo_euler(t0, X0, h, n)
    res_euler$Metodo <- "Euler"
    
    res_rk4 <- runge_kutta4(t0, X0, h, n)
    res_rk4$Metodo <- "RK4"
    
    res_rk45 <- rk45(f_cafeina, t0, X0, h, n)
    res_rk45$Metodo <- "RK45"
    
    res_combined <- rbind(
      data.frame(t = res_euler$t, X = res_euler$X, Metodo = res_euler$Metodo),
      data.frame(t = res_rk4$t, X = res_rk4$X, Metodo = res_rk4$Metodo),
      data.frame(t = res_rk45$Tiempo, X = res_rk45$Cafeina, Metodo = res_rk45$Metodo)
    )
    
    output$plot_result_avanzadas <- renderPlot({
      ggplot(res_combined, aes(x = t, y = X, color = Metodo)) +
        geom_line(size = 1.2) +
        labs(title = "Comparación de Métodos Numéricos", x = "Tiempo (hr)", y = "Cafeína (mg)", color = "Método") +
        theme_minimal(base_size = 15) +
        scale_color_manual(values = c("Euler" = "blue", "RK4" = "red", "RK45" = "green"))
    })
  })
}

shinyApp(ui = ui, server = server)