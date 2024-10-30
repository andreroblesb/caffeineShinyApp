library(shiny)
library(deSolve)
library(ggplot2)
library(bslib)


f_cafeina <- function(t, X, a) {
  return(-0.15 * X + a)
}

f_farmacos1 <- function(t, X) {
  a <- 0.03
  k <- 0.02
  
  Xg <- X[1] 
  Xb <- X[2]  

  dxdt_Xg <- -a * Xg
  dxdt_Xb <- a * Xg - k * Xb

  return(c(dxdt_Xg, dxdt_Xb))
}

metodo_euler_farmaco <- function(t0, Xg0, Xb0, h, n, a, k) {
  t <- numeric(n + 1)
  Xg <- numeric(n + 1)
  Xb <- numeric(n + 1)
  
  t[1] <- t0
  Xg[1] <- Xg0
  Xb[1] <- Xb0
  
  for (i in 1:n) {
    Xg[i + 1] <- Xg[i] + h * (a * Xg[i] - k * Xb[i])
    Xb[i + 1] <- Xb[i] + h * (k * Xb[i] - a * Xg[i])
    t[i + 1] <- t[i] + h
  }
  
  return(data.frame(time = t, Xg = Xg, Xb = Xb))
}

# metodo Euler
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
  
  return(data.frame(Tiempo = t, Cafeina = y))
}

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


ui <- navbarPage(
  title = "Ecuaciones diferenciales",
  theme = bslib::bs_theme(bootswatch = "solar"),
  
  tabPanel("Home",
           h1("Simulación de ecuaciones diferenciales para el consumo de cafeína y fármacos. ", style = "color: beige;"),
           p("Con el uso de ecuaciones diferenciales, en esta app nos dimos la tarea de observar cómo se comporta el metabolismo de la cafeína y observamos la absorción y liberación de un fármaco en el hígado. 
             Para cada una de las modelaciones, tenemos una ecuación en el caso de la cafeína, y dos ecuaciones en el caso del fármaco, que representan el comportamiento dentro del cuerpo humano respecto al tiempo.", style = "color: beige;"),
           tags$ul(
             tags$li("Cafeína", style = "color: beige;"),
             p("X = -0.15x + a", style = "color: beige;"),
             tags$ul(
               tags$li("Donde:", style = "color: beige;"),
               tags$li("X es la cantidad de café en el cuerpo", style = "color: beige;"),
               tags$li("0.15 representa el porcentaje que elimina la cafeína por cada hora", style = "color: beige;"),
               tags$li("a cantidad de consumo de café constante", style = "color: beige;")
             ),
             tags$li("Fármaco", style = "color: beige;"),
             p("Absorción", style = "color: beige;"), 
             p("Xb = aXg", style = "color: beige;"),
             tags$ul(
               tags$li("Donde:", style = "color: beige;"),
               tags$li("Xg cantidad de fármaco en el sistema", style = "color: beige;"),
               tags$li("a es el porcentaje de eliminación por minuto", style = "color: beige;")),
             p("Liberación ", style = "color: beige;"),
             p("Xb = -aXg", style = "color: beige;"),
             tags$ul(
               tags$li("Donde:", style = "color: beige;"),
               tags$li("Xg cantidad de fármaco en el sistema", style = "color: beige;"),
               tags$li("a es el porcentaje de eliminación por minuto", style = "color: beige;")
             ),
           ),
           p("Este simulador compara diferentes métodos numéricos para modelar la eliminación de cafeína en el cuerpo humano.", style = "color: beige;"),
           p("Métodos numéricos:", style = "color: beige;"),
           tags$ul(
             tags$li("Método de Euler", style = "color: beige;"),
             tags$li("Método de Runge-Kutta de Orden 4 (RK4)", style = "color: beige;"),
             tags$li("Método de Runge-Kutta-Fehlberg 4(5) (RK45)", style = "color: beige;"),
             tags$li("Método de Runge-Kutta-Fehlberg adaptativo (RK45 Fehlberg)", style = "color: beige;")
           ),
  ),
  
  tabPanel("Simulador",
           sidebarLayout(
             sidebarPanel(
               h3("Parámetros de Entrada"),
               sliderInput("h", "Paso de tiempo (h):", 0.1, min = 0.01, max = 1, step = 0.01),
               sliderInput("n", "Número de pasos (n):", 1, min = 1, max = 200),
               actionButton("solve_euler", "Euler - farmaco", icon = icon("chart-line"), class = "btn-primary", width = '215px'),
               actionButton("solve_rk4", "Resolver con RK4", icon = icon("chart-line"), class = "btn-danger", width = '215px'),
               actionButton("solve_rk45", "Resolver con RK45", icon = icon("chart-line"), class = "btn-success", width = '215px'),
               actionButton("solve_rkf45", "Resolver con RK45 Fehlberg", icon = icon("chart-line"), class = "btn-secundary", width = '215px'),
               actionButton("solve_comparar", "Comparar Métodos", icon = icon("balance-scale"), class = "btn-info", width = '215px')
             ),
             
             mainPanel(
               plotOutput("plot_result_simulador", height = "600px")
             )
           )
  ),
  
  tabPanel("Opciones avanzadas",
           sidebarLayout(
             sidebarPanel(
               h3("Parámetros de Entrada"),
               sliderInput("h", "Paso de tiempo (h):", 0.1, min = 0.01, max = 1, step = 0.01),
               sliderInput("n", "Número de pasos (n):", 1, min = 1, max = 200),
               sliderInput("X0", "Cantidad inicial de cafeína (mg):", 90, min = 1, max = 200),
               sliderInput("a", "Consumo adicional (mg/h):", 1, min = 1, max = 200),
               actionButton("solve_euler2", "Resolver con Euler", icon = icon("chart-line"), class = "btn-primary", width = '215px'),
               actionButton("solve_rk42", "Resolver con RK4", icon = icon("chart-line"), class = "btn-danger", width = '215px'),
               actionButton("solve_rk452", "Resolver con RK45", icon = icon("chart-line"), class = "btn-success", width = '215px'),
               actionButton("solve_rkf452", "Resolver con RK45 Fehlberg", icon = icon("chart-line"), class = "btn-secundary", width = '215px'),
               actionButton("solve_comparar2", "Comparar Métodos", icon = icon("balance-scale"), class = "btn-info", width = '215px')
             ),
             
             mainPanel(
               plotOutput("plot_result_avanzadas", height = "600px")
             )
           )
  )
)


server <- function(input, output) {
  observeEvent(input$solve_euler, {
    h <- as.numeric(input$h)
    n <- as.numeric(input$n)
    Xg0 <- input$a
    Xb0 <- 0 
    k <- 0.02
    a <- 0.03
    
    res <- metodo_euler_farmaco(0, Xg0, Xb0, h, n, a, k)
    
    output$plot_result_simulador <- renderPlot({
      ggplot(res, aes(x = time)) +
        geom_line(aes(y = Xg, color = "Xg"), size = 1.2) +
        geom_point(aes(y = Xg, color = "Xg"), size = 2) +
        geom_line(aes(y = Xb, color = "Xb"), size = 1.2) +
        geom_point(aes(y = Xb, color = "Xb"), size = 2) +
        labs(title = "Método de Euler", x = "Tiempo (hr)", y = "Concentración (mg)") +
        theme_minimal(base_size = 15) +
        scale_color_manual(values = c("Xg" = "blue", "Xb" = "red")) +
        guides(color = guide_legend(title = "Componente"))
    })
  })
  
  observeEvent(input$solve_rk4, {
    h <- input$h
    n <- input$n
    t0 <- 0
    X0 <- 90
    a <- 0
    
    res <- runge_kutta4(t0, X0, h, n, a)
    print(length(res$X))
    
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
    a <- 0
    
    res <- rk45(t0, y0, h, n, a)
    print(length(res$X))
    
    output$plot_result_simulador <- renderPlot({
      ggplot(res, aes(x = Tiempo, y = Cafeina)) +
        geom_line(color = "green", size = 1.2) +
        geom_point(color = "green", size = 2) +
        labs(title = "Método de RK45", x = "Tiempo (hr)", y = "Cafeína (mg)") +
        theme_minimal(base_size = 15)
    })
  })
  
 
  observeEvent(input$solve_rkf45, {
    h <- input$h
    n <- input$n
    t0 <- 0
    X0 <- 90
    a <- 0
    
    res <- rkf45(t0, X0, h, n, a)
    print(length(res$X))
    
    output$plot_result_simulador <- renderPlot({
      ggplot(res, aes(x = t, y = X)) +
        geom_line(color = "purple", size = 1.2) +
        geom_point(color = "purple", size = 2) +
        labs(title = "Método de RK45 Fehlberg", x = "Tiempo (hr)", y = "Cafeína (mg)") +
        theme_minimal(base_size = 15)
    })
  })
  
  observeEvent(input$solve_comparar, {
    h <- input$h
    n <- input$n
    t0 <- 0
    X0 <- 90
    a <- 0
    
   
    res_euler <- metodo_euler(t0, X0, h, n, a)
    res_euler$Metodo <- "Euler"
    
   
    res_rk4 <- runge_kutta4(t0, X0, h, n, a)
    res_rk4$Metodo <- "RK4"
    
    
    res_rk45 <- rk45(t0, X0, h, n, a)
    res_rk45$Metodo <- "RK45"
    
    
    res_rkf45 <- rkf45(t0, X0, h, n, a)
    res_rkf45$Metodo <- "RK45 Fehlberg"
    

    res_combined <- rbind(
      data.frame(t = res_euler$t, X = res_euler$X, Metodo = res_euler$Metodo),
      data.frame(t = res_rk4$t, X = res_rk4$X, Metodo = res_rk4$Metodo),
      data.frame(t = res_rk45$Tiempo, X = res_rk45$Cafeina, Metodo = res_rk45$Metodo),
      data.frame(t = res_rkf45$t[1:length(res_rk45$t)], X = res_rkf45$X[1:length(res_rk45$t)], Metodo = res_rkf45$Metodo)
    )
    
    output$plot_result_simulador <- renderPlot({
      ggplot(res_combined, aes(x = t, y = X, color = Metodo)) +
        geom_line(size = 1.2) +
        labs(title = "Comparación de Métodos Numéricos", x = "Tiempo (hr)", y = "Cafeína (mg)", color = "Método") +
        theme_minimal(base_size = 15) +
        scale_color_manual(values = c("Euler" = "blue", "RK4" = "red", "RK45" = "green", "RK45 Fehlberg" = "purple"))
    })
  })
  
  # Opciones avanzadas - Euler
  observeEvent(input$solve_euler2, {
    h <- input$h
    n <- input$n
    t0 <- 0
    X0 <- input$X0
    a = input$a
    
    res <- metodo_euler(t0, X0, h, n, a)
    
    output$plot_result_avanzadas <- renderPlot({
      ggplot(res, aes(x = t, y = X)) +
        geom_line(color = "blue", size = 1.2) +
        geom_point(color = "blue", size = 2) +
        labs(title = "Método de Euler", x = "Tiempo (hr)", y = "Cafeína (mg)") +
        theme_minimal(base_size = 15)
    })
  })
  
  # Opciones avanzadas - RK4
  observeEvent(input$solve_rk42, {
    h <- input$h
    n <- input$n
    t0 <- 0
    X0 <- input$X0 
    a <- input$a
    
    res <- runge_kutta4(t0, X0, h, n, a)
    
    output$plot_result_avanzadas <- renderPlot({
      ggplot(res, aes(x = t, y = X)) +
        geom_line(color = "red", size = 1.2) +
        geom_point(color = "red", size = 2) +
        labs(title = "Método de RK4", x = "Tiempo (hr)", y = "Cafeína (mg)") +
        theme_minimal(base_size = 15)
    })
  })
  
  # Opciones avanzadas - RK45
  observeEvent(input$solve_rk452, {
    h <- input$h
    n <- input$n
    t0 <- 0
    y0 <- input$X0 
    a <- input$a
    
    res <- rk45(t0, y0, h, n, a)
    
    output$plot_result_avanzadas <- renderPlot({
      ggplot(res, aes(x = Tiempo, y = Cafeina)) +
        geom_line(color = "green", size = 1.2) +
        geom_point(color = "green", size = 2) +
        labs(title = "Método de RK45", x = "Tiempo (hr)", y = "Cafeína (mg)") +
        theme_minimal(base_size = 15)
    })
  })
  observeEvent(input$solve_rkf452, {
    h <- input$h
    n <- input$n
    t0 <- 0
    X0 <- input$X0 
    a <- input$a
    res <- rkf45(t0, X0, h, n, a)
    
    output$plot_result_avanzadas <- renderPlot({
      ggplot(res, aes(x = t, y = X)) +
        geom_line(color = "purple", size = 1.2) +
        geom_point(color = "purple", size = 2) +
        labs(title = "Método de RK45 Fehlberg", x = "Tiempo (hr)", y = "Cafeína (mg)") +
        theme_minimal(base_size = 15)
    })
  })
  
  # Comparar los tres métodos - Opciones avanzadas
  observeEvent(input$solve_comparar2, {
    h <- input$h
    n <- input$n
    t0 <- 0
    X0 <- input$X0 
    a = input$a

    res_euler <- metodo_euler(t0, X0, h, n, a)
    res_euler$Metodo <- "Euler"
    
    
    res_rk4 <- runge_kutta4(t0, X0, h, n, a)
    res_rk4$Metodo <- "RK4"
    
    
    res_rk45 <- rk45(t0, X0, h, n, a)

    res_rk45$Metodo <- "RK45"
    
    
    res_rkf45 <- rkf45(t0, X0, h, n, a)
    res_rkf45$Metodo <- "RK45 Fehlberg"
   
    res_combined <- rbind(
      data.frame(t = res_euler$t, X = res_euler$X, Metodo = res_euler$Metodo),
      data.frame(t = res_rk4$t, X = res_rk4$X, Metodo = res_rk4$Metodo),
      data.frame(t = res_rk45$Tiempo, X = res_rk45$Cafeina, Metodo = res_rk45$Metodo),
      data.frame(t = res_rkf45$t, X = res_rkf45$X, Metodo = res_rkf45$Metodo)
    )
    
    output$plot_result_avanzadas <- renderPlot({
      ggplot(res_combined, aes(x = t, y = X, color = Metodo)) +
        geom_line(size = 1.2) +
        labs(title = "Comparación de Métodos Numéricos", x = "Tiempo (hr)", y = "Cafeína (mg)", color = "Método") +
        theme_minimal(base_size = 15) +
        scale_color_manual(values = c("Euler" = "blue", "RK4" = "red", "RK45" = "green", "RK45 Fehlberg" = "grey"))
    })
  })
}

shinyApp(ui = ui, server = server)
