library(shiny)

# Define UI for random distribution app ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Tabsets"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      # Input: Slider for the number of observations to generate ----
      sliderInput("m",
                  "m:",
                  value = 1,
                  min = 1,
                  max = 8),      
      sliderInput("beta",
                  "EPD: beta",
                  value = 1,
                  min = 1,
                  max = 4),
      
      sliderInput("mfun",
                  "function for m",
                  value = 1,
                  min = 1,
                  max = 5)
      ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Tabset w/ plot, summary, and table ----
      tabsetPanel(type = "tabs",
                  tabPanel("Plot", plotOutput("plot"))
                  )
      )
      
    )
  )


# Define server logic for random distribution app ----
server <- function(input, output) {
  l<-levels(b$Estimator)
  ms<-8
  # Reactive expression to generate the requested distribution ----
  # This is called whenever the inputs change. The output functions
  # defined below then use the value computed from this expression
  d <- reactive({
    dist <- switch(input$dist,
                   norm = rnorm,
                   unif = runif,
                   lnorm = rlnorm,
                   exp = rexp,
                   rnorm)
    
    dist(input$n)
  })
  
  # Generate a plot of the data ----
  # Also uses the inputs to build the plot label. Note that the
  # dependencies on the inputs and the data reactive expression are
  # both tracked, and all expressions are called in the sequence
  # implied by the dependency graph.
  output$plot <- renderPlot({
    dist <- input$dist
    beta <- input$beta
    m <- input$m
    mfun <- input$mfun
    
    plot(a[,1],main="Estimatates of the EVI",ylab=bquote(xi),xlab="K",ylim=c(0,1.5),type="l",lwd=3)
    lines(a[,m],lwd=4,col="purple")
    lines(a[,ms+mfun],lwd=4,col="red4")
    lines(a[,ms+5+beta],lwd=4,col="steelblue")
    abline(h=0.5)
    legend("bottom", bty="n",col=c(1,"purple","red4","steelblue"),lwd=4,
           sapply(c(bquote(H[k]),bquote(H[xi]^Pareto~(m==.(m))),bquote(H[xi]^Pareto~(m=="M*(k/n)^"~.(mfun))),bquote(EPD~(beta=="-"~.(beta)))), as.expression))
    
    })
  
  }

# Create Shiny app ----
shinyApp(ui, server)