library(shiny)
load("data500.RData")
values <- list(A = c('k/log(k)', 'M*(k/n)^1', 'M*(k/n)^2','M*(k/n)^3','M*(k/n)^4'))
# Define UI for random distribution app ----
ui <- fluidPage(
  # Use the Google webfont "Source Sans Pro"
  tags$link(
    href=paste0("http://fonts.googleapis.com/css?",
                "family=Source+Sans+Pro:300,600,300italic"),
    rel="stylesheet", type="text/css"),
  tags$style(type="text/css",
             "body {font-family: 'Source Sans Pro'}"
  ),
  
  h2("Exper for promoters"),
  br(),
  h3("Inputs"),
  code("please feel free to experiment with these sliders."),
  br(),
  code("(ignore the selection dropbox, I don't know why I put it there...)")
  ,
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
  
      selectInput('selection', 'selection', 'A'),
      uiOutput('selectUI'),
      sliderInput(inputId = "target", label = "m function",
                  min = 1, max = length(values$A),
                  step = 1,
                  value = length(values$A)),
      
      sliderInput("m",
                  "m:",
                  value = 1,
                  min = 1,
                  max = 52),      
      sliderInput("beta",
                  "EPD: beta",
                  value = 2,
                  min = 2,
                  max = 4)
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
server <- function(input, output,session) {
  dc <- c("#3366cc", "#dc3912", "#ff9900", "#109618", "#990099", "#0099c6", "#dd4477")
  ms<-52
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
  output$selectUI <- renderUI({
    sel_values <- paste(paste0('"', values[[input$selection]], '"'), collapse = ',')
    print(sel_values)
    list(
      (HTML(
        sprintf('
                <script type="text/javascript">
                $(document).ready(function() {
                var vals = [%s];
                $(\'#target\').data(\'ionRangeSlider\').update(
                {values:vals,
                min: 0,
                max: %s,
                from:%s})
                })
                </script>
                ', sel_values, 
                length(values[[input$selection]]),
                length(values[[input$selection]]))))
    )}
    
    
  )
    
  # Generate a plot of the data ----
  # Also uses the inputs to build the plot label. Note that the
  # dependencies on the inputs and the data reactive expression are
  # both tracked, and all expressions are called in the sequence
  # implied by the dependency graph.
  output$plot <- renderPlot({
    dist <- input$dist
    beta <- input$beta
    m <- input$m
    mfun <- input$target
    
    plot(a[,1],main="Estimatates of the EVI",ylab=bquote(xi),xlab="K",ylim=c(0,1.5),type="l",lwd=3)
    lines(a[,m],lwd=4,col=dc[1])
    lines(a[,ms+mfun+1],lwd=4,col=dc[2])
    lines(a[,ms+5+beta],lwd=4,col=dc[3])
    abline(h=0.5)
    legend("bottom", bty="n",col=c(1,dc[1:3]),lwd=4,ncol=2,
           sapply(c(bquote(H[k]),bquote(H[xi]^Pareto~(m==.(m))),bquote(H[xi]^Pareto~(m=="M*(k/n)^"~.(mfun))),bquote(EPD~(beta=="-"~.(beta)))), as.expression))
    
    })
  
  
  }

# Create Shiny app ----
shinyApp(ui, server)
