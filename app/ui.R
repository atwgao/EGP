fluidPage(
  # Use the Google webfont "Source Sans Pro"
  tags$link(
    href=paste0("http://fonts.googleapis.com/css?",
                "family=Source+Sans+Pro:300,600,300italic"),
    rel="stylesheet", type="text/css"),
  tags$style(type="text/css",
             "body {font-family: 'Source Sans Pro'}"
  ),
  
  h2("Experiments for promoters"),
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

