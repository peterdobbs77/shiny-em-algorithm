library(shiny)
library(ggplot2)

ui <- fluidPage(
  titlePanel("Shiny-EM-Algorithm",windowTitle = "Shiny-EM-Algorithm"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file1", "Choose CSV File",
                accept = c("text/csv",
                           "text/comma-separated-values",
                           ".csv")),
      checkboxInput("header", "Header", TRUE),
      # radioButtons("separator","Separator",
      #              choices=c("Comma","Semicolon","Tab"),
      #              selected="Comma"),
      numericInput("numModes", "How many modes are in mixture model?",2)
    ),
    mainPanel(
      textOutput("description"),
      tableOutput("contents"),
      plotOutput("distPlot"),
      
    )
  )
)

server <- function(input, output) {
  output$description <- renderText({
    inFile <- input$file1
    if(is.null(inFile))
      return("EXAMPLE")
  })
  
  output$contents <- renderTable({
    inFile <- input$file1
    if (is.null(inFile))
      return(head(faithful[1]))
    
    data <- read.csv(inFile$datapath, header=input$header)
    head(data)
  })
  
  output$distPlot <- renderPlot({
    inFile <- input$file1
    if (is.null(inFile))
      return(ggplot(faithful) + 
               geom_histogram(aes(x=eruptions,y=..density..),
                              binwidth=0.20, fill="cyan", color="black")+
               geom_density(aes(x=eruptions), color="red"))
    
    # read the selected data file 
    data <- read.csv(inFile$datapath, header=input$header)
    
    # generate bins based on input$bins from ui.R
    x    <- data[, 1]
    bins <- seq(min(x), max(x), length.out = 10)
    
    # draw the histogram with the specified number of bins
    hist(x, xlab="X")
    
    #
    # TODO: apply EM-Algorithm here
    #
    #
    
  })
  
}

shinyApp(ui, server)