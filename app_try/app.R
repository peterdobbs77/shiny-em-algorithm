library(shiny)
library(ggplot2)

init <- function(x, nModes){
  # initialize model parameters
  pi <- rep(1/nModes,times=nModes)
  mu <- seq(from=min(x),to=max(x),by=(max(x)-min(x))/(nModes-1))
  sigma <- rep(1,nModes)
  
  x.j <- matrix(NA,nrow=nModes,ncol=3)
  x.j[,1] <- pi
  x.j[,2] <- mu
  x.j[,3] <- sigma
  
  theta <- list("pi" = x.j[,1],
                "mu" = x.j[,2],
                "sigma" = x.j[,3])
  theta
}

e_step <- function(x,nModes,theta){
  post.proc <- matrix(0,nrow=length(data),ncol=nModes)
  post.proc.j <- matrix(0,nrow=length(data),ncol=nModes)
  sum.probs <- 0
  for (i in 1:nModes){
    post.proc[,i] <- theta$pi[i] * dnorm(x,theta$mu[i],theta$sigma[i])
    sum.probs <- sum.probs + post.proc[,i]
  }
  
  for (i in 1:nModes){
    post.proc.j[,i] <- post.proc[,i]/sum.probs
  }
  
  comp.ln <- log(sum.probs)
  comp.ln.sum <- sum(comp.ln)
  logLike <- comp.ln.sum
  list("logLike"=logLike,
       "post"=post.proc.j)
}

m_step <- function(x,nModes,post){
  theta <- matrix(0, nrow = nModes, ncol = 3)
  colnames(theta) <- paste(c('pi','mu','sigma'))
  pi <- mu <- sigma <- temp <- rep(0,nModes)
  
  for (i in 1:nModes){
    temp[i] <- sum(post[,i])
    pi[i] <- temp[i]/length(x)
    mu[i] <- (1/temp[i])*sum(post[,i] * x)
    sigma[i] <- sqrt(sum(post[,i]*(x-mu[i])^2) / temp[i])
    
    theta[i] <- c(pi[i],mu[i],sigma[i])
  }
  theta
}

iterate <- function(x,nModes,times){
  theta <- init(x,nModes)
  
  # repeat E and M steps until convergence
  for(i in 1:times){
    if(i==1){
      e.step <- e_step(x,nModes,theta)
      m.step <- m_step(x,nModes,e.step$post)
      logLike.c <- e.step$logLike
      logLike.i <- e.step$logLike
      i <- i + 1
    }
    theta$pi <- m.step$pi
    theta$mu <- m.step$mu
    theta$sigma <- m.step$sigma
    e.step <- e_step(x,nModes,theta)
    m.step <- m_step(x,nModes,e.step$post)
    
    #
    logLike.c <- c(logLike.c, e.step$logLike)
    
    # check convergence
    logLike.err <- abs(logLike.i - e.step$logLike)
    if(logLike.err < 1e-4)
      break
    
    logLike.i <- e.step$logLike
  }
  
  list("result"=m.step,
       "logLike"=logLike.c)
}

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
      numericInput("dataColumn", "Column #:",1),
      numericInput("numModes", "How many modes are in mixture model?",2),
      numericInput("times", "EM-Step:",0),
      checkboxInput("fixTheta","Fix expression(mu)"),
      checkboxGroupInput("displayOptions","Display:",
                         choiceNames = 
                           list("initial","final","legend")),
      sliderInput("bins","Number of bins:",min=1,max=50,value=30)
    ),
    mainPanel(
      textOutput("description"),
      tableOutput("contents"),
      plotOutput("distPlot")
    )
  )
)

server <- function(input, output) {
  
  
  data <- reactive({
    # try to read the selected data file 
    inFile <- input$file1
    if(is.null(inFile)){return()} 
    read.csv(file=inFile$datapath, sep=",", header = input$header)
  })
  
  times <- reactive({input$times})
  
  nModes <- reactive({input$nModes})
  
  output$description <- renderText({
    if(is.null(data()))
      return("EXAMPLE")
    return(NULL)
  })
  
  output$contents <- renderTable({
    if (is.null(data()))
      return(summary(faithful[1]))
    summary(data)
  })
  
  output$distPlot <- renderPlot({
    if (is.null(data()))
      return(ggplot(faithful) + 
               geom_histogram(aes(x=eruptions,y=..density..),
                              bins=input$bins, fill="cyan", color="black")+
               geom_density(aes(x=eruptions), color="red"))
    
    x <- data.frame(x=data()[, input$dataColumn])
    
    ggplot(x) + 
      geom_histogram(aes(x=x,y=..density..),
                     bins=input$bins, fill="cyan", color="black")+
      geom_density(aes(x=x), color="red")
  })
  
  #
  # TODO: apply EM-Algorithm here to create pdf
  #
  #
  output$steps <- renderPrint({
    if(is.null(data()))
      return(NULL)
    k <- iterate(x,nModes(),times())
    print(k$result)
    print(k$logLike)
  })
  
}

shinyApp(ui, server)