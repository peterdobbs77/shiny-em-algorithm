library(shiny)
library(ggplot2)

init <- function(data, nModes){
  # initialize model parameters
  pi <- 1
  mu <- mean(data)
  sigma <- 1
  if(nModes>1){
    pi <- rep(1/nModes,nModes)
    mu <- seq(from=min(data),to=max(data),by=(max(data)-min(data))/(nModes-1))
    sigma <- rep(1,nModes)
  }
  
  data.j <- matrix(NA,nrow=nModes,ncol=3)
  data.j[,1] <- pi
  data.j[,2] <- mu
  data.j[,3] <- sigma
  
  list("pi" = data.j[,1],
       "mu" = data.j[,2],
       "sigma" = data.j[,3])
}

e_step <- function(x,nModes,theta){
  post.proc <- matrix(0,nrow=length(x),ncol=nModes)
  post.proc.j <- matrix(0,nrow=length(x),ncol=nModes)
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
  pi <- mu <- sigma <- temp <- rep(0,nModes)
  
  for (i in 1:nModes){
    temp[i] <- sum(post[,i])
    pi[i] <- temp[i]/length(x)
    mu[i] <- (1/temp[i])*sum(post[,i] * x)
    sigma[i] <- sqrt(sum(post[,i]*(x-mu[i])^2) / temp[i])
    
    theta[i,] <- c(pi[i],mu[i],sigma[i])
  }
  list("pi"=theta[,1],
       "mu"=theta[,2],
       "sigma"=theta[,3])
}

iterate <- function(x,nModes,times){
  #browser()
  theta <- init(x,nModes)
  
  # initialize locals
  iter <- 0
  cur.logLike <- 0
  vec.logLike <- 0
  
  # repeat E and M steps until convergence
  for(i in 1:times){
    if(i == 1){
      e.step <- e_step(x,nModes,theta)
      m.step <- m_step(x,nModes,e.step$post)
      vec.logLike <- e.step$logLike
      cur.logLike <- e.step$logLike
    } else {
      theta$pi <- m.step$pi
      theta$mu <- m.step$mu
      theta$sigma <- m.step$sigma
      e.step <- e_step(x,nModes,theta)
      m.step <- m_step(x,nModes,e.step$post)
      
      #
      vec.logLike <- c(vec.logLike, e.step$logLike)
      
      # check convergence
      err.logLike <- abs((cur.logLike - e.step$logLike))
      if(err.logLike < 1e-7){
        break
      } else {
        cur.logLike <- e.step$logLike
      }
    }
    iter <- iter+1
  }
  
  list("result"=m.step,
       "logLike"=vec.logLike,
       "iterations"=iter)
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
      radioButtons("separator","Separator",
                   choices=c("Comma","Semicolon","Tab"),
                   selected="Comma"),
      numericInput("dataColumn", "Column #:",1,min=1),
      numericInput("numModes", "How many modes are in mixture model?",2,min=1),
      numericInput("times", "EM-Step:",0,min=0),
      checkboxInput("fixTheta","Fix mu"),
      # checkboxGroupInput("displayOptions","Display:",
      #                    choiceNames = 
      #                      list("initial","final","legend")),
      sliderInput("bins","Number of bins:",min=1,max=50,value=30)
    ),
    mainPanel(
      tableOutput("contents"),
      plotOutput("fitCriteriaPlot"),
      plotOutput("distPlot")
    )
  )
)

server <- function(input, output) {
  
  data <- reactive({
    # try to read the selected data file 
    inFile <- input$file1
    if(is.null(inFile))
      return()
    read.csv(file=inFile$datapath,
             sep=",",
             header=input$header)
  })
  
  
  output$contents <- renderTable({
    if (is.null(data()))
      return(summary(faithful[input$dataColumn]))
    x <- data.frame(x=data()[, input$dataColumn])
    summary(x)
  })
  
  output$fitCriteriaPlot <- renderPlot({
    x<-0
    if(is.null(data()))
      x <- faithful$eruptions
    else
      x <- data()[,input$dataColumn]
    
    m <- 2000
    nummode <- 1:6
    cur.logL <- iterate(x,1,m)
    vec.logL <- cur.logL$logLike[length(cur.logL)]
    cat("nModes:",1L,"\n")
    cat("iterations:",cur.logL$iterations,"\n\n")
    for (i in 2:max(nummode)){
      cur.logL <- iterate(x,i,m)
      vec.logL <- c(vec.logL, cur.logL$logLike[length(cur.logL)])
      cat("nModes:",i,"\n")
      cat("iterations:",cur.logL$iterations,"\n\n")
    }
    
    k <- 3*nummode-1
    n <- length(x)
    
    aic <- -2*vec.logL + 2*k + 2*k*(k+1)/(n-k-1)
    AIC.df <- data.frame(modes=nummode,
                         values=aic)
    
    bic <- -2*vec.logL + k*log(n)
    BIC.df <- data.frame(modes=nummode,
                         values=bic)
    
    y <- data.frame(modes=nummode,
                    values=c(aic, bic),
                    criteria=c(rep("AIC",max(nummode)),
                             rep("BIC",max(nummode))))
    
    minAIC <- AIC.df[which.min(AIC.df$values),]$modes
    
    minBIC <- BIC.df[which.min(BIC.df$values),]$modes
    
    ggplot(y, aes(x=modes,y=values,group=criteria,color=criteria)) + geom_line() +
      geom_point(aes(x=minAIC,min(aic)),size=2,color="red") +
      geom_point(aes(x=minBIC,min(bic)),size=2,color="cyan") +
      ggtitle("Information Criterion") + xlab("# of Modes") + ylab("AIC/BIC") +
      theme(plot.title = element_text(size=20,face="bold",hjust=0.5))
    
    
  })
  
  output$distPlot <- renderPlot({
    if (is.null(data()))
      return(ggplot(faithful) + 
               geom_histogram(aes(x=eruptions,y=..density..),
                              bins=input$bins, fill="cyan", color="black") +
               geom_density(aes(x=eruptions), color="red"))
    
    x <- data.frame(x=data()[, input$dataColumn])
    
    ggplot(x) + 
      geom_histogram(aes(x=x,y=..density..),
                     bins=input$bins, fill="cyan", color="black")+
      geom_density(aes(x=x), color="red")
  })
  
}

shinyApp(ui, server)