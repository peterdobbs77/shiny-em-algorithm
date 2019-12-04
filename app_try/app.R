library(shiny)
library(ggplot2)

init <- function(data, nModes){
    # initialize model parameters
    pi <- 1
    mu <- mean(data)
    sd <- 1
    if(nModes>1){
        pi <- rep(1/nModes,nModes)
        mu <- seq(from=min(data),
                  to=max(data),
                  by=(max(data)-min(data))/(nModes-1))
        sd <- runif(nModes,0,1)
    }
    
    data.j <- matrix(NA,nrow=nModes,ncol=3)
    data.j[,1] <- pi
    data.j[,2] <- mu
    data.j[,3] <- sd
    
    list("pi" = data.j[,1],
         "mu" = data.j[,2],
         "sd" = data.j[,3])
}

e_step <- function(x,nModes,theta){
    post.proc <- matrix(0,nrow=length(x),ncol=nModes)
    post.proc.j <- matrix(0,nrow=length(x),ncol=nModes)
    sum.probs <- 0
    for (k in 1:nModes){
        post.proc[,k] <- theta$pi[k] * dnorm(x,theta$mu[k],theta$sd[k])
        sum.probs <- sum.probs + post.proc[,k]
    }
    
    for (k in 1:nModes){
        post.proc.j[,k] <- post.proc[,k]/sum.probs
    }
    
    comp.ln <- log(sum.probs)
    comp.ln.sum <- sum(comp.ln)
    logLike <- comp.ln.sum
    list("logLike"=logLike,
         "post"=post.proc.j)
}

m_step <- function(x,nModes,post){
    theta <- matrix(0, nrow = nModes, ncol = 3)
    pi <- mu <- sd <- temp <- rep(0,nModes)
    
    for (i in 1:nModes){
        temp[i] <- sum(post[,i])
        pi[i] <- temp[i]/length(x)
        mu[i] <- (1/temp[i])*sum(post[,i] * x)
        sd[i] <- sqrt(sum(post[,i]*(x-mu[i])^2) / temp[i])
        
        theta[i,] <- c(pi[i],mu[i],sd[i])
    }
    list("pi"=theta[,1],
         "mu"=theta[,2],
         "sd"=theta[,3])
}

iterate <- function(x,nModes,times){
    x <- as.double(x)
    theta <- init(x,nModes)
    
    df.pi <- data.frame(pi=t(theta$pi))
    df.mu <- data.frame(mu=t(theta$mu))
    df.sd <- data.frame(sd=t(theta$sd))
    
    # initialize locals
    e.step <- 0
    m.step <- 0
    cur.logLike <- 0
    vec.logLike <- 0
    
    # browser()
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
            theta$sd <- m.step$sd
            e.step <- e_step(x,nModes,theta)
            m.step <- m_step(x,nModes,e.step$post)
            
            # add results to running lists
            df.pi <- rbind(df.pi,theta$pi)
            df.mu <- rbind(df.mu,theta$mu)
            df.sd <- rbind(df.sd,theta$sd)
            vec.logLike <- c(vec.logLike, e.step$logLike)
            
            # check convergence
            err.logLike <- as.double(abs((cur.logLike - e.step$logLike)))
            if(err.logLike < 1e-4){
                break
            } else {
                cur.logLike <- e.step$logLike
            }
        }
    }
    
    df.theta <- cbind(df.pi,df.mu,df.sd)
    
    list("result"=m.step,
         "logLike"=vec.logLike,
         "theta"=df.theta,
         "responsibility"=e.step$post)
}

ui <- fluidPage(
    titlePanel("Shiny-EM-Algorithm", windowTitle = "Shiny-EM-Algorithm"),
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
            numericInput("numModes", "How many modes are in mixture model?",1,min=1,max=10),
            numericInput("numSteps", "EM-Step:",0,min=0),
            # checkboxInput("fixTheta","Fix mu"),
            # checkboxGroupInput("displayOptions","Display:",
            #                    choiceNames = 
            #                      list("initial","final","legend")),
            sliderInput("bins","Number of bins:",min=1,max=50,value=30),
            # hr(),
            # tagList("Url link:",a("Github Page",href="https://github.com/peterdobbs77/shiny-em-algorithm"))
        ),
        mainPanel(
            tabsetPanel(type="tabs",
                        tabPanel("Data", tableOutput("contents")),
                        tabPanel("Model Selection", plotOutput("fitCriteriaPlot")),
                        tabPanel("Plot", plotOutput("distPlot")),
                        tabPanel("Summary",tableOutput("theta"))
            )
        )
    )
)

server <- function(input, output, session) {
    
    data <- reactive({
        # try to read the selected data file 
        inFile <- input$file1
        if(is.null(inFile))
            return(NULL)
        read.csv(file=inFile$datapath,
                 sep=",",
                 header=input$header)
    })
    
    output$contents <- renderTable({ data() })
    
    output$fitCriteriaPlot <- renderPlot({
        if(is.null(data())) return(NULL)
        x <- as.numeric(data()[,input$dataColumn])
        
        m <- 2000
        nummode <- 1:8
        cur.logL <- iterate(x,1,m)
        vec.logL <- cur.logL$logLike[length(cur.logL$logLike)]
        vec.iter <- nrow(cur.logL$theta)
        for (i in 2:max(nummode)){
            cur.logL <- iterate(x,i,m)
            vec.logL <- c(vec.logL, cur.logL$logLike[length(cur.logL$logLike)])
            vec.iter <- c(vec.iter, nrow(cur.logL$theta))
        }
        
        k <- 3*nummode-1
        n <- length(x)
        
        aic <- -2*vec.logL + 2*k
        AIC.df <- data.frame(modes=nummode,
                             values=aic)
        
        bic <- -2*vec.logL + k*log(n)
        BIC.df <- data.frame(modes=nummode,
                             values=bic)
        
        y <- data.frame(modes=nummode,
                        iters=vec.iter,
                        values=c(aic, bic),
                        criteria=c(rep("AIC",max(nummode)),
                                   rep("BIC",max(nummode))))
        # print(y)
        
        minAIC <- AIC.df[which.min(AIC.df$values),]$modes
        
        minBIC <- BIC.df[which.min(BIC.df$values),]$modes
        
        
        updateNumericInput(session, "numModes", value = y[which.min(y$values),]$modes)
        updateNumericInput(session, "numSteps", value = y[which.min(y$values),]$iters)
        
        ggplot(y, aes(x=modes,y=values,group=criteria,color=criteria)) + geom_line() +
            geom_point(aes(x=minAIC,min(aic)),size=2,color="red") +
            geom_point(aes(x=minBIC,min(bic)),size=2,color="cyan") +
            ggtitle("Information Criterion") + xlab("# of Modes") + ylab("AIC/BIC") +
            theme(plot.title = element_text(size=16,face="bold",hjust=0.5))
    })
    
    output$distPlot <- renderPlot({
        if(is.null(data())) return(NULL)
        d <- data()[, input$dataColumn]
        nModes <- input$numModes
        nSteps <- input$numSteps
        res <- iterate(d,nModes,nSteps)
        
        m<-10000
        comp <- sample(1:nModes,prob=res$result$pi,size=m,replace=TRUE)
        dist <- rnorm(n=m,mean=res$result$mu[comp],sd=res$result$sd[comp])
        
        bins <- seq(from = min(d),
                    to = max(d),
                    length.out = input$bins+1)
        #browser()
        hist(d,breaks=bins,xlab="Data",prob=T)
        lines(density(dist,adjust=0.5),lwd=2,col="red")
    })
    
    output$theta <- renderTable({
        if(is.null(data())) return(NULL)
        x <- data()[,input$dataColumn]
        
        nModes <- input$numModes
        nSteps <- input$numSteps
        
        res <- iterate(x,nModes,nSteps)
        res$theta
    })
    
}

shinyApp(ui, server)
