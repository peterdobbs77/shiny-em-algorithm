library(ggplot2)

df <- data.frame(x=faithful$waiting,y=faithful$eruptions)
ggplot(df, aes(x=x, y=y)) + geom_point()

p <- ggplot(df) + 
  geom_histogram(aes(x=y, y=..density..),
                 binwidth=0.25, fill="cyan", color="black")
p + geom_density(data=faithful,aes(x=eruptions), color="red")

df <- data.frame(x=faithful$eruptions)

theta <- list(
  pi=0.5,
  mu=c(2,5),
  sig=c(1,1)
)

#logL<-NULL; mus<-seq(min(df$x),max(df$x),length=100); mu<- 1

#logLike <- function(mu,data=df$x){sum(dnorm(df$x,mu))}
#max.logL <- optim(1,logLike,df$x,control=list(fnscale=-1))

#for(i in 1:length(mus)){logL[i]<-logLike(mus[i],df$x)}
