library(ggplot2)

ggplot(faithful, aes(x=waiting, y=eruptions)) + geom_point()

p <- ggplot(faithful) + 
  geom_histogram(aes(x=eruptions, y=..density..),
                 bins=30, fill="cyan", color="black")
p + geom_density(data=faithful,aes(x=eruptions), color="red")

x <- faithful$eruptions

# number of modes
nModes <- 5
# TODO: use AIC to compute optimal number of modes


# initialize model parameters
pi <- rep(1/nModes,times=nModes)
mu <- seq(from=min(x),to=max(x),by=(max(x)-min(x))/(nModes-1))
s2 <- rep(1,times=nModes)

x.j <- matrix(NA,nrow=nModes,ncol=3)
x.j[,1] <- pi
x.j[,2] <- mu
x.j[,3] <- s2

theta <- list("pi" = x.j[,1],
              "mu" = x.j[,2],
              "s2" = x.j[,3])


logLike.mu <- function(mu,data){
  sum(dnorm(data$x,mu))
}
#max.logL <- optim(1,logLike,x,control=list(fnscale=-1))

#for(i in 1:70){logL[i]<-logLike(mus[i],df$x)}

# f(x; pi, mu, s2) = sum(pi * fi(x; mu_i, sqrt(s2_i))
