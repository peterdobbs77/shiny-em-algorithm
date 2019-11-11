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
div <- split(sort(x), ceiling(seq_along(x)/(length(x)/nModes)))

pi <- rep(1/nModes,times=nModes)
mu <- unlist(lapply(div,mean))
s2 <- unlist(lapply(div,var))

logLike.mu <- function(mu,data){
  sum(dnorm(data$x,mu))
}
#max.logL <- optim(1,logLike,x,control=list(fnscale=-1))

#for(i in 1:70){logL[i]<-logLike(mus[i],df$x)}

# f(x; pi, mu, s2) = sum(pi * fi(x; mu_i, sqrt(s2_i))
