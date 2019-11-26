library(ggplot2)


init <- function(x, nModes){
  # initialize model parameters
  pi <- rep(1/nModes,times=nModes)
  mu <- seq(from=min(x),to=max(x),by=(max(x)-min(x))/(nModes-1))
  sd <- rep(1,nModes)
  
  x.j <- matrix(NA,nrow=nModes,ncol=3)
  x.j[,1] <- pi
  x.j[,2] <- mu
  x.j[,3] <- sd
  
  theta <- list("pi" = x.j[,1],
                "mu" = x.j[,2],
                "sd" = x.j[,3])
  theta
}

e_step <- function(x,nModes,theta){
  post.proc <- matrix(0,nrow=length(x),ncol=nModes)
  post.proc.j <- matrix(0,nrow=length(x),ncol=nModes)
  sum.probs <- 0
  for (i in 1:nModes){
    post.proc[,i] <- theta$pi[i] * dnorm(x,theta$mu[i],theta$sd[i])
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
  theta <- init(x,nModes)
  df.pi <- data.frame(pi=t(theta$pi))
  df.mu <- data.frame(mu=t(theta$mu))
  df.sd <- data.frame(sd=t(theta$sd))
  # repeat E and M steps until convergence
  for(i in 1:times){
    if(i==1){
      e.step <- e_step(x,nModes,theta)
      m.step <- m_step(x,nModes,e.step$post)
      logLike.c <- data.frame(logL=e.step$logLike)
      logLike.i <- e.step$logLike
      i <- i + 1
    }
    theta$pi <- m.step$pi
    theta$mu <- m.step$mu
    theta$sd <- m.step$sd
    e.step <- e_step(x,nModes,theta)
    m.step <- m_step(x,nModes,e.step$post)
    
    # add results to running lists
    df.pi <- rbind(df.pi,theta$pi)
    df.mu <- rbind(df.mu,theta$mu)
    df.sd <- rbind(df.sd,theta$sd)
    logLike.c <- rbind(logLike.c, e.step$logLike)
    
    # check convergence
    logLike.err <- abs(logLike.i - e.step$logLike)
    if(logLike.err < 1e-2)
      break
    
    logLike.i <- e.step$logLike
  }
  
  df.theta <- cbind(df.pi,df.mu,df.sd)
  
  list("result"=m.step,
       "logLike"=logLike.c,
       "theta"=df.theta)
}


ggplot(faithful, aes(x=waiting, y=eruptions)) + geom_point()

p <- ggplot(faithful) + 
  geom_histogram(aes(x=eruptions, y=..density..),
                 bins=30, fill="cyan", color="black")
p + geom_density(data=faithful,aes(x=eruptions), color="red")

x <- faithful$eruptions

# number of modes
nModes <- 5
# TODO: use AIC to compute optimal number of modes

theta <- init(x,nModes)
print(theta)

m <- 2000

k <- iterate(x,nModes,m)

print(k$result)
print(k$logLike)
print(k$theta)

df <- k$theta

#logLike.mu <- sum(dnorm(x,theta$mu))
#max.logL <- optim(1,logLike,x,control=list(fnscale=-1))

#for(i in 1:70){logL[i]<-logLike(mus[i],df$x)}

# f(x; pi, mu, s2) = sum(pi * fi(x; mu_i, sqrt(s2_i))


