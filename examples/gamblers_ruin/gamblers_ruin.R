source('../signaledTranspPSystems.R')

gambler <- function(p,a,b){
  q <- 1-p
  r <- p/q
  n <- a+b
  if (r!=1){
    return((1-(r^a))/(1-(r^n)))
  } else {
    return(a/n)
  }
}

prepareSystem <- function(a,b,p,nveh){
  # number of membranes
  n <- a+b
  # create membranes
  mlist <- list()
  for (i in 0:n){
    if (i == 0 || i == n){
      # endpoints - no rules
      m <- createMembrane(i
                          ,objects=new('multisets',objects=list())
                          ,rules=new('rules',rules=list()))
      mlist[[length(mlist)+1]] <- m
    } else {
      # midpoints - send object a in both directions
      lrs <- createRLS(label=i,object=createMultiset('a'))
      rrs1 <- createRRS(label=i+1,object=createMultiset('a'))
      rrs2 <- createRRS(label=i-1,object=createMultiset('a'))
      r1 <- createRule(lrs,rrs1,t=1,p=p)
      r2 <- createRule(lrs,rrs2,t=1,p=1-p)
      m <- createMembrane(i
                          ,objects=new('multisets',objects=list())
                          ,rules=new('rules',rules=list(r1,r2)))
      mlist[[length(mlist)+1]] <- m
    }
  }
  # input objects into membrane
  olist <- list()
  for (j in 1:nveh){
    olist[[length(olist)+1]] <- createMultiset('a')
  }
  mlist[[a+1]]@objects@objects <- olist
  # create membrane system
  membSystem <- createMembraneStructure(membranes=mlist
                                        ,synapses=list()
                                        ,pendingRules=list()
                                        ,pendingSynapses=list()
                                        ,evts=list())
  # return
  return(membSystem)
}

# calculate distribution of length of games
gameLength <- function(m){
  lengths <- c()
  currentN <- 0
  for(i in 1:nrow(m)){
    N <- m[i,1]+m[i,ncol(m)]
    if (currentN != N){
      lengths <- c(lengths,rep(i,N-currentN))
      currentN <- N
    }
  }
  return(lengths)
}

# a - initial money of player 1, n - totalmoney
expectedGameLength <- function(a,n){
  return(a*(n-a))
}

# different probabilities
res <- list()
probs =  c(.1,.2,.3,.4,.45,.49,.5,.55,.6,.7)
for (pi in probs){
  print(pi)
  ms1 <- prepareSystem(5,10,pi,1000)
  sim <- simulate(ms1)
  M <- t(simplify2array(lapply(sim,function(y){
    unlist(lapply(y,function(x){return(length(x[[2]]))}))
  })))
  res[[length(res)+1]] <- M[,1]
}

plot(1,1,type='n',xlim=c(0,400),ylim=c(0,1000),xlab='krok',ylab='hráčů ve hře')
for (i in c(4,5,6,8,9,10)){
  lines(c(0:length(res[[i]])),c(1000,1000-res[[i]]),lwd=2,lty=2)
  text(350,min(1000-res[[i]]),paste("p=",sprintf("%.2f", probs[i]),sep=""))
}
lines(c(0:length(res[[7]])),c(1000,1000-res[[7]]),lwd=2,lty=2)
text(350,375,paste("p=0.50",sep=""))
title('ruinování hráče, 1000 běhů simulace')



# different amount
res2 <- list()
for (i in 1:5){
  ms1 <- prepareSystem(1*i,2*i,0.5,100)
  sim <- simulate(ms1)
  M <- t(simplify2array(lapply(sim,function(y){
    unlist(lapply(y,function(x){return(length(x[[2]]))}))
  })))
  res2[[length(res2)+1]] <- M[,1]
}

plot(1,1,type='n',xlim=c(0,300),ylim=c(0,1000),xlab='krok',ylab='hráčů ve hře')
for (i in 1:5){
  lines(c(0:length(res2[[i]])),c(1000,1000-res2[[i]]),lwd=2,lty=2)
}
text(20,270,"a=1")
text(70,270,"a=2")
text(105,270,"a=3")
text(210,270,"a=4")
text(290,270,"a=5")
title('rychlost ruinování hráče')


# length of simulation

# different amount
res3 <- list()
for (i in 1:5){
  ms1 <- prepareSystem(1*i,2*i,0.5,1000)
  sim <- simulate(ms1)
  M <- t(simplify2array(lapply(sim,function(y){
    unlist(lapply(y,function(x){return(length(x[[2]]))}))
  })))
  res3[[length(res3)+1]] <- gameLength(M)
}

# test of expected game duration with
ms_expected <- prepareSystem(1,99,0.5,1000)
sim <- simulate(ms_expected)
M <- t(simplify2array(lapply(sim,function(y){
  unlist(lapply(y,function(x){return(length(x[[2]]))}))
})))
gl <- gameLength(M)
