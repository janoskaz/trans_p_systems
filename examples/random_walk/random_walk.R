### RANDOM WALK ###

source('../signaledTranspPSystems.R')

# probabilities in both directions
p <- 0.5
q <- 1-p

# membrane 1
lrs1 <- createRLS(label=1,object=createMultiset('a'))
rrs1a <- createRRS(label=2,object=createMultiset('a'))
r1a <- createRule(lrs=lrs1,rrs=rrs1a,t=1,p=p) # probability p in one direction
rrs1b <- createRRS(label=3,object=createMultiset('a'))
r1b <- createRule(lrs=lrs1,rrs=rrs1b,t=1,p=q) # probability q in other direction

m1 <- createMembrane(memblabel=1,
                     objects=new('multisets',objects=list(createMultiset('a'))),
                     rules=new('rules',rules=list(r1a,r1b)))

# membrane 2
lrs2 <- createRLS(label=2,object=createMultiset('a'))
r2a <- createRule(lrs=lrs2,rrs=rrs1b,t=1,p=p) # probability p in one direction
rrs2b <- createRRS(label=1,object=createMultiset('a'))
r2b <- createRule(lrs=lrs2,rrs=rrs2b,t=1,p=q) # probability q in other direction

m2 <- createMembrane(memblabel=2,
                     objects=new('multisets',objects=list()),
                     rules=new('rules',rules=list(r2a,r2b)))

# membrane 3
lrs3 <- createRLS(label=3,object=createMultiset('a'))
r3a <- createRule(lrs=lrs3,rrs=rrs2b,t=1,p=p) # probability p in one direction
r3b <- createRule(lrs=lrs3,rrs=rrs1a,t=1,p=q) # probability q in other direction

m3 <- createMembrane(memblabel=3,
                     objects=new('multisets',objects=list()),
                     rules=new('rules',rules=list(r3a,r3b)))

synapses = list(createSynapse(from=1,to=2),
                createSynapse(from=1,to=3),
                createSynapse(from=2,to=1),
                createSynapse(from=2,to=3),
                createSynapse(from=3,to=1),
                createSynapse(from=3,to=2))

# membrane system
membraneSystem <- createMembraneStructure(
  membranes=list(m1,m2,m3),
  synapses=synapses,
  evts=list(),
  pendingRules=list())


# simulate
sim <- simulate(membraneSystem,20)

# repeated simulation
results <- list()

for (j in 1:10){
  # simulation
  sim <- simulate(membraneSystem,100)
  # results as matrix
  M <- t(simplify2array(lapply(sim,function(y){
    unlist(lapply(y,function(x){return(length(x[[2]]))}))
  })))
  #calculate path
  
  currentPosition <- 1
  currentNumber <- 0
  res <- c(0,0)
  for (i in 1:nrow(M)){
    w <- which(M[i,] == 1)
    if (w == 1){
      if (currentPosition == 2){
        currentNumber <- currentNumber-1
      } else {
        currentNumber <- currentNumber+1
      }
    } else if (w == 2){
      if (currentPosition == 3){
        currentNumber <- currentNumber-1
      } else {
        currentNumber <- currentNumber+1
      }
    } else {
      if (currentPosition == 1){
        currentNumber <- currentNumber-1
      } else {
        currentNumber <- currentNumber+1
      }
    }
    currentPosition <- w
    res <-rbind(res,c(i,currentNumber))
  }
  results[[length(results)+1]] <- res
}


plot(1,1,type='n',xlim=c(0,100),ylim=c(-35,35),xlab='krok',ylab='pozice')
for (i in 1:10){
  lines(results[[i]],col=colors()[sample(1:658,size=1)])
}
lines(0:100,rep(0,101),lwd=2,lty=2)
title('10 náhodných procházek, p=0.50')
