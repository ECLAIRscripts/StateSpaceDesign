 

library(parallel)
library(MASS)

starts <- rep(100,40)

fx <- function(nstart) kmeans(Boston, 4,nstart=nstart)
numCores <- detectCores()
numCores

system.time(
  results <- lapply(starts, fx)
)

system.time(
  results <- mclapply(starts, fx, mc.cores = numCores)
)


library(foreach)
library(iterators)
library(doParallel)

registerDoParallel(numCores)
foreach (i=1:3, .combine=rbind) %dopar% {
  sqrt(i)
}

stopImplicitCluster()