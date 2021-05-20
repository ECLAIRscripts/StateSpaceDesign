

sample.ini(120,3)
 
  sample.int(120,3)
 
  hist(sample.int(120,3), breaks=12)
 
 hist(sample.int(1200,30), breaks=12)
 
  hist(sample.int(120000,3000),breaks=12)
  hist(sample.int(1200000,30000),breaks=12)
  
  
  hist(sample.int(365,9), breaks=365)
  
  hist(sample.int(36500,900),breaks=365)
  hist(sample.int(36500000,90000),breaks=365)
  
  par(mfrow=c(2,2)) 
  hist(sample.int(5898405,150000),breaks=12, main='1mon')  
  hist(sample.int(5898405,150000),breaks=365, main='1day')  
  hist(sample.int(5898405,150000),breaks=365*24, main='24h')  
  hist(sample.int(5898405,150000),breaks=365*12, main='12h')
  dev.new()
  dev.off
  plot.new()
  
  
  par(mfrow=c(1,2)) 
  hist(sample.int(5898405,150000),breaks=365*6, main='6h') 
  hist(sample.int(length(SAMPLES),num_samples), breaks=num_samples,main='1 space')
  
  2136*0.25
  
  hist(sample.int(2136,534),breaks=12, main='6h') 
  
  
  A=matrix(runif(3*5),3,5)
  
  BSP <- vector(mode = "list");
  BSP[[1]] <- (1:nrow(A));
  #print(BSP[[1]])
  for(i in 1:ncol(A)){
    nBSP <- vector(mode = "list");
    for(j in 1:length(BSP))
    {
      print(A[BSP[[j]]])
      print(c(j,i))
      print(A[BSP[[j]],i])
      nBSP[[2*(j-1) + 1]] <- as.vector(BSP[[j]][which(IN_MAT[BSP[[j]],i] < median(A[BSP[[j]],i]))]);
      nBSP[[2*(j-1) + 2]] <- as.vector(BSP[[j]][which(IN_MAT[BSP[[j]],i] > median(A[BSP[[j]],i]))]);
    }
    BSP <- nBSP;
  }
  
  
 