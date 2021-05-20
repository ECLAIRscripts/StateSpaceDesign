#parameters

num_samples <- 150;

#--------------------------


allpoints <-read.csv("BSP_design.csv",header=TRUE)
allpoints <- allpoints[-1]



tpot_inv <- as.vector(allpoints[,"tpot_inv"])         
tpot_pbl <- as.vector(allpoints[,"tpot_pbl"])
q_inv <- as.vector(allpoints[,"q_inv"])         
lwp <- as.vector(allpoints[,"lwp"])         
pbl <- as.vector(allpoints[,"pblh"])         
cdnc <- as.vector(allpoints[,"cdnc"])  

IN_MAT <- as.matrix(cbind(q_inv,tpot_inv,lwp,tpot_pbl,pbl,cdnc))


# create medians vector (BSP)
BSP <- vector(mode = "list");
BSP[[1]] <- (1:nrow(IN_MAT));
#print(BSP[[1]])
for(i in 1:ncol(IN_MAT)){
  nBSP <- vector(mode = "list");
  for(j in 1:length(BSP))
  {
    #print(IN_MAT[BSP[[j]],i])
    nBSP[[2*(j-1) + 1]] <- as.vector(BSP[[j]][which(IN_MAT[BSP[[j]],i] < median(IN_MAT[BSP[[j]],i]))]);
    nBSP[[2*(j-1) + 2]] <- as.vector(BSP[[j]][which(IN_MAT[BSP[[j]],i] > median(IN_MAT[BSP[[j]],i]))]);
  }
  BSP <- nBSP;
}

# now sample design

samp_count <- 0;
SAMPLES <- vector(mode = "list");
for(i in 1:length(BSP)){
  ind <- sample(BSP[[i]],size = 1);
  #print(ind)
  #print(IN_MAT[ind,])
  SAMPLES[[i]] <- IN_MAT[ind,];
  
  #print(IN_MAT[BSP[[i]][ind],])
}



DESIGN <- matrix(nrow = num_samples,ncol = ncol(IN_MAT));
if(length(SAMPLES)>=num_samples){
  sel_ind <- sample(1:length(SAMPLES),size = num_samples);
  DES_SAMPLES <- SAMPLES[sel_ind];
  DESIGN <- do.call(rbind,DES_SAMPLES);
}else{
  add_ind <- sample(1:length(BSP),size = (num_samples - length(SAMPLES)),replace = TRUE);
  base_ind <- length(SAMPLES);
  j <- 1;
  for(i in add_ind){
    ind <- sample(BSP[[i]],size = 1);
    #print(ind)
    #print(IN_MAT[ind,])
    SAMPLES[[base_ind + j]] <- IN_MAT[ind,];
    j <- j + 1;
    #print(IN_MAT[BSP[[i]][ind],])
  }
  DESIGN <- do.call(rbind,SAMPLES);
}




for(i in 1:ncol(IN_MAT)){
  ind <- which.min(DESIGN[i:nrow(DESIGN),i]);
  ind <- ind + i - 1;
  tmp <- DESIGN[i,];
  DESIGN[i,] <- DESIGN[ind,]
  DESIGN[ind,] <- tmp
}

for(i in 1:ncol(IN_MAT)){
  ind <- which.max(DESIGN[(i+ncol(IN_MAT)):nrow(DESIGN),i]);
  ind <- ind + i - 1 + ncol(IN_MAT);
  tmp <- DESIGN[i+ncol(IN_MAT),];
  DESIGN[i+ncol(IN_MAT),] <- DESIGN[ind,]
  DESIGN[ind,] <- tmp
}

#c("q_inv","tpot_inv","clw_max","tpot_pbl","pblh","num_pbl")
#colnames(DESIGN) <- c("pblh","tpot_inv","clw_max","tpot_pbl","num_pbl","q_inv");

#IN_MAT <- as.matrix(cbind(q_inv,tpot_inv,lwp,tpot_pbl,pbl,cdnc));

colnames(DESIGN) <- c("q_inv","tpot_inv","lwp","tpot_pbl","pblh","cdnc")

write.csv(DESIGN,file = "BSP_design_subsample.csv")

install.packages('animation')
library(animation)
imgs <- list.files(pattern="*.png")
saveVideo({
  for(img in imgs){
    im <- magick::image_read(img)
    plot(as.raster(im))
  }  
})

