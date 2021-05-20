# extend design   

#generate 120 points of design, for the 90 points find out 90 points which are close to them and then replace them with the original.

# This code is for adding new 30/n points to the existing points, note that the existing points are from the orignial design. 

small_design <- read.csv("design_small.csv",header=TRUE)
small_design <- small_design[,-1]

big_design <- read.csv("BSP_design.csv",header=TRUE)
big_design <- big_design[,-1]

weights <- apply(big_design,2,max) - apply(big_design,2,min)
weights <- 1.0 / weights

euc.dist <- function(x1, x2) sqrt(sum((weights * (x1 - x2)) ^ 2))

extended_design <- big_design

assigned_ind <- c()

for (i in 1:nrow(small_design))
{
  vec = small_design[i,]
  dist <- apply(big_design,1,euc.dist,x2=vec)
  ind <- which.min(dist)
  extended_design[i,]<-small_design[i,]
  big_design[ind,] = rep(-999,9)
  assigned_ind <- c(assigned_ind,ind)
}
extended_design[91:135,] <- big_design[-assigned_ind,]  # new 30 points adding into the design

write.csv(extended_design,file = "BSP_design_extended.csv")