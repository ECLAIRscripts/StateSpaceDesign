design <- read.csv("design.csv",header=TRUE)
design <- big_design[,-1]

old <- design[1:75,]
new <- design[76:100,]



ind1=5
ind2=9

plot(old[,ind1],old[,ind2],type = 'p',xlim = c(min(design[,ind1]),max(design[,ind1])),ylim=c(min(design[,ind2]),max(design[,ind2])))
points(new[,ind1],new[,ind2],col='red',type='p')
