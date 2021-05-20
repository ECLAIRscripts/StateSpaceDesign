install.packages("ncdf4")
library("ncdf4")
setwd("/home/liujia/Documents/Main_project_FMI/ECLAIR_files/State_space_design")
#parameters
num_samples <- 90;

#--------------------------


nc <- nc_open("final_revised_data.nc", write=FALSE)

pblh <- as.vector(ncvar_get( nc, 'pblh', verbose=FALSE ))         # same..
tinv <- as.vector(ncvar_get( nc, 't_inv', verbose=FALSE ))        # tpot_inv
clw_max <- as.vector(ncvar_get( nc, 'clw_max', verbose=FALSE ))   # same..
tmean <- as.vector(ncvar_get( nc, 't_mean', verbose=FALSE ))      # tpot_pbl
nummean <- as.vector(ncvar_get( nc, 'num_mean', verbose=FALSE ))  # num_pbl
h2o_inv <- as.vector(ncvar_get( nc, 'h2o_inv', verbose=FALSE ))   # q_inv
h2o_mean <- as.vector(ncvar_get( nc, 'h2o_mean', verbose=FALSE ))   # q_pbl

#  select approp. subset

aclc_flat <- as.vector(ncvar_get( nc, 'aclc_flat', verbose=FALSE ))

# 
#aclc_flat2<- aclc_flat [seq(1, length(aclc_flat), 20)]
#plot(aclc_flat2)

sel_ind <- which(aclc_flat >= 0.1 & nummean >= 10 * 10^6)
#length(which(aclc_flat >= 0.1 & nummean >= 10 * 10^6))/length(nummean)
#18.55% data points are taken into account?

pblh <- pblh[sel_ind]
tinv <- tinv[sel_ind]
clw_max <- clw_max[sel_ind]
tmean <- tmean[sel_ind]
nummean <- nummean[sel_ind]
h2o_inv <- h2o_inv[sel_ind]
h2o_mean <- h2o_mean[sel_ind]

jpeg("pblh_hist.jpg");
  hist(1000*pblh,breaks = 40,xlab = "pblh in (m)");

dev.off();

    IN_MAT <- as.matrix(cbind(pblh,tinv,clw_max,tmean,nummean,h2o_inv));

IN_MAT[,1] <- IN_MAT[,1] * 1000;
IN_MAT[,5] <- IN_MAT[,5] * 10^(-6);

#  filter for constraints

calc_psat_w <- function(T){
  # Function calculates the saturation vapor pressure (Pa) of liquid water as a function of temperature (K)
  #
  # thrm.f90:  real function rslf(p,t)
  c0<-0.6105851e+03
  c1<-0.4440316e+02
  c2<-0.1430341e+01
  c3<-0.2641412e-01
  c4<-0.2995057e-03
  c5<-0.2031998e-05
  c6<-0.6936113e-08
  c7<-0.2564861e-11
  c8<- -0.3704404e-13
  #
  x<-max(-80.,T-273.16)
  calc_psat_w <- c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))
}

calc_rh <- function(rw,T,press){
  # Calculate RH (%) from water vapor mixing ratio rw (r=m_w/m_air [kg/kg]), temperature (K) and pressure (Pa)
  #
  # r=m_w//m_air=pw/Rm/(pair/R)=pw/(p-pw)*R/Rm => pw=p*r/(R/Rm+r)
  #
  R  <- 287.04	# Specific gas constant for dry air (R_specific=R/M), J/kg/K
  Rm <- 461.5	# Specific gas constant for water
  ep <- R/Rm
  #
  psat    <- calc_psat_w(T)
  calc_rh <- press*rw/(ep+rw)/psat*100
}

calc_sat_mixr <- function(p,T){
  # Function calculates saturation mixing ratio for water (kg/kg)
  #
  # thrm.f90: real function rslf(p,t)
  #
  # r=m_w//m_air
  # R/Rm=287.04/461.5=.622
  #
  esl<-calc_psat_w(T)
  calc_sat_mixr <- 0.622*esl/(p-esl)
}



calc_lwc_altitude <- function(p_surf,theta,rw,zz){
  # Calculate cloud water mixing ratio at a given altitude z (m) when liquid water potential 
  # temperature (theta [k]) and water vapor mixing ratio (rw [kg/kg]) are constants. 
  # Surface pressure p_surf is given in Pa.
  #
  # Constants
  R<-287.04	# Specific gas constant for dry air (R_specific=R/M), J/kg/K
  Rm<-461.5	# -||- for water
  ep2<-Rm/R-1.0 #M_air/M_water-1
  cp<-1005.0	# Specific heat for a constant pressure
  rcp<-R/cp
  cpr<-cp/R
  g<-9.8
  p00<-1.0e+05
  alvl <- 2.5e+06 #  ! latent heat of vaporization
  #
  # a) Integrate to cloud base altitude
  dz<-1.			# 1 m resolution
  z<-0.				# The first altitude
  press<-p_surf	# Start from surface
  RH<-0
  while (z<zz){
    # Temperature (K) 
    tavg<-theta*(press/p00)**rcp
    #
    # Current RH (%)
    RH<-calc_rh(rw,tavg,press)
    if (RH>100) break
    #
    # From z to z+dz
    z <- z +dz
    # Virtual temperature: T_virtual=T*(1+ep2*rl)
    xsi<-(1+ep2*rw)
    # Pressure (Pa)
    press<- press - g*dz*press/(R*tavg*xsi)
  }
  # No cloud or cloud water
  if (RH<100) {
    calc_lwc_altitude <- 0
    return(0)
  }
  #
  # b) Integrate up to given altitude
  while (z<zz){
    # From z to z+dz
    z<- z+dz
    #
    # Moist adiabatic lapse rate
    q_sat <- calc_sat_mixr(press,tavg)
    tavg<- tavg - g*(1+alvl*q_sat/(R*tavg))/(cp+alvl**2*q_sat/(Rm*tavg**2))*dz
    #
    # New pressure
    xsi<-(1+ep2*q_sat)
    press<- press-g*dz*press/(R*tavg*xsi)
  }
  #
  # Return cloud water mixing ratio = totol - vapor
  calc_lwc_altitude <- rw-q_sat
}

q_pbl <- function(p_surf,theta,lwc,zz){
  # Solve total water mixing ratio (rw, kg/kg) from surface pressure (p_surf, Pa), liquid water potential
  # temperature (theta, K) and liquid water mixing ratio (lwc) at altitude zz (m)
  #
  # Constants
  R<-287.04	# Specific gas constant for dry air (R_specific=R/M), J/kg/K
  Rm<-461.5	# -||- for water
  ep2<-Rm/R-1.0 #M_air/M_water-1
  cp<-1005.0	# Specific heat for a constant pressure
  rcp<-R/cp
  cpr<-cp/R
  g<-9.8
  p00<-1.0e+05
  alvl <- 2.5e+06 #  ! latent heat of vaporization
  #
  # Mimimum water vapor mixing ratio is at least lwc
  q_min<-lwc
  #
  # Maximum water vapor mixing ratio is unlimited, but should be smaller
  # than that for a cloud which base is at surface
  t_surf<-theta*(p_surf/p00)**rcp
  q_max<-calc_sat_mixr(p_surf,t_surf)
  #
  k<-0
  while (k<100){
    q_new<-(q_min+q_max)/2
    
    lwc_calc<-calc_lwc_altitude(p_surf,theta,q_new,zz)
    #print(lwc_calc)
    #print(lwc)
    #	
    if (abs(lwc-lwc_calc)<1e-7)
      break
    if (lwc<lwc_calc)
      q_max<-q_new
    else
      q_min<-q_new
    k<- k+1
    # Failed
    if (k==50) q_new <- -999
  }
  q_pbl <- q_new
}

p_surf=101780. # Surface pressure (Pa)

comp_qpb <- function(tpotpbl,lcw,pblh){
  
  
  comp_qpb <- 1000*q_pbl(p_surf,tpotpbl,0.001*lcw,pblh)
  
}

calc_cloud_base <- function(p_surf,theta,rw){
  # Calulate cloud base heigh when liquid water potential temperature (theta [kK) and water
  # vapor mixing ratio (rw [kg/kg]) are constants. Surface pressure p_surf is given in Pa.
  # For more information, see "lifted condensation level" (LCL).
  #
  # Constants
  R <- 287.04	# Specific gas constant for dry air (R_specific=R/M), J/kg/K
  Rm <- 461.5	# -||- for water
  ep2 <- Rm/R - 1.0 #M_air/M_water-1
  cp<-1005.0	# Specific heat for a constant pressure
  rcp<-R/cp
  cpr<-cp/R
  g<-9.8
  p00<-1.0e+05
  #
  # Integrate to cloud base altitude
  dz<-1.			# 1 m resolution
  z<-0.				# The first altitude
  press<-p_surf	# Start from surface
  RH<-0
  while (RH<100 & z<10000){
    # Temperature (K)
    tavg <- theta*(press/p00)**rcp
    #
    # Current RH (%)
    RH<-calc_rh(rw,tavg,press)
    if (RH>100) break
    #
    # From z to z+dz
    z <- z + dz
    # Virtual temperature: T_virtual=T*(1+ep2*rl)
    xsi<-(1+ep2*rw)
    # Pressure (Pa)
    press<- press - g*dz*press/(R*tavg*xsi)
  }
  #
  # No cloud
  if (RH<100) return(-999)
  #
  # Return cloud base altitude
  calc_cloud_base <- 0.001*z
  
}

is_cloudy <- function(q_pblv,t_pbl,pblh){
  p_surf <- 101780.0
  if((pblh - 0.05 > calc_cloud_base(p_surf,t_pbl,q_pblv*1e-3)) & calc_cloud_base(p_surf,t_pbl,q_pblv*1e-3) >0.03) is_cloudy <- TRUE
  else is_cloudy <- FALSE
}


# to speed up the computation, e.g., by parallel!  Large samples, origial is 2000, increase upto 10,000 due the the acceptance rate, currently
# is 20%

library(foreach)
library(iterators)
library(doParallel)

 




CONS_LIST <- vector(mode = "list");
j <- 1;
for(i in 1:nrow(IN_MAT)){
  print(paste0("Now checking row ",i));
  val <- IN_MAT[i,];
  #print(val)
  qpbl <- comp_qpb(tpotpbl = val[4],lcw = val[3],pblh = val[1]);
  permissable <- TRUE;
  
  if(!(1 < val[6] & val[6] < qpbl)) permissable <- FALSE;
  if(!(1 < val[2] & val[2] < 15))                         permissable <- FALSE;
  if(!(80 < val[1] & val[1] < 3000))                      permissable <- FALSE;
  if(!(val[3] > 0.1))                                     permissable <- FALSE;
  if(!(is_cloudy(qpbl,val[4],val[1])))                  permissable <- FALSE;
  if(permissable){
    CONS_LIST[[j]] <- val;
    j <- j + 1;
  }
  print(paste0("Permissable: ",permissable));
}

IN_MAT <- do.call(rbind,CONS_LIST);

#---save data in R s
save.image("myrun.RData")
# -------------------------


# create medians vector (BSP) binary space partition. ROI hyper cube
BSP <- vector(mode = "list");
BSP[[1]] <- (1:nrow(IN_MAT));
#print(BSP[[1]])
for(i in 1:6){
  nBSP <- vector(mode = "list");
  for(j in 1:length(BSP))
  {
    #print(IN_MAT[BSP[[j]],i])
    nBSP[[2*(j-1) + 1]] <- as.vector(BSP[[j]][which(IN_MAT[BSP[[j]],i] < median(IN_MAT[BSP[[j]],i]))]);
    nBSP[[2*(j-1) + 2]] <- as.vector(BSP[[j]][which(IN_MAT[BSP[[j]],i] > median(IN_MAT[BSP[[j]],i]))]);
  }
  BSP <- nBSP;
}

# now sample design , uniform sampling from the partitions.

samp_count <- 0;
SAMPLES <- vector(mode = "list");
for(i in 1:length(BSP)){
  ind <- sample(BSP[[i]],size = 1);
  #print(ind)
  #print(IN_MAT[ind,])
  SAMPLES[[i]] <- IN_MAT[ind,];
  
  #print(IN_MAT[BSP[[i]][ind],])
}



DESIGN <- matrix(nrow = num_samples,ncol = 6);
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


# reorder columns
TMPMAT <- DESIGN;
DESIGN[,1] <- TMPMAT[,6]  
DESIGN[,2] <- TMPMAT[,2]
DESIGN[,3] <- TMPMAT[,3]
DESIGN[,4] <- TMPMAT[,4]
DESIGN[,5] <- TMPMAT[,1]
DESIGN[,6] <- TMPMAT[,5]

#reordering of the rows, first min, then max and then random.
for(i in 1:6){
  ind <- which.min(DESIGN[i:nrow(DESIGN),i]);
  ind <- ind + i - 1;
  tmp <- DESIGN[i,];
  DESIGN[i,] <- DESIGN[ind,]
  DESIGN[ind,] <- tmp
}

for(i in 1:6){
  ind <- which.max(DESIGN[(i+6):nrow(DESIGN),i]);
  ind <- ind + i - 1 + 6;
  tmp <- DESIGN[i+6,];
  DESIGN[i+6,] <- DESIGN[ind,]
  DESIGN[ind,] <- tmp
}

#c("q_inv","tpot_inv","clw_max","tpot_pbl","pblh","num_pbl")
  #colnames(DESIGN) <- c("pblh","tpot_inv","clw_max","tpot_pbl","num_pbl","q_inv");


  colnames(DESIGN) <- c("q_inv","tpot_inv","clw_max","tpot_pbl","pblh","num_pbl")

write.csv(DESIGN,file = "BSP_design.csv")


