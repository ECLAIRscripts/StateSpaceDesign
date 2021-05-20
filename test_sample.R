  library("ncdf4")

#parameters
num_samples <- 200;

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

sel_ind <- which(aclc_flat >= 0.1 & nummean >= 10 * 10^6)  #constraints given by other colleagues,

sel_ind <- sample(sel_ind,2000); #otherwise the constraint testing takes too long...

pblh <- pblh[sel_ind]
tinv <- tinv[sel_ind]
clw_max <- clw_max[sel_ind]
tmean <- tmean[sel_ind]
nummean <- nummean[sel_ind]
h2o_inv <- h2o_inv[sel_ind]
h2o_mean <- h2o_mean[sel_ind]



#"q_inv","tpot_inv","clw_max","tpot_pbl","pblh","num_pbl"
IN_MAT <- as.matrix(cbind(h2o_inv,tinv,clw_max,tmean,pblh,nummean));



IN_MAT[,5] <- IN_MAT[,5] * 1000;
IN_MAT[,6] <- IN_MAT[,6] * 10^(-6);

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
  c8<--0.3704404e-13
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

CONS_LIST <- vector(mode = "list");
j <- 1;
for(i in 1:nrow(IN_MAT)){
  print(paste0("Now checking row ",i));
  val <- IN_MAT[i,];
  #print(val)
  qpbl <- comp_qpb(tpotpbl = val[4],lcw = val[3],pblh = val[5]);
  permissable <- TRUE;
  
  if(!(1 < val[1] & val[1] < qpbl)) permissable <- FALSE;
  if(!(1 < val[2] & val[2] < 15))                         permissable <- FALSE;
  if(!(80 < val[5] & val[5] < 3000))                      permissable <- FALSE;
  if(!(val[3] > 0.1))                                     permissable <- FALSE;
  if(!(is_cloudy(qpbl,val[4],val[5])))                  permissable <- FALSE;
  if(permissable){
    CONS_LIST[[j]] <- val;
    j <- j + 1;
  }
  print(paste0("Permissable: ",permissable));
}


IN_MAT <- do.call(rbind,CONS_LIST);

reg_index <- which(IN_MAT[,3] <= 0.2 & IN_MAT[,6] <= 75)

min(IN_MAT[,1])

reg_count <- length(reg_index);
tot_count <- length(IN_MAT[,1]);

reg_ratio <- reg_count / tot_count;

reg_sample_count <- num_samples * reg_ratio ;

reg_sample_ind <- sample(1:reg_count,reg_sample_count);

IN_REG <- IN_MAT[reg_index,];

min(IN_REG[,1])

reg_samples <-IN_REG[reg_sample_ind,];

IN_OTH <- IN_MAT[-reg_index,];

other_sample_ind <- sample(1:nrow(IN_OTH),num_samples - reg_sample_count);


min(IN_OTH[,1])

other_samples <- IN_OTH[other_sample_ind,]
min(other_samples[,1])


#the satisfied points from sampling.
all_samples <- rbind(reg_samples,other_samples)

colnames(all_samples) <- c("q_inv","tpot_inv","clw_max","tpot_pbl","pblh","num_pbl")

for(i in 1:6){
  ind <- which.min(all_samples[i:nrow(all_samples),i]);
  ind <- ind + i - 1;
  tmp <- all_samples[i,];
  all_samples[i,] <- all_samples[ind,]
  all_samples[ind,] <- tmp
}

for(i in 1:6){
  ind <- which.max(all_samples[(i+6):nrow(all_samples),i]);
  ind <- ind + i - 1 + 6;
  tmp <- all_samples[i+6,];
  all_samples[i+6,] <- all_samples[ind,]
  all_samples[ind,] <- tmp
}

hist(h2o_inv)
hist(all_samples[,1])

min(IN_MAT[,1])

min(all_samples[,1])

write.csv(all_samples,file = "validation.csv")



