rm(list=ls()) # clear all

## Call the library :
  library(readr)
  library(ggplot2)
  library(tidyverse)
  library(dplyr)
  library(signal)
#=======================================================================================
  ##This function is used to extract features of ECG
  ## Input:
    ##   data:   ECG signal
  ##   fs:     sample frequency
  ##Output:
    ##   R_value, R_loc: value and location of R peak
  ##   Q_value, Q_loc: value and location of Q peak
  ##   S_value, S_loc: value and location of S peak
  ##   J_value, J_loc: value and location of J peak
  ##   T_value, T_loc: value and location of T peak
  ##   P_value, P_loc: value and location of P peak
  ##   RR:             RR interval
  ##   PR:             PR segment
  ##   QT:             QT segment
  ##   BPM:            Beat per min
  ##   tqrs:           duration of QRS complex
  ##   trr:            duration of RR interval
  ##   tpr:            duration of PR segment
  ##   tqt:            duration of QT segment

#=======================================================================================
vus<- read.csv("D:/THESIS/vus.csv",header = FALSE,col.names=c("V0","V1"))
dataname <- "vus"
channel <- "V0"
fs <- 128
ecg_feature_extraction(dataname,channel,fs)

#============================= filter baseline =========================================
ecg_feature_extraction <- function(dataname,channel,fs) { 
  d <- get(dataname)
  d <- d[[channel]]
  t <-(1:length(d)-1)/fs
    # bf_1<- butter(5,c(0.5,30)/(fs/2),type=c("pass")) # Trial : Using bandpass filter with butterworth filter 5th, range(0.5,30)
    # c_1<-filtfilt(bf_1,d)
  
  bf<- butter(3,0.5/(fs/2),type=c("low"))  # filter the basline
  base <- filtfilt(bf,d) 
  d <- d-base                             # delete the baseline in the signal
  
  bf<- butter(6,30/(fs/2),type=c("low"))   # filter noise
  c<- filtfilt(bf,d)
  c<- c/max(c)
  c<- c - mean(c)

#============================ Make impulse response =================================== 

h <- c(-1,-2,0,2,1)/8 
 ## Cach 1:
    int_c <- (5-1)/(fs * 1/40)
    b <- interp1(1:5,c(1,2,0,-2,-1)*(1/8)*fs,seq(1,5,by=int_c))
    filt <- Arma(b,1)
    x4<- filtfilt(filt,c)
    x4 <- x4/max(x4)   # normalize the result
    
 ## Cach 2: 
    #x4 <- conv(c,h)
    #x4 <- x4[(1: length(d))+2]
    #x4 <- x4/max(x4)   # normalize the result
#=============================== Squaring ============================================
  x5 <- x4**2
  x5 <- x5/max( abs(x5))  
  
#======================= Make impulse response ========================== 
  x6 <- conv(x5 ,seq(1,1,length.out=0.15*fs)/round(0.150*fs));
  ## Apply filter
  delay<- round((0.15*fs)/2) 
  x6<- x6[(1:length(d))+delay]
  x6  <- x6 / max(x6)
  
  thresh <- mean(x6)
  poss_reg <- as.numeric(x6>thresh)
  
  left <- which(diff(c(0,poss_reg))==1)
  right <-  which(diff(c(poss_reg,0))==-1)
  
  R_value <- seq(0,0,length.out=length(left))
  R_loc <-  seq(0,0,length.out=length(left))
  need_to_remove <-  seq(0,0,length.out=length(left))
  j <- 1
  
## ============= R peak detection and rejection the fail peak =============
    # ======= detect R peak ======
 for ( i in 1: length(left)){
   R_value[i]<- max(c[left[i]:right[i]])
   R_loc[i] <-  which(c[left[i]:right[i]]==  R_value[i])
   
   #R_value[i] <- max(d_raw[left[i]:right[i]]) # The amplitude in the raw-data
   #R_loc[i] <-  which(c[left[i]:right[i]]==  R_value[i]) 
   R_loc[i] <- R_loc[i] -1 + left[i]  ## add offset
   # R_loc_time[i] <- t[R_loc[i]]
 }
    # ====== reject fail peak ======
  for ( i in 5:length(R_loc)){
    value1 <- (R_loc[i-2]-R_loc[i-3])/(R_loc[i-3]-R_loc[i-4])
    value2 <- (R_loc[i-1]-R_loc[i-3])/(R_loc[i-3]-R_loc[i-4])
    value3 <- (R_loc[i]-R_loc[i-1])/(R_loc[i-3]-R_loc[i-4])
    if(value1 < 0.5 & abs(value2-value3)<0.75) {
      need_to_remove[j] <- i
      j<-j+1
    }
  }
  
  need_to_remove <- need_to_remove[need_to_remove != 0]
   if ( length(need_to_remove) != 0 ){
     R_value<- R_value[-need_to_remove]
     R_loc<- R_loc[-need_to_remove]
     left<- left[-need_to_remove]
     right<-right[-need_to_remove]
   }
##=================== Detect other peaks and interval ====================  
   # ====== Detect peaks: Q,S,J,K,P,T ========= #
  Q_value<- seq(0,0,length.out=length(left))
  Q_loc<- seq(0,0,length.out=length(left))
  
  S_value<- seq(0,0,length.out=length(left))
  S_loc<- seq(0,0,length.out=length(left))
  
  J_value<- seq(0,0,length.out=length(left))
  J_loc<- seq(0,0,length.out=length(left))
  
  K_value<- seq(0,0,length.out=length(left))
  K_loc<- seq(0,0,length.out=length(left))
  
  P_value<- seq(0,0,length.out=length(left)-1)
  P_loc<- seq(0,0,length.out=length(left)-1)
  
  T_value<- seq(0,0,length.out=length(left)-1)
  T_loc<- seq(0,0,length.out=length(left)-1)
  
  # ====== Determined other peaks: Q,S,J,K,P,T ====== #
  RR <- seq(0,0,length.out=length(left)-1)
  trr<- seq(0,0,length.out=length(left)-1)
  HRV <- seq(0,0,length.out=length(left)-1)
  QR <- seq(0,0,length.out=length(left)-1)
  QR_loc<- seq(0,0,length.out=length(left)-1)
  tqr<- seq(0,0,length.out=length(left)-1)
  
  for ( i in 1: length(R_loc)){
    
    # Q_peak: which the value is min in the left_site of R_peak
    Q_value[i]<- min(c[left[i]:R_loc[i]]) 
    Q_loc[i] <-  which(c[left[i]:R_loc[i]] ==  min(c[left[i]:R_loc[i] ]))
    Q_loc[i] <- Q_loc[i] -1 + left[i]  ## add offset 

    # S_peak: which the value is min in the right_site of R_peak
    S_value[i]<- min(c[R_loc[i]: right[i]])
    S_loc[i] <-  which(c[R_loc[i]:right[i]]==S_value[i])
    S_loc[i] <- S_loc[i] + R_loc[i] -1  ## add offset 
    
    # J point: is the off set of QRS_complex
    J_loc[i] <- right[i]
    J_value[i] <-c[J_loc[i]]
    
    # K_point: 
    K_loc[i] <- left[i]
    K_value[i] <-c[K_loc[i]]
    # QR amplitude :
    QR[i] <-abs(R_value[i] - Q_value[i])
    
  if ( i != 1){
    # ====== RR interval ====== 
    RR[i-1]<- R_loc[i] -R_loc[i-1]
    trr[i-1] <- RR[i-1]/fs
    
    # ====== BPM (vent rate) ======
    HRV[i-1] <-  60/trr[i-1]
    
    # ====== T peak ====== 
    intv <-floor(R_loc[i-1] + (0.15*RR[i-1])) : floor(R_loc[i-1] + (0.55*RR[i-1])) 
    T_value[i-1] <- max(c[intv ])
    T_loc[i-1] <- which(c[intv] ==T_value[i-1] )
    T_loc[i-1] <- T_loc[i-1]+ R_loc[i-1] + floor(0.15*RR[i-1]) - 1
    
    # ====== P peak ====== 
    intv <- floor(left[i] - 0.15*RR[i-1]): Q_loc[i]
    P_value[i-1] <- max(c[intv])
    P_loc[i-1] <-  which(c[intv]==P_value[i-1])
    P_loc[i-1] <- P_loc[i-1] + round(left[i] - 0.15*RR[i-1]) - 1  # add offset
    # ====== QR_interval:P wave region ====== 
    QR_loc[i-1]<- R_loc[i] - Q_loc[i-1]
    tqr[i-1]<- QR_loc[i-1]/fs
  }
    
  }
  data_feature <- list(R_loc,R_value,Q_loc,Q_value,S_loc,S_value,J_loc,J_value,K_loc,K_value,P_loc,P_value,T_loc,T_value,QR,RR,trr,HRV,QR_loc,tqr)
  names(data_feature) <- c("R_loc","R_value","Q_loc","Q_value","S_loc","S_value","J_loc","J_value","K_loc","K_value",
                           "P_loc","P_value","T_loc","T_value","QR","RR","trr","HRV","QR_loc","tqr")
  assign(paste(dataname,"_",channel,"_features",sep=""),data_feature,envir =  .GlobalEnv)
  assign(paste(dataname,"_",channel,"_pro",sep=""),c,envir =  .GlobalEnv)
  assign("t",t,envir =.GlobalEnv)
  }

  
  ##=============Visuallize data with marked R-peak ================== ##
 ## Take the sample from 1:500 with fs = 128 Hz
    R_loc <- vus_V0_features[["R_loc"]]
    R_loc_time<- R_loc[R_loc<500]
    
    x <- t[1:500]
    y <- vus_V0_pro[1:500]
    
    x_1 <- t[R_loc_time]
    y_1<- vus_V0_pro[ R_loc_time]
    
    
    P_loc <- vus_V0_features[["P_loc"]]
    P_loc_time<- P_loc[P_loc< 500 ]
    x_2 <- t[P_loc[P_loc<500]]
    y_2 <- vus_V0_pro[P_loc_time]
    
    T_loc <- vus_V0_features[["T_loc"]]
    T_loc_time<- T_loc[T_loc< 500 ]
    x_3 <- t[T_loc[T_loc<500]]
    y_3 <- vus_V0_pro[T_loc_time]
    
    
    plot(x,y,type="l") 
    points(x_1,y_1,col="green")  # R-peak
    points(x_2,y_2,col="red")    # P-peak
    points(x_3,y_3,col="yellow")


##=========================== P-wave variation ======================================= #
    # P_fearure <- function(data){
  #======== P-onset & P_offset ========#
    P_loc<- vus_V0_features[["P_loc"]]  # Take P_location from the list feature extraction
    d <- round(fs * 0.075)              # Create the window with 75ms in 2 sides of P location
    P_bef<- P_loc - d   # before P-peak a interval d
    P_aft <- P_loc + d  # after P-peak a interval d
    P_on <- seq(0,0,length.out=length(P_loc))
    P_off <- seq(0,0,length.out=length(P_loc))
   
    T_bef<- T_loc - 10  # before P-peak a interval d
    T_aft <- T_loc + 10
    
    
    
     w<- c()
     for ( i in 1: length(P_loc)){
       PP_1  <- vus_V0_pro[P_loc[i]:P_bef[i]]
       PP_2  <- vus_V0_pro[P_loc[i]:P_aft[i]]
       loc_on <- which(diff(sign(diff(PP_1)))== 2) + 1
       loc_off<- which(diff(sign(diff(PP_2)))== 2) + 1
       P_on[i] <- P_bef[i] + loc_on      # The location of P_onset
       P_off[i] <- P_loc[i]+ loc_off     # The location of P_offset
       #names(P) <- c("y-on",paste("y-",as.character(2:(length(P)-1)),sep=""),"y-off")
       #w <- rbind(w,P)
     }
    
    plot((P_bef[7]):(P_aft[7]),vus_V0_pro[(P_bef[7]):(P_aft[7])],type="l")
    points(P_bef[7],vus_V0_pro[P_bef[7]],col="green")
    points(P_aft[7],vus_V0_pro[P_aft[7]],col="green")
    points(P_loc[7],vus_V0_pro[P_loc[7]],col="red")
    
    plot((T_bef[4]-10):(T_aft[4]+10),vus_V0_pro[(T_bef[4]-10):(T_aft[4]+10)],type="l")
    points(T_bef[4],vus_V0_pro[T_bef[4]],col="green")
    points(T_aft[4],vus_V0_pro[T_aft[4]],col="green")
    points(T_loc[4],vus_V0_pro[T_loc[4]],col="red")
    
    
    #P_duration <- P_offset - P_onset
    #P_duration_t <- P_duration*(1000/128)   # time in ms unit
  #}
  
  #plot(P_onset[2]:P_offset[2],vus[P_onset[2]:P_offset[2],1],type="l")
  #P_loc[2],vus[P_loc[2],1],col="green")
    
    
  