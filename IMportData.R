rm(list=ls()) # clear all

  ## Call the library :
  library(readr)
  library(ggplot2)
  library(tidyverse)
  library(dplyr)
  library(signal)
  

data.direct <- file.path("D:/THESIS/Data/CSV/Raw_ECG")

#================ Import multi-CSV file from the identified direct ==============#
  data.filenames <- list.files(data.direct,pattern= "*.csv")
  data.filenames <-str_replace(data.filenames,".csv","") 
  for (i in 1:length(data.filenames))assign(data.filenames[i], 
                                            read.csv(paste(data.direct,"/",data.filenames[i],".csv",sep=""),skip=1,header = TRUE,nrow=-2,col.names=c("EslapeTime","ECG0","ECG1")))
  
  rm(i)
#=================== PREPROCESSING the imported signal =====================#
 patient <- "n01" # Pick up Specific data
 data <- get(patient)
 bf<- butter(5,c(0.005,0.3),type=c("pass")) # Create the butterword band-pass filter order 5th in the range [0.5,30]/100 Hz
 
 # Apply the filter in to the raw signal
   data[,"ECG0"]<- filtfilt(bf,data[[2]])
   data[,"ECG1"] <- filtfilt(bf,data[[3]])
   
 # The new data after processing with pattern = data_pro   
 assign(paste(patient,"_pro",sep=""),data)
 