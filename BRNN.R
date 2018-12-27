## Call the packages:
library(tidyverse)
library(Formula)
library(brnn)

## Load the AWSUM_feature
AWSUM_feature <- read.csv("D:/THESIS/Data/AWSUM_feature.csv", header=FALSE)


newdata <- c() 
## Create node
  for ( i in 1: length(AWSUM_feature[,1])) {
    for ( j in 1:(length(AWSUM_feature[1,])-5)){
      add <- AWSUM_feature[i,j:(j+5)]
      names(add)<- c("y_1","y_2","y_3","y_4","y_5","y_t")
      newdata<- rbind(newdata,add)
    }
  }
  rm(i,j)
data<- newdata
## INput and Output
  X<- as.matrix(newdata[1:1000,1:5])
  Y <- newdata[1:1000,6]
## Apllied the bayesian network
 
   bn<-brnn(X,Y,neurons=1,normal_size=TRUE,epoch=1000)

## Prediction: 
  Xtest <- as.matrix(newdata[1:16,1:5])
  Ytest <- newdata[1:16,6]
  y_predict <- predict.brnn(bn,Xtest)
  plot(Ytest[350:550],type="l")
  lines(y_predict[350:550],col="red")
  err<- Ytest - y_predict
  
  
## [28112018] - Learning Neural network
# Call library
  library(neuralnet)
  library(tidyverse)
  library(dplyr)
  
  data <- read.csv("C:/Users/Asus/Downloads/cereal.csv",header=T)
  data <- select(data,calories,protein,fat,sodium,fiber,rating)
 
  
  samplesize <- 0.60 * nrow(data)
# Random sampling
  set.seed(80)
  index <- sample( seq_len ( nrow ( data ) ), size = samplesize )
# Create training and test set
  datatrain <- data[ index, ]
  datatest <- data[ -index, ]
  
  datatest <- data[ 85:(85+17),]
  

## Scale data for neural network
  max <- apply(data , 2 , max)
  min <- apply(data, 2 , min)
  scaled <- as.data.frame(scale(data, center = min, scale = max - min))
  
  trainNN <- scaled[index , ]
  
#testNN <- scaled[-index , ]
 testNN <- scaled[85:(85+17),]
 
 
  set.seed(2)
  NN = neuralnet(y_t ~ y_1 + y_2  +y_3 + y_4 + y_5, trainNN, hidden = 3 , act.fct="tanh",linear.output = T )
  plot(NN)
## Prediction using neural network
  predict_testNN <- neuralnet::compute(NN, testNN[,c(1:5)])
  predict_testNN <- (predict_testNN$net.result * (max(data$y_t) - min(data$y_t))) + min(data$y_t)
  
  plot(predict_testNN,col="red",type="p")
  points(datatest$y_t, col='blue')
    
  
## Calculate Root Mean Square Error (RMSE)
  RMSE.NN <- (sum((datatest$y_t - predict_testNN)^2) / nrow(datatest)) ^ 0.5  