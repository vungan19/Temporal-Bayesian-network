library(bnlearn)
  #data(coronary)
  #bn_df<- data.frame(coronary)
  #res <- hc(bn_df)
  #res$arcs <- res$arcs[-which((res$arcs[,'from'] == "M..Work" & res$arcs[,'to'] == "Family")),]

newdata <- c() 

## Create node
for ( i in 1: length(AWSUM_feature[,1])) {
  for ( j in 1:(length(AWSUM_feature[1,])-5)){
    add <- AWSUM_feature[i,j:(j+5)]
    names(add)<- c("y_t-5","y_t-4","y_t-3","y_t-2","y_t-1","y_t")
    newdata<- rbind(newdata,add)
  }
}
data <- newdata
samplesize <- 0.60 * nrow(data)
# Random sampling
set.seed(80)
index <- sample( seq_len ( nrow ( data ) ), size = samplesize )
# Create training and test set
datatrain <- data[ index, ]
datatest <- data[ -index, ]

res<- hc(datatrain)  # Learning the structure
plot(res)          ## Plot the learning structure

add_arc <- data.frame("from"=c( "y_t-5","y_t-5","y_t-4","y_t-4","y_t-3","y_t-2","y_t-2"),"to"= c("y_t-3","y_t-2","y_t-1","y_t","y_t","y_t-1","y_t"))
plot(res)

fittedbn <- bn.fit(res, data = datatrain)
coef_t<- fittedbn[["y_t"]][["coefficients"]]
y_t_train <- coef_t[1] + coef_t[2]*datatest[["y_t-5"]]+coef_t[3]*datatest[["y_t-3"]] +coef_t[4]*datatest[["y_t-1"]]

summary(fittedbn)

plot(datatest[["y_t"]],y_t_train,col="blue",pch=16,y_lab=" Predicted value at next window",xlab="real value")
abline(0,1)


plot(datatest[["y_t"]],type="l")
lines(y_t_train,col="red")

plot(newdata[["y_t"]],type="p")
points(y_t_new,col="red")
##=====================================================================================================================##
## Learning Bayesian Stucture [12/12/2018]
 # Call library:
library(bnstruct)
data <- AWSUM_feature
## Create a BNDataset
dataset<- BNDataset(data = data,
                               discreteness = rep('c',6),
                               variables = c("y_t-5","y_t-4","y_t-3","y_t-2","y_t-1","y_t"),
                               node.sizes = c(5,5,5,5,5,5))
net.1 <- learn.network(dataset,initial.network = "random.chain")
net.2 <- learn.network(dataset,initial.network = net.1)
net.3 <- learn.network(dataset,algo="mmhc",initial.network = net.2)

net.1 <- learn.network(dataset,
                       algo = "sem",
                       scoring.func = "AIC")
