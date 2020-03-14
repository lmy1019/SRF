# features picked by rf
#X = X[1:n,c(2,4,6,43,52,53)], L = 9(0.174), fold = 10
#X = X[1:n,c(2,4,6,43,52,53)], L = 7(0.391), fold = 10
#X = X[1:n,c(2,4,6,43,52,53)], L = 5(0.588), fold = 10
#X = X[1:n,c(2,4,6,43,52,53)], L = 3(0.752), fold = 10

# features picked by cox
#X = X[1:n,id.features[1:5]], L = 9(0.174), fold = 10
#X = X[1:n,id.features[1:5]], L = 7(0.391), fold = 10
#X = X[1:n,id.features[1:5]], L = 5(0.588), fold = 10
#X = X[1:n,id.features[1:5]], L = 3(0.752), fold = 10

# features picked by cox
#X = X[1:n,id.features[1:10]], L = 9(0.174), fold = 10
#X = X[1:n,id.features[1:10]], L = 7(0.391), fold = 10
#X = X[1:n,id.features[1:10]], L = 5(0.588), fold = 10
#X = X[1:n,id.features[1:10]], L = 3(0.752), fold = 10

#X = X[1:n,id.features[1:20]], L = 9(0.174), fold = 10

#X = X[1:n,id.features[1:30]], L = 9(0.174), fold = 10

#X = X[1:n,id.features[1:40]], L = 9(0.174), fold = 10

# 47 43 25  7 46


rm(list=ls())
load("/Users/lmy/Dropbox/Mingyang-HLi/projects/survival_random_forest/realdata/x_submit.Rdata")
library(grf)
setwd("/Users/lmy/Dropbox/Mingyang-HLi/projects/survival_random_forest/code/rscript")
source("naive.estimators.R")
source("util.R")
source("lu.R")
source("wang.R")
library(dummies)
library(emil)





helper<-function(pred, truth, std, weights){
  return( 
    c( mean( (truth<=pred+1.96*std) & (truth>=pred-1.96*std) ), 
       mean(abs(pred-truth)*weights),
       sqrt( mean( (pred-truth)^2 *weights ) ) 
    ) 
  )
}


num_features = 69#69
L=11


data = x_nomiss 
data = data[data[,dim(data)[2]]>=1e-10,]
data = data[complete.cases(data),]
#n = dim(data)[1]
n = dim(data)[1]
p = dim(data)[2]-2
X = cbind(  data[, 2:42] ,  dummy(data$gender, sep = "_")[,-1], data[, 44],dummy(data$nyha, sep = "_")[,-1], data[, 46],dummy(data$race, sep = "_")[,-1],data[, 48:49], dummy(data$ischemic, sep = "_")[,-1],data[, 51])
Z = data[1:n,p+2]+rnorm(n, 0, 1e-5)
delta = data[1:n,p+1]; 

fold = 10;  unit = floor(n/fold);   base = (fold-1)*unit
output=NULL;

# features picked by cox
# surv_object <- Surv(time = Z, event = delta)
# pv = c()
# for( i in 1:dim(X)[2]){
#   cox <- coxph(surv_object ~ ., data = data.frame(X[,i]) )
#   pv = c(pv, summary(cox)$waldtest["pvalue"])
# }
# id.features=order(pv)

#X = X[1:n,c(2,4,6,43,52,53)]
#X = X[1:n,id.features[1:num_features]]


#numeric_feature.idx = setdiff(2:p, c(1, p-1, p-4, p-6, p-8)  )

for( i in 1:fold ){
  

  
  #training
  print(i)
  tr.idx = ( 0:(base-1) + (i-1)*unit ) %% n +1
  X.tr = X[tr.idx, ];  #X.tr <- scale(X.tr)
  
  Z.tr = Z[tr.idx];  delta.tr = delta[tr.idx]
  pv = c()
  surv_object <- Surv(time = Z.tr, event = delta.tr)
  for( l in 1:dim(X)[2]){
    cox <- coxph(surv_object ~ ., data = data.frame(X.tr[,l]) )
    pv = c(pv, summary(cox)$waldtest["pvalue"])
  }
  id.features=order(pv)
  id.features= id.features[1:num_features]
  X.tr = X.tr[,id.features]
  
  
  
  id1=order(Z.tr);   Z.tr=Z.tr[id1];  X.tr=X.tr[id1,];   delta.tr = delta.tr[id1]
  surv_object <- Surv(time = Z.tr, event = 1-delta.tr)
  G = Ghat.new(surv_object, data.frame (X.tr), sum(Z.tr<=L)) 
  Z.tr.res =pmin(Z.tr, L)
  delta.tr.res =delta.tr;  delta.tr.res [Z.tr>=L]=1
  s.forest = custom_forest(X.tr, Z.tr.res , delta.tr.res , G )
  #testing
  te.idx = setdiff(1:n, tr.idx)
  X.test = X[te.idx,  id.features];  #X.test <- scale(X.test)
  Z.test = Z[te.idx];  
  delta.test = delta[te.idx];   
  id2=order(Z.test);   Z.test=Z.test[id2];  X.test=X.test[id2,];   delta.test = delta.test[id2]
  Z.test.res = pmin(Z.test, L)
  delta.test.res = delta.test;   delta.test.res[Z.test>=L]=1;
  #truth.idx = (delta.test.res==1);  n.test = sum(truth.idx)
  #X.test = X.test[truth.idx , ]
  truth = Z.test.res
  s.result  = predict(s.forest, X.test,estimate.variance = TRUE)
  #summary
  pred = s.result$predictions
  #print(pred)
  std = sqrt(s.result$variance.estimates)
  naives = naive.estimators( X.tr, Z.tr, delta.tr, X.test,  truth, L  )
  lu = lu.method(X.tr, Z.tr , Z.tr.res, delta.tr, delta.tr.res, X.test, truth, L)
  wang = wang.method(X.tr, Z.tr , Z.tr.res, delta.tr, delta.tr.res, X.test, truth, L, delta.tr.res/(1-G))
  
  #we should use \hatG from training...
  #surv_object <- Surv(time = Z.test, event = 1-delta.test)
  #weights = Ghat.new(surv_object, data.frame (X.test), sum(Z.test<=L))
  
  weights = Ghat.old(surv_object, data.frame (X.tr),  Z.test.res, data.frame (X.test))
  
  n.test = dim(X.test)[1]
  if(i == 1){
             output = rbind( helper(pred, truth, std,weights),
                           helper(naives$km.predictions,truth,rep(0, n.test),weights),
                           helper(naives$cox.predictions,truth,rep(0, n.test),weights),
                           helper(lu$lu.id.predictions,truth,rep(0, n.test),weights),
                           helper(lu$lu.exp.predictions,truth,rep(0, n.test),weights),
                           helper(wang$wang.id.predictions,truth,rep(0, n.test),weights),
                           helper(wang$wang.exp.predictions,truth,rep(0, n.test),weights)
                           )
  }else{
    
    output = output + rbind( helper(pred, truth, std,weights),
                             helper(naives$km.predictions,truth,rep(0, n.test),weights),
                             helper(naives$cox.predictions,truth,rep(0, n.test),weights),
                             helper(lu$lu.id.predictions,truth,rep(0, n.test),weights),
                             helper(lu$lu.exp.predictions,truth,rep(0, n.test),weights),
                             helper(wang$wang.id.predictions,truth,rep(0, n.test),weights),
                             helper(wang$wang.exp.predictions,truth,rep(0, n.test),weights)
                          )
    
  }
  
  
  
  
  
  feature = truth  #put interested feature here for plot
  
  x = sort(feature)
  
  y = truth[order(feature)]
  
  yhat = pred; 
  yhat.km = naives$km.predictions;
  yhat.c = naives$cox.predictions;
  yhat.lu.id = lu$lu.id.predictions;
  yhat.lu.exp = lu$lu.exp.predictions;
  yhat.wang.id = wang$wang.id.predictions;
  yhat.wang.exp = wang$wang.exp.predictions;
  stdhat = std[order(feature), ];
  
  plot(x,y,'l',col="red", xlim=c( min( feature )-0.01 , max( feature )+0.01 ), ylim=c(min( truth )-0.01 , max( truth )+0.01), ylab="RMST", xlab="Truth")
  
  lines(x , yhat, type="p", col="green"  ,  cex=0.3)
  
  #lines(x , yhat.r,type="p", col="76",  cex=0.3)
  
  lines(x , yhat.c, type="p", col="orange"  ,  cex=0.3)
  
  lines(x , yhat.km, type="p", col="yellow"  ,  cex=0.3)
  
  lines(x , yhat.lu.id, type="p", col="magenta"  ,  cex=0.3)
  
  lines(x , yhat.lu.exp, type="p", col="cyan"  ,  cex=0.3)
  
  lines(x , yhat.wang.id, type="p", col="black"  ,  cex=0.3)
  
  lines(x , yhat.wang.exp, type="p", col="blue"  ,  cex=0.3)
  
  legend("bottomright", 
         legend=c("srf","naive.cox", "naive.km", "lu.id", "lu.exp","wang.id", "wang.exp"),
         col=c("green", "orange", "yellow","magenta" , "cyan", "black", "blue"), 
         lty=1:2, cex=0.8)
  
  #lines(x , yhat + 1.96*stdhat,type="l", col="ivory3",  cex=0.03, lty=5)
  
  #lines(x , yhat - 1.96*stdhat,type="l", col="ivory3",  cex=0.03, lty=5)  
  
  
   
}

output = output / fold

row.names(output) <- c('custom.forest',
                       'naive.km.estimator', 
                       'naive.cox.estimator',
                       'lu.id.estimator',
                       'lu.exp.estimator',
                       'wang.id.estimator',
                       'wang.exp.estimator')
colnames(output)  <- c('coverage.prob', 'mae','mse')

cat("L: ");cat(L);cat(", ")
cat("num_features: ");cat(num_features);cat(", \n")
print(output)









