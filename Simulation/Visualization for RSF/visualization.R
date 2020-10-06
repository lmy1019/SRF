rm(list=ls())
library(grf)
setwd("/Users/lmy/Dropbox/Mingyang-HLi/projects/survival_random_forest/code/rscript")
#setwd("/Users//mingyal/Dropbox/Mingyang-HLi/projects/survival_random_forest/code/rscript")
source("naive.estimators.R")
source("util.R")
source("lu.R")
source("wang.R")


n=3000; p = 5; L=190; hazardC=2.6*1e-3; SNR=0.3; func='exp';  times=1; model = 1
n=3000; p = 5; L=190; hazardC=2.6*1e-3; SNR=0.3; func='exp';  times=1; model = 2
n=3000; p = 5; L=5.3; hazardC=0.08; SNR=0.3; func='id';  times=1; model = 1
n=3000; p = 5; L=5.3; hazardC=0.08; SNR=0.3; func='id';  times=1; model = 2
n=3000; p = 5; L=5.3; hazardC=0.08; SNR=0.3; func='logexp';  times=1; model = 1
n=3000; p = 5; L=5.3; hazardC=0.08; SNR=0.3; func='logexp';  times=1; model = 2




#---------------------------



intercept=5;

betaT = rep(0,p); betaT[1]=0.25; betaT[2]=0.25;

betaC = rep(0,p); betaC[1]=log(2);

sigma =  noise.std(SNR, func, intercept, betaT);

cat('sigma: ');cat(sigma); cat(', ')
cat('func: ');cat(func);cat(', ')
cat('p: ');cat(p);cat(', ')
cat('n: ');cat(n);cat(',\n')
cat('model: ');cat(model);cat(', ')
cat('L: '); cat(L);cat(', ')
cat('hazardC: ');cat(hazardC);cat(', ')
cat('SNR: '); cat(SNR);cat(',\n')

#generating simulated data

#for( i in 1:times ){

i = times

cat(i); cat('...')

#test data

set.seed(i)

X.test =   matrix(runif(n*p,-1,1),n,p) 

pred = NULL; variance = NULL; output=NULL; truth = NULL

error.term = rnorm(n,0,sigma)

truth = pmin( do.call( func,   list(  (X.test^model)%*%betaT + intercept) ) , L )

survtime.test = pmin( do.call( func,   list(  (X.test^model)%*%betaT + intercept) ) , L )+error.term

censtime.test = rexp( n,exp(X.test%*%betaC)*hazardC )

Z.test =  pmin(survtime.test, censtime.test)

Z.test.res =pmin(Z.test, L)

delta.test = (censtime.test>survtime.test)+0; 

delta.test.res =delta.test;  delta.test.res[Z.test>=L]=1


#training data

X.tr = matrix(runif(n*p,-1,1),n,p)

error.term = rnorm(n,0,sigma)

survtime = NULL

survtime = do.call(func,  list(  (X.tr^model)%*%betaT + intercept )   ) + error.term

censtime = rexp( n,exp(X.tr%*%betaC)*hazardC )

Z.tr =  pmin(survtime, censtime)

##sort X.tr, Z.tr, Z.tr.res, delta.tr, delta.tr.res to speed up computing IPCW
id1=order(Z.tr);   Z.tr=Z.tr[id1];  X.tr=X.tr[id1,];   survtime = survtime[id1];   censtime = censtime[id1]

Z.tr.res =pmin(Z.tr, L)

delta.tr = (censtime>survtime)+0; 

delta.tr.res =delta.tr;  delta.tr.res [Z.tr>=L]=1

surv_object <- Surv(time = Z.tr, event = 1-delta.tr)

data <- data.frame(x = X.tr)

#---------------------------
#training

#compute Ghat
cox <- coxph(surv_object ~ ., data = data)
G =c()
for(i in 1:n){
  #if(i%%1000==0){print(i)}
  cox.summary = summary(survfit(cox, data[i,]))
  z=Z.tr[i]
  count = which(sort(c(cox.summary$time, z+1e-10))==z+1e-10)
  count = count-1
  if(count <1){count = 1}
  G =c(G , 1-cox.summary$surv[count])
}
G[is.na(G )] <- 0
G[G==1]<-1-1e-10

s.forest = custom_forest(X.tr, Z.tr.res , delta.tr.res , G , num.trees = 10000 )  #training survival forest
#r.forest = regression_forest(X.tr, Z.tr.res )               #training regression forest for comparision when low censoring 

#---------------------------------------
#prediction

s.result  = predict(s.forest, X.test,estimate.variance = TRUE)
#r.result = predict(r.forest, X.test,estimate.variance = TRUE)
naives = naive.estimators( X.tr, Z.tr, delta.tr, X.test,  truth, L  )
lu = lu.method(X.tr, Z.tr , Z.tr.res, delta.tr, delta.tr.res, X.test, truth, L)
wang = wang.method(X.tr, Z.tr , Z.tr.res, delta.tr, delta.tr.res, X.test, truth, L, delta.tr.res/(1-G))



#---------------------------------------
#compute Ghat
data.test <- data.frame(x = X.test)
G.test =c()
for(i in 1:n){
  #if(i%%1000==0){print(i)}
  cox.summary = summary(survfit(cox, data.test[i,]))
  z=Z.test[i]
  count = which(sort(c(cox.summary$time, z+1e-10))==z+1e-10)
  count = count-1
  if(count <1){count = 1}
  G.test =c(G.test , 1-cox.summary$surv[count])
}
G.test[is.na(G.test )] <- 0
G.test[G.test==1]<-1-1e-10
#---------------------------------------
#summarize result

pred = bind_helper(pred, cbind(s.result$predictions,
                               #r.result$predictions, 
                               naives$cox.predictions, 
                               naives$km.predictions,
                               lu$lu.id.predictions,
                               lu$lu.exp.predictions,
                               wang$wang.id.predictions,
                               wang$wang.exp.predictions) 
) 


cat('done...')

#pred = pred;  



output = summerize(pred, NULL, Z.test.res,  mean(delta.tr),  sum(survtime>=L)/n, delta.test.res/(1-G.test) )

row.names(output) <- c('custom.forest' ,
                       #'regression.forest', 
                       'naive.cox.estimator', 
                       'naive.km.estimator',
                       'lu.id.estimator',
                       'lu.exp.estimator',
                       'wang.id.estimator',
                       'wang.exp.estimator'
)
colnames(output)  <- c('coverage.prob', 'mae','mse', 'uncensored.rate', 'truncated.rate')





#plot the true RMST, Predicted RMST and CI with respect to the first feature
#--------------------------------------
#output[['predictions']]=s.result
#plot(sort(X.test[,1]),truth[order(X.test[,1])],'l',xlim=c(0,1), ylim=c(100,700),col="red",  xlab="Feature 1", ylab="RMST", main="Visualization of Performance of SRF", sub = paste0("Coverage Prob: ", coverage.s, ",RMSE: ", round(sqrt(mse.s),digits=3), ', MAE: ', round((mae.s),digits=3)))
plot(sort(truth),truth[order(truth)],'l',xlim=c(min(truth)-(median(truth)-min(truth))*0.05, (max(truth)+(-median(truth)+max(truth))*0.05)), ylim=c(min(truth)-(median(truth)-min(truth))*0.05, (max(truth)+(-median(truth)+max(truth))*0.05)),col="red",  xlab="Truth", ylab="Pediction", main=paste("Link:", func, ", model:", model))
lines(sort(truth),s.result$predictions[order(truth)],'l',col="green")
lines(sort(truth),naives$cox.predictions[order(truth)],'l',col="black")
lines(sort(truth),wang$wang.id.predictions[order(truth)],type = "l",col="yellow")
lines(sort(truth),wang$wang.exp.predictions[order(truth)],type = "l",col="blue")
legend('topleft', legend=c("True RMST", "Predicted RMST", "Naive.Cox", "Wang method(id link)", "Wang method(exp link)"), col=c("red", "green", "black", "yellow", "blue"), lty=1:2, cex=0.8)
print(output)

write.csv(output, paste("/Users/lmy/Dropbox/Mingyang-HLi/projects/survival_random_forest/Revision_material/code for Table2/Visualization for RSF/", func, '_', 'model', model, '_', p, 'test.csv', sep = ""), row.names = FALSE)
