#########################################################################################################
##test srf predict
library(survival)
library(survminer)
library(dplyr)
library(grf)

test_custom_forest <- function(n = 500, p = 2, L=1e17, hazardC=1, sigma=0.01, times=5){
  
  n=250; p = 8; L=5.5; hazardC=0.03; sigma=0.01;times = 5
  X=matrix(runif(n*p,0,1),n,p); X.test = matrix(runif(n*p,0,1),n,p)
  betaT = rep(0,p); betaT[1]=2;
  betaC = rep(0,p); betaC[1]=1;
  survtime = log( exp( X%*%betaT+5 )+1 )+rnorm(n,0,sigma); Y.test = log( exp( X.test%*%betaT+5 ) + 1 )+rnorm(n,0,sigma); truth = pmin( log( exp( X.test%*%betaT+5 ) + 1 ),L); Y.test.res.L=pmin(Y.test, L)
  censtime = rexp( n,exp(X%*%betaC)*hazardC )
  Z =  pmin(survtime, censtime); Z.res.L=pmin(Z, L); 
  delta = (censtime>survtime)+0; delta.res.L=delta;  delta.res.L[Z>=L]=1
  surv_object <- Surv(time = Z, event = 1-delta)
  data <- data.frame (x = X)
  
  output={}
  output[['predictions']]=NULL
  output[['sd']]=NULL
  output[['truth']]=truth

  #---------------------------
  for(i in 1:times){
    cat(i);cat('...')
    cox <- coxph(surv_object ~ ., data = data)
    G =c()
    for(i in 1:n){
      if(i%%1000==0){print(i)}
      cox.summary = summary(survfit(cox, data[i,]))
      z=Z[i]
      count = which(sort(c(cox.summary$time, z+1e-10))==z+1e-10)
      count = count-1
      if(count <1){count = 1}
      G =c(G , 1-cox.summary$surv[count])
    }
    G[is.na(G )] <- 0
    G[G==1]<-1-1e-10
    #r.forest = regression_forest(X, Z.res.L)
    s.forest = custom_forest(X, Z.res.L, delta.res.L, G, num.trees = 1000 )

    # r.pred = predict(r.forest, X.test,estimate.variance = TRUE)
    # r.prediction=r.pred$predictions
    # r.var = r.pred$variance.estimates
    # coverage.r = sum(Y.test.res.L<(r.prediction+1.96*sqrt(r.var)) & Y.test.res.L>(r.prediction-1.96*sqrt(r.var))  )/n
    # mse.r = sqrt(sum(( r.prediction-Y.test.res.L)^2))/n
    # mae.r = (sum(abs( r.prediction-Y.test.res.L)))/n
    
    s.result  = predict(s.forest, X.test,estimate.variance = TRUE)
    s.predictions =s.result $predictions
    s.var  = s.result $variance.estimates
    #coverage.s = sum(truth<(s.predictions +1.96*sqrt(s.var )) & truth>(s.predictions -1.96*sqrt(s.var ))  )/n
    mse.s = sqrt(sum(( s.predictions -Y.test.res.L)^2))/n
    mae.s = (sum(abs( s.predictions -Y.test.res.L)))/n
    
    #--------------------------------------
  
    if(is.null(output[['predictions']])){
      output[['predictions']]=s.predictions
      output[['sd']]=sqrt(s.var)
    }
    else{
      output[['predictions']]=cbind(output[['predictions']],s.predictions  )
      output[['sd']]=cbind(output[['sd']], sqrt(s.var)  )
    }
  }
  #plot
  plot(sort(X.test[,1]),truth[order(X.test[,1])],'l',col="red", xlim=c(0,1), ylim=c(4,7))
  p = rowMeans(output[['predictions']])
  lines(sort(X.test[,1]),p[order(X.test[,1])],'l',col="green")
  #lines(sort(X.test[,1]),r.prediction[order(X.test[,1])],'l',col="blue")
  std = rowMeans(output[['sd']])
  lines(sort(X.test[,1]),p[order(X.test[,1])]+1.96*std,type="l", col="black",  cex=0.3)
  lines(sort(X.test[,1]),p[order(X.test[,1])]-1.96*std,type="l", col="black",  cex=0.3)
  
  output[['summary of each x']]=data.frame(coverage=c(), mae=c(mae.s), mse=c(mse.s), uncensored=sum(delta)/n, trate=sum(survtime>=L)/n)

  
  return(output)
}


  
#data.frame(coverage=c(coverage.s), mae=c(mae.s), mse=c(mse.s), uncensored=sum(delta)/n, trate=sum(survtime>=L)/n)
  



# for(i in 1:3){
#   cat(i);cat('...')
#   
#     output=test_custom_forest(1000,8,5.5,0.03,0.05)
#   }else{
#     output=output+test_custom_forest(1000,8,5.5,0.03,0.05)
#   }
# }
# output=output/3
# output