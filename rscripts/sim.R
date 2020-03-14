rm(list=ls())
library(grf)
setwd("/Users/lmy/Dropbox/Mingyang-HLi/projects/survival_random_forest/code/rscript")
#setwd("/Users//mingyal/Dropbox/Mingyang-HLi/projects/survival_random_forest/code/rscript")
source("naive.estimators.R")
source("util.R")
source("lu.R")
source("wang.R")


sim_srf<- function(n = 500, p = 2, L = 1e17, hazardC = 0.2, SNR = 2, func='id',  times=5){
#print(sim_srf(2000,2,226,2.6*1e-3,snr,'exp',5))
  #n=10000;p = 8; L=5.4; hazardC=0.065; SNR=1; func='id';  times=5
  
  #---------------------------
  #generating simulated data
  #test data
  
  X.test =   matrix(runif(n*p,-1,1),n,p) 
  
  intercept=5;
  
  betaT = rep(0,p); betaT[1]=0.25; betaT[2]=0.25;
  
  betaC = rep(0,p); betaC[1]=log(2);
  
  pred = NULL; variance = NULL; output=NULL;
  
  if(func == 'id' | func == 'exp' | func=='logexp' ){
    
    truth = pmin( do.call( func,   list(  (X.test^2)%*%betaT + intercept) ) , L )
    
  } else if(func == 'friedman'){
    
    if(p<5){
      
      stop("p is smaller than 5")
      
    }
    
    truth = pmin( do.call( func, list(X.test) ), L )
    
  }
  
  sigma =  noise.std(SNR, func, intercept, betaT);
  
  
  cat('sigma: ');cat(sigma); cat(', ')
  cat('func: ');cat(func);cat(', ')
  cat('p: ');cat(p);cat(', ')
  cat('n: ');cat(n);cat(',\n')
  
  for( i in 1:times ){
    
    cat(i); cat('...')
  
    #training data
    
    X.tr = matrix(runif(n*p,-1,1),n,p)
    
    error.term = rnorm(n,0,sigma)
    
    survtime = NULL
    
    if(func == 'id' | func == 'exp' | func=='logexp' ){
      
      survtime = do.call(func,  list(  (X.tr^2)%*%betaT + intercept )   ) + error.term
      
    } else if(func == 'friedman'){
      
      survtime = do.call( func,   list(X.test) ) + error.term 
      
    }

    
    #Y.test = func(X.test%*%betaT)+error.term;  Y.test.res =pmin( Y.test , L ) ;
    
    censtime = rexp( n,exp(X.tr%*%betaC)*hazardC )
    
    Z.tr =  pmin(survtime, censtime)
    
    ##sort X.tr, Z.tr, Z.tr.res, delta.tr, delta.tr.res to speed up computing IPCW
    id1=order(Z.tr);   Z.tr=Z.tr[id1];  X.tr=X.tr[id1,];   survtime = survtime[id1];   censtime = censtime[id1]
    
    Z.tr.res =pmin(Z.tr, L)
    
    delta.tr = (censtime>survtime)+0; 
    
    delta.tr.res =delta.tr;  delta.tr.res [Z.tr>=L]=1
    
    surv_object <- Surv(time = Z.tr, event = 1-delta.tr)
    
    data <- data.frame (x = X.tr)
    
    #---------------------------
    #training
    #G.old = Ghat.old(surv_object, data, Z.tr)  #compute Ghat
    G = Ghat.new(surv_object, data, sum(Z.tr<=L))         #compute Ghat
    
    s.forest = custom_forest(X.tr, Z.tr.res , delta.tr.res , G )  #training survival forest
    #r.forest = regression_forest(X.tr, Z.tr.res )               #training regression forest for comparision when low censoring 
    
    #---------------------------------------
    #prediction
    
    s.result  = predict(s.forest, X.test,estimate.variance = TRUE)
    # r.result = predict(r.forest, X.test,estimate.variance = TRUE)
    # naives = naive.estimators( X.tr, Z.tr, delta.tr, X.test,  truth, L  )
    # lu = lu.method(X.tr, Z.tr , Z.tr.res, delta.tr, delta.tr.res, X.test, truth, L)
    # wang = wang.method(X.tr, Z.tr , Z.tr.res, delta.tr, delta.tr.res, X.test, truth, L, delta.tr.res/(1-G))
    
    #---------------------------------------
    #summarize result
    
    pred = bind_helper(pred, cbind(s.result$predictions 
                                   #r.result$predictions, 
                                   #naives$cox.predictions, 
                                   #naives$km.predictions,
                                   # lu$lu.id.predictions,
                                   # lu$lu.exp.predictions,
                                   # wang$wang.id.predictions,
                                   # wang$wang.exp.predictions
                                   ) 
                       ) 
    
    variance = bind_helper(variance, cbind(s.result$variance.estimates 
                                           #r.result$variance.estimates, 
                                           #naives$cox.variance.estimates, 
                                           #naives$km.variance.estimates,
                                           #lu$lu.id.variance.estimates,
                                           #lu$lu.exp.variance.estimates,
                                           #wang$wang.id.variance.estimates,
                                           #wang$wang.exp.variance.estimates
                                           ) 
                           ) 
    
  }
  
  cat('done...')
  
  pred = pred / times;  std=sqrt(variance / times) ; 
  
  output = summerize(pred, std, truth,  mean(delta.tr),  sum(survtime>=L)/n )
  
  row.names(output) <- c('custom.forest' 
                         #'regression.forest', 
                         #'naive.cox.estimator', 
                         #'naive.km.estimator',
                         #'lu.id.estimator',
                         #'lu.exp.estimator',
                         #'wang.id.estimator',
                         #'wang.exp.estimator'
                         )
  colnames(output)  <- c('coverage.prob', 'mae','mse', 'uncensored.rate', 'truncated.rate')
  
  #---------------------------------------
  #plot
  
  feature = truth  #put interested feature here for plot
  
  x = sort(feature)
  
  y = truth[order(feature)]
  
  #yhat = pred[order(feature),][,1]; 
  yhat = pred[order(feature)]; 
  
  # yhat.r = pred[order(feature), ][,2] ; 
  # yhat.c = pred[order(feature), ][,3];
  # yhat.km = pred[order(feature), ][,4];
  # yhat.lu.id = pred[order(feature), ][,5];
  # yhat.lu.exp = pred[order(feature), ][,6];
  # yhat.wang.id = pred[order(feature), ][,5];
  # yhat.wang.exp = pred[order(feature), ][,6];
  
  #stdhat = std[order(feature), ][,1 ]; 
  stdhat = std[order(feature)]; 
  
  #stdhat.r = std[order(feature), ][,2 ]
  
  plot(x,y,'l',col="red", xlim=c( min( feature )-0.01 , max( feature )+0.01 ), ylim=c(min( truth-1.96*stdhat )-0.01 , max( truth+1.96*stdhat )+0.01), ylab="RMST", xlab="Truth")
  
  lines(x , yhat, type="p", col="green"  ,  cex=0.3)
  
  #lines(x , yhat.r+1.96*stdhat.r,type="l", col="black",  cex=0.3)
  #lines(x , yhat.r-1.96*stdhat.r,type="l", col="black",  cex=0.3)
  
  # lines(x , yhat.c, type="p", col="orange"  ,  cex=0.3)
  # 
  # lines(x , yhat.km, type="p", col="yellow"  ,  cex=0.3)
  # 
  # lines(x , yhat.lu.id, type="p", col="magenta"  ,  cex=0.3)
  # 
  # lines(x , yhat.lu.exp, type="p", col="cyan"  ,  cex=0.3)
  # 
  # lines(x , yhat.wang.id, type="p", col="black"  ,  cex=0.3)
  # 
  # lines(x , yhat.wang.exp, type="p", col="blue"  ,  cex=0.3)
  
  # legend("topleft", 
  #        legend=c("srf","naive.cox", "naive.km", "lu.id", "lu.exp","wang.id", "wang.exp"),
  #        col=c("green", "orange", "yellow","magenta" , "cyan", "black", "blue"), 
  #        lty=1:2, cex=0.8)
  
  lines(x , yhat + 1.96*stdhat,type="l", col="ivory3",  cex=0.03, lty=5)
  
  lines(x , yhat - 1.96*stdhat,type="l", col="ivory3",  cex=0.03, lty=5)
  
  #lines(x, y, ,type="l", col="red")
  
  #lines(x , yhat + 1.96*stdhat.r,type="l", col="black",  cex=0.1)
  
  #lines(x , yhat - 1.96*stdhat.r,type="l", col="black",  cex=0.1)
  
  #print(sum(abs(pred[,1]-pred[,2]))/n)
  
  #print(sum(abs(std[,1]-std[,2]))/n)
  
  
  
  
  return( output )

}  
  
  
#Note1: increase SNR may not increase coverage.prob
#sim_srf(1000, 3, 1e17, 1e-17, 1,logexp ,10)
#SNR 
#0.01 2.9669, 
#custom.forest                 1 0.09727412 0.01822500 
#regression.forest             1 0.10246104 0.01960876
#0.1  0.9382 
#custom.forest             0.999 0.04425544 0.003484341      
#regression.forest         1.000 0.04307084 0.003234238
#1    0.2966
#custom.forest             0.986 0.02950164 0.001821891
#regression.forest         0.993 0.02665284 0.001455189
#10   0.0938, 
#custom.forest             0.915 0.02245890 0.0011491006
#regression.forest         0.933 0.01959888 0.0008818021
#100  0.0296
#custom.forest             0.799 0.02257604 0.001326754
#regression.forest         0.835 0.01903960 0.001001441 








# for(  snr in c( 0.8,1, 2 ) ){
#     
#   for(n in c(1000, 2000, 5000)) {
#   
#   
#      for(  i in c(2,4,6, 8)  ){
#        
#        print(c('id',snr, n, i))
#   
#        print( sim_srf(n,i,5.435,0.1,snr,id,5)       )
#   
#      }
#   }
# }

# for(  snr in c( 0.8,1, 2 ) ){
#   
#   for(n in c(1000, 2000, 5000)) {
#     
#     
#     for(  i in c(2,4,6, 8)  ){
#       
#       print(c('logexp',snr, n, i))
#       
#       print( sim_srf(n,i,5.435,0.1,snr,logexp,5)       )
#       
#     }
#   }
# }







start_time <- Sys.time()
print( start_time   )

print(sim_srf(10000,8,226,2.6*1e-3,0.8,'exp',5))

end_time <- Sys.time()
print(  end_time - start_time   )