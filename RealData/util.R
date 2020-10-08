library(survival)
library(survminer)
library(grf)
library(dplyr)
library(cubature) 




Ghat.new<-function(surv_object, data, thres){
  
  #weights1=delta.tr.res/rep(fitc1$surv, table(Z.tr.res))
  
  n = dim(data)[1]
  
  cox <- coxph(surv_object ~ ., data = data)
    
  cox.surv = survfit(cox, data)$surv
  
  if(thres==0){
    
    warnings(" All test points are truncated ")
    
  }
  
  if( n> thres ){
    
    G = 1 - c(   diag(cox.surv)[1:thres] ,  cox.surv[(thres+1):n,thres]  )
    
  }else{
    
    G = 1 - c(   diag(cox.surv) )
  }
  
  G[is.na(G )] <- 0  ;  G[G==1]<-1-1e-10
  
  return(G)
}



Ghat.old<-function(surv_object, data, newZ, newdata){

  cox <- coxph(surv_object ~ ., data = data)

  G =c()

  n = dim(newdata)[1]

  for(i in 1:n){
    print(i)

    #if(i%%1000==0){cat(i); cat(', ')}

    cox.summary = summary(survfit(cox, newdata[i,]))

    z=newZ[i]

    count = which(sort(c(cox.summary$time, z+1e-10))==z+1e-10)-1

    if(count <1){count = 1}

    G =c(G , 1-cox.summary$surv[count])

  }

  #cat("\n")

  G[is.na(G )] <- 0  ;  G[G==1]<-1-1e-10

  return(G)
}

mean.sq.error<-function(a){
  
  l = length(a)
  
  return(  sum(a^2)/l   )
  
}

mean.abs.error<-function(a){
  
  l = length(a)
  
  return( sum( abs( a ) )/l )
  
}



summerize<-function(pred, std, truth, uncensor.rate, t.rate, weights){
  
  df <- data.frame(matrix(ncol = 5, nrow = 0))
  
  num.col = dim(pred)[2]
  
  n = dim(pred)[1]
  
  for(i in 1:num.col){
      coverage.prob = 0
      if(is.null(std)){
          coverage.prob = sum( (truth < ( pred[,i] + 1.96*std[,i] ) & truth > ( pred[,i] -1.96*std[,i] ) ) + 0 )/n
      }
      mse =  mean(((pred[,i] -truth)**2) * weights)

      mae = mean((abs(pred[,i] -truth)) * weights)

      df = rbind( df, c(  coverage.prob, mae, mse , uncensor.rate, t.rate  ) )

  }
  
  #colnames(df) <- c('coverage.prob', 'mae', 'mse', 'uncensor.rate', 't.rate')

  return( df  )
  
}



id <- function(x){  return(x) }



logexp <- function(x){ return( log( exp(x)+1 ) )   }

friedman<-function(x){
  
  #return( sin(pi*x[,1]*x[,2]) + (x[,3])^2+x[,4]+x[,5]+3      )
  
  return( x[,1]*x[,2]  + 2      )
  
}


bind_helper<- function(a, b){
  if( is.null(a) ){  
    a = b
  }else{
    a = a + b 
  }
  return(a)
}



noise.std<-function(SNR, func, intercept, betaT){
  signal.var=0
  
  d = length(betaT)
  
  #exp
  if( func == 'friedman'){
    
    int<-function(x){  10*sin(pi*x[1]*x[2]) + 20*(x[3]-0.5)^2+10*x[4]+5*x[5]     }
    
  } else if( func=='exp' ){
    
      int<-function(x){  exp( sum( betaT*x )+intercept )    }
    
  #id    
  } else if( func=='id' ){
      
      int<-function(x){  sum( betaT*x ) + intercept     }
  
  #logexp    
  } else if( func== 'logexp' ){
    
      int<-function(x){  log( exp( sum( betaT*x ) + intercept  ) +1  )    }
    
  } 
  
  int.sq <- function(x){  int(x)^2     }
  
  signal.var = ( adaptIntegrate( int.sq, lowerLimit = rep(0,d), upperLimit = rep(1,d) )$integral 
                 - ( adaptIntegrate(int, lowerLimit = rep(0,d), upperLimit = rep(1,d))$integral )^2 
               )
  
  return( sqrt(signal.var/SNR)  )
  
}
  
  # lines(x , y ,'p',col="green",cex=0.3)
  # 
  # sd = std[ order(X.test.feature) ]
  # 
  # lines(x , y + 1.96*sd,type="l", col="black",  cex=0.3)
  # 
  # lines(x , y - 1.96*sd,type="l", col="black",  cex=0.3)

test.Ghat.new<-function(){
  
  n=1000;p = 6; L=5.435; hazardC=0.1; SNR=1e17; func='exp';  times=1;
  
  X.tr = matrix(runif(n*p,0,1),n,p)
  
  error.term = rnorm(n,0,sigma)
  
  survtime = NULL
  
  if(func == 'id' | func == 'exp' | func=='logexp' ){
    
    survtime = do.call(func,  list(  (X.tr^2)%*%betaT + intercept+error.term )   )
    
  } else if(flag == 'friedman'){
    
    #survtime = func(  X.test ) + error.term
    
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
  
  start_time <- Sys.time()
  G.old = Ghat.old(surv_object, data, Z.tr)  #compute Ghat
  end_time <- Sys.time()
  print(  end_time - start_time   )
  start_time <- Sys.time()
  G.new = Ghat.new(surv_object, data)
  end_time <- Sys.time()
  print(  end_time - start_time   )
  print(sum(abs(G.old-G.new))   )
  
}