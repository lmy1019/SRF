library(survival)
library(survRM2)


#testdata set for cox.rmst
# D=rmst2.sample.data()
# 
# D=D[D$arm==0,]
# 
# time=D$time
# 
# status=D$status
# 
# X = D[,4:8]
# 
# n = dim(x)[1]
# 
# tr.set=1:50
# 
# te.set = 51:154
# 
# X.tr = X[tr.set,];   Z.tr = time[tr.set];   delta.tr = status[tr.set]
# 
# X.te = X[te.set, ][status[te.set]==1,];    truth = time[te.set ][status[te.set]==1];
# 
# tau = 9;


#X.tr = X; Z.tr = Z; delta.tr = delta; X.te = X.test; tau = L;
















naive.estimators<-function(X.tr, Z.tr , delta.tr, X.te, truth, tau){

############# cox estimator  
    
  cox=  coxph(Surv(Z.tr, delta.tr) ~ ., data=data.frame(X.tr)) 
  
  n = dim(X.te)[1]
  
  time.period = summary( survfit(cox, data.frame(X.te))  )$time
  
  idx=time.period<=tau
  
  wk.time=c(time.period[idx],tau)
  
  pred = c();  var = c();
  
  for(i in 1:n){
    
    ft.cox = summary(survfit(  cox, data.frame(X.te)[i,]  ) )
    
    wk.surv=ft.cox$surv[idx]
    
    wk.n.risk =ft.cox$n.risk[idx]
    
    wk.n.event=ft.cox$n.event[idx]
    
    #estimator of rmst by $int_0^tau \hatS(t) dt$
    time.diff <- diff(c(0, wk.time))
    
    areas <- time.diff * c(1, wk.surv)
    
    rmst = sum(areas)
    
    pred=c(pred,rmst)
    
    wk.var <- ifelse((wk.n.risk-wk.n.event)==0, 0, wk.n.event /(wk.n.risk *(wk.n.risk - wk.n.event)))
    
    wk.var =c(wk.var,0)
    
    rmst.var = sum( cumsum(rev(areas[-1]))^2 * rev(wk.var)[-1])
    
    var = c(var, rmst.var)
    
  }
  
  df = cbind(pred, var)
  
############# km estimator
  
  ft.km = summary(survfit(Surv(Z.tr, delta.tr) ~ 1, data=data.frame(X.tr)) )
  
  time.period = ft.km$time
  
  idx=time.period<=tau
  
  wk.time=c(time.period[idx],tau)
  
  wk.surv=ft.km$surv[idx]
  
  time.diff <- diff(c(0, wk.time))
  
  areas <- time.diff * c(1, wk.surv)
  
  rmst = sum(areas)
  
############## combine results 
  
  df = cbind(df, rep(rmst, n), rep(0,n))
  
  colnames(df)  <- c('cox.predictions', 'cox.variance.estimates', 'km.predictions', 'km.variance.estimates')
  
  df=data.frame(df)
  
  return (  df  )

}

#cox.rmst(X.tr, Z.tr , delta.tr, X.te, truth, tau)