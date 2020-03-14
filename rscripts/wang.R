wang.method<-function(X.tr, Z.tr , Z.tr.res, delta.tr, delta.tr.res, X.te, truth, tau, weights)
{
  
  
  
  pred = NULL
  
  n = dim(X.te)[1]
  
  fitt=lm( Z.tr.res~as.matrix(X.tr), weights=weights)
  beta0=fitt$coef
  beta0[is.na(beta0)]=0
  pred = cbind(pred, (as.matrix(cbind(rep(1,n), X.te))%*%beta0)  )
  pred = cbind(pred, rep(0, dim(pred)[1]) )
  
  fitt=glm(Z.tr.res~as.matrix(X.tr), family="quasipoisson", weights=weights)
  beta0=fitt$coef
  beta0[is.na(beta0)]=0
  pred = cbind(pred,  exp(as.matrix(cbind(rep(1,n), X.te))%*%beta0) )
  pred = cbind(pred, rep(0, dim(pred)[1]) )
  
  colnames(pred)  <- c('wang.id.predictions', 'wang.id.variance.estimates', 'wang.exp.predictions', 'wang.exp.variance.estimates')
  
  pred=data.frame(pred)
  
  return(pred)
}