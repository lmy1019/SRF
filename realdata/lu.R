


lu.method<-function(X.tr, Z.tr , Z.tr.res, delta.tr, delta.tr.res, X.te, truth, tau, weights)
  {
    
    
    fit=survfit(  Surv(Z.tr, 1-delta.tr) ~1  )  
    
    surv.prob = fit$surv
    
    surv.prob = c(  surv.prob[1:sum(Z.tr<=tau)], rep(  surv.prob[sum(Z.tr<=tau)],  sum(Z.tr>tau) )       )
    
    weights=delta.tr.res/rep(surv.prob, table(Z.tr )   )
    
    weights[is.na(weights )] <- 1e10
    weights[weights >1e10] <- 1e10
    
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
    
    colnames(pred)  <- c('lu.id.predictions', 'lu.id.variance.estimates', 'lu.exp.predictions', 'lu.exp.variance.estimates')
    
    pred=data.frame(pred)
    
    return(pred)
}