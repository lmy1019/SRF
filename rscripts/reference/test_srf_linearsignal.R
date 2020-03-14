n = 50; p = 10
X = matrix(rnorm(n*p), n, p)
Y = X[,1] * rnorm(n)
r.forest = regression_forest(X, Y)
X.test = matrix(0, 101, p)
X.test[,1] = seq(-2, 2, length.out = 101)
r.pred = predict(r.forest, X.test,estimate.variance = TRUE)








##test srf precompute_prediction_values
n = 40; p = 3
X = matrix(c( 0.7971644,  0.1643566,  1.9003346, -0.7435434, -0.8592456, -0.4087193, 
              -0.1475978,  0.9578696,  0.4545381,-1.2659856 ,-0.6239932,  0.6410597 , 
              0.8782860 , 1.2274549 ,-1.1981016, 6.197086e-01,1.672072e+00,-1.358019e+00,-1.513019e-01,1.664467e+00,-1.320852e+00,-1.956744e+00,
              -1.061079e+00,-1.262412e+00,6.049343e-01,2.847008e-01,6.149476e-01,6.512440e-01,-1.539953e+00,
              -2.613098e+00,-6.837740e-01,1.175292e+00,9.444710e-02,1.071358e+00,9.933585e-01,1.204634e-01,
              -6.635685e-02,-8.220334e-01,-3.751672e-01,1.183517e-01,-8.792098e-02,-6.357123e-01,-1.817490e+00,
              -8.382318e-01,5.611334e-01,7.557604e-01,1.054660e+00,-2.061587e+00,1.779892e-01,9.381476e-01,
              9.554074e-01,7.954354e-01,8.326882e-02,-9.757865e-02,1.399318e-01,-5.348283e-01,-1.024704e-01,
              9.605577e-01,-2.538372e-01,5.934649e-01,-2.286059e-01,8.543835e-01,1.243647e+00,5.815926e-02,
              -1.368670e+00,-2.464866e+00,-1.799905e+00,1.078503e+00,5.305150e-01,-6.467312e-01,-3.013599e-01,
              5.897217e-01,1.611932e-01,-1.459800e+00,-1.895871e+00,1.486191e+00,-2.329366e-01,-1.328581e+00,
              -7.923263e-01,-5.259562e-01,1.370358e-01,1.753682e-01,1.393375e-01,-5.200221e-01,1.579518e+00,
              -8.347396e-01,7.335344e-01,1.474749e+00,-7.451170e-01,9.995486e-02,-3.179263e-01,-1.735930e-01,
              1.706422e-01,1.627429e-01,-7.446587e-01,-8.747157e-01,1.164555e+00,1.688395e-02,2.156049e+00,
              -3.558398e-01,-7.148459e-01,-1.303688e+00,4.793875e-05,4.467090e-01,5.633411e-01,-6.279641e-01,
              -1.915379e+00,-1.794203e-01,-7.269782e-01,-1.024835e+00,-1.086034e+00,-1.988938e-01,1.673888e+00,
              -4.049213e-01,1.044266e-01,-1.284444e+00,1.010600e+00,-2.812033e-01,2.565989e-01,3.128000e-01), ncol = p)
betaT = c(0.36741111,0.61058509, 0.09764042)
betaC = c(0.1427991, 0.8416534, 0.7830785)
#survtime = rexp( X%*%betaT )
survtime = c(0.16984996,0.23491326,1.23799147,1.53806081,3.81319649,0.46062594,0.21549024,0.44181379,2.22807835,
             3.75809026,0.04025734,0.07906310,0.07313295,1.70790569,5.19917266,2.12111427,0.55573856,0.14198803,
             1.16559030,1.04683733,1.70724537,0.05047355,0.60899928,0.29062628,0.08461173,2.27773343,2.01192825,
             0.22720283,0.93904951,0.64970171,0.32100514,1.37751772,0.46318687,2.40557938,0.09603730,0.78967110,
             1.29171458,1.11538525,1.10726933,0.05055736)
#censtime = rexp( X%*%betaC )
censtime = c(0.16944709,0.44992185,0.22385162,2.56972802,0.24446657,0.51919091,0.55554215,0.10730781,0.23380832,
             0.44265684,1.23570742,1.79921638,0.01986540,0.52752181,0.32540804,0.85904828,1.45510586,0.08212552,
             0.71671124,0.50287087,0.09257617,3.23769139,0.59339833,1.27986229,2.85855721,0.35579162,0.98529894,
             0.02511905,0.44762430,2.39700182,0.57965442,0.20966340,0.38044764,1.20496874,1.44420241,1.37005885,
             0.03425792,0.07034222,1.64248583,0.68067947)
delta = (censtime>survtime)+0
Z =  pmin(survtime, censtime)
G = rep(0.5, n)
s.forest= custom_forest(X, Z, delta, G, honesty=FALSE)
options(digits=12)
s .result = predict(s.forest, X)
pred = s .result[seq(length(s .result)-40+2, length(s .result), 2)]
point = s .result[seq(length(s .result)-40+1, length(s .result), 2)]

unique_pred = unique(pred)
for(u in unique_pred){
  pc = point[abs(pred-u)<0.00001 & ! is.na(pred)]
  #print(pc)
  print(sum(Z[pc+1]*delta[pc+1]/0.5)/sum(delta[pc+1]/0.5))
}




#########################################################################################################
##test srf predict
library(survival)
library(survminer)
library(dplyr)

test_custom_forest <- function(n = 500, p = 2, L=1e17, hazardC=1, sigma=1){
  
  n=5000; p = 2; L=5.5; hazardC=0.2;sigma=0.01;
  X=matrix(rbinom(n*p,1,0.5),n,p); X.test = matrix(rbinom(n*p,1,0.5),n,p)
  #betaT = rep(0.5/p,p);
  betaT = rep(0.5/p,p); betaT[1]=1;
  betaC = rep(0,p); betaC[1]=1;
  #linear+noise_free
  survtime = X%*%betaT+rnorm(n,0,sigma)+5; Y.test = X.test%*%betaT+rnorm(n,0,sigma)+5; truth = X.test%*%betaT+5; Y.test.res.L=pmin(Y.test, L)
  censtime = rexp( n,exp(-X%*%betaC)*hazardC )
  Z =  pmin(survtime, censtime); Z.res.L=pmin(Z, L); 
  delta = (censtime>survtime)+0; delta.res.L=delta;  delta.res.L[Z>=L]=1
  surv_object <- Surv(time = Z, event = 1-delta)
  data <- data.frame (x = X)
  
  
  #---------------------------
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
  r.forest = regression_forest(X, Z.res.L)
  s.forest = custom_forest(X, Z.res.L, delta.res.L, G )
  
  #---------------------------------------
  
  s.result  = predict(s.forest, X.test,estimate.variance = TRUE)
  s.predictions =s.result $predictions
  s.var  = s.result $variance.estimates
  coverage.s = sum(Y.test.res.L<(s.predictions +1.96*sqrt(s.var )) & Y.test.res.L>(s.predictions -1.96*sqrt(s.var ))  )/n
  mse.s = sqrt(sum(( s.predictions -Y.test.res.L)^2))/n
  mae.s = (sum(abs( s.predictions -Y.test.res.L)))/n
  
  r.pred = predict(r.forest, X.test,estimate.variance = TRUE)
  r.prediction=r.pred$predictions
  r.var = r.pred$variance.estimates
  coverage.r = sum(Y.test.res.L<(r.prediction+1.96*sqrt(r.var)) & Y.test.res.L>(r.prediction-1.96*sqrt(r.var))  )/n
  mse.r = sqrt(sum(( r.prediction-Y.test.res.L)^2))/n
  mae.r = (sum(abs( r.prediction-Y.test.res.L)))/n
  
  output = data.frame(coverage=c(coverage.s,coverage.r), mae=c(mae.s, mae.r), mse=c(mse.s, mse.r), median=c( median(sqrt(s.var )),median(sqrt(r.var)) ) )
  row.names(output) <- c('custom.forest', 'regression.forest')
  #plot
  
  
  #plot(sort(X.test[,1]),truth[order(X.test[,1])],'l',col="red", xlim=c(-1,1), ylim=c(4.8,6))
  #lines(sort(X.test[,1]),s.predictions[order(X.test[,1])],'l',col="green")
  #lines(sort(X.test[,1]),r.prediction[order(X.test[,1])],'l',col="blue")
  
  
  #output
  return(output)
}


#senatity check1
#when censoring is all 0,  the performance of custom_forest should be increasing as regression forest.
#hazardC=1e-17
#n,  coverage        mae         mse     median
#100 0.82 0.28250809 0.054263900 0.15145394    (0.81 0.29409789 0.055014073 0.14437461)
#500 0.938 0.074341904 0.014166848 0.041643742 (0.938 0.075069320 0.014148393 0.045553829)
#1000 0.965 0.038343593 0.0052161126 0.025809047 (0.965 0.038559176 0.0051935989 0.024422171)
#2000 0.9905 0.016483333 0.0017158993 0.016648921 (0.9930 0.016344596 0.0017257976 0.016524812)
#5000 0.9932 0.0073805255 0.00052281312 0.0087114774 (0.9950 0.0071793641 0.00051817691 0.0090110131)

