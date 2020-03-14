rm(list=ls())
library(TCGA2STAT)
library(plyr)
library(stringr)

library(grf)
setwd("/Users/lmy/Dropbox/Mingyang-HLi/projects/survival_random_forest/code/rscript")
source("naive.estimators.R")
source("util.R")
source("lu.R")
source("wang.R")
library(dummies)
library(emil)


tcga_performance <- function(L=3, num_features=-1){
  helper<-function(pred, truth, std, weights){
    return( 
      c( #mean( (truth<=pred+1.96*std) & (truth>=pred-1.96*std) ), 
         mean(abs(pred-truth)*weights),
         sqrt( mean( (pred-truth)^2 *weights ) ) 
      ) 
    )
  }
  
  cat("L: ");cat(L);cat(", ")
  cat("num_features: ");cat(num_features);cat(", \n")
  
  # data pre-processing: remove multiple measurements for one patient,
  # and find intersection of parients of protein and clinic
  #------------------------------------------------------------------------------------
  protein <- read.csv("/Users/lmy/Dropbox/Mingyang-HLi/projects/survival_random_forest/realdata/TCGA_BRCA_PROTEIN/protein.csv")
  protein$bcr_patient_barcode<-tolower(substr( str_trim( unlist(strsplit(toString(protein$Sample_REF), split="[,]")) ), 1, 12 ) )
  #there is a piece after bcr_patient_barcode "01A" "11A" "01B" "06A" "11B"
  protein$meas_idx <- substr( str_trim( unlist(strsplit(toString(protein$Sample_REF), split="[,]")) ), 14, 16 ) 
  protein=protein[protein$meas_idx=="01A",]
  clinic <- read.csv("/Users/lmy/Dropbox/Mingyang-HLi/projects/survival_random_forest/realdata/TCGA_BRCA_PROTEIN/clinic.csv")
  #clinic = clinic[!is.na(clinic$years_to_birth),]
  clinic$bcr_patient_barcode<-clinic$Hybridization.REF
  clinic= clinic[clinic$bcr_patient_barcode %in% protein$bcr_patient_barcode, ]
  clinic = clinic[  !( is.na(clinic$days_to_death) & is.na(clinic$days_to_last_followup) ),    ]
  clinic$days_to_death[is.na(clinic$days_to_death)] <- -1
  clinic$days_to_last_followup[is.na(clinic$days_to_last_followup)] <- -1
  clinic = clinic[ !( clinic$days_to_last_followup<=0 & clinic$days_to_death<=0 ),  ]
  clinic$observed.time=apply(clinic[,c(5,6)],1,max,na.rm=TRUE) 
  protein= protein[protein$bcr_patient_barcode %in% clinic$bcr_patient_barcode, ]
  
  #data pre-processing: rm feature with missing data, and with same value 
  #------------------------------------------------------------------------------------
  protein = protein[,!apply(is.na(protein), 2, any)]
  protein = protein[vapply(protein, function(x) length(unique(x)) > 1, logical(1L))]
  clinic = clinic[,!apply(is.na(clinic), 2, any)]
  clinic = clinic[vapply(clinic, function(x) length(unique(x)) > 1, logical(1L))]
  tcga = merge(protein, clinic, by = "bcr_patient_barcode")
  
  tcga = tcga[sample(nrow(tcga)),]

  
  non.feature.list = c("bcr_patient_barcode" ,
                       "Sample_REF", 
                       "Hybridization.REF",
                       "vital_status",
                       "days_to_death",
                       "days_to_last_followup",
                       "observed.time")
  
  #data pre-processing: create design-matrix X, observed-time Z, delta
  #------------------------------------------------------------------------------------
  
  X = tcga[, ! (colnames(tcga) %in% non.feature.list) ]
  n = dim(X)[1]
  Z = ( tcga$observed.time+rnorm(n, 0, 1e-5) )
  delta = tcga$vital_status
  X = cbind(  X[, 1:216]
              #dummy(X$pathology_T_stage, sep = "_")[,-1],
              #dummy(X$pathology_N_stage, sep = "_")[,-1],
              #dummy(X$pathology_M_stage, sep = "_")[,-1],
              #dummy(X$gender, sep = "_")[,-1],
              #dummy(X$histological_type, sep = "_")[,-1],
              #X$years_to_birth
          )
  X=scale(X)
  
  #num_features = p
  fold = 10;  unit = floor(n/fold);   base = (fold-1)*unit
  output=NULL;
  #L=1000
  
  
  
  
  for( i in 1:fold ){
    
    
    tr.idx = ( 0:(base-1) + (i-1)*unit ) %% n +1
    X.tr = X[tr.idx, ];  #X.tr <- scale(X.tr)
    Z.tr = Z[tr.idx];  delta.tr = delta[tr.idx]
    
    
    #feature selection
    if(num_features>=1){
      pv = c()
      surv_object <- Surv(time = Z.tr, event = delta.tr)
      for( l in 1:dim(X)[2]){
        cox <- coxph(surv_object ~ ., data = data.frame(X.tr[,l]) )
        pv = c(pv, summary(cox)$waldtest["pvalue"])
      }
      id.features=order(pv)
      id.features= id.features[1:num_features]
    }else{
      id.features = 1:dim(X)[2]
    }
    X.tr = X.tr[,id.features]
    
    
    
    #training
    cat(i); cat("...");

    id1=order(Z.tr);   Z.tr=Z.tr[id1];  X.tr=X.tr[id1,];   delta.tr = delta.tr[id1]
    surv_object <- Surv(time = Z.tr, event = 1-delta.tr)
    G = Ghat.new(surv_object, data.frame (X.tr), sum(Z.tr<=L)) 
    G[G>=1-1e-8]<-1-1e-8
    Z.tr.res =pmin(Z.tr, L)
    delta.tr.res =delta.tr;  delta.tr.res [Z.tr>=L]=1
    s.forest = custom_forest(X.tr, Z.tr.res , delta.tr.res , G , num.trees = 10000)
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
    
    #print(output)
    
    
    
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
    
    plot(x,y,'l',col="red", xlim=c( min( feature )-1 , max( feature )+1 ), ylim=c(min( truth )-1 , max( truth )+5), ylab="RMST", xlab="Truth")
    
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
  colnames(output)  <- c('mae','rmse')
  
  
  #print(output)
  return(output)

}
# for(l in c(1,3,5)){
#   cat("L: "); cat(l); cat('\n')
#   for( i in 1:10){
#     output = tcga_performance(l ,-1)
#     print(output)
#   }
# }

# 365*5, -1 works 1 time
for( i in 1:3){
  output = tcga_performance(365*5 ,-1)
  print(output)
}
  

