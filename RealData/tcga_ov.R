rm(list=ls())
library(TCGA2STAT)
library(plyr)
library(stringr)
library(grf)
#change to the path of realdata
setwd("/Users/lmy/Dropbox/Mingyang-HLi/projects/survival_random_forest/code/realdata")
source("naive.estimators.R")
source("util.R")
source("lu.R")
source("wang.R")
library(dummies)
library(emil)
set.seed(1)




# data pre-processing: remove multiple measurements for one patient,
# and find intersection of parients of protein and clinic
#------------------------------------------------------------------------------------
protein <- read.csv("./TCGA_OV_PROTEIN/protein.csv")
protein$bcr_patient_barcode<-tolower(substr( str_trim( unlist(strsplit(toString(protein$Sample.REF), split="[,]")) ), 1, 12 ) )
#there is a piece after bcr_patient_barcode "01A" "11A" "01B" "06A" "11B"
protein$meas_idx <- substr( str_trim( unlist(strsplit(toString(protein$Sample.REF), split="[,]")) ), 14, 16 ) 
protein=protein[protein$meas_idx=="01A",]
clinic <- read.csv("./TCGA_OV_PROTEIN/clinic.csv")
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
#sort(tcga[tcga$vital_status=="alive",]$observed.time)

# non.feature.list = c("bcr_patient_barcode" ,
#                      "Sample_REF", 
#                      "bcr_patient_uuid",
#                      "patient_id",
#                      "vital_status",
#                      "days_to_death",
#                      "days_to_last_followup",
#                      "observed.time")

#data pre-processing: create design-matrix X, observed-time Z, delta
#------------------------------------------------------------------------------------
X = cbind( tcga[, 3:197], tcga$age_at_initial_pathologic_diagnosis)
n = dim(X)[1]
Z = ( tcga$observed.time+rnorm(n, 0, 1e-5) )/365
delta = as.integer(tcga$vital_status=="dead")
#num_features = p
fold = 10;  unit = floor(n/fold);   base = (fold-1)*unit;
output=NULL
pvalue.cutoff=c()

MAE = NULL; 
RMSE = NULL; 
num_features=5


#Change Restriction time here:
##############################
L =3     
##############################
print(L)
  
for( k in 1:fold ){
  
    set.seed(1)
  
    tr.idx = ( 0:(base-1) + (k-1)*unit ) %% n +1
    
    #print(c(dim(X), tr.idx[1], tr.idx[length(tr.idx)]))
    
    X.tr = X[tr.idx, ];  X.tr <- scale(X.tr)
    Z.tr = Z[tr.idx];  delta.tr = delta[tr.idx]
    
    
    #feature selection
    if(num_features>=1){
      pv = c()
      surv_object <- Surv(time = Z.tr, event = delta.tr)
      for( l in 1:dim(X)[2]){
        cox <- coxph(surv_object ~ ., data = data.frame(X.tr[,l]) )
        pv = c(pv, summary(cox)$waldtest["pvalue"])
      }
      pvalue.cutoff=c(pvalue.cutoff, sort(pv)[10])
      id.features=order(pv)
      id.features= id.features[1:num_features]
    }else{
      id.features = 1:dim(X)[2]
    }
    id.features =sort(id.features)
    X.tr = X.tr[,id.features]
    
    #print(id.features)
    
    #training
    cat('Fold: '); cat(k); cat("...");if(k==fold){cat("\n")}
    id1=order(Z.tr);   Z.tr=Z.tr[id1];  X.tr=X.tr[id1,];   delta.tr = delta.tr[id1]
    Z.tr.res =pmin(Z.tr, L)
    delta.tr.res =delta.tr;  delta.tr.res [Z.tr>=L]=1
    
    surv_object <- Surv(time = Z.tr, event = 1-delta.tr)
    data <- data.frame (x = X.tr)
    cox <- coxph(surv_object ~ ., data = data)
    G =c()
    for(i in 1:dim(X.tr)[1]){
      #print(i)
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
    s.forest = custom_forest(X.tr, Z.tr.res, delta.tr.res, G, num.trees = 2000 )
    
    #test dataset
    te.idx = setdiff(1:n, tr.idx)
    X.test = X[te.idx,  id.features];  X.test <- scale(X.test)
    Z.test = Z[te.idx];  
    delta.test = delta[te.idx];   
    id2=order(Z.test);   Z.test=Z.test[id2];  X.test=X.test[id2,];   delta.test = delta.test[id2]
    Z.test.res = pmin(Z.test, L)
    delta.test.res = delta.test;   delta.test.res[Z.test>=L]=1;
    
    #---------------------------------------
    #prediction
    s.result  = predict(s.forest, X.test,estimate.variance = TRUE)
    naives = naive.estimators( X.tr, Z.tr, delta.tr, X.test,  truth, L  )
    lu = lu.method(X.tr, Z.tr , Z.tr.res, delta.tr, delta.tr.res, X.test, truth, L)
    wang = wang.method(X.tr, Z.tr , Z.tr.res, delta.tr, delta.tr.res, X.test, truth, L, delta.tr.res/(1-G))
    
    #---------------------------------------
    #compute Ghat
    data.test <- data.frame(x = X.test)
    G.test =c()
    for(i in 1:dim(X.test)[1]){
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
    
    
    prediction_matrix =     cbind(s.result$predictions,
                                  #r.result$predictions, 
                                  naives$cox.predictions, 
                                  naives$km.predictions,
                                  lu$lu.id.predictions,
                                  lu$lu.exp.predictions,
                                  wang$wang.id.predictions,
                                  wang$wang.exp.predictions)
    
    curRMSE =     sqrt( colMeans((( prediction_matrix -c(Z.test.res))^2)*c(delta.test.res/(1-G.test))  )  )
    
    curMAE =      colMeans((abs(prediction_matrix -c(Z.test.res)))*c(delta.test.res/(1-G.test))  ) 
    
    
    RMSE = cbind(RMSE,curRMSE)
    MAE = cbind(MAE,curMAE)
}    
    
    
output = cbind(rowMeans(MAE),rowMeans(RMSE))    

row.names(output) <- c('custom.forest' ,
                       #'regression.forest', 
                       'naive.cox.estimator', 
                       'naive.km.estimator',
                       'lu.id.estimator',
                       'lu.exp.estimator',
                       'wang.id.estimator',
                       'wang.exp.estimator'
)
colnames(output)  <- c('mae','rmse')    

print(output)

write.csv(output, paste("/Users/lmy/Dropbox/Mingyang-HLi/projects/survival_random_forest/Revision_material/RealData/", 'L=',L, '.csv', sep = ""), row.names = FALSE)


labels = c('custom.forest' ,'naive.cox.estimator', 'naive.km.estimator','lu.id.estimator','lu.exp.estimator','wang.id.estimator','wang.exp.estimator')
boxplot(  data.frame(t(RMSE)) , ylim = c(quantile(RMSE,probs = c(0,1))), xaxt = "n",  xlab = "",ylab=paste('RMSE(L=', L, ')') )
axis(1, labels = FALSE)
text(x =  seq_along(labels), y = par("usr")[3] - (par("usr")[4] - par("usr")[3])/30, srt = 45, adj = 1,labels = labels, xpd = TRUE)


labels = c('custom.forest' ,'naive.cox.estimator', 'naive.km.estimator','lu.id.estimator','lu.exp.estimator','wang.id.estimator','wang.exp.estimator')
boxplot(  data.frame(t(MAE)) , ylim = c(quantile(MAE,probs = c(0, 1))), xaxt = "n",  xlab = "",ylab=paste('MAE(L=', L, ')') )
axis(1, labels = FALSE)
text(x =  seq_along(labels), y = par("usr")[3] - (par("usr")[4] - par("usr")[3])/30, srt = 45, adj = 1,labels = labels, xpd = TRUE)

