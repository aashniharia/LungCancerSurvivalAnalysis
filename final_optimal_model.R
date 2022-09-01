###
#Lung Cancer prediction of overall Survival -for section 4.3
####

#########
# Finding Optimal Model for Lung Cancer Dataset
#########
rm(list = ls())
library(survival)
library(caret)
library(plotrix)
library(survminer)
library(glmnet)
library(gbm)
library(randomForestSRC)
###############################################################################
## data pre-processing
###############################################################################


#reading the data
dat = read.csv("C:\\MSDSA\\SEM3\\lung.csv", stringsAsFactors = TRUE)


#removing variables of no interest
dat <- subset(dat, select = c(surtim,surind,stg,run.entropy_GLSZM,gradb1,max.seg,sex,het1,qcod_HIST,energy))



#converting stage to factor
dat$stg <- revalue(dat$stg, c("1"=1,"1A"=1,"1B"=1,"2A"=2,"2B"=2,"3A"=2,"3B"=2,"4"=3))
dat$stg=as.factor(dat$stg)
#class(dat$stg)
#levels(dat$stg)

#imputing the age
mm <- median(dat$age[dat$sex==1], na.rm = TRUE)
mf <- median(dat$age[dat$sex==0] , na.rm = TRUE)
dat$age[is.na(dat$age) & dat$sex==1] <- mm
dat$age[is.na(dat$age) & dat$sex==0] <- mf



data_after_prefiltering=dat
data_after_prefiltering$stg=as.numeric(data_after_prefiltering$stg)

###############################################################################
## Function for plotting KM
###############################################################################

split.kms <- function(zos,os,oevn,NQ=100,zqmin=.05,zqmax=.95){
  qmin = quantile(zos,zqmin,na.rm=TRUE)
  qmax = quantile(zos,zqmax,na.rm=TRUE)
  p0q = numeric(NQ)
  med0s = seq(qmin,qmax,length=NQ)
  izs = matrix(NA,ncol=NQ,nrow=length(zos))
  for(iq in 1:NQ){
    IZ = as.integer(zos<med0s[iq])
    p0q[iq] = 1-pchisq(survdiff(Surv(os,oevn)~factor(IZ))$chisq,1)
    izs[,iq] = IZ
  }
  best.med0 = med0s[which.min(p0q)]
  IZ = as.integer(zos<best.med0)
  return(list(iz=IZ,izs=izs,best.t=best.med0,
              ps=p0q,best.p=p0q[which.min(p0q)]))
}

###############################################################################
# Model without prefiltering
###############################################################################

repetitions = 10
folds = 5

iterationCount = folds*repetitions
threshold = (50*iterationCount)/100
subset_cols = ncol(data_after_prefiltering)-2
subset_rows = nrow(data_after_prefiltering)
subset_names = colnames(subset(data_after_prefiltering, select = -c(surtim, surind)))

stepCox_features_selected = matrix(0,nrow=iterationCount,ncol=subset_cols)
colnames(stepCox_features_selected) = subset_names
stepCox_predictions_matrix = matrix(0,nrow=iterationCount,ncol=subset_rows)



final.stepCox.train.conc  = vector("list",iterationCount)
final.stepCox.test.conc = vector("list",iterationCount)
final.stepCox.models = vector("list",iterationCount)


lasso.train.conc = lasso.test.conc = numeric(iterationCount)
lasso.features.selected = matrix(0, nrow=iterationCount, ncol=subset_cols)


colnames(lasso.features.selected) = subset_names
lasso_prediction_matrix  = matrix(0, nrow=iterationCount, ncol=nrow(data_after_prefiltering))


ridge.train.conc = ridge.test.conc =  numeric(iterationCount)
ridge_coefficients_matrix = matrix(0, nrow=iterationCount, ncol=subset_cols)
colnames(ridge_coefficients_matrix) = subset_names
ridge_predection_matrix  = matrix(0, nrow=iterationCount, ncol=nrow(data_after_prefiltering))

gbm_Model_Summary = vector("list",iterationCount) 
gbm.train.conc = gbm.test.conc = vector("double",iterationCount) 
gbm_coefffiecients_matrix = matrix(0, nrow=iterationCount, ncol= subset_cols)
colnames(gbm_coefffiecients_matrix) = subset_names
gbm_prediction_matrix = matrix(0,nrow=iterationCount,ncol=subset_rows)

rsf_model_fits = rsf_predictions = vector("list",iterationCount)
rsf.train.conc = rsf.test.conc= vector("list",iterationCount)
rsf.test.conc = vector("list",iterationCount)
rsf.test.conc = vector("list",iterationCount)
rsf.features.selected = matrix(0,nrow=iterationCount,ncol=subset_cols)
colnames(rsf.features.selected) = subset_names


n=nrow(data_after_prefiltering)
model_counter = 1
set.seed(2)
for(eachRepetition in 1:repetitions){
  shuffled_data = data_after_prefiltering[sample(1:nrow(data_after_prefiltering),nrow(data_after_prefiltering)), ]
  print(dim(shuffled_data))
  cvIndex = createFolds(factor(shuffled_data$surind), folds,returnTrain = T)
   for(eachFold in 1:length(cvIndex)){
    
     train_data = shuffled_data[cvIndex[[eachFold]],]
     test_data = shuffled_data[-cvIndex[[eachFold]],]

    
     x.train = subset(train_data, select = -c(surtim, surind))
     y.train = subset(train_data, select = c(surtim, surind))

    surv.train = Surv(y.train$surtim,y.train$surind)
    x.train.m = model.matrix(surv.train~.+0,data=x.train)

    y.test = subset(test_data, select = c(surtim, surind))
    surv.test = Surv(y.test$surtim,y.test$surind)
    x.test = subset(test_data, select = -c(surtim, surind))
    x.test.m = model.matrix(surv.test ~.+0,data=x.test)
	
	

    ####
	# Cox Proportional Hazard
	#####
    step_cox_fit =  coxph(Surv(surtim,surind)~.,data=train_data)
    stp_cox_predictions = predict(step_cox_fit,newdata = test_data,type="lp")
    final.stepCox.train.conc[model_counter]<-concordance(step_cox_fit, newdata = train_data)$concordance
    final.stepCox.test.conc[model_counter] <- concordance(step_cox_fit, newdata = test_data)$concordance
	
	####
	#Ridge Cox
	#####
	
    ridge_optimal_lambda = cv.glmnet(x.train.m,  surv.train, family = "cox",alpha=0,type.measure ="C", folds=folds,maxit=100000)
    ridge_model_fit = glmnet(x.train.m, surv.train,family ="cox", alpha=0,lambda = ridge_optimal_lambda$lambda.min)
    ridge_train_predictions = predict(ridge_model_fit,newx=x.train.m,type="link")[,1] 
    ridge_test_predictions = predict(ridge_model_fit,newx=x.test.m,type="link")[,1] 
    for(eachRidgePrediction in 1:length(ridge_test_predictions)){
      ridge_temp_var = ridge_test_predictions[eachRidgePrediction]
      ridge_index = as.numeric(names(ridge_temp_var))
      ridge_predection_matrix[model_counter,ridge_index] = as.numeric(ridge_temp_var)[1]
    }
    ridge.train.conc[model_counter] =  Cindex(ridge_train_predictions,y= surv.train)
    ridge.test.conc[model_counter] =  Cindex(ridge_test_predictions,y= surv.test)
    ridge.coefs = coef(ridge_model_fit)[,1]
    for(eachRidgeCoef in 1:length(ridge.coefs)){
      ridge_coefficients_matrix[model_counter,eachRidgeCoef] = abs(as.numeric(ridge.coefs[eachRidgeCoef]))
    }
    
	####
	#LASSO Cox
	#####

    lasso_optimal_lambda = cv.glmnet(x.train.m,  surv.train, family = "cox",type.measure ="C", folds=folds,maxit=100000)
    lasso_model_fit = glmnet(x.train.m, surv.train,family ="cox", lambda = lasso_optimal_lambda$lambda.min)
    lasso_train_predictions = predict(lasso_model_fit,newx=x.train.m,type="link")[,1]
    lasso_test_predictions = predict(lasso_model_fit,newx=x.test.m,type="link")[,1] # the fitted relative-risk for "cox";
    for(eachLassoPrediction in 1:length(lasso_test_predictions)){
      lasso_var_temp = lasso_test_predictions[eachLassoPrediction]
      lasso_index = as.numeric(names(lasso_var_temp))
      lasso_prediction_matrix[model_counter,lasso_index] = as.numeric(lasso_var_temp)[1]
    }
    lasso.train.conc[model_counter] =  Cindex(lasso_train_predictions,y= surv.train)
    lasso.test.conc[model_counter] =  Cindex(lasso_test_predictions,y= surv.test)
    lasso.features.selected[model_counter,] = as.numeric(coef(lasso_model_fit)!=0)


   

	####
	#Random Survival Forest
	#####
    rsf_model = rfsrc(Surv(surtim, surind) ~ .,data = train_data, importance = TRUE)
    rsf_model_predictions = predict.rfsrc(rsf_model,newdata = test_data,type="lp")
	rsf_predictions[[model_counter]] = rsf_model_predictions
    rsf_model_fits[[model_counter]] = rsf_model
    
    #
	rsf.test.conc[model_counter] = Cindex(rsf_model_predictions$predicted,y=surv.test)
    rsf.train.conc[model_counter] = Cindex(rsf_model$predicted.oob,y=surv.train)
    
    #
    rsf_var_imp = subsample(rsf_model)
    rsf_features = head(subset_names[order(rsf_model$importance, decreasing=TRUE)],10)
    for(eachRsfFeature in 1:length(rsf_features)){
      rsf.features.selected[model_counter,rsf_features[eachRsfFeature]] = 1
    }
    ####
	#Gradient Boosting Machine
	#####
    gbm_model_fit = gbm(surv.train~., data=x.train, distribution="coxph")
    gbm_Model_Summary[[model_counter]] = summary(gbm_model_fit)
    gbm_train_predictions = predict(gbm_model_fit, newdata=train_data, type="link")
    gbm_test_predictions = predict(gbm_model_fit, newdata=test_data, type="link")
    gbm.train.conc[model_counter] = Cindex(gbm_train_predictions, y=surv.train)
    gbm.test.conc[model_counter] = Cindex(gbm_test_predictions, y=surv.test)
    gbm_indexs = rownames(test_data)
    for(eachGbmIndex in 1:length(gbm_indexs)){
      gbm.index.pos = as.numeric(gbm_indexs[eachGbmIndex])
      gbm_prediction_matrix[model_counter,gbm.index.pos] = gbm_train_predictions[eachGbmIndex]
    }

    for(row_var in 1:iterationCount ){
      for(col_var in 1:subset_cols){
        gbm_index_var=which(subset_names==gbm_Model_Summary[[row_var]][col_var,]$var)
        gbm_coefffiecients_matrix[row_var,gbm_index_var]=gbm_Model_Summary[[row_var]][col_var,]$rel.inf
      }
    }

    model_counter = model_counter + 1
  }
}



#Calculating concordances for Cox Proportional Hazards Model
round(median(unlist(final.stepCox.train.conc)),4)
round(std.error(unlist(final.stepCox.train.conc)),4)
round(median(unlist(final.stepCox.test.conc)),4)
round(std.error(unlist(final.stepCox.test.conc)),4)

#Calculating concordances for LASSO Cox
round(median(lasso.train.conc),4)
round(std.error(lasso.train.conc),4)
round(median(lasso.test.conc),4)
round(std.error(lasso.test.conc),4)

#Calculating concordances for Ridge Cox
round(median(ridge.train.conc),4)
round(std.error(ridge.train.conc),4)
round(median(ridge.test.conc),4)
round(std.error(ridge.test.conc),4)

#Calculating concordances for Gradient Boosting Machine
round(median(gbm.train.conc),4)
round(std.error(gbm.train.conc),4)
round(median(gbm.test.conc),4)
round(std.error(gbm.test.conc),4)

#Calculating concordances for Random Survival Forest
round(median(unlist(rsf.train.conc)),4)
round(std.error(unlist(rsf.train.conc)),4)
round(median(unlist(rsf.test.conc)),4)
round(std.error(unlist(rsf.test.conc)),4)


##
#Visualisation of concordances
##

par(mar = c(2, 2, 2, 4), xpd = TRUE)
boxplot( unlist(final.stepCox.train.conc),unlist(final.stepCox.test.conc),lasso.train.conc, lasso.test.conc,
         ridge.train.conc,ridge.test.conc,gbm.train.conc,gbm.test.conc, unlist(rsf.train.conc),
         unlist(rsf.test.conc),  col= c(0,"#FFCCE5",0,"#FFCCE5",0,"#FFCCE5",0,"#FFCCE5",0,"#FFCCE5",0,"#FFCCE5") ,
         names=c("Cox","","LASSO","","Ridge","", "GBM","","RSF",""),
         main="Train and Test concordances")
legend("topright",c("train concordance","test concordance"),fill=c(0,"#FFCCE5"),inset = c(-0.1,-0.13))


##################################
## Kaplan Meier Curves
#################################

Cox_mean_predictions = as.array(apply(stepCox_predictions_matrix,2,mean))
stp_kmo = split.kms(Cox_mean_predictions, os=data_after_prefiltering$surtim, oevn=data_after_prefiltering$surind)
stp_temp_var = stp_kmo$iz
stp_km.split = survfit(Surv(data_after_prefiltering$surtim, data_after_prefiltering$surind)~stp_temp_var, data=data_after_prefiltering)
ggsurvplot(stp_km.split,title="Overall Survival For Cox Model",pval = TRUE,risk.table=TRUE,conf.int = TRUE,legend.labs =
             c("high risk", "low risk"))


lasso_mean_predictions = as.array(apply(lasso_prediction_matrix,2,mean))
lasso_kmo = split.kms(lasso_mean_predictions, os=data_after_prefiltering$surtim, oevn=data_after_prefiltering$surind)
lasso_temp_var = lasso_kmo$iz
lasso_km.split = survfit(Surv(data_after_prefiltering$surtim, data_after_prefiltering$surind)~lasso_temp_var, data=data_after_prefiltering)
ggsurvplot(lasso_km.split,title="Overall Survival for Penalised Cox Model(LASSO)",pval = TRUE,risk.table=TRUE,conf.int = TRUE,legend.labs =
             c("high risk", "low risk"))


ridge_mean_predictions = as.array(apply(ridge_predection_matrix,2,mean))
ridge_kmo = split.kms(ridge_mean_predictions, os=data_after_prefiltering$surtim, oevn=data_after_prefiltering$surind)
ridge_temp_var = ridge_kmo$iz
ridge_km.split = survfit(Surv(data_after_prefiltering$surtim, data_after_prefiltering$surind)~ridge_temp_var, data=data_after_prefiltering)
ggsurvplot(ridge_km.split,title="Overall Survival for Penalised Cox Model(Ridge)",pval = TRUE,risk.table=TRUE,conf.int = TRUE,legend.labs =
             c("high risk", "low risk"))

gbm_mean_predictions = as.array(apply(gbm_prediction_matrix,2,mean))
gbm_kmo = split.kms(gbm_mean_predictions, os=data_after_prefiltering$surtim, oevn=data_after_prefiltering$surind)
gbm_temp_var = gbm_kmo$iz
gbm_km.split = survfit(Surv(data_after_prefiltering$surtim, data_after_prefiltering$surind)~gbm_temp_var, data=data_after_prefiltering)
ggsurvplot(gbm_km.split,title="Overall Survival for Boosting Cox Model",pval = TRUE,conf.int = TRUE,risk.table=TRUE,legend.labs =
             c("high risk", "low risk"))
