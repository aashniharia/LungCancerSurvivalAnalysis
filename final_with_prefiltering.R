###
#Lung Cancer prediction of overall Survival -for section 4.2
####

#########
# Using Prefiltering
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
dat <- subset(dat, select = -c(inj, wt, dimx, dimy, dimz,ptid,X,id,xid,range_HIST,max_HIST,min_HIST))



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



#convert to numeric
dat$sex =( as.numeric(dat$sex))
dat$stg =( as.numeric(dat$stg))






###############################################################################
## removing correlated variables
###############################################################################

#Removing highly correlated variables

cor_dat=subset(dat, select = -c(surtim,surind,stg))
cor_list=list()
x=cor(cor_dat)
cor_list=findCorrelation(
  x,
  cutoff = 0.9,
  verbose = FALSE,
  names = FALSE,
  # exact = ncol(x) < 100
)

#eliminating the correlated variables

remaining_dat=subset(cor_dat,select=-c(cor_list))
all_remaining_predictors=subset(dat,select=c(stg))

remaining_dat_stg=subset(dat,select=c(surtim,surind,stg))


#all remaining predictors excluding surtim and surind
all_predictors=cbind(all_remaining_predictors,remaining_dat)

#all the remaining features/variables
remaining_filtered_data=cbind(remaining_dat_stg,remaining_dat)



data_after_prefiltering=remaining_filtered_data
data_after_prefiltering$stg=as.numeric(data_after_prefiltering$stg)
data_after_prefiltering$sex=as.numeric(data_after_prefiltering$sex)
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
subset_rows = nrow(data_after_prefiltering)
subset_cols = ncol(data_after_prefiltering)-2
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
	#Stepwise Selection + Cox Proportional Hazard
	#####
	start_cox = coxph(Surv(surtim,surind) ~ 1, data = train_data)
    full_cox = coxph(Surv(surtim,surind) ~ ., data = train_data)
    fit_step = step(start_cox, direction = "both", scope = full_cox$formula)
    stpCox_names = names(fit_step$coefficients)
    for(eachStepName in 1:length(stpCox_names)){
      stepCox_features_selected[model_counter,stpCox_names[eachStepName]] = 1
    }
   
    full_form = fit_step$formula[3]
    full_form_formula = formula(paste("Surv(surtim,surind)~",full_form))

    #fitting cox model on features selected by stepwise
    step_cox_fit = coxph(full_form_formula,data=train_data)
    stp_cox_predictions = predict(step_cox_fit,newdata = test_data,type="lp")

    for(eachPred in 1:length(stp_cox_predictions)){
      stp_temp_var = stp_cox_predictions[eachPred]
      step_indx = as.numeric(names(stp_temp_var))
      stepCox_predictions_matrix[model_counter,step_indx] = as.numeric(stp_temp_var)[1]
    }

    final.stepCox.models[model_counter] <- full_form
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
    lasso_test_predictions = predict(lasso_model_fit,newx=x.test.m,type="link")[,1] 
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
####
#Recursive Feature Elimination with Random Survival Forest
#####
set.seed(2)
dim(data_after_prefiltering)
rfe.train.conc = vector("integer",subset_cols)
rfe.test.conc = vector("integer",subset_cols)
chosen_variables = list()

model_used = dat
rfe_counter = 1
while(length(model_used)>2){
  shuffled_orig_dat = model_used[sample(1:nrow(model_used),nrow(model_used)), ]
  cvIndex_rfe = createFolds(factor(shuffled_orig_dat$surind), folds,returnTrain = T)
  init_train_data = shuffled_orig_dat[cvIndex_rfe[[eachFold]],]
  init_test_data = shuffled_orig_dat[-cvIndex_rfe[[eachFold]],]
  init.surv.train = Surv(init_train_data$surtim,init_train_data$surind)
  init.surv.test = Surv(init_test_data$surtim,init_test_data$surind)
  RFSRC_mod_fit=rfsrc(Surv(surtim,surind) ~ .,data = init_train_data, importance = TRUE)
  pred = predict.rfsrc(RFSRC_mod_fit,newdata = init_test_data,type="lp")
  rfe.train.conc[rfe_counter] = Cindex(RFSRC_mod_fit$predicted.oob, y = init.surv.train)
  rfe.test.conc[rfe_counter] = Cindex(pred$predicted, y = init.surv.test)
  response_subset=subset(model_used, select = c(surtim,surind))
  var_imp=vimp.rfsrc(RFSRC_mod_fit)$importance
  sorted_vars= sort(var_imp,decreasing=TRUE)
  selected_vars = head(sorted_vars,-1)
  chosen_variables[[rfe_counter]] = toString(names(selected_vars))
  new_vars_subset=subset(model_used,select=c(names(selected_vars)))
  next_iteration_subset=cbind(new_vars_subset,response_subset)
  model_used=next_iteration_subset
  rfe_counter = rfe_counter + 1
}

# boxplot(rfe.train.conc, rfe.test.conc)
##
chosen_variables = unlist(chosen_variables)
max_test_conc = which.max(rfe.test.conc)
final_rfe_features = chosen_variables[max_test_conc]

####
#most predective features of models
####

final_lasso_features = which(colSums(lasso.features.selected)>threshold) == TRUE

ridge_colsum = colSums(ridge_coefficients_matrix)

final_ridge_features = names(head(ridge_colsum[order(-ridge_colsum)],6))

gbm_colsum = colSums(gbm_coefffiecients_matrix)
final_gbm_features = names(head(gbm_colsum[order(-gbm_colsum)],6))

rsf_colsum=colSums(rsf.features.selected)
final_rsf_features= names(head(rsf_colsum[order(-rsf_colsum)],6))

#Calculating concordances for Stepwise selectiction + Cox Proportional Hazards Model
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

#Calculating concordances for Recursive Feature Elimination + RSF
round(median(rfe.train.conc),4)
round(std.error(rfe.train.conc),4)
round(median(rfe.test.conc),4)
round(std.error(rfe.test.conc),4)


##
#Visualisation of concordances
##
par(mar = c(2, 2, 2, 4), xpd = TRUE)
boxplot( unlist(final.stepCox.train.conc),unlist(final.stepCox.test.conc),lasso.train.conc, lasso.test.conc,
         ridge.train.conc,ridge.test.conc,gbm.train.conc,gbm.test.conc, unlist(rsf.train.conc),
         unlist(rsf.test.conc),rfe.train.conc,rfe.test.conc ,  col= c(0,7,0,7,0,7,0,7,0,7,0,7) ,
         names=c("Cox","","LASSO","","Ridge","", "GBM","","RSF","","RFE+RSF",""),
         main="Train and Test concordances")
legend("topright",c("train concordance","test concordance"),fill=c(0,7),inset = c(-0.1,-0.1))

##################################
## Kaplan Meier Curves
#################################

stepCox_mean_predictions = as.array(apply(stepCox_predictions_matrix,2,mean))
stp_kmo = split.kms(stepCox_mean_predictions, os=data_after_prefiltering$surtim, oevn=data_after_prefiltering$surind)
stp_temp_var = stp_kmo$iz
stp_km.split = survfit(Surv(data_after_prefiltering$surtim, data_after_prefiltering$surind)~stp_temp_var, data=data_after_prefiltering)
ggsurvplot(stp_km.split,title="Overall Survival For Cox Model",pval = TRUE,conf.int = TRUE,surv.median.line="hv",risk.table=TRUE,break.x.by = 100,legend.labs =
             c("high risk", "low risk"))


lasso_mean_predictions = as.array(apply(lasso_prediction_matrix,2,mean))
lasso_kmo = split.kms(lasso_mean_predictions, os=data_after_prefiltering$surtim, oevn=data_after_prefiltering$surind)
lasso_temp_var = lasso_kmo$iz
lasso_km.split = survfit(Surv(data_after_prefiltering$surtim, data_after_prefiltering$surind)~lasso_temp_var, data=data_after_prefiltering)
ggsurvplot(lasso_km.split,title="Overall Survival for Penalised Cox Model(LASSO)",pval = TRUE,surv.median.line="hv",break.x.by = 100,risk.table=TRUE,conf.int = TRUE,legend.labs =
             c("high risk", "low risk"))


ridge_mean_predictions = as.array(apply(ridge_predection_matrix,2,mean))
ridge_kmo = split.kms(ridge_mean_predictions, os=data_after_prefiltering$surtim, oevn=data_after_prefiltering$surind)
ridge_temp_var = ridge_kmo$iz
ridge_km.split = survfit(Surv(data_after_prefiltering$surtim, data_after_prefiltering$surind)~ridge_temp_var, data=data_after_prefiltering)
ggsurvplot(ridge_km.split,title="Overall Survival for Penalised Cox Model(Ridge)",pval = TRUE,risk.table=TRUE,conf.int = TRUE,surv.median.line="hv",break.x.by = 100,legend.labs =
             c("high risk", "low risk"))

gbm_mean_predictions = as.array(apply(gbm_prediction_matrix,2,mean))
gbm_kmo = split.kms(gbm_mean_predictions, os=data_after_prefiltering$surtim, oevn=data_after_prefiltering$surind)
gbm_temp_var = gbm_kmo$iz
gbm_km.split = survfit(Surv(data_after_prefiltering$surtim, data_after_prefiltering$surind)~gbm_temp_var, data=data_after_prefiltering)
ggsurvplot(gbm_km.split,title="Overall Survival for Boosting Cox Model",pval = TRUE,surv.median.line="hv",break.x.by = 100,conf.int = TRUE,risk.table=TRUE,legend.labs =
             c("high risk", "low risk"))


