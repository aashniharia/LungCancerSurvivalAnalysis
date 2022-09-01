###
#Lung Cancer prediction of overall Survival -for section 4.1
####

#########
# Without using prefiltering
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
## data visualization
###############################################################################

dat$sex=as.factor(data_without_prefiltering$sex)
dat$surind=as.factor(dat$surind)
class(dat$stg)
class(dat$sex)
class(dat$surind)

#visualizing survival time based on stage
graph_1 = ggplot(dat, aes(x = stg, y = surtim)) +  geom_boxplot(fill = 'purple')
graph_1 + ggtitle("Survival Time Based On Stage") +
  xlab("Cancer Stage") + ylab("Survival Time(Days)")

#visualizing survival time based on sex
ggplot(dat, aes(x = sex, y = surtim, fill = sex)) +  geom_boxplot() + ggtitle("Survival Time Based On Gender") +
  xlab("Patient Gender") + ylab("Survival Time(Days)") + scale_fill_discrete(name = "Sex", labels = c("Female", "Male"))




#visualizing stage based on sex
ggplot(dat, aes(x = stg,fill = sex)) +
  geom_bar(position = position_dodge(preserve = "single")) + scale_fill_discrete(name = "Sex", labels = c("Female", "Male")) +
  ggtitle("Count of Stage Based on Gender") +
  xlab("Cancer Stage") + ylab("Count")



#visualizing the number of patients in each stage after imputatiom
ggplot(dat, aes(x = stg)) +   geom_bar(fill = "cornflowerblue", color = "black") +
  labs(x = 'Cancer Stage', y = "Count", title  = "Patients By Stage")


#visualizing the number of patients based on gender 
ggplot(dat, aes(x = sex, fill = sex)) +
  geom_bar() +  labs(x = 'Patient Gender', y = "Count", title  = "Patients By Gender") +
  scale_fill_discrete(name = "Sex", labels = c("Female", "Male"))


#visualizing the number of patients based on surind 
ggplot(dat, aes(x = surind,fill=surind)) +
  geom_bar() +  labs(x = 'Survival Index', y = "Count", title  = "Survival Index") +
  scale_fill_discrete(labels = c("Alive", "Dead"))+
  scale_fill_brewer(palette = "Accent")


#distribution of survival time
# ggplot(dat, aes(x = surtim,fill=surtim))+geom_histogram(binwidth=3,color="black", fill="lavender")

hist(dat$surtim,col="lavender",xlab="Survival time (Days)",main="Distribution of Survival Time(Days) ")


hist(dat$age[dat$sex==0],main="Histogram for Age of Females Patients", 
     xlab="Age", 
     border="black", 
     col="pink",)
hist(dat$age[dat$sex==1],main="Histogram for Age of Male Patients", 
     xlab="Age", 
     border="black", 
     col="cornflowerblue",)
#inference: both the histograms are skewed and thus we would impute age using median values


# Kaplan Meier Curves
response_y = Surv(time = dat$surtim, event = dat$surind)
length(response_y)
is.Surv(response_y)
# estimating survival without using any x variables
km.mod = survfit(response_y ~ 1,data=dat)

#model fit
km.mod
#n -> total number of individuals
#events -> 58 deaths
#median-> median survival time
# 95% ci for median lower and upper

length(km.mod)
#summary of the model
summary(km.mod)


# plotting the kaplan meier curve without and X variables
ggsurvplot(
  km.mod,
  conf.int = T,
  xlab = "Time(Days)",
  ylab = " % Alive= S(t)",
  title = "Overall Survival Curve",
  surv.median.line="hv"
  
)
median(dat$surtim)
class(dat$stg)
km_curve1=survfit(response_y~stg,data=dat)
ggsurvplot(km_curve1,pval=TRUE,title="Survival by Cancer Stage",xlab = "Time(Days)",conf.int = TRUE,risk.table=TRUE,
           legend.labs = c("Stage 1", "Stage 2", "Stage 3"),)


class(dat$sex)
dat$sex=as.factor(dat$sex)
km_curve2=survfit(response_y~sex,data=dat)
ggsurvplot(km_curve2,pval=TRUE,conf.int = TRUE,risk.table=TRUE,title="Survival by Gender",legend.labs = c("Female", "Male"),xlab = "Time(Days)")

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
# Models without prefiltering
###############################################################################
data_without_prefiltering=dat
repetitions = 10
folds = 5

iterationCount = folds*repetitions
threshold = (50*iterationCount)/100

subset_cols = ncol(data_without_prefiltering)-2
subset_rows = nrow(data_without_prefiltering)
subset_names = colnames(subset(data_without_prefiltering, select = -c(surtim, surind)))


lasso.train.conc = lasso.test.conc = numeric(iterationCount)
lasso.features.selected = matrix(0, nrow=iterationCount, ncol=subset_cols)
colnames(lasso.features.selected) = subset_names
lasso_prediction_matrix  = matrix(0, nrow=iterationCount, ncol=nrow(data_without_prefiltering))


ridge.train.conc = ridge.test.conc =  numeric(iterationCount)
ridge_coefficients_matrix = matrix(0, nrow=iterationCount, ncol=subset_cols)
colnames(ridge_coefficients_matrix) = subset_names
ridge_predection_matrix  = matrix(0, nrow=iterationCount, ncol=nrow(data_without_prefiltering))

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


n=nrow(data_without_prefiltering)
model_counter = 1
set.seed(2)
for(eachRepetition in 1:repetitions){
  shuffled_data = data_without_prefiltering[sample(1:nrow(data_without_prefiltering),nrow(data_without_prefiltering)), ]
  #print(dim(shuffled_data))
  cvIndex = createFolds(factor(shuffled_data$surind), folds,returnTrain = T)
   for(eachFold in 1:length(cvIndex)){
    
    train_data = shuffled_data[cvIndex[[eachFold]],]
    test_data = shuffled_data[-cvIndex[[eachFold]],]

    x.train = subset(train_data, select = -c(surtim, surind))
    y.train = subset(train_data, select = c(surtim, surind))

    surv.train = Surv(y.train$surtim,y.train$surind)
    x.train.m = model.matrix(surv.train~.+0,data=x.train)
 
	x.test = subset(test_data, select = -c(surtim, surind))
    x.test.m = model.matrix(surv.test ~.+0,data=x.test)
 
    y.test = subset(test_data, select = c(surtim, surind))
    surv.test = Surv(y.test$surtim,y.test$surind)
    
    ####
	#Ridge Cox
	#####
	
    ridge_optimal_lambda = cv.glmnet(x.train.m,  surv.train, family = "cox",alpha=0,type.measure ="C", folds=folds,maxit=100000)
    ridge_model_fit = glmnet(x.train.m, surv.train,family ="cox", alpha=0,lambda = ridge_optimal_lambda$lambda.min)
    ridge_train_predictions = predict(ridge_model_fit,newx=x.train.m,type="link")[,1] # the fitted relative-risk for "cox";
    ridge_test_predictions = predict(ridge_model_fit,newx=x.test.m,type="link")[,1] # the fitted relative-risk for "cox";
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
	#Gradient Boosting Machine
	#####
    gbm_model_fit = gbm(surv.train~., data=x.train, distribution="coxph")
    gbm_Model_Summary[[model_counter]] = summary(gbm_model_fit)
    gbm_train_predictions = predict(gbm_model_fit, newdata=train_data, type="link")
    gbm_test_predictions = predict(gbm_model_fit, newdata=test_data, type="link")
	gbm_indexs = rownames(test_data)
	gbm.test.conc[model_counter] = Cindex(gbm_test_predictions, y=surv.test)
    gbm.train.conc[model_counter] = Cindex(gbm_train_predictions, y=surv.train)
    #creating the predictions matrix for gbm
    for(eachGbmIndex in 1:length(gbm_indexs)){
      gbm.index.pos = as.numeric(gbm_indexs[eachGbmIndex])
      gbm_prediction_matrix[model_counter,gbm.index.pos] = gbm_train_predictions[eachGbmIndex]
    }
	#storing the relative influence
    for(row_var in 1:iterationCount ){
      for(col_var in 1:subset_cols){
        gbm_index_var=which(subset_names==gbm_Model_Summary[[row_var]][col_var,]$var)
        gbm_coefffiecients_matrix[row_var,gbm_index_var]=gbm_Model_Summary[[row_var]][col_var,]$rel.inf
      }
    }

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

    model_counter = model_counter + 1
  }
}

####
#Recursive Feature Elimination with Random Survival Forest
#####
set.seed(2)
dim(data_without_prefiltering)
rfe.train.conc = vector("integer",subset_cols)
rfe.test.conc = vector("integer",subset_cols)
chosen_variables = list()

model_used = data_without_prefiltering
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
boxplot( lasso.train.conc, lasso.test.conc,
         ridge.train.conc,ridge.test.conc,gbm.train.conc,gbm.test.conc, unlist(rsf.train.conc),
         unlist(rsf.test.conc),rfe.train.conc,rfe.test.conc ,  col= c(0,6,0,6,0,6,0,6,0,6,0,6) ,
         names=c("Lasso","","Ridge","", "GBM","","RSF","","RFE+RSF",""),
         main="Train and Test concordances")
legend("topright",c("train concordance","test concordance"),fill=c(0,6),inset = c(-0.1,-0.1))
