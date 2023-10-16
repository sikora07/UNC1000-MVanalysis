library(readr)
library(ggplot2)
library(mice)
library(pROC)
library(caret)

yesno_to_01<- function(sel_vector){
  sel_vector[which(sel_vector=='Yes')] = 1
  sel_vector[which(sel_vector == 'No')] = 0
  sel_vector = factor(sel_vector,levels = c(0,1))
  return(sel_vector)
}

####################################
## import dataset
raw_data = read_csv("MV_dataclean_0906.csv")
raw_data = raw_data[which( raw_data$dur_MV_days_0_included >= 1 ),] ## population condition
dim(raw_data)

## split dataset
set.seed(123)
train_idx  = sample(1:nrow(raw_data),nrow(raw_data)*0.9, replace=F)


outcome_name = 'prolong_mv'

## basic variables that Andrea selected 
basic_variables = c( 'Dev_duration', 'fluid_overload_MVday0',
                      'HCO3_lt20', 
                     'max_FiO2', 'max_PEEP', 'PaFi', 
                     'PH_lt725Day1', 'Pulse_gt110',
                     'max_PaCO2','max_minute_vent',
                     'max_vent_PaCO2')

true_basic_var = basic_variables

## MRC concerned variables 
raw_data$MRCday0
raw_data$MRCday0_ge10
raw_data$MRCday1
raw_data$MRCday1_ge10
raw_data$score_24h 
raw_data$MRC_ICU24_ge10

## compared score
raw_data$APACHE_24h_score
raw_data$SOFA_24h_score


###################################################
## extract dataset with only basic variables & format transformation
basic_raw = data.frame( raw_data[,c(outcome_name, basic_variables)] )

basic_raw$prolong_mv = as.factor(basic_raw$prolong_mv)
class(basic_raw$Dev_duration)
class( basic_raw$fluid_overload_MVday0 )
basic_raw$fluid_overload_MVday0  = yesno_to_01( basic_raw$fluid_overload_MVday0 )
basic_raw$fluid_overload_MVday1 = yesno_to_01(basic_raw$fluid_overload_MVday1)
basic_raw$HCO3_lt20 = yesno_to_01( basic_raw$HCO3_lt20 )
basic_raw$max_FiO2 
basic_raw$max_PEEP
basic_raw$PaFi
basic_raw$PH_lt725Day1 = yesno_to_01( basic_raw$PH_lt725Day1 )
basic_raw$Pulse_gt110 = yesno_to_01( basic_raw$Pulse_gt110 )
basic_raw$resp  ## 24 hours
basic_raw$min_vent_PaCO2 

str(basic_raw)
colnames(basic_raw)

## normalization
for(i in 1:ncol(basic_raw)){
  if( class(basic_raw[,i]) == 'numeric'){
    basic_raw[,i] = as.numeric(scale(basic_raw[,i]))
  }
}

tra_std = basic_raw[train_idx,]
test_std = basic_raw[-train_idx,]
#######################################################
## univariate p value
all( colnames(basic_raw)[-1] == basic_variables )

p_univ = c()
for(idx in 1:length(basic_variables)){
  xname = basic_variables[idx]
  if( class( basic_raw[,outcome_name] )=='factor' ){
    uni_model = summary(glm(formula( paste0(paste0(outcome_name, '~'), xname) ),
                    family = binomial(), data = basic_raw))
    p_univ = c( p_univ, uni_model$coefficients[-1,4] )
  }
}

univ_df = data.frame(basic_var = basic_variables, pvalue = p_univ)
univ_df$sig_var = ifelse(univ_df$pvalue<0.05, 1, 0 )

write.csv(univ_df, 'results/basic/univ_basic_var_pvalue.csv')

#######################################################
## include all basic variables, with no variable selection step
## impute dataset five times
## for each imputation, extract a small part as test dataset, leaving the rest as training with the data class being 'mids'.
## modeling using multiple-imputation with rubin's rule
## test on five extracted test subsets, calculate mean of AUC and Acc.
nfold = 5
#n_test = floor(nrow(basic_raw)*0.1)

set.seed(123)
imp_df = mice(basic_raw, m=nfold, printFlag = F)
comple_df = complete(imp_df, action = 'all', include = F) ## imputed data stacked df 


## extract testing 
train_list = {}
test_list = {}
for(i in 1:length(comple_df)){
  # if(i==nfold){
  #   idx_test = ((i-1) * n_test + 1): nrow(df_cv)
  # }else{
  #   idx_test = ((i-1) * n_test + 1): (i*n_test)
  # }
  #idx_test = ((i-1) * n_test + 1): (i*n_test)
  #print(c(idx_test[1], idx_test[length(idx_test)]))
  #train_list[[i]] = (comple_df[[i]])[-idx_test,]
  #test_list[[i]] = (comple_df[[i]])[idx_test,]
  train_list[[i]] = (comple_df[[i]])[train_idx,]
  test_list[[i]] = (comple_df[[i]])[-train_idx,]
  
  train_list[[i]]$.imp = i
  train_list[[i]]$.id = 1:nrow( train_list[[i]])
}

## transform to mice: original dataset and imputation dataset
raw_temp = basic_raw[train_idx,]
raw_temp$.imp = 0
raw_temp$.id = 1:nrow(raw_temp)
train_temp = rbind(raw_temp, do.call('rbind',train_list)) ## back to comple_df 
imp_sampled <- as.mids(train_temp,.imp = '.imp',.id = '.id')

## formula
basic_formu = paste0(outcome_name, '~')
for(name in basic_variables){
  basic_formu = paste0(paste0(basic_formu, '+'),name)
}

with_df = with(imp_sampled, expr=glm(formula(basic_formu), family = binomial()))

summ_df = summary(pool(with_df))

## aOR and CI of coefficients
summ_df$aOR<- exp(summ_df$estimate)
summ_df$low_ci <- exp(summ_df$estimate+summ_df$std.error*qnorm(0.025))
summ_df$up_ci <- exp(summ_df$estimate + summ_df$std.error*qnorm(0.975))

write.csv(summ_df, "results_of_mv_before/basic/coef_all_basic_imp_with_test_out.csv")

## test_list
## tempdf = test_list[[1]]
##colnames( (model.matrix(formula(paste0(outcome_name,'~.')), data = tempdf) ))==summ_df$term
thres_type='default'

perf_test_nfold = list()
for(i in 1:length(test_list)){
  #print(i)
  test_fold = (model.matrix(formula(paste0(outcome_name,'~.')), data = test_list[[i]]))
  eta_test = test_fold %*% summ_df$estimate
  prob_test = exp(eta_test) / ( 1 + exp(eta_test) ) ## probability 
  roc_test = roc(test_list[[i]][,outcome_name], prob_test, positive=1, ci=T)
  if(thres_type=='default'){
    prob_thres = 0.5
  }else{
    idx_youden = which.max(roc_test$sensitivities + roc_test$specificities - 1)
    prob_thres = roc_test$thresholds[idx_youden]
  }
  y_pred = as.factor(ifelse(prob_test> prob_thres,1, 0))
  cv_confusion =  confusionMatrix(y_pred,test_list[[i]][,outcome_name],positive = "1")
  
  perf_test_nfold[[i]] = c(roc_test$ci[2], roc_test$ci[1] ,roc_test$ci[3],
                           cv_confusion$overall[1], cv_confusion$overall[3], 
                           cv_confusion$overall[4], prob_thres)
}

perf_test_df = data.frame( do.call('rbind', perf_test_nfold))
colnames(perf_test_df)[c(1:3,7)] = c('AUC', 'AUCLower', 'AUCUpper', 'pred_threshold')

perf_test_df = rbind(perf_test_df, colMeans(perf_test_df))  
rownames(perf_test_df)[nfold+1] = 'mean'    ## the last row is the cv result

write.csv(perf_test_df, "results_of_mv_before/basic/perf_basic_test.csv")


###################################################################
###################################################################
###################################################################
## for each MRC 
sum(is.na(raw_data$MRCday0))
sum(is.na(raw_data$MRCday1))
sum(is.na(raw_data$score_24h))
sum(is.na(raw_data$MRC_ICU24_ge10))
raw_data$MRCday0_ge10 = as.factor(raw_data$MRCday0_ge10)
raw_data$MRCday1_ge10 = as.factor(raw_data$MRCday1_ge10)
raw_data$MRC_ICU24_ge10 = as.factor(raw_data$MRC_ICU24_ge10)

mrc_concerned_variables = c('MRCday0', 'MRCday0_ge10')#,
                            #'MRCday1', 'MRCday1_ge10',
                            #'score_24h', 'MRC_ICU24_ge10')

#dim(tra_std)[1] + dim(test_std)[1]

univ_mrc_coef = data.frame()
perf_mrc_tab = NA
for(name in mrc_concerned_variables){
  formu_expr = ( paste0(paste0(outcome_name,'~'), name) )
  glm_fit = glm(formu_expr, family = binomial(), data = raw_data[train_idx,])
  sum_fit = summary(glm_fit)
  univ_mrc_coef = rbind(univ_mrc_coef, t(data.frame(sum_fit$coefficients[-1,])))
  
  pred_fit = predict(glm_fit, newdata =  raw_data[-train_idx,], type = 'response')
  roc_test = roc(factor((raw_data[-train_idx,]$prolong_mv)), pred_fit, positive=1, ci=T)
  if(thres_type=='default'){
    prob_thres = 0.5
  }else{
    idx_youden = which.max(roc_test$sensitivities + roc_test$specificities - 1)
    prob_thres = roc_test$thresholds[idx_youden]
  }
  
  y_pred = as.factor(ifelse(pred_fit> prob_thres,1, 0))
  cv_confusion =  confusionMatrix(y_pred,factor((raw_data[-train_idx,]$prolong_mv)),positive = "1")
  if( is.na(perf_mrc_tab) ){
    perf_mrc_tab = c(roc_test$ci[2], roc_test$ci[1] ,roc_test$ci[3],
                     cv_confusion$overall[1], cv_confusion$overall[3], 
                     cv_confusion$overall[4], prob_thres)
  }else{
    perf_mrc_tab = rbind(perf_mrc_tab, 
                         c(roc_test$ci[2], roc_test$ci[1] ,roc_test$ci[3],
                           cv_confusion$overall[1], cv_confusion$overall[3], 
                           cv_confusion$overall[4], prob_thres))
  }
}

rownames(univ_mrc_coef) = mrc_concerned_variables
rownames(perf_mrc_tab)= mrc_concerned_variables
colnames(perf_mrc_tab)[c(1:3,7)] = c('AUC', 'AUCLower', 'AUCUpper', 'pred_threshold')

write.csv(univ_mrc_coef, 'results_of_mv_before/mrc_only/mrc_univ_coef.csv')
write.csv(perf_mrc_tab, 'results_of_mv_before/mrc_only/perf_mrc_univ_on_test_once_split.csv')

# summary(glm(prolong_mv~MRCday0_ge10 + MRCday1_ge10 + MRC_ICU24_ge10, family = binomial(),
#             data = raw_data))
# 
# summary(glm(prolong_mv~MRCday0 + MRCday1 + score_24h, family = binomial(),
#             data = raw_data))
##############################################################################################
##############################################################################################
## basic + mrc 
#true_basic_var= basic_variables

mrc_concerned_variables
raw_data$APACHE_24h_score
raw_data$SOFA_24h_score
compared_score = c('APACHE_24h_score', 'SOFA_24h_score')

basic_variables = c(compared_score[2], true_basic_var)
basic_raw = data.frame( raw_data[,c(outcome_name, basic_variables)] )

basic_raw$prolong_mv = as.factor(basic_raw$prolong_mv)
class(basic_raw$Dev_duration)
class( basic_raw$fluid_overload_MVday0 )
basic_raw$fluid_overload_MVday0  = yesno_to_01( basic_raw$fluid_overload_MVday0 )
basic_raw$fluid_overload_MVday1 = yesno_to_01(basic_raw$fluid_overload_MVday1)
basic_raw$HCO3_lt20 = yesno_to_01( basic_raw$HCO3_lt20 )
basic_raw$max_FiO2 
basic_raw$max_PEEP
basic_raw$PaFi
basic_raw$PH_lt725Day1 = yesno_to_01( basic_raw$PH_lt725Day1 )
basic_raw$Pulse_gt110 = yesno_to_01( basic_raw$Pulse_gt110 )
basic_raw$resp  ## 24 hours
basic_raw$min_vent_PaCO2 

str(basic_raw)
colnames(basic_raw)

## normalization
for(i in 1:ncol(basic_raw)){
  if( class(basic_raw[,i]) == 'numeric'){
    basic_raw[,i] = as.numeric(scale(basic_raw[,i]))
  }
}

tra_std = basic_raw[train_idx,]
test_std = basic_raw[-train_idx,]

###########################################################
nfold = 5
#n_test = floor(nrow(basic_raw)*0.1)

set.seed(123)
imp_df = mice(basic_raw, m=nfold, printFlag = F)
comple_df = complete(imp_df, action = 'all', include = F)

train_list = {}
test_list = {}
for(i in 1:length(comple_df)){
  # if(i==nfold){
  #   idx_test = ((i-1) * n_test + 1): nrow(df_cv)
  # }else{
  #   idx_test = ((i-1) * n_test + 1): (i*n_test)
  # }
  #idx_test = ((i-1) * n_test + 1): (i*n_test)
  #print(c(idx_test[1], idx_test[length(idx_test)]))
  #train_list[[i]] = (comple_df[[i]])[-idx_test,]
  #test_list[[i]] = (comple_df[[i]])[idx_test,]
  train_list[[i]] = (comple_df[[i]])[train_idx,]
  test_list[[i]] = (comple_df[[i]])[-train_idx,]
  
  train_list[[i]]$.imp = i
  train_list[[i]]$.id = 1:nrow( train_list[[i]])
}

## transform to mice: original dataset and imputation dataset
raw_temp = basic_raw[train_idx,]
raw_temp$.imp = 0
raw_temp$.id = 1:nrow(raw_temp)
train_temp = rbind(raw_temp, do.call('rbind',train_list))
imp_sampled <- as.mids(train_temp,.imp = '.imp',.id = '.id')


basic_formu = paste0(outcome_name, '~')
for(name in basic_variables){
  basic_formu = paste0(paste0(basic_formu, '+'),name)
}

with_df = with(imp_sampled, expr=glm(formula(basic_formu), family = binomial()))

summ_df = summary(pool(with_df))

## aOR and CI
summ_df$aOR<- exp(summ_df$estimate)
summ_df$low_ci <- exp(summ_df$estimate+summ_df$std.error*qnorm(0.025))
summ_df$up_ci <- exp(summ_df$estimate + summ_df$std.error*qnorm(0.975))

write.csv(summ_df, "results/compared_score/sofa/coef_all_basic_imp_with_test_out.csv")

## test_list
## tempdf = test_list[[1]]
##colnames( (model.matrix(formula(paste0(outcome_name,'~.')), data = tempdf) ))==summ_df$term
thres_type='default'

perf_test_nfold = list()
for(i in 1:length(test_list)){
  #print(i)
  test_fold = (model.matrix(formula(paste0(outcome_name,'~.')), data = test_list[[i]]))
  eta_test = test_fold %*% summ_df$estimate
  prob_test = exp(eta_test) / ( 1 + exp(eta_test) )
  roc_test = roc(test_list[[i]][,outcome_name], as.numeric(prob_test), positive=1, ci=T)
  if(thres_type=='default'){
    prob_thres = 0.5
  }else{
    idx_youden = which.max(roc_test$sensitivities + roc_test$specificities - 1)
    prob_thres = roc_test$thresholds[idx_youden]
  }
  y_pred = as.factor(ifelse(prob_test> prob_thres,1, 0))
  cv_confusion =  confusionMatrix(y_pred,test_list[[i]][,outcome_name],positive = "1")
  
  perf_test_nfold[[i]] = c(roc_test$ci[2], roc_test$ci[1] ,roc_test$ci[3],
                           cv_confusion$overall[1], cv_confusion$overall[3], 
                           cv_confusion$overall[4], prob_thres)
}

perf_test_df = data.frame( do.call('rbind', perf_test_nfold))
colnames(perf_test_df)[c(1:3,7)] = c('AUC', 'AUCLower', 'AUCUpper', 'pred_threshold')

perf_test_df = rbind(perf_test_df, colMeans(perf_test_df))  
rownames(perf_test_df)[nfold+1] = 'mean'    ## the last row is the cv result

write.csv(perf_test_df, "results/compared_score/sofa/perf_basic_test.csv")

#################################################################
#################################################################

