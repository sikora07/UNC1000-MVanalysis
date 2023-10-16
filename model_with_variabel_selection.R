## this is based on the organized basic_raw

library(readr)
library(ggplot2)
library(mice)
library(pROC)
library(caret)
library(readxl)
library(glmnet)

yesno_to_01<- function(sel_vector){
  sel_vector[which(sel_vector=='Yes')] = 1
  sel_vector[which(sel_vector == 'No')] = 0
  sel_vector = factor(sel_vector,levels = c(0,1))
  return(sel_vector)
}

## This function can achieve variable selection, estimation, performance evaluation, and automatically save the results to local path.
mv_estimate <- function(input_df, outcome_name, save_path, gt_sel_freq, nfold=5, thres_type='default'){
  ## step2 imputation
  set.seed(123)
  imp_df = mice(input_df, m=nfold, printFlag = F)
  comple_df = complete(imp_df, action = 'all', include = F) ## imputed data stacked df
  
  ## select variables for each fold of imputaions via CV
  formula_str = paste0(outcome_name, ' ~.')
  variable_cv_list = c()
  for(i in 1:length(comple_df)){
    train_current = comple_df[[i]]
    train_cv = model.matrix(formula(formula_str), data = train_current)[,-1]
    
    set.seed(123)
    cv_lasso = cv.glmnet(x=train_cv, y= train_current[,outcome_name], family='binomial')
    glm_lasso = glmnet(x=train_cv, y= train_current[,outcome_name], family='binomial', lambda = cv_lasso$lambda.1se)
    ########################################################################
    ########################################################################
    ########################################################################
     
    sel_incc_dir = unlist(rownames(glm_lasso$beta)[which(glm_lasso$beta!=0)])
    variable_cv_list = c(variable_cv_list, sel_incc_dir)
    
  }
  
  var_seldf = data.frame(table(variable_cv_list))
  write.csv(var_seldf, paste0(save_path,'freq_of_selvar.csv'))
  
  inter_sel_names = as.character(var_seldf$variable_cv_list[which( var_seldf$Freq > gt_sel_freq)]) ## > (nfold/2)
  sel_names = c()
  for(i in 1:length(inter_sel_names)){
    # print(i)
    sel_names = c(sel_names, colnames(train_current)[ grep(substr(inter_sel_names[i],1,nchar(inter_sel_names[i])-1), colnames(train_current)) ])
  }
  
  ## Ending of variable selection, next step is to estimate and predict using CV dataset
  ## extract testing 
  set.seed(123)
  train_idx  = sample(1:nrow(raw_data),nrow(raw_data)*0.9, replace=F)
  
  train_list = {}
  test_list = {}
  for(i in 1:length(comple_df)){
    train_list[[i]] = (comple_df[[i]])[train_idx,]
    test_list[[i]] = (comple_df[[i]])[-train_idx,]
    
    train_list[[i]]$.imp = i
    train_list[[i]]$.id = 1:nrow( train_list[[i]])
  }
  
  ## transform to mice format: original dataset and imputation dataset
  raw_temp = input_df[train_idx,]
  raw_temp$.imp = 0
  raw_temp$.id = 1:nrow(raw_temp)
  train_temp = rbind(raw_temp, do.call('rbind',train_list)) ## back to comple_df 
  imp_sampled <- as.mids(train_temp,.imp = '.imp',.id = '.id')
  
  ## formula
  basic_formu = paste0(outcome_name, '~')
  for(name in sel_names){
    basic_formu = paste0(paste0(basic_formu, '+'),name)
  }
  
  with_df = with(imp_sampled, expr=glm(formula(basic_formu), family = binomial()))
  
  summ_df = summary(pool(with_df))
  
  ## extimates OR: aOR和CI
  summ_df$aOR<- exp(summ_df$estimate)
  summ_df$low_ci <- exp(summ_df$estimate+summ_df$std.error*qnorm(0.025))
  summ_df$up_ci <- exp(summ_df$estimate + summ_df$std.error*qnorm(0.975))
  
  write.csv(summ_df, paste0(save_path,'coef_of_selvar.csv'))
  
  perf_test_nfold = list()
  for(i in 1:length(test_list)){
    #print(i)
    test_fold = (model.matrix(formula(paste0(outcome_name,'~.')), data = test_list[[i]]))
    eta_test = test_fold[,as.character(summ_df$term)] %*% summ_df$estimate
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
  
  write.csv(perf_test_df,  paste0(save_path,'perf_of_selvar.csv'))
  return(list(estimation = summ_df, evaluation = perf_test_df))
}



basic_results = mv_estimate(basic_raw, outcome_name, 
                            save_path='/Users/Desktop/MV_project/oldmvdata_with_vs/',
                            gt_sel_freq=4, thres_type='default')
basic_results$evaluation
basic_results$estimation

