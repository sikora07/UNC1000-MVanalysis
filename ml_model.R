library(randomForest) ## rf
library(mice)
library(xgboost) ## xgb
library(pROC)
library(caret)
library(mlr)
library(glmnet)
library(e1071) ## svm
library(randomForest)
library(janitor)


## modeling function
## meghod = 'rf' /  'svm'  /  'xgb' 
ml_func<- function(train_raw, test_raw, y_name, variable_name,method = 'rf', thres_type='default'){
  ## random forest 
  data_train = train_raw[, c(y_name, variable_name)]
  data_test = test_raw[, c(y_name, variable_name)]
  
  if(method == 'svm'){
    fit_svm = svm(formula(paste0(y_name,'~.')), data = data_train,type = 'C-classification', kernel = "radial",probability=T)
    pred_result = predict(fit_svm, newdata = data_test, probability = T)
    prob_pre <-attr(pred_result, "probabilities")[,1]
    
  }else if(method=='xgb'){
    set.seed(123)
    fit_xgb = xgboost(data = data.matrix(data_train[,variable_name]), 
            label = as.numeric(as.character(data_train[,y_name])), 
            eta = 0.1, max_depth = 15, 
            nround=25,  subsample = 0.5,
            colsample_bytree = 0.5,
            eval_metric = "error")
    
    prob_pre <- predict(fit_xgb, data.matrix(data_test[,variable_name]) )
    
  }else if(method=='rf'){
    ## binary y
    ## find best mtry
    set.seed(123)
    rf_tuning_mtry <- tuneRF(x=data_train[,variable_name],y=data_train[, y_name],
                             stepFactor = 1.1, trace=TRUE,plot=TRUE,
                             ntreeTry = 1000)
    #rf_tuning_mtry[1]
    
    ## find most used node size 
    set.seed(123)
    rf_tuning_node <- randomForest(formula(paste0(y_name,'~.')),data=data_train,
                                   importance=TRUE,localImp=TRUE,proximity=TRUE,
                                   ntree = 1000,mtry=rf_tuning_mtry[1])
    best_node_size =  as.numeric(names(which.max( table(treesize(rf_tuning_node)) )))
    
    ## build the tree
    set.seed(123)
    rf_nontuning <- randomForest(formula(paste0(y_name,'~.')),data=data_train,
                                 importance=TRUE,localImp=TRUE,proximity=TRUE,
                                 mtry=rf_tuning_mtry[1], nodesize = best_node_size)
    
    prob_pre <-predict(rf_nontuning,data_test, type='prob')[,2]
  }
  

  
  ####################################################################################################
  pred_roc <- roc(response=data_test[,y_name],predictor = prob_pre,ci=T,positive=1)
  
  #pred_roc
  plot(pred_roc,print.thres=F,print.auc=TRUE,col='darkorange')
  
  if(thres_type=='default'){
    pred_y  = factor(ifelse(prob_pre > 0.5 , 1, 0))
    pred_thres = 0.5
  }else{
    max_youden = which.max(pred_roc$sensitivities + pred_roc$specificities -1)
    pred_y  = factor(ifelse(prob_pre >pred_roc$thresholds[max_youden] , 1, 0))
    pred_thres = pred_roc$thresholds[max_youden]
  }
  
  con_matrix <- confusionMatrix(pred_y,data_test[,y_name],positive = "1")
  
  output <- rbind(as.matrix(con_matrix$byClass),
                  as.matrix(con_matrix$overall),
                  AUC=pred_roc[["auc"]], 
                  AUC_lowCI = pred_roc[['ci']][1],
                  AUC_upCI = pred_roc[['ci']][3],
                  mtry = rf_tuning_mtry[,1],
                  nodesize = best_node_size, 
                  OOB = rf_tuning_mtry[,2],
                  pred_thres = pred_thres)
  
  if(method == 'rf'){
    MDSplot(rf_nontuning, test_raw$mv_in5)
    title('train classification')
    MDSplot(rf_nontuning, test_raw$mv_in5)
    title('test classification')
    
    
    ## importance 
    import_df<- data.frame(importance(rf_nontuning), check.names = FALSE)
    import_df <- import_df[order(import_df$MeanDecreaseGini,decreasing = T),]
    #gini_name = rownames(import_df)[order(import_df$MeanDecreaseGini,decreasing = T)]
    #acc_name =  rownames(import_df)[order(import_df$MeanDecreaseAccuracy,decreasing = T)]
    
    rf_model = rf_nontuning
    varImpPlot(rf_nontuning)
    
    
    return(list(rf_model = rf_nontuning, 
                roc_rf = pred_roc, 
                output = output, 
                import_df=import_df))
  }else if(method=='svm'){
    return(list(
      model = fit_svm,
      roc_result = pred_roc,
      output = output)
    )
    }
  else if(method=='xgb'){
    return(list(
      model = fit_xgb,
      roc_result = pred_roc,
      output = output)
    )
  }
}
