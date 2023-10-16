## variable list which you need to evaluate
merge_basic_x = unique( c(base_variable, day1_variables, icu_24_variables, basic_icu24) )

mv_df = rawdf[-which(is.na(rawdf$dur_MV_days)),c('dur_MV_days', merge_basic_x)] ## complete outcome

mv_df$mv_in5 = as.factor(ifelse(mv_df$dur_MV_days >=5 ,1,0))

mv_df = (mv_df[which(mv_df$dur_MV_days>1),])


#### univariate analysis 
outcome_str =   'prolong_mv' ##  'dur_MV_days' ## 
x_all_name =merge_basic_x
p_value_x_all = c()
#names(p_value_x_all) = x_all_name
for(i in 1:length(x_all_name)){
  #print(i)
  mv_one = mv_df[, c(outcome_str, x_all_name[i])]
  
  if(class(mv_one[,outcome_str]) == 'factor'){ ## this one is for mv_beyond_5days as predictor
    print('factor outcome')
    fit_sum = summary(glm(prolong_mv~., data = mv_one, family = binomial()))
    pvalue = pchisq(fit_sum$null.deviance - fit_sum$deviance, 
                    fit_sum$df.null-fit_sum$df.residual,
                    lower.tail = F)
  }else{
    print('numeric outcome') ## this one is for mv_duration_days as predictor
    table_smm = summary(lm(log(dur_MV_days)~., data= mv_one))
    pvalue =table_smm$coefficients[,4][-1]
  }
  names(pvalue)= x_all_name[i]
  p_value_x_all = c(p_value_x_all,pvalue)
}

p_value_x_all_df = data.frame(names(p_value_x_all), p_value_x_all)
colnames(p_value_x_all_df)[1] = 'x_all_name'
p_value_x_all_df$merge_namesignif = 0
p_value_x_all_df$merge_namesignif[(p_value_x_all_df$p_value_x_all<0.05)] = 1


# numer_y_p_value = p_value_x_all_df
cat_y_p_value = p_value_x_all_df


colnames(mv_df)

## 
cat_basic_imp = mice(mv_df[,c('prolong_mv', x_all_name)], 5)
num_basic_imp = mice(mv_df[,c('dur_MV_days', x_all_name)], 5)

name = x_all_name[1]
for(i in 2:length(x_all_name)){
  tempname = x_all_name[i]
  name = paste0(paste0(name,'+'),tempname)
}

cat_with = with(data= cat_basic_imp, 
                expr = glm(formula( paste0('prolong_mv~', name)), family = binomial()))

cat_imp5 = summary(pool(cat_with))


num_with = with(data= num_basic_imp, 
                expr = lm(formula( paste0('dur_MV_days~', name))))

num_imp5 = summary(pool(num_with))


dim(cat_imp5)
dim(cat_y_p_value)

dim(num_imp5)
dim(numer_y_p_value)

colnames(num_imp5)
num_imp5$p.value

as.character(cat_imp5$term[-1])[1:3]
cat_y_p_value$x_all_name
cat_y_p_value$x_all_name[1:3]

pvalue_merge=data.frame(cbind(cat_y_p_value$x_all_name,cat_y_p_value$p_value_x_all, cat_imp5$p.value[-1],
                              numer_y_p_value$p_value_x_all, num_imp5$p.value[-1]))
colnames(pvalue_merge) = c('predictor', 'univ_p_prolong_mv', 'multi_p_prolong_mv',
                           'univ_p_durmv', 'multiv_p_durmv')

View(pvalue_merge)

write.csv(pvalue_merge, '/Users/Desktop/MV_project/pvalue_of_basic_predictors.csv')



library(gtsummary)
tempdesc <- tbl_summary(mv_df,missing = 'ifany')
tempdesc

summary(mv_df$Dev_duration)

