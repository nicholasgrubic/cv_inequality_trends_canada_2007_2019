##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
# CV INDICATOR DESCRIPTIVE TABLES (STRATIFIED BY AGE CATEGORY AND SEX)            
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################

##############################################################################################################################
##############################################################################################################################
# OBESITY  
##############################################################################################################################
##############################################################################################################################

##############################################################################################################################
# STRATIFIED BY AGE CATEGORY
##############################################################################################################################

##############################################################################################################################
# Predicted outcome probabilities (predictive margins)
# Adjusted for age + sex + white + immigrant
##############################################################################################################################
#Define possible values for age_cat
age_cats <- c("6-17 years", "18-39 years", "40-64 years")
strata_combinations <- expand.grid(age_cat = age_cats)
print(strata_combinations)

#Initialize an empty list to store results
results_list <- list()

#Iterate over each combination of strata
for (i in seq_len(nrow(strata_combinations))) {
  age_cat_i <- strata_combinations$age_cat[i]
  
  #Model for calculation of predicted probabilities (age-stratified with adjustment for age, sex, white, and immigration status)
  pred_prob_model <- svyglm(obesity ~ age + sex + white + immigrant,
                            design = svy_design,
                            family = quasibinomial(link="logit"),
                            subset = age_cat == age_cat_i)
  
  #Calculate predicted probabilities
  preds <- as.data.frame(svypredmeans(adjustmodel=pred_prob_model,groupfactor=~cycle)) %>%
    #t-distribution Wald-based logit transformation for CI calculation
    mutate(prob_logit = log(mean/(1 - mean)),
           SE_logit = SE/(mean*(1 - mean)),
           lower_CI_logit = prob_logit - qt(p=0.975, df=68) * SE_logit, 
           upper_CI_logit = prob_logit + qt(p=0.975, df=68) * SE_logit,
           lower_CI = plogis(lower_CI_logit),
           upper_CI = plogis(upper_CI_logit)) %>%
    rownames_to_column(var = "values") %>% mutate(cycle = sub(".*\\.", "", values)) %>% rename(prob = mean) %>%
    mutate(age_cat = age_cat_i) %>%
    select(age_cat, cycle, prob, SE, lower_CI, upper_CI) %>%
    arrange(cycle)
  
  results_list[[i]] <- preds
  rm(preds,pred_prob_model)
  
}

pred_probs_age_obesity <- do.call(rbind, results_list)

remove(results_list)

##############################################################################################################################
# Logistic regression for linear trends
# Adjusted for age + sex + white + immigrant
##############################################################################################################################
#Define possible values for age_cat
age_cats <- c("6-17 years", "18-39 years", "40-64 years")
strata_combinations <- expand.grid(age_cat = age_cats)
print(strata_combinations)  

#Compute p for trend
prev_trend_results_age_obesity <- data.frame(
  age_cat = character(),
  cycle_coef = numeric(),
  p_trend = numeric(),
  stringsAsFactors = FALSE
)

results_list <- list()
for (i in seq_len(nrow(strata_combinations))) {
  #Extract current values
  age_cat_i <- strata_combinations$age_cat[i]
  
  #Run prevalence trend models
  prev_trend_model_age <- svyglm(obesity ~ cycle_cont + age + sex + white + immigrant,
                                 design = svy_design,
                                 family = quasibinomial(link="logit"),
                                 subset = age_cat == age_cat_i)
  
  #Get the p trend
  cycle_coef <- coef(prev_trend_model_age)["cycle_cont"]
  cycle_pval <- summary(prev_trend_model_age)$coefficients["cycle_cont", "Pr(>|t|)"]
  
  #Append the results to the data frame
  prev_trend_results_age_obesity <- rbind(prev_trend_results_age_obesity, data.frame(
    age_cat = age_cat_i,
    cycle_coef = cycle_coef,
    p_trend = cycle_pval,
    stringsAsFactors = FALSE
  ))
}

prev_trend_results_age_obesity <- prev_trend_results_age_obesity %>% arrange(age_cat)

rm(results_list,prev_trend_model_age)

##############################################################################################################################
##############################################################################################################################
# STRATIFIED BY AGE CATEGORY AND SEX
##############################################################################################################################
##############################################################################################################################

##############################################################################################################################
# Predicted outcome probabilities (predictive margins)
# Adjusted for age + white + immigrant
##############################################################################################################################
#Define possible values for age_cat and sex
age_cats <- c("6-17 years", "18-39 years", "40-64 years")
sexes <- c("M", "F")
strata_combinations <- expand.grid(age_cat = age_cats, sex = sexes)
print(strata_combinations)

#Initialize an empty list to store results
results_list <- list()

#Iterate over each combination of strata
for (i in seq_len(nrow(strata_combinations))) {
  age_cat_i <- strata_combinations$age_cat[i]
  sex_i <- strata_combinations$sex[i]
  
  #Model for calculation of predicted probabilities (age- and sex-stratified with adjustment for age, white, and immigration status)
  pred_prob_model <- svyglm(obesity ~ age + white + immigrant,
                            design = svy_design,
                            family = quasibinomial(link="logit"),
                            subset = age_cat == age_cat_i & sex == sex_i)
  
  #Calculate predicted probabilities
  preds <- as.data.frame(svypredmeans(adjustmodel=pred_prob_model,groupfactor=~cycle)) %>%
    #t-distribution Wald-based logit transformation for CI calculation
    mutate(prob_logit = log(mean/(1 - mean)),
           SE_logit = SE/(mean*(1 - mean)),
           lower_CI_logit = prob_logit - qt(p=0.975, df=68) * SE_logit, 
           upper_CI_logit = prob_logit + qt(p=0.975, df=68) * SE_logit,
           lower_CI = plogis(lower_CI_logit),
           upper_CI = plogis(upper_CI_logit)) %>%
    rownames_to_column(var = "values") %>% mutate(cycle = sub(".*\\.", "", values)) %>% rename(prob = mean) %>%
    mutate(age_cat = age_cat_i, sex = sex_i) %>%
    select(age_cat, sex, cycle, prob, SE, lower_CI, upper_CI) %>%
    arrange(sex, cycle)
  
  results_list[[i]] <- preds
  rm(preds,pred_prob_model)
  
}

pred_probs_age_sex_obesity <- do.call(rbind, results_list)

remove(results_list)

##############################################################################################################################
# Logistic regression for linear trends
# Adjusted for age + white + immigrant
##############################################################################################################################
#Define possible values for age_cat and sex
age_cats <- c("6-17 years", "18-39 years", "40-64 years")
sexes <- c("M", "F")
strata_combinations <- expand.grid(age_cat = age_cats, sex = sexes)
print(strata_combinations)  

#Compute p for trend
prev_trend_results_age_sex_obesity <- data.frame(
  age_cat = character(),
  sex = character(),
  cycle_coef = numeric(),
  p_trend = numeric(),
  stringsAsFactors = FALSE
)

results_list <- list()
for (i in seq_len(nrow(strata_combinations))) {
  #Extract current values
  sex_i <- strata_combinations$sex[i]
  age_cat_i <- strata_combinations$age_cat[i]
  
  #Run prevalence trend models
  prev_trend_model_age_sex <- svyglm(obesity ~ cycle_cont + age + white + immigrant,
                                     design = svy_design,
                                     family = quasibinomial(link="logit"),
                                     subset = sex == sex_i & age_cat == age_cat_i)
  
  #Get the p trend
  cycle_coef <- coef(prev_trend_model_age_sex)["cycle_cont"]
  cycle_pval <- summary(prev_trend_model_age_sex)$coefficients["cycle_cont", "Pr(>|t|)"]
  
  #Append the results to the data frame
  prev_trend_results_age_sex_obesity <- rbind(prev_trend_results_age_sex_obesity, data.frame(
    age_cat = age_cat_i,
    sex = sex_i,
    cycle_coef = cycle_coef,
    p_trend = cycle_pval,
    stringsAsFactors = FALSE
  ))
}
prev_trend_results_age_sex_obesity <- prev_trend_results_age_sex_obesity %>% arrange(age_cat, sex)

remove(results_list,prev_trend_model_age_sex)

##############################################################################################################################
##############################################################################################################################
# HBP
##############################################################################################################################
##############################################################################################################################

##############################################################################################################################
# STRATIFIED BY AGE CATEGORY
##############################################################################################################################

##############################################################################################################################
# Predicted outcome probabilities (predictive margins)
# Adjusted for age + sex + white + immigrant
##############################################################################################################################
#Define possible values for age_cat
age_cats <- c("6-17 years", "18-39 years", "40-64 years")
strata_combinations <- expand.grid(age_cat = age_cats)
print(strata_combinations)

#Initialize an empty list to store results
results_list <- list()

#Iterate over each combination of strata
for (i in seq_len(nrow(strata_combinations))) {
  age_cat_i <- strata_combinations$age_cat[i]
  
  #Model for calculation of predicted probabilities (age-stratified with adjustment for age, sex, white, and immigration status)
  pred_prob_model <- svyglm(hbp ~ age + sex + white + immigrant,
                            design = svy_design,
                            family = quasibinomial(link="logit"),
                            subset = age_cat == age_cat_i)
  
  #Calculate predicted probabilities
  preds <- as.data.frame(svypredmeans(adjustmodel=pred_prob_model,groupfactor=~cycle)) %>%
    #t-distribution Wald-based logit transformation for CI calculation
    mutate(prob_logit = log(mean/(1 - mean)),
           SE_logit = SE/(mean*(1 - mean)),
           lower_CI_logit = prob_logit - qt(p=0.975, df=68) * SE_logit, 
           upper_CI_logit = prob_logit + qt(p=0.975, df=68) * SE_logit,
           lower_CI = plogis(lower_CI_logit),
           upper_CI = plogis(upper_CI_logit)) %>%
    rownames_to_column(var = "values") %>% mutate(cycle = sub(".*\\.", "", values)) %>% rename(prob = mean) %>%
    mutate(age_cat = age_cat_i) %>%
    select(age_cat, cycle, prob, SE, lower_CI, upper_CI) %>%
    arrange(cycle)
  
  results_list[[i]] <- preds
  rm(preds,pred_prob_model)
  
}

pred_probs_age_hbp <- do.call(rbind, results_list)

remove(results_list)

##############################################################################################################################
# Logistic regression for linear trends
# Adjusted for age + sex + white + immigrant
##############################################################################################################################
#Define possible values age_cat
age_cats <- c("6-17 years", "18-39 years", "40-64 years")
strata_combinations <- expand.grid(age_cat = age_cats)
print(strata_combinations)  

#Compute p for trend
prev_trend_results_age_hbp <- data.frame(
  age_cat = character(),
  cycle_coef = numeric(),
  p_trend = numeric(),
  stringsAsFactors = FALSE
)

results_list <- list()
for (i in seq_len(nrow(strata_combinations))) {
  #Extract current values
  age_cat_i <- strata_combinations$age_cat[i]
  
  #Run prevalence trend models
  prev_trend_model_age <- svyglm(hbp ~ cycle_cont + age + sex + white + immigrant,
                                 design = svy_design,
                                 family = quasibinomial(link="logit"),
                                 subset = age_cat == age_cat_i)
  
  #Get the p trend
  cycle_coef <- coef(prev_trend_model_age)["cycle_cont"]
  cycle_pval <- summary(prev_trend_model_age)$coefficients["cycle_cont", "Pr(>|t|)"]
  
  #Append the results to the data frame
  prev_trend_results_age_hbp <- rbind(prev_trend_results_age_hbp, data.frame(
    age_cat = age_cat_i,
    cycle_coef = cycle_coef,
    p_trend = cycle_pval,
    stringsAsFactors = FALSE
  ))
}

prev_trend_results_age_hbp <- prev_trend_results_age_hbp %>% arrange(age_cat)

rm(results_list,prev_trend_model_age)

##############################################################################################################################
##############################################################################################################################
# STRATIFIED BY AGE CATEGORY AND SEX
##############################################################################################################################
##############################################################################################################################

##############################################################################################################################
# Predicted outcome probabilities (predictive margins)
# Adjusted for age + white + immigrant
##############################################################################################################################
#Define possible values for age_cat and sex
age_cats <- c("6-17 years", "18-39 years", "40-64 years")
sexes <- c("M", "F")
strata_combinations <- expand.grid(age_cat = age_cats, sex = sexes)
print(strata_combinations)

#Initialize an empty list to store results
results_list <- list()

#Iterate over each combination of strata
for (i in seq_len(nrow(strata_combinations))) {
  age_cat_i <- strata_combinations$age_cat[i]
  sex_i <- strata_combinations$sex[i]
  
  #Model for calculation of predicted probabilities (age- and sex-stratified with adjustment for age, white, and immigration status)
  pred_prob_model <- svyglm(hbp ~ age + white + immigrant,
                            design = svy_design,
                            family = quasibinomial(link="logit"),
                            subset = age_cat == age_cat_i & sex == sex_i)
  
  #Calculate predicted probabilities
  preds <- as.data.frame(svypredmeans(adjustmodel=pred_prob_model,groupfactor=~cycle)) %>%
    #t-distribution Wald-based logit transformation for CI calculation
    mutate(prob_logit = log(mean/(1 - mean)),
           SE_logit = SE/(mean*(1 - mean)),
           lower_CI_logit = prob_logit - qt(p=0.975, df=68) * SE_logit, 
           upper_CI_logit = prob_logit + qt(p=0.975, df=68) * SE_logit,
           lower_CI = plogis(lower_CI_logit),
           upper_CI = plogis(upper_CI_logit)) %>%
    rownames_to_column(var = "values") %>% mutate(cycle = sub(".*\\.", "", values)) %>% rename(prob = mean) %>%
    mutate(age_cat = age_cat_i, sex = sex_i) %>%
    select(age_cat, sex, cycle, prob, SE, lower_CI, upper_CI) %>%
    arrange(sex, cycle)
  
  results_list[[i]] <- preds
  rm(preds,pred_prob_model)
  
}

pred_probs_age_sex_hbp <- do.call(rbind, results_list)

remove(results_list)

##############################################################################################################################
# Logistic regression for linear trends
# Adjusted for age + white + immigrant
##############################################################################################################################
#Define possible values for age_cat and sex
age_cats <- c("6-17 years", "18-39 years", "40-64 years")
sexes <- c("M", "F")
strata_combinations <- expand.grid(age_cat = age_cats, sex = sexes)
print(strata_combinations)  

#Compute p for trend
prev_trend_results_age_sex_hbp <- data.frame(
  age_cat = character(),
  sex = character(),
  cycle_coef = numeric(),
  p_trend = numeric(),
  stringsAsFactors = FALSE
)

results_list <- list()
for (i in seq_len(nrow(strata_combinations))) {
  #Extract current values
  sex_i <- strata_combinations$sex[i]
  age_cat_i <- strata_combinations$age_cat[i]
  
  #Run prevalence trend models
  prev_trend_model_age_sex <- svyglm(hbp ~ cycle_cont + age + white + immigrant,
                                     design = svy_design,
                                     family = quasibinomial(link="logit"),
                                     subset = sex == sex_i & age_cat == age_cat_i)
  
  #Get the p trend
  cycle_coef <- coef(prev_trend_model_age_sex)["cycle_cont"]
  cycle_pval <- summary(prev_trend_model_age_sex)$coefficients["cycle_cont", "Pr(>|t|)"]
  
  #Append the results to the data frame
  prev_trend_results_age_sex_hbp <- rbind(prev_trend_results_age_sex_hbp, data.frame(
    age_cat = age_cat_i,
    sex = sex_i,
    cycle_coef = cycle_coef,
    p_trend = cycle_pval,
    stringsAsFactors = FALSE
  ))
}
prev_trend_results_age_sex_hbp <- prev_trend_results_age_sex_hbp %>% arrange(age_cat, sex)

remove(results_list,prev_trend_model_age_sex)

##############################################################################################################################
##############################################################################################################################
# DYSLIPIDEMIA
##############################################################################################################################
##############################################################################################################################

##############################################################################################################################
# STRATIFIED BY AGE CATEGORY
##############################################################################################################################

##############################################################################################################################
# Predicted outcome probabilities (predictive margins)
# Adjusted for age + sex + white + immigrant
##############################################################################################################################
#Define possible values for age_cat
age_cats <- c("6-17 years", "18-39 years", "40-64 years")
strata_combinations <- expand.grid(age_cat = age_cats)
print(strata_combinations)

#Initialize an empty list to store results
results_list <- list()

#Iterate over each combination of strata
for (i in seq_len(nrow(strata_combinations))) {
  age_cat_i <- strata_combinations$age_cat[i]
  
  #Model for calculation of predicted probabilities (age-stratified with adjustment for age, sex, white, and immigration status)
  pred_prob_model <- svyglm(dyslipidemia_abnormal ~ age + sex + white + immigrant,
                            design = svy_design,
                            family = quasibinomial(link="logit"),
                            subset = age_cat == age_cat_i)
  
  #Calculate predicted probabilities
  preds <- as.data.frame(svypredmeans(adjustmodel=pred_prob_model,groupfactor=~cycle)) %>%
    #t-distribution Wald-based logit transformation for CI calculation
    mutate(prob_logit = log(mean/(1 - mean)),
           SE_logit = SE/(mean*(1 - mean)),
           lower_CI_logit = prob_logit - qt(p=0.975, df=68) * SE_logit, 
           upper_CI_logit = prob_logit + qt(p=0.975, df=68) * SE_logit,
           lower_CI = plogis(lower_CI_logit),
           upper_CI = plogis(upper_CI_logit)) %>%
    rownames_to_column(var = "values") %>% mutate(cycle = sub(".*\\.", "", values)) %>% rename(prob = mean) %>%
    mutate(age_cat = age_cat_i) %>%
    select(age_cat, cycle, prob, SE, lower_CI, upper_CI) %>%
    arrange(cycle)
  
  results_list[[i]] <- preds
  rm(preds,pred_prob_model)
  
}

pred_probs_age_dyslipidemia_abnormal <- do.call(rbind, results_list)

remove(results_list)

##############################################################################################################################
# Logistic regression for linear trends
# Adjusted for age + sex + white + immigrant
##############################################################################################################################
#Define possible values for age_cat
age_cats <- c("6-17 years", "18-39 years", "40-64 years")
strata_combinations <- expand.grid(age_cat = age_cats)
print(strata_combinations)  

#Compute p for trend
prev_trend_results_age_dyslipidemia_abnormal <- data.frame(
  age_cat = character(),
  cycle_coef = numeric(),
  p_trend = numeric(),
  stringsAsFactors = FALSE
)

results_list <- list()
for (i in seq_len(nrow(strata_combinations))) {
  #Extract current values
  age_cat_i <- strata_combinations$age_cat[i]
  
  #Run prevalence trend models
  prev_trend_model_age <- svyglm(dyslipidemia_abnormal ~ cycle_cont + age + sex + white + immigrant,
                                 design = svy_design,
                                 family = quasibinomial(link="logit"),
                                 subset = age_cat == age_cat_i)
  
  #Get the p trend
  cycle_coef <- coef(prev_trend_model_age)["cycle_cont"]
  cycle_pval <- summary(prev_trend_model_age)$coefficients["cycle_cont", "Pr(>|t|)"]
  
  #Append the results to the data frame
  prev_trend_results_age_dyslipidemia_abnormal <- rbind(prev_trend_results_age_dyslipidemia_abnormal, data.frame(
    age_cat = age_cat_i,
    cycle_coef = cycle_coef,
    p_trend = cycle_pval,
    stringsAsFactors = FALSE
  ))
}

prev_trend_results_age_dyslipidemia_abnormal <- prev_trend_results_age_dyslipidemia_abnormal %>% arrange(age_cat)

rm(results_list,prev_trend_model_age)

##############################################################################################################################
##############################################################################################################################
# STRATIFIED BY AGE CATEGORY AND SEX
##############################################################################################################################
##############################################################################################################################

##############################################################################################################################
# Predicted outcome probabilities (predictive margins)
# Adjusted for age + white + immigrant
##############################################################################################################################
#Define possible values for age_cat and sex
age_cats <- c("6-17 years", "18-39 years", "40-64 years")
sexes <- c("M", "F")
strata_combinations <- expand.grid(age_cat = age_cats, sex = sexes)
print(strata_combinations)

#Initialize an empty list to store results
results_list <- list()

#Iterate over each combination of strata
for (i in seq_len(nrow(strata_combinations))) {
  age_cat_i <- strata_combinations$age_cat[i]
  sex_i <- strata_combinations$sex[i]
  
  #Model for calculation of predicted probabilities (age- and sex-stratified with adjustment for age, white, and immigration status)
  pred_prob_model <- svyglm(dyslipidemia_abnormal ~ age + white + immigrant,
                            design = svy_design,
                            family = quasibinomial(link="logit"),
                            subset = age_cat == age_cat_i & sex == sex_i)
  
  #Calculate predicted probabilities
  preds <- as.data.frame(svypredmeans(adjustmodel=pred_prob_model,groupfactor=~cycle)) %>%
    #t-distribution Wald-based logit transformation for CI calculation
    mutate(prob_logit = log(mean/(1 - mean)),
           SE_logit = SE/(mean*(1 - mean)),
           lower_CI_logit = prob_logit - qt(p=0.975, df=68) * SE_logit, 
           upper_CI_logit = prob_logit + qt(p=0.975, df=68) * SE_logit,
           lower_CI = plogis(lower_CI_logit),
           upper_CI = plogis(upper_CI_logit)) %>%
    rownames_to_column(var = "values") %>% mutate(cycle = sub(".*\\.", "", values)) %>% rename(prob = mean) %>%
    mutate(age_cat = age_cat_i, sex = sex_i) %>%
    select(age_cat, sex, cycle, prob, SE, lower_CI, upper_CI) %>%
    arrange(sex, cycle)
  
  results_list[[i]] <- preds
  rm(preds,pred_prob_model)
  
}

pred_probs_age_sex_dyslipidemia_abnormal <- do.call(rbind, results_list)

remove(results_list)

##############################################################################################################################
# Logistic regression for linear trends
# Adjusted for age + white + immigrant
##############################################################################################################################
#Define possible values for age_cat and sex
age_cats <- c("6-17 years", "18-39 years", "40-64 years")
sexes <- c("M", "F")
strata_combinations <- expand.grid(age_cat = age_cats, sex = sexes)
print(strata_combinations)  

#Compute p for trend
prev_trend_results_age_sex_dyslipidemia_abnormal <- data.frame(
  age_cat = character(),
  sex = character(),
  cycle_coef = numeric(),
  p_trend = numeric(),
  stringsAsFactors = FALSE
)

results_list <- list()
for (i in seq_len(nrow(strata_combinations))) {
  #Extract current values
  sex_i <- strata_combinations$sex[i]
  age_cat_i <- strata_combinations$age_cat[i]
  
  #Run prevalence trend models
  prev_trend_model_age_sex <- svyglm(dyslipidemia_abnormal ~ cycle_cont + age + white + immigrant,
                                     design = svy_design,
                                     family = quasibinomial(link="logit"),
                                     subset = sex == sex_i & age_cat == age_cat_i)
  
  #Get the p trend
  cycle_coef <- coef(prev_trend_model_age_sex)["cycle_cont"]
  cycle_pval <- summary(prev_trend_model_age_sex)$coefficients["cycle_cont", "Pr(>|t|)"]
  
  #Append the results to the data frame
  prev_trend_results_age_sex_dyslipidemia_abnormal <- rbind(prev_trend_results_age_sex_dyslipidemia_abnormal, data.frame(
    age_cat = age_cat_i,
    sex = sex_i,
    cycle_coef = cycle_coef,
    p_trend = cycle_pval,
    stringsAsFactors = FALSE
  ))
}
prev_trend_results_age_sex_dyslipidemia_abnormal <- prev_trend_results_age_sex_dyslipidemia_abnormal %>% arrange(age_cat, sex)

remove(results_list,prev_trend_model_age_sex)

##############################################################################################################################
##############################################################################################################################
# DIABETES
##############################################################################################################################
##############################################################################################################################

##############################################################################################################################
# STRATIFIED BY AGE CATEGORY
##############################################################################################################################

##############################################################################################################################
# Predicted outcome probabilities (predictive margins)
# Adjusted for age + sex + white + immigrant
##############################################################################################################################
#Define possible values for age_cat
age_cats <- c("6-17 years", "18-39 years", "40-64 years")
strata_combinations <- expand.grid(age_cat = age_cats)
print(strata_combinations)

#Initialize an empty list to store results
results_list <- list()

#Iterate over each combination of strata
for (i in seq_len(nrow(strata_combinations))) {
  age_cat_i <- strata_combinations$age_cat[i]
  
  #Model for calculation of predicted probabilities (age-stratified with adjustment for age, sex, white, and immigration status)
  pred_prob_model <- svyglm(diabetes ~ age + sex + white + immigrant,
                            design = svy_design,
                            family = quasibinomial(link="logit"),
                            subset = age_cat == age_cat_i)
  
  #Calculate predicted probabilities
  preds <- as.data.frame(svypredmeans(adjustmodel=pred_prob_model,groupfactor=~cycle)) %>%
    #t-distribution Wald-based logit transformation for CI calculation
    mutate(prob_logit = log(mean/(1 - mean)),
           SE_logit = SE/(mean*(1 - mean)),
           lower_CI_logit = prob_logit - qt(p=0.975, df=68) * SE_logit, 
           upper_CI_logit = prob_logit + qt(p=0.975, df=68) * SE_logit,
           lower_CI = plogis(lower_CI_logit),
           upper_CI = plogis(upper_CI_logit)) %>%
    rownames_to_column(var = "values") %>% mutate(cycle = sub(".*\\.", "", values)) %>% rename(prob = mean) %>%
    mutate(age_cat = age_cat_i) %>%
    select(age_cat, cycle, prob, SE, lower_CI, upper_CI) %>%
    arrange(cycle)
  
  results_list[[i]] <- preds
  rm(preds,pred_prob_model)
  
}

pred_probs_age_diabetes <- do.call(rbind, results_list)

remove(results_list)

##############################################################################################################################
# Logistic regression for linear trends
# Adjusted for age + sex + white + immigrant
##############################################################################################################################
#Define possible values for age_cat and sex
age_cats <- c("6-17 years", "18-39 years", "40-64 years")
strata_combinations <- expand.grid(age_cat = age_cats)
print(strata_combinations)  

#Compute p for trend
prev_trend_results_age_diabetes <- data.frame(
  age_cat = character(),
  cycle_coef = numeric(),
  p_trend = numeric(),
  stringsAsFactors = FALSE
)

results_list <- list()
for (i in seq_len(nrow(strata_combinations))) {
  #Extract current values
  age_cat_i <- strata_combinations$age_cat[i]
  
  #Run prevalence trend models
  prev_trend_model_age <- svyglm(diabetes ~ cycle_cont + age + sex + white + immigrant,
                                 design = svy_design,
                                 family = quasibinomial(link="logit"),
                                 subset = age_cat == age_cat_i)
  
  #Get the p trend
  cycle_coef <- coef(prev_trend_model_age)["cycle_cont"]
  cycle_pval <- summary(prev_trend_model_age)$coefficients["cycle_cont", "Pr(>|t|)"]
  
  #Append the results to the data frame
  prev_trend_results_age_diabetes <- rbind(prev_trend_results_age_diabetes, data.frame(
    age_cat = age_cat_i,
    cycle_coef = cycle_coef,
    p_trend = cycle_pval,
    stringsAsFactors = FALSE
  ))
}

prev_trend_results_age_diabetes <- prev_trend_results_age_diabetes %>% arrange(age_cat)

rm(results_list,prev_trend_model_age)

##############################################################################################################################
##############################################################################################################################
# STRATIFIED BY AGE CATEGORY AND SEX
##############################################################################################################################
##############################################################################################################################

##############################################################################################################################
# Predicted outcome probabilities (predictive margins)
# Adjusted for age + white + immigrant
##############################################################################################################################
#Define possible values for age_cat and sex
age_cats <- c("6-17 years", "18-39 years", "40-64 years")
sexes <- c("M", "F")
strata_combinations <- expand.grid(age_cat = age_cats, sex = sexes)
print(strata_combinations)

#Initialize an empty list to store results
results_list <- list()

#Iterate over each combination of strata
for (i in seq_len(nrow(strata_combinations))) {
  age_cat_i <- strata_combinations$age_cat[i]
  sex_i <- strata_combinations$sex[i]
  
  #Model for calculation of predicted probabilities (age- and sex-stratified with adjustment for age, white, and immigration status)
  pred_prob_model <- svyglm(diabetes ~ age + white + immigrant,
                            design = svy_design,
                            family = quasibinomial(link="logit"),
                            subset = age_cat == age_cat_i & sex == sex_i)
  
  #Calculate predicted probabilities
  preds <- as.data.frame(svypredmeans(adjustmodel=pred_prob_model,groupfactor=~cycle)) %>%
    #t-distribution Wald-based logit transformation for CI calculation
    mutate(prob_logit = log(mean/(1 - mean)),
           SE_logit = SE/(mean*(1 - mean)),
           lower_CI_logit = prob_logit - qt(p=0.975, df=68) * SE_logit, 
           upper_CI_logit = prob_logit + qt(p=0.975, df=68) * SE_logit,
           lower_CI = plogis(lower_CI_logit),
           upper_CI = plogis(upper_CI_logit)) %>%
    rownames_to_column(var = "values") %>% mutate(cycle = sub(".*\\.", "", values)) %>% rename(prob = mean) %>%
    mutate(age_cat = age_cat_i, sex = sex_i) %>%
    select(age_cat, sex, cycle, prob, SE, lower_CI, upper_CI) %>%
    arrange(sex, cycle)
  
  results_list[[i]] <- preds
  rm(preds,pred_prob_model)
  
}

pred_probs_age_sex_diabetes <- do.call(rbind, results_list)

remove(results_list)

##############################################################################################################################
# Logistic regression for linear trends
# Adjusted for age + white + immigrant
##############################################################################################################################
#Define possible values for age_cat and sex
age_cats <- c("6-17 years", "18-39 years", "40-64 years")
sexes <- c("M", "F")
strata_combinations <- expand.grid(age_cat = age_cats, sex = sexes)
print(strata_combinations)  

#Compute p for trend
prev_trend_results_age_sex_diabetes <- data.frame(
  age_cat = character(),
  sex = character(),
  cycle_coef = numeric(),
  p_trend = numeric(),
  stringsAsFactors = FALSE
)

results_list <- list()
for (i in seq_len(nrow(strata_combinations))) {
  #Extract current values
  sex_i <- strata_combinations$sex[i]
  age_cat_i <- strata_combinations$age_cat[i]
  
  #Run prevalence trend models
  prev_trend_model_age_sex <- svyglm(diabetes ~ cycle_cont + age + white + immigrant,
                                     design = svy_design,
                                     family = quasibinomial(link="logit"),
                                     subset = sex == sex_i & age_cat == age_cat_i)
  
  #Get the p trend
  cycle_coef <- coef(prev_trend_model_age_sex)["cycle_cont"]
  cycle_pval <- summary(prev_trend_model_age_sex)$coefficients["cycle_cont", "Pr(>|t|)"]
  
  #Append the results to the data frame
  prev_trend_results_age_sex_diabetes <- rbind(prev_trend_results_age_sex_diabetes, data.frame(
    age_cat = age_cat_i,
    sex = sex_i,
    cycle_coef = cycle_coef,
    p_trend = cycle_pval,
    stringsAsFactors = FALSE
  ))
}
prev_trend_results_age_sex_diabetes <- prev_trend_results_age_sex_diabetes %>% arrange(age_cat, sex)

remove(results_list,prev_trend_model_age_sex)