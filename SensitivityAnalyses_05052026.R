##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
# SENSITIVITY ANALYSES               
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################

##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
# MULTIPLE IMPUTATION              
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################

##############################################################################################################################
# MI PREPARATION
##############################################################################################################################
#p value calculation: https://stackoverflow.com/questions/49078764/get-p-values-from-results-of-svyglm-when-using-multiple-imputations-in-r

#Call MI packages
library(mitools)

#Import MI data
obj1_mi_data <- read_sas("J:/Grubic_11125/Data/CHMS Linked Source Files/obj1_cohort_mi.sas7bdat")

#Factorize character variables
vars_to_factor <- c("cycle","immigrant","age_cat","sex","race_ethnic","white","ei_quartile_moecd") #character variables

obj1_mi_data[vars_to_factor] <- lapply(obj1_mi_data[vars_to_factor], function(x) factor(x, exclude=NA))
obj1_mi_data$cycle_cont <- as.numeric(obj1_mi_data$cycle)
str(obj1_mi_data)

#Split dataset into list of imputations
imps_list <- split(obj1_mi_data, obj1_mi_data$`_Imputation_`)

#Define the survey design using replicate weights
mi_svy_designs <- lapply(imps_list, function(d) {
  svrepdesign(
    degf = 68, #Recommended by Statistics Canada for national level analyses of cycles 1-6 
    #(file:///S:/CHMS%20Cycle%206/CHMS_ECMS_C6_W6_v1/documentation/english/user_guide/CHMS_c1_c6_w2_comb_f1_T15.1_v1_EN.pdf)
    weights = ~WGT_FULL, 
    repweights = d[, grep("^BSW", names(d))],
    type = "bootstrap",
    data = d
)
})

##############################################################################################################################
# Descriptive characters of respondents with and without complete data
##############################################################################################################################

#Output complete cohort from imputation data
obj1_missing_data <- subset(obj1_mi_data,obj1_mi_data$`_Imputation_` == 1)

#Define the survey design using replicate weights
svy_design_missing_data <- svrepdesign(data = obj1_missing_data, 
                                degf = 68, #Recommended by Statistics Canada for national level analyses of cycles 1-6 
                                #(file:///S:/CHMS%20Cycle%206/CHMS_ECMS_C6_W6_v1/documentation/english/user_guide/CHMS_c1_c6_w2_comb_f1_T15.1_v1_EN.pdf)
                                 weights = ~WGT_FULL, 
                                 repweights = obj1_missing_data[, grep("BSW", names(obj1_missing_data))], 
                                 type = "bootstrap") 

#Generate descriptive characteristics
descriptive_chars_missing <- svyby(~ age + age_cat + sex + cycle + ei_quartile_moecd + eq_moecd_hhi,
                                by = ~missing_data,
                                design = svy_design_missing_data,
                                vartype = c("ci", "cvpct"),
                                FUN = svymean,
                                NA.rm = TRUE) #exclude missing

#Sociodemographic characteristics
svy_stats_age <- descriptive_chars_missing %>% 
  mutate(age = paste0(round(`age`,1)," (",round(`ci_l.age`,1),"-",round(`ci_u.age`,1),"), CV%= ",round(`cv%.age`,1))) %>%
  select(missing_data,age) %>%
  pivot_longer(cols = -missing_data, names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = missing_data, values_from = value)

svy_stats_age_cat <- descriptive_chars_missing %>% 
  mutate(age_cat_6_17_years = paste0(round(`age_cat6-17 years`*100,1)," (",round(`ci_l.age_cat6-17 years`*100,1),"-",round(`ci_u.age_cat6-17 years`*100,1),"), CV%= ",round(`cv%.age_cat6-17 years`,2))) %>%
  mutate(age_cat_18_39_years = paste0(round(`age_cat18-39 years`*100,1)," (",round(`ci_l.age_cat18-39 years`*100,1),"-",round(`ci_u.age_cat18-39 years`*100,1),"), CV%= ",round(`cv%.age_cat18-39 years`,2))) %>%
  mutate(age_cat_40_64_years = paste0(round(`age_cat40-64 years`*100,1)," (",round(`ci_l.age_cat40-64 years`*100,1),"-",round(`ci_u.age_cat40-64 years`*100,1),"), CV%= ",round(`cv%.age_cat40-64 years`,2))) %>%
  select(missing_data,age_cat_6_17_years,age_cat_18_39_years,age_cat_40_64_years) %>%
  pivot_longer(cols = -missing_data, names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = missing_data, values_from = value)

svy_stats_sex <- descriptive_chars_missing %>% 
  mutate(sexF = paste0(round(`sexF`*100,1)," (",round(`ci_l.sexF`*100,1),"-",round(`ci_u.sexF`*100,1),"), CV%= ",round(`cv%.sexF`,2))) %>%
  mutate(sexM = paste0(round(`sexM`*100,1)," (",round(`ci_l.sexM`*100,1),"-",round(`ci_u.sexM`*100,1),"), CV%= ",round(`cv%.sexM`,2))) %>%
  select(missing_data,sexF,sexM) %>%
  pivot_longer(cols = -missing_data, names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = missing_data, values_from = value)

svy_stats_cycle <- descriptive_chars_missing %>% 
  mutate(cycle1 = paste0(round(`cycle1`*100,1)," (",round(`ci_l.cycle1`*100,1),"-",round(`ci_u.cycle1`*100,1),"), CV%= ",round(`cv%.cycle1`,2))) %>%
  mutate(cycle2 = paste0(round(`cycle2`*100,1)," (",round(`ci_l.cycle2`*100,1),"-",round(`ci_u.cycle2`*100,1),"), CV%= ",round(`cv%.cycle2`,2))) %>%
  mutate(cycle3 = paste0(round(`cycle3`*100,1)," (",round(`ci_l.cycle3`*100,1),"-",round(`ci_u.cycle3`*100,1),"), CV%= ",round(`cv%.cycle3`,2))) %>%
  mutate(cycle4 = paste0(round(`cycle4`*100,1)," (",round(`ci_l.cycle4`*100,1),"-",round(`ci_u.cycle4`*100,1),"), CV%= ",round(`cv%.cycle4`,2))) %>%
  mutate(cycle5 = paste0(round(`cycle5`*100,1)," (",round(`ci_l.cycle5`*100,1),"-",round(`ci_u.cycle5`*100,1),"), CV%= ",round(`cv%.cycle5`,2))) %>%
  mutate(cycle6 = paste0(round(`cycle6`*100,1)," (",round(`ci_l.cycle6`*100,1),"-",round(`ci_u.cycle6`*100,1),"), CV%= ",round(`cv%.cycle6`,2))) %>%
  select(missing_data,cycle1,cycle2,cycle3,cycle4,cycle5,cycle6) %>%
  pivot_longer(cols = -missing_data, names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = missing_data, values_from = value)

svy_stats_ei_quartile_moecd <- descriptive_chars_missing %>% 
  mutate(ei_quartile_moecd1 = paste0(round(`ei_quartile_moecd1`*100,1)," (",round(`ci_l.ei_quartile_moecd1`*100,1),"-",round(`ci_u.ei_quartile_moecd1`*100,1),"), CV%= ",round(`cv%.ei_quartile_moecd1`,2))) %>%
  mutate(ei_quartile_moecd2 = paste0(round(`ei_quartile_moecd2`*100,1)," (",round(`ci_l.ei_quartile_moecd2`*100,1),"-",round(`ci_u.ei_quartile_moecd2`*100,1),"), CV%= ",round(`cv%.ei_quartile_moecd2`,2))) %>%
  mutate(ei_quartile_moecd3 = paste0(round(`ei_quartile_moecd3`*100,1)," (",round(`ci_l.ei_quartile_moecd3`*100,1),"-",round(`ci_u.ei_quartile_moecd3`*100,1),"), CV%= ",round(`cv%.ei_quartile_moecd3`,2))) %>%
  mutate(ei_quartile_moecd4 = paste0(round(`ei_quartile_moecd4`*100,1)," (",round(`ci_l.ei_quartile_moecd4`*100,1),"-",round(`ci_u.ei_quartile_moecd4`*100,1),"), CV%= ",round(`cv%.ei_quartile_moecd4`,2))) %>%
  select(missing_data,ei_quartile_moecd1,ei_quartile_moecd2,ei_quartile_moecd3,ei_quartile_moecd4) %>%
  pivot_longer(cols = -missing_data, names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = missing_data, values_from = value)

svy_stats_eq_moecd_hhi <- descriptive_chars_missing %>% 
  mutate(eq_moecd_hhi = paste0(round(`eq_moecd_hhi`,1)," (",round(`ci_l.eq_moecd_hhi`,1),"-",round(`ci_u.eq_moecd_hhi`,1),"), CV%= ",round(`cv%.eq_moecd_hhi`,1))) %>%
  select(missing_data,eq_moecd_hhi) %>%
  pivot_longer(cols = -missing_data, names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = missing_data, values_from = value)

descriptive_chars_missing_table <- rbind(svy_stats_age,svy_stats_age_cat,svy_stats_sex,svy_stats_cycle,
                                         svy_stats_ei_quartile_moecd,svy_stats_eq_moecd_hhi)

remove(svy_stats_age,svy_stats_age_cat,svy_stats_sex,svy_stats_cycle,
       svy_stats_ei_quartile_moecd,svy_stats_eq_moecd_hhi)

##############################################################################################################################
# Health inequalities - stratified by age category (with interactions to parse out cycle-specific effects) - w/o svycontrast
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
  
  #Run RII and SII models on MI data and combine estimates and SEs
  RII_model <- function(design) {
    svyglm(r_poor_cv_health_1 ~ ridit_eqi_moecd*cycle + age + sex + white + immigrant,
           design = design,
           family = quasipoisson(link="log"),
           subset = age_cat == age_cat_i)
  }
  RII_model <- lapply(mi_svy_designs, RII_model)
  RII_model_pooled <- MIcombine(MIextract(RII_model, fun = coef), MIextract(RII_model, fun = vcov))
  summary(RII_model_pooled)
  
  RII_trend_model <- function(design) {
    svyglm(r_poor_cv_health_1 ~ ridit_eqi_moecd*cycle_cont + age + sex + white + immigrant,
           design = design,
           family = quasipoisson(link="log"),
           subset = age_cat == age_cat_i)
  }
  RII_trend_model <- lapply(mi_svy_designs, RII_trend_model)
  RII_trend_model_pooled <- MIcombine(MIextract(RII_trend_model, fun = coef), MIextract(RII_trend_model, fun = vcov))
  summary(RII_trend_model_pooled)
  
  SII_model <- function(design) {
    svyglm(r_poor_cv_health_1 ~ ridit_eqi_moecd*cycle + age + sex + white + immigrant,
           design = design,
           family = gaussian(link="identity"),
           subset = age_cat == age_cat_i)
  }
  SII_model <- lapply(mi_svy_designs, SII_model)
  SII_model_pooled <- MIcombine(MIextract(SII_model, fun = coef), MIextract(SII_model, fun = vcov))
  summary(SII_model_pooled)
  
  SII_trend_model <- function(design) {
    svyglm(r_poor_cv_health_1 ~ ridit_eqi_moecd*cycle_cont + age + sex + white + immigrant,
           design = design,
           family = gaussian(link="identity"),
           subset = age_cat == age_cat_i)
  }
  SII_trend_model <- lapply(mi_svy_designs, SII_trend_model)
  SII_trend_model_pooled <- MIcombine(MIextract(SII_trend_model, fun = coef), MIextract(SII_trend_model, fun = vcov))
  summary(SII_trend_model_pooled)
  
  #Define contrast vectors
  r_c1 <- c("ridit_eqi_moecd" = 1)
  r_c2 <- c("ridit_eqi_moecd" = 1, "ridit_eqi_moecd:cycle2" = 1)
  r_c3 <- c("ridit_eqi_moecd" = 1, "ridit_eqi_moecd:cycle3" = 1)
  r_c4 <- c("ridit_eqi_moecd" = 1, "ridit_eqi_moecd:cycle4" = 1)
  r_c5 <- c("ridit_eqi_moecd" = 1, "ridit_eqi_moecd:cycle5" = 1)
  r_c6 <- c("ridit_eqi_moecd" = 1, "ridit_eqi_moecd:cycle6" = 1)
  
  #Fill in zeros for other coefficients and compute contrast estimates, SE, RII/SII, and CIs
  #RII
  r_c1_full <- rep(0, length(coef(RII_model_pooled)))
  names(r_c1_full) <- names(coef(RII_model_pooled))
  r_c1_full[names(r_c1)] <- r_c1
  est <- sum(r_c1_full * coef(RII_model_pooled))
  se <- sqrt(t(r_c1_full) %*% vcov(RII_model_pooled) %*% r_c1_full)
  RII_parms_r_c1 <- data.frame(Est = est, SE = se, RII = exp(est), LCL = exp(est - 1.96 * se), UCL = exp(est + 1.96 * se))
  
  r_c2_full <- rep(0, length(coef(RII_model_pooled)))
  names(r_c2_full) <- names(coef(RII_model_pooled))
  r_c2_full[names(r_c2)] <- r_c2
  est <- sum(r_c2_full * coef(RII_model_pooled))
  se <- sqrt(t(r_c2_full) %*% vcov(RII_model_pooled) %*% r_c2_full)
  RII_parms_r_c2 <- data.frame(Est = est, SE = se, RII = exp(est), LCL = exp(est - 1.96 * se), UCL = exp(est + 1.96 * se))
  
  r_c3_full <- rep(0, length(coef(RII_model_pooled)))
  names(r_c3_full) <- names(coef(RII_model_pooled))
  r_c3_full[names(r_c3)] <- r_c3
  est <- sum(r_c3_full * coef(RII_model_pooled))
  se <- sqrt(t(r_c3_full) %*% vcov(RII_model_pooled) %*% r_c3_full)
  RII_parms_r_c3 <- data.frame(Est = est, SE = se, RII = exp(est), LCL = exp(est - 1.96 * se), UCL = exp(est + 1.96 * se))
  
  r_c4_full <- rep(0, length(coef(RII_model_pooled)))
  names(r_c4_full) <- names(coef(RII_model_pooled))
  r_c4_full[names(r_c4)] <- r_c4
  est <- sum(r_c4_full * coef(RII_model_pooled))
  se <- sqrt(t(r_c4_full) %*% vcov(RII_model_pooled) %*% r_c4_full)
  RII_parms_r_c4 <- data.frame(Est = est, SE = se, RII = exp(est), LCL = exp(est - 1.96 * se), UCL = exp(est + 1.96 * se))
  
  r_c5_full <- rep(0, length(coef(RII_model_pooled)))
  names(r_c5_full) <- names(coef(RII_model_pooled))
  r_c5_full[names(r_c5)] <- r_c5
  est <- sum(r_c5_full * coef(RII_model_pooled))
  se <- sqrt(t(r_c5_full) %*% vcov(RII_model_pooled) %*% r_c5_full)
  RII_parms_r_c5 <- data.frame(Est = est, SE = se, RII = exp(est), LCL = exp(est - 1.96 * se), UCL = exp(est + 1.96 * se))
  
  r_c6_full <- rep(0, length(coef(RII_model_pooled)))
  names(r_c6_full) <- names(coef(RII_model_pooled))
  r_c6_full[names(r_c6)] <- r_c6
  est <- sum(r_c6_full * coef(RII_model_pooled))
  se <- sqrt(t(r_c6_full) %*% vcov(RII_model_pooled) %*% r_c6_full)
  RII_parms_r_c6 <- data.frame(Est = est, SE = se, RII = exp(est), LCL = exp(est - 1.96 * se), UCL = exp(est + 1.96 * se))
  
  RII_parms <- rbind(RII_parms_r_c1,RII_parms_r_c2,RII_parms_r_c3,RII_parms_r_c4,RII_parms_r_c5,RII_parms_r_c6)
  
  #SII
  est <- sum(r_c1_full * coef(SII_model_pooled))
  se <- sqrt(t(r_c1_full) %*% vcov(SII_model_pooled) %*% r_c1_full)
  SII_parms_r_c1 <- data.frame(Est = est, SE = se, SII = est, LCL = est - 1.96 * se, UCL = est + 1.96 * se)
  
  est <- sum(r_c2_full * coef(SII_model_pooled))
  se <- sqrt(t(r_c2_full) %*% vcov(SII_model_pooled) %*% r_c2_full)
  SII_parms_r_c2 <- data.frame(Est = est, SE = se, SII = est, LCL = est - 1.96 * se, UCL = est + 1.96 * se)
  
  est <- sum(r_c3_full * coef(SII_model_pooled))
  se <- sqrt(t(r_c3_full) %*% vcov(SII_model_pooled) %*% r_c3_full)
  SII_parms_r_c3 <- data.frame(Est = est, SE = se, SII = est, LCL = est - 1.96 * se, UCL = est + 1.96 * se)
  
  est <- sum(r_c4_full * coef(SII_model_pooled))
  se <- sqrt(t(r_c4_full) %*% vcov(SII_model_pooled) %*% r_c4_full)
  SII_parms_r_c4 <- data.frame(Est = est, SE = se, SII = est, LCL = est - 1.96 * se, UCL = est + 1.96 * se)
  
  est <- sum(r_c5_full * coef(SII_model_pooled))
  se <- sqrt(t(r_c5_full) %*% vcov(SII_model_pooled) %*% r_c5_full)
  SII_parms_r_c5 <- data.frame(Est = est, SE = se, SII = est, LCL = est - 1.96 * se, UCL = est + 1.96 * se)
  
  est <- sum(r_c6_full * coef(SII_model_pooled))
  se <- sqrt(t(r_c6_full) %*% vcov(SII_model_pooled) %*% r_c6_full)
  SII_parms_r_c6 <- data.frame(Est = est, SE = se, SII = est, LCL = est - 1.96 * se, UCL = est + 1.96 * se)
  
  SII_parms <- rbind(SII_parms_r_c1,SII_parms_r_c2,SII_parms_r_c3,SII_parms_r_c4,SII_parms_r_c5,SII_parms_r_c6)
  
  #Get RII and SII parameters
  RII_coefs <- as.numeric(RII_parms$Est)  
  RII_SEs <- as.numeric(RII_parms$SE)             
  RII_values <- as.numeric(RII_parms$RII)   
  RII_CIs <- RII_parms[4:5]   
  RII_trend_coef <- as.numeric(coef(RII_trend_model_pooled)["ridit_eqi_moecd:cycle_cont"])
  RII_trend_SE <- as.numeric(SE(RII_trend_model_pooled)["ridit_eqi_moecd:cycle_cont"])
  RII_trend_pvalue <- pt(abs(RII_trend_coef/RII_trend_SE),df=RII_trend_model[[1]]$df.residual,lower.tail=FALSE)*2
  
  SII_coefs <- as.numeric(SII_parms$Est)     
  SII_SEs <- as.numeric(SII_parms$SE) 
  SII_values <- as.numeric(SII_parms$SII)   
  SII_CIs <- SII_parms[4:5] 
  SII_trend_coef <- as.numeric(coef(SII_trend_model_pooled)["ridit_eqi_moecd:cycle_cont"])
  SII_trend_SE <- as.numeric(SE(SII_trend_model_pooled)["ridit_eqi_moecd:cycle_cont"])
  SII_trend_pvalue <- pt(abs(SII_trend_coef/SII_trend_SE),df=SII_trend_model[[1]]$df.residual,lower.tail=FALSE)*2
  
  RII <- data.frame(
    Age_cat = age_cat_i,
    Metric = 'RII',
    Cycle = c(1:6),
    Coef = RII_coefs,
    SE = RII_SEs,
    Value = RII_values,
    Lower_95_CI = as.numeric(RII_CIs$LCL),
    Upper_95_CI = as.numeric(RII_CIs$UCL),
    Trend_coef = RII_trend_coef,
    Trend_SE = RII_trend_SE,
    Trend_pvalue = RII_trend_pvalue 
  )
  
  SII <- data.frame(
    Age_cat = age_cat_i,
    Metric = 'SII',
    Cycle = c(1:6),
    Coef = SII_coefs,
    SE = SII_SEs,
    Value = SII_values,
    Lower_95_CI = as.numeric(SII_CIs$LCL),
    Upper_95_CI = as.numeric(SII_CIs$UCL),
    Trend_coef = SII_trend_coef,
    Trend_SE = SII_trend_SE,
    Trend_pvalue = SII_trend_pvalue 
  )
  
  results_list[[i]] <- rbind(RII,SII)
  rm(RII_model,RII_trend_model,SII_model,SII_trend_model,
     RII_model_pooled,RII_trend_model_pooled,SII_model_pooled,SII_trend_model_pooled,
     r_c1,r_c2,r_c3,r_c4,r_c5,r_c6,r_c1_full,r_c2_full,r_c3_full,r_c4_full,r_c5_full,r_c6_full,est,se,
     RII_parms_r_c1,RII_parms_r_c2,RII_parms_r_c3,RII_parms_r_c4,RII_parms_r_c5,RII_parms_r_c6,RII_parms,
     SII_parms_r_c1,SII_parms_r_c2,SII_parms_r_c3,SII_parms_r_c4,SII_parms_r_c5,SII_parms_r_c6,SII_parms,
     RII_coefs,RII_SEs,RII_values,RII_CIs,RII_trend_coef,RII_trend_SE,RII_trend_pvalue,RII,
     SII_coefs,SII_SEs,SII_values,SII_CIs,SII_trend_coef,SII_trend_SE,SII_trend_pvalue,SII)
}

health_inequalities_age_mi<- do.call(rbind, results_list)
rownames(health_inequalities_age_mi) <- NULL

##############################################################################################################################
# Health inequalities - stratified by age category and sex (with interactions to parse out cycle-specific effects) - w/o svycontrast
##############################################################################################################################
#Define possible values for age_cat
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
  
  #Run RII and SII models on MI data and combine estimates and SEs
  RII_model <- function(design) {
    svyglm(r_poor_cv_health_1 ~ ridit_eqi_moecd*cycle + age + white + immigrant,
           design = design,
           family = quasipoisson(link="log"),
           subset = age_cat == age_cat_i & sex == sex_i)
  }
  RII_model <- lapply(mi_svy_designs, RII_model)
  RII_model_pooled <- MIcombine(MIextract(RII_model, fun = coef), MIextract(RII_model, fun = vcov))
  summary(RII_model_pooled)
  
  RII_trend_model <- function(design) {  
    svyglm(r_poor_cv_health_1 ~ ridit_eqi_moecd*cycle_cont + age + white + immigrant,
           design = design,
           family = quasipoisson(link="log"),
           subset = age_cat == age_cat_i & sex == sex_i)
  }
  RII_trend_model <- lapply(mi_svy_designs, RII_trend_model)
  RII_trend_model_pooled <- MIcombine(MIextract(RII_trend_model, fun = coef), MIextract(RII_trend_model, fun = vcov))
  summary(RII_trend_model_pooled) 
  
  RII_trend_sex_model <- function(design) {
    svyglm(r_poor_cv_health_1 ~ ridit_eqi_moecd*cycle_cont*sex + age + white + immigrant,
           design = design,
           family = quasipoisson(link="log"),
           subset = age_cat == age_cat_i)
  }
  RII_trend_sex_model <- lapply(mi_svy_designs, RII_trend_sex_model)
  RII_trend_sex_model_pooled <- MIcombine(MIextract(RII_trend_sex_model, fun = coef), MIextract(RII_trend_sex_model, fun = vcov))
  summary(RII_trend_sex_model_pooled) 
  
  SII_model <- function(design) {
    svyglm(r_poor_cv_health_1 ~ ridit_eqi_moecd*cycle + age + white + immigrant,
           design = design,
           family = gaussian(link="identity"),
           subset = age_cat == age_cat_i & sex == sex_i)
  }
  SII_model <- lapply(mi_svy_designs, SII_model)
  SII_model_pooled <- MIcombine(MIextract(SII_model, fun = coef), MIextract(SII_model, fun = vcov))
  summary(SII_model_pooled)
  
  SII_trend_model <- function(design) {  
    svyglm(r_poor_cv_health_1 ~ ridit_eqi_moecd*cycle_cont + age + white + immigrant,
           design = design,
           family = gaussian(link="identity"),
           subset = age_cat == age_cat_i & sex == sex_i)
  }
  SII_trend_model <- lapply(mi_svy_designs, SII_trend_model)
  SII_trend_model_pooled <- MIcombine(MIextract(SII_trend_model, fun = coef), MIextract(SII_trend_model, fun = vcov))
  summary(SII_trend_model_pooled) 
  
  SII_trend_sex_model <- function(design) {
    svyglm(r_poor_cv_health_1 ~ ridit_eqi_moecd*cycle_cont*sex + age + white + immigrant,
           design = design,
           family = gaussian(link="identity"),
           subset = age_cat == age_cat_i)
  }
  SII_trend_sex_model <- lapply(mi_svy_designs, SII_trend_sex_model)
  SII_trend_sex_model_pooled <- MIcombine(MIextract(SII_trend_sex_model, fun = coef), MIextract(SII_trend_sex_model, fun = vcov))
  summary(SII_trend_sex_model_pooled) 
  
  #Define contrast vectors
  r_c1 <- c("ridit_eqi_moecd" = 1)
  r_c2 <- c("ridit_eqi_moecd" = 1, "ridit_eqi_moecd:cycle2" = 1)
  r_c3 <- c("ridit_eqi_moecd" = 1, "ridit_eqi_moecd:cycle3" = 1)
  r_c4 <- c("ridit_eqi_moecd" = 1, "ridit_eqi_moecd:cycle4" = 1)
  r_c5 <- c("ridit_eqi_moecd" = 1, "ridit_eqi_moecd:cycle5" = 1)
  r_c6 <- c("ridit_eqi_moecd" = 1, "ridit_eqi_moecd:cycle6" = 1)
  
  #Fill in zeros for other coefficients and compute contrast estimates, SE, RII/SII, and CIs
  #RII
  r_c1_full <- rep(0, length(coef(RII_model_pooled)))
  names(r_c1_full) <- names(coef(RII_model_pooled))
  r_c1_full[names(r_c1)] <- r_c1
  est <- sum(r_c1_full * coef(RII_model_pooled))
  se <- sqrt(t(r_c1_full) %*% vcov(RII_model_pooled) %*% r_c1_full)
  RII_parms_r_c1 <- data.frame(Est = est, SE = se, RII = exp(est), LCL = exp(est - 1.96 * se), UCL = exp(est + 1.96 * se))
  
  r_c2_full <- rep(0, length(coef(RII_model_pooled)))
  names(r_c2_full) <- names(coef(RII_model_pooled))
  r_c2_full[names(r_c2)] <- r_c2
  est <- sum(r_c2_full * coef(RII_model_pooled))
  se <- sqrt(t(r_c2_full) %*% vcov(RII_model_pooled) %*% r_c2_full)
  RII_parms_r_c2 <- data.frame(Est = est, SE = se, RII = exp(est), LCL = exp(est - 1.96 * se), UCL = exp(est + 1.96 * se))
  
  r_c3_full <- rep(0, length(coef(RII_model_pooled)))
  names(r_c3_full) <- names(coef(RII_model_pooled))
  r_c3_full[names(r_c3)] <- r_c3
  est <- sum(r_c3_full * coef(RII_model_pooled))
  se <- sqrt(t(r_c3_full) %*% vcov(RII_model_pooled) %*% r_c3_full)
  RII_parms_r_c3 <- data.frame(Est = est, SE = se, RII = exp(est), LCL = exp(est - 1.96 * se), UCL = exp(est + 1.96 * se))
  
  r_c4_full <- rep(0, length(coef(RII_model_pooled)))
  names(r_c4_full) <- names(coef(RII_model_pooled))
  r_c4_full[names(r_c4)] <- r_c4
  est <- sum(r_c4_full * coef(RII_model_pooled))
  se <- sqrt(t(r_c4_full) %*% vcov(RII_model_pooled) %*% r_c4_full)
  RII_parms_r_c4 <- data.frame(Est = est, SE = se, RII = exp(est), LCL = exp(est - 1.96 * se), UCL = exp(est + 1.96 * se))
  
  r_c5_full <- rep(0, length(coef(RII_model_pooled)))
  names(r_c5_full) <- names(coef(RII_model_pooled))
  r_c5_full[names(r_c5)] <- r_c5
  est <- sum(r_c5_full * coef(RII_model_pooled))
  se <- sqrt(t(r_c5_full) %*% vcov(RII_model_pooled) %*% r_c5_full)
  RII_parms_r_c5 <- data.frame(Est = est, SE = se, RII = exp(est), LCL = exp(est - 1.96 * se), UCL = exp(est + 1.96 * se))
  
  r_c6_full <- rep(0, length(coef(RII_model_pooled)))
  names(r_c6_full) <- names(coef(RII_model_pooled))
  r_c6_full[names(r_c6)] <- r_c6
  est <- sum(r_c6_full * coef(RII_model_pooled))
  se <- sqrt(t(r_c6_full) %*% vcov(RII_model_pooled) %*% r_c6_full)
  RII_parms_r_c6 <- data.frame(Est = est, SE = se, RII = exp(est), LCL = exp(est - 1.96 * se), UCL = exp(est + 1.96 * se))
  
  RII_parms <- rbind(RII_parms_r_c1,RII_parms_r_c2,RII_parms_r_c3,RII_parms_r_c4,RII_parms_r_c5,RII_parms_r_c6)
  
  #SII
  est <- sum(r_c1_full * coef(SII_model_pooled))
  se <- sqrt(t(r_c1_full) %*% vcov(SII_model_pooled) %*% r_c1_full)
  SII_parms_r_c1 <- data.frame(Est = est, SE = se, SII = est, LCL = est - 1.96 * se, UCL = est + 1.96 * se)
  
  est <- sum(r_c2_full * coef(SII_model_pooled))
  se <- sqrt(t(r_c2_full) %*% vcov(SII_model_pooled) %*% r_c2_full)
  SII_parms_r_c2 <- data.frame(Est = est, SE = se, SII = est, LCL = est - 1.96 * se, UCL = est + 1.96 * se)
  
  est <- sum(r_c3_full * coef(SII_model_pooled))
  se <- sqrt(t(r_c3_full) %*% vcov(SII_model_pooled) %*% r_c3_full)
  SII_parms_r_c3 <- data.frame(Est = est, SE = se, SII = est, LCL = est - 1.96 * se, UCL = est + 1.96 * se)
  
  est <- sum(r_c4_full * coef(SII_model_pooled))
  se <- sqrt(t(r_c4_full) %*% vcov(SII_model_pooled) %*% r_c4_full)
  SII_parms_r_c4 <- data.frame(Est = est, SE = se, SII = est, LCL = est - 1.96 * se, UCL = est + 1.96 * se)
  
  est <- sum(r_c5_full * coef(SII_model_pooled))
  se <- sqrt(t(r_c5_full) %*% vcov(SII_model_pooled) %*% r_c5_full)
  SII_parms_r_c5 <- data.frame(Est = est, SE = se, SII = est, LCL = est - 1.96 * se, UCL = est + 1.96 * se)
  
  est <- sum(r_c6_full * coef(SII_model_pooled))
  se <- sqrt(t(r_c6_full) %*% vcov(SII_model_pooled) %*% r_c6_full)
  SII_parms_r_c6 <- data.frame(Est = est, SE = se, SII = est, LCL = est - 1.96 * se, UCL = est + 1.96 * se)
  
  SII_parms <- rbind(SII_parms_r_c1,SII_parms_r_c2,SII_parms_r_c3,SII_parms_r_c4,SII_parms_r_c5,SII_parms_r_c6)
  
  #Get RII and SII parameters
  RII_coefs <- as.numeric(RII_parms$Est)       
  RII_SEs <- as.numeric(RII_parms$SE)            
  RII_values <- as.numeric(RII_parms$RII)      
  RII_CIs <- RII_parms[4:5]  
  RII_trend_coef <- as.numeric(coef(RII_trend_model_pooled)["ridit_eqi_moecd:cycle_cont"])
  RII_trend_SE <- as.numeric(SE(RII_trend_model_pooled)["ridit_eqi_moecd:cycle_cont"])
  RII_trend_pvalue <- pt(abs(RII_trend_coef/RII_trend_SE),df=RII_trend_model[[1]]$df.residual,lower.tail=FALSE)*2
  
  RII_trend_sex_coef <- as.numeric(coef(RII_trend_sex_model_pooled)["ridit_eqi_moecd:cycle_cont:sexM"])
  RII_trend_sex_SE <- as.numeric(SE(RII_trend_sex_model_pooled)["ridit_eqi_moecd:cycle_cont:sexM"])
  RII_trend_sex_pvalue <- pt(abs(RII_trend_sex_coef/RII_trend_sex_SE),df=RII_trend_sex_model[[1]]$df.residual,lower.tail=FALSE)*2  
  
  SII_coefs <- as.numeric(SII_parms$Est)        
  SII_SEs <- as.numeric(SII_parms$SE)   
  SII_values <- as.numeric(SII_parms$SII)   
  SII_CIs <- SII_parms[4:5]
  SII_trend_coef <- as.numeric(coef(SII_trend_model_pooled)["ridit_eqi_moecd:cycle_cont"])
  SII_trend_SE <- as.numeric(SE(SII_trend_model_pooled)["ridit_eqi_moecd:cycle_cont"])
  SII_trend_pvalue <- pt(abs(SII_trend_coef/SII_trend_SE),df=SII_trend_model[[1]]$df.residual,lower.tail=FALSE)*2
  
  SII_trend_sex_coef <- as.numeric(coef(SII_trend_sex_model_pooled)["ridit_eqi_moecd:cycle_cont:sexM"])
  SII_trend_sex_SE <- as.numeric(SE(SII_trend_sex_model_pooled)["ridit_eqi_moecd:cycle_cont:sexM"])
  SII_trend_sex_pvalue <- pt(abs(SII_trend_sex_coef/SII_trend_sex_SE),df=SII_trend_sex_model[[1]]$df.residual,lower.tail=FALSE)*2        
  
  RII <- data.frame(
    Age_cat = age_cat_i,
    Sex = sex_i,
    Metric = 'RII',
    Cycle = c(1:6),
    Coef = RII_coefs,
    SE = RII_SEs,
    Value = RII_values,
    Lower_95_CI = as.numeric(RII_CIs$LCL),
    Upper_95_CI = as.numeric(RII_CIs$UCL),
    Trend_coef = RII_trend_coef,
    Trend_SE = RII_trend_SE,
    Trend_pvalue = RII_trend_pvalue,
    Trend_sex_coef = RII_trend_sex_coef,
    Trend_sex_SE = RII_trend_sex_SE,
    Trend_sex_pvalue = RII_trend_sex_pvalue
  )
  
  SII <- data.frame(
    Age_cat = age_cat_i,
    Sex = sex_i,
    Metric = 'SII',
    Cycle = c(1:6),
    Coef = SII_coefs,
    SE = SII_SEs,
    Value = SII_values,
    Lower_95_CI = as.numeric(SII_CIs$LCL),
    Upper_95_CI = as.numeric(SII_CIs$UCL),
    Trend_coef = SII_trend_coef,
    Trend_SE = SII_trend_SE,
    Trend_pvalue = SII_trend_pvalue,
    Trend_sex_coef = SII_trend_sex_coef,
    Trend_sex_SE = SII_trend_sex_SE,
    Trend_sex_pvalue = SII_trend_sex_pvalue
  )
  
  results_list[[i]] <- rbind(RII,SII)
  
  rm(RII_model,RII_trend_model,RII_trend_sex_model,SII_model,SII_trend_model,SII_trend_sex_model,
     RII_model_pooled,RII_trend_model_pooled,RII_trend_sex_model_pooled,SII_model_pooled,SII_trend_model_pooled,SII_trend_sex_model_pooled,
     r_c1,r_c2,r_c3,r_c4,r_c5,r_c6,r_c1_full,r_c2_full,r_c3_full,r_c4_full,r_c5_full,r_c6_full,est,se,
     RII_parms_r_c1,RII_parms_r_c2,RII_parms_r_c3,RII_parms_r_c4,RII_parms_r_c5,RII_parms_r_c6,RII_parms,
     SII_parms_r_c1,SII_parms_r_c2,SII_parms_r_c3,SII_parms_r_c4,SII_parms_r_c5,SII_parms_r_c6,SII_parms,
     RII_coefs,RII_SEs,RII_values,RII_CIs,RII_trend_coef,RII_trend_SE,RII_trend_pvalue,RII_trend_sex_coef,RII_trend_sex_SE,RII_trend_sex_pvalue,RII,
     SII_coefs,SII_SEs,SII_values,SII_CIs,SII_trend_coef,SII_trend_SE,SII_trend_pvalue,SII_trend_sex_coef,SII_trend_sex_SE,SII_trend_sex_pvalue,SII)
}

health_inequalities_age_sex_mi <- do.call(rbind, results_list)
rownames(health_inequalities_age_sex_mi) <- NULL

##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
# OOECD AND SQRT EQUIVALENCE SCALES              
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################

##############################################################################################################################
##############################################################################################################################
# Original OECD equivalence scale              
##############################################################################################################################
##############################################################################################################################

##############################################################################################################################
# Stratified by age category (with interactions to parse out cycle-specific effects) - w/o svycontrast
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
  
  #Run RII and SII models
  RII_model <- svyglm(r_poor_cv_health_1 ~ ridit_eqi_ooecd*cycle + age + sex + white + immigrant,
                      design = svy_design,
                      family = quasipoisson(link="log"),
                      subset = age_cat == age_cat_i)
  RII_trend_model <- svyglm(r_poor_cv_health_1 ~ ridit_eqi_ooecd*cycle_cont + age + sex + white + immigrant,
                            design = svy_design,
                            family = quasipoisson(link="log"),
                            subset = age_cat == age_cat_i)
  SII_model <- svyglm(r_poor_cv_health_1 ~ ridit_eqi_ooecd*cycle + age + sex + white + immigrant,
                      design = svy_design,
                      family = gaussian(link="identity"),
                      subset = age_cat == age_cat_i)
  SII_trend_model <- svyglm(r_poor_cv_health_1 ~ ridit_eqi_ooecd*cycle_cont + age + sex + white + immigrant,
                            design = svy_design,
                            family = gaussian(link="identity"),
                            subset = age_cat == age_cat_i)
  
  #Define contrast vectors
  r_c1 <- c("ridit_eqi_ooecd" = 1)
  r_c2 <- c("ridit_eqi_ooecd" = 1, "ridit_eqi_ooecd:cycle2" = 1)
  r_c3 <- c("ridit_eqi_ooecd" = 1, "ridit_eqi_ooecd:cycle3" = 1)
  r_c4 <- c("ridit_eqi_ooecd" = 1, "ridit_eqi_ooecd:cycle4" = 1)
  r_c5 <- c("ridit_eqi_ooecd" = 1, "ridit_eqi_ooecd:cycle5" = 1)
  r_c6 <- c("ridit_eqi_ooecd" = 1, "ridit_eqi_ooecd:cycle6" = 1)
  
  #Fill in zeros for other coefficients and compute contrast estimates, SE, RII/SII, and CIs
  #RII
  r_c1_full <- rep(0, length(coef(RII_model)))
  names(r_c1_full) <- names(coef(RII_model))
  r_c1_full[names(r_c1)] <- r_c1
  est <- sum(r_c1_full * coef(RII_model))
  se <- sqrt(t(r_c1_full) %*% vcov(RII_model) %*% r_c1_full)
  RII_parms_r_c1 <- data.frame(Est = est, SE = se, RII = exp(est), LCL = exp(est - 1.96 * se), UCL = exp(est + 1.96 * se))
  
  r_c2_full <- rep(0, length(coef(RII_model)))
  names(r_c2_full) <- names(coef(RII_model))
  r_c2_full[names(r_c2)] <- r_c2
  est <- sum(r_c2_full * coef(RII_model))
  se <- sqrt(t(r_c2_full) %*% vcov(RII_model) %*% r_c2_full)
  RII_parms_r_c2 <- data.frame(Est = est, SE = se, RII = exp(est), LCL = exp(est - 1.96 * se), UCL = exp(est + 1.96 * se))
  
  r_c3_full <- rep(0, length(coef(RII_model)))
  names(r_c3_full) <- names(coef(RII_model))
  r_c3_full[names(r_c3)] <- r_c3
  est <- sum(r_c3_full * coef(RII_model))
  se <- sqrt(t(r_c3_full) %*% vcov(RII_model) %*% r_c3_full)
  RII_parms_r_c3 <- data.frame(Est = est, SE = se, RII = exp(est), LCL = exp(est - 1.96 * se), UCL = exp(est + 1.96 * se))
  
  r_c4_full <- rep(0, length(coef(RII_model)))
  names(r_c4_full) <- names(coef(RII_model))
  r_c4_full[names(r_c4)] <- r_c4
  est <- sum(r_c4_full * coef(RII_model))
  se <- sqrt(t(r_c4_full) %*% vcov(RII_model) %*% r_c4_full)
  RII_parms_r_c4 <- data.frame(Est = est, SE = se, RII = exp(est), LCL = exp(est - 1.96 * se), UCL = exp(est + 1.96 * se))
  
  r_c5_full <- rep(0, length(coef(RII_model)))
  names(r_c5_full) <- names(coef(RII_model))
  r_c5_full[names(r_c5)] <- r_c5
  est <- sum(r_c5_full * coef(RII_model))
  se <- sqrt(t(r_c5_full) %*% vcov(RII_model) %*% r_c5_full)
  RII_parms_r_c5 <- data.frame(Est = est, SE = se, RII = exp(est), LCL = exp(est - 1.96 * se), UCL = exp(est + 1.96 * se))
  
  r_c6_full <- rep(0, length(coef(RII_model)))
  names(r_c6_full) <- names(coef(RII_model))
  r_c6_full[names(r_c6)] <- r_c6
  est <- sum(r_c6_full * coef(RII_model))
  se <- sqrt(t(r_c6_full) %*% vcov(RII_model) %*% r_c6_full)
  RII_parms_r_c6 <- data.frame(Est = est, SE = se, RII = exp(est), LCL = exp(est - 1.96 * se), UCL = exp(est + 1.96 * se))
  
  RII_parms <- rbind(RII_parms_r_c1,RII_parms_r_c2,RII_parms_r_c3,RII_parms_r_c4,RII_parms_r_c5,RII_parms_r_c6)
  
  #SII
  est <- sum(r_c1_full * coef(SII_model))
  se <- sqrt(t(r_c1_full) %*% vcov(SII_model) %*% r_c1_full)
  SII_parms_r_c1 <- data.frame(Est = est, SE = se, SII = est, LCL = est - 1.96 * se, UCL = est + 1.96 * se)
  
  est <- sum(r_c2_full * coef(SII_model))
  se <- sqrt(t(r_c2_full) %*% vcov(SII_model) %*% r_c2_full)
  SII_parms_r_c2 <- data.frame(Est = est, SE = se, SII = est, LCL = est - 1.96 * se, UCL = est + 1.96 * se)
  
  est <- sum(r_c3_full * coef(SII_model))
  se <- sqrt(t(r_c3_full) %*% vcov(SII_model) %*% r_c3_full)
  SII_parms_r_c3 <- data.frame(Est = est, SE = se, SII = est, LCL = est - 1.96 * se, UCL = est + 1.96 * se)
  
  est <- sum(r_c4_full * coef(SII_model))
  se <- sqrt(t(r_c4_full) %*% vcov(SII_model) %*% r_c4_full)
  SII_parms_r_c4 <- data.frame(Est = est, SE = se, SII = est, LCL = est - 1.96 * se, UCL = est + 1.96 * se)
  
  est <- sum(r_c5_full * coef(SII_model))
  se <- sqrt(t(r_c5_full) %*% vcov(SII_model) %*% r_c5_full)
  SII_parms_r_c5 <- data.frame(Est = est, SE = se, SII = est, LCL = est - 1.96 * se, UCL = est + 1.96 * se)
  
  est <- sum(r_c6_full * coef(SII_model))
  se <- sqrt(t(r_c6_full) %*% vcov(SII_model) %*% r_c6_full)
  SII_parms_r_c6 <- data.frame(Est = est, SE = se, SII = est, LCL = est - 1.96 * se, UCL = est + 1.96 * se)
  
  SII_parms <- rbind(SII_parms_r_c1,SII_parms_r_c2,SII_parms_r_c3,SII_parms_r_c4,SII_parms_r_c5,SII_parms_r_c6)
  
  #Get RII and SII parameters
  RII_coefs <- as.numeric(RII_parms$Est)  
  RII_SEs <- as.numeric(RII_parms$SE)             
  RII_values <- as.numeric(RII_parms$RII)   
  RII_CIs <- RII_parms[4:5]   
  RII_trend_coef <- as.numeric(coef(RII_trend_model)["ridit_eqi_ooecd:cycle_cont"])
  RII_trend_SE <- as.numeric(SE(RII_trend_model)["ridit_eqi_ooecd:cycle_cont"])
  RII_trend_pvalue <- summary(RII_trend_model)$coefficients["ridit_eqi_ooecd:cycle_cont", "Pr(>|t|)"]
  
  SII_coefs <- as.numeric(SII_parms$Est)     
  SII_SEs <- as.numeric(SII_parms$SE) 
  SII_values <- as.numeric(SII_parms$SII)   
  SII_CIs <- SII_parms[4:5] 
  SII_trend_coef <- as.numeric(coef(SII_trend_model)["ridit_eqi_ooecd:cycle_cont"])
  SII_trend_SE <- as.numeric(SE(SII_trend_model)["ridit_eqi_ooecd:cycle_cont"])
  SII_trend_pvalue <- summary(SII_trend_model)$coefficients["ridit_eqi_ooecd:cycle_cont", "Pr(>|t|)"]
  
  RII <- data.frame(
    Age_cat = age_cat_i,
    Metric = 'RII',
    Cycle = c(1:6),
    Coef = RII_coefs,
    SE = RII_SEs,
    Value = RII_values,
    Lower_95_CI = as.numeric(RII_CIs$LCL),
    Upper_95_CI = as.numeric(RII_CIs$UCL),
    Trend_coef = RII_trend_coef,
    Trend_SE = RII_trend_SE,
    Trend_pvalue = RII_trend_pvalue 
  )
  
  SII <- data.frame(
    Age_cat = age_cat_i,
    Metric = 'SII',
    Cycle = c(1:6),
    Coef = SII_coefs,
    SE = SII_SEs,
    Value = SII_values,
    Lower_95_CI = as.numeric(SII_CIs$LCL),
    Upper_95_CI = as.numeric(SII_CIs$UCL),
    Trend_coef = SII_trend_coef,
    Trend_SE = SII_trend_SE,
    Trend_pvalue = SII_trend_pvalue 
  )
  
  results_list[[i]] <- rbind(RII,SII)
  rm(RII_model,RII_trend_model,SII_model,SII_trend_model,
     r_c1,r_c2,r_c3,r_c4,r_c5,r_c6,r_c1_full,r_c2_full,r_c3_full,r_c4_full,r_c5_full,r_c6_full,est,se,
     RII_parms_r_c1,RII_parms_r_c2,RII_parms_r_c3,RII_parms_r_c4,RII_parms_r_c5,RII_parms_r_c6,RII_parms,
     SII_parms_r_c1,SII_parms_r_c2,SII_parms_r_c3,SII_parms_r_c4,SII_parms_r_c5,SII_parms_r_c6,SII_parms,
     RII_coefs,RII_SEs,RII_values,RII_CIs,RII_trend_coef,RII_trend_SE,RII_trend_pvalue,RII,
     SII_coefs,SII_SEs,SII_values,SII_CIs,SII_trend_coef,SII_trend_SE,SII_trend_pvalue,SII)
}

health_inequalities_age_ooecd <- do.call(rbind, results_list)
rownames(health_inequalities_age_ooecd) <- NULL
rm(results_list)

##############################################################################################################################
# Stratified by age category and sex (with interactions to parse out cycle-specific effects) - w/o svycontrast
##############################################################################################################################
#Define possible values for age_cat
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
  
  #Run RII and SII models, stratified by sex
  RII_model <- svyglm(r_poor_cv_health_1 ~ ridit_eqi_ooecd_sex*cycle + age + white + immigrant,
                      design = svy_design,
                      family = quasipoisson(link="log"),
                      subset = age_cat == age_cat_i & sex == sex_i)
  RII_trend_model <- svyglm(r_poor_cv_health_1 ~ ridit_eqi_ooecd_sex*cycle_cont + age + white + immigrant,
                            design = svy_design,
                            family = quasipoisson(link="log"),
                            subset = age_cat == age_cat_i & sex == sex_i)
  RII_trend_sex_model <- svyglm(r_poor_cv_health_1 ~ ridit_eqi_ooecd_sex*cycle_cont*sex + age + white + immigrant,
                                design = svy_design,
                                family = quasipoisson(link="log"),
                                subset = age_cat == age_cat_i)
  SII_model <- svyglm(r_poor_cv_health_1 ~ ridit_eqi_ooecd_sex*cycle + age + white + immigrant,
                      design = svy_design,
                      family = gaussian(link="identity"),
                      subset = age_cat == age_cat_i & sex == sex_i)
  SII_trend_model <- svyglm(r_poor_cv_health_1 ~ ridit_eqi_ooecd_sex*cycle_cont + age + white + immigrant,
                            design = svy_design,
                            family = gaussian(link="identity"),
                            subset = age_cat == age_cat_i & sex == sex_i)
  SII_trend_sex_model <- svyglm(r_poor_cv_health_1 ~ ridit_eqi_ooecd_sex*cycle_cont*sex + age + white + immigrant,
                                design = svy_design,
                                family = gaussian(link="identity"),
                                subset = age_cat == age_cat_i)
  
  #Define contrast vectors
  r_c1 <- c("ridit_eqi_ooecd_sex" = 1)
  r_c2 <- c("ridit_eqi_ooecd_sex" = 1, "ridit_eqi_ooecd_sex:cycle2" = 1)
  r_c3 <- c("ridit_eqi_ooecd_sex" = 1, "ridit_eqi_ooecd_sex:cycle3" = 1)
  r_c4 <- c("ridit_eqi_ooecd_sex" = 1, "ridit_eqi_ooecd_sex:cycle4" = 1)
  r_c5 <- c("ridit_eqi_ooecd_sex" = 1, "ridit_eqi_ooecd_sex:cycle5" = 1)
  r_c6 <- c("ridit_eqi_ooecd_sex" = 1, "ridit_eqi_ooecd_sex:cycle6" = 1)
  
  #Fill in zeros for other coefficients and compute contrast estimates, SE, RII/SII, and CIs
  #RII
  r_c1_full <- rep(0, length(coef(RII_model)))
  names(r_c1_full) <- names(coef(RII_model))
  r_c1_full[names(r_c1)] <- r_c1
  est <- sum(r_c1_full * coef(RII_model))
  se <- sqrt(t(r_c1_full) %*% vcov(RII_model) %*% r_c1_full)
  RII_parms_r_c1 <- data.frame(Est = est, SE = se, RII = exp(est), LCL = exp(est - 1.96 * se), UCL = exp(est + 1.96 * se))
  
  r_c2_full <- rep(0, length(coef(RII_model)))
  names(r_c2_full) <- names(coef(RII_model))
  r_c2_full[names(r_c2)] <- r_c2
  est <- sum(r_c2_full * coef(RII_model))
  se <- sqrt(t(r_c2_full) %*% vcov(RII_model) %*% r_c2_full)
  RII_parms_r_c2 <- data.frame(Est = est, SE = se, RII = exp(est), LCL = exp(est - 1.96 * se), UCL = exp(est + 1.96 * se))
  
  r_c3_full <- rep(0, length(coef(RII_model)))
  names(r_c3_full) <- names(coef(RII_model))
  r_c3_full[names(r_c3)] <- r_c3
  est <- sum(r_c3_full * coef(RII_model))
  se <- sqrt(t(r_c3_full) %*% vcov(RII_model) %*% r_c3_full)
  RII_parms_r_c3 <- data.frame(Est = est, SE = se, RII = exp(est), LCL = exp(est - 1.96 * se), UCL = exp(est + 1.96 * se))
  
  r_c4_full <- rep(0, length(coef(RII_model)))
  names(r_c4_full) <- names(coef(RII_model))
  r_c4_full[names(r_c4)] <- r_c4
  est <- sum(r_c4_full * coef(RII_model))
  se <- sqrt(t(r_c4_full) %*% vcov(RII_model) %*% r_c4_full)
  RII_parms_r_c4 <- data.frame(Est = est, SE = se, RII = exp(est), LCL = exp(est - 1.96 * se), UCL = exp(est + 1.96 * se))
  
  r_c5_full <- rep(0, length(coef(RII_model)))
  names(r_c5_full) <- names(coef(RII_model))
  r_c5_full[names(r_c5)] <- r_c5
  est <- sum(r_c5_full * coef(RII_model))
  se <- sqrt(t(r_c5_full) %*% vcov(RII_model) %*% r_c5_full)
  RII_parms_r_c5 <- data.frame(Est = est, SE = se, RII = exp(est), LCL = exp(est - 1.96 * se), UCL = exp(est + 1.96 * se))
  
  r_c6_full <- rep(0, length(coef(RII_model)))
  names(r_c6_full) <- names(coef(RII_model))
  r_c6_full[names(r_c6)] <- r_c6
  est <- sum(r_c6_full * coef(RII_model))
  se <- sqrt(t(r_c6_full) %*% vcov(RII_model) %*% r_c6_full)
  RII_parms_r_c6 <- data.frame(Est = est, SE = se, RII = exp(est), LCL = exp(est - 1.96 * se), UCL = exp(est + 1.96 * se))
  
  RII_parms <- rbind(RII_parms_r_c1,RII_parms_r_c2,RII_parms_r_c3,RII_parms_r_c4,RII_parms_r_c5,RII_parms_r_c6)
  
  #SII
  est <- sum(r_c1_full * coef(SII_model))
  se <- sqrt(t(r_c1_full) %*% vcov(SII_model) %*% r_c1_full)
  SII_parms_r_c1 <- data.frame(Est = est, SE = se, SII = est, LCL = est - 1.96 * se, UCL = est + 1.96 * se)
  
  est <- sum(r_c2_full * coef(SII_model))
  se <- sqrt(t(r_c2_full) %*% vcov(SII_model) %*% r_c2_full)
  SII_parms_r_c2 <- data.frame(Est = est, SE = se, SII = est, LCL = est - 1.96 * se, UCL = est + 1.96 * se)
  
  est <- sum(r_c3_full * coef(SII_model))
  se <- sqrt(t(r_c3_full) %*% vcov(SII_model) %*% r_c3_full)
  SII_parms_r_c3 <- data.frame(Est = est, SE = se, SII = est, LCL = est - 1.96 * se, UCL = est + 1.96 * se)
  
  est <- sum(r_c4_full * coef(SII_model))
  se <- sqrt(t(r_c4_full) %*% vcov(SII_model) %*% r_c4_full)
  SII_parms_r_c4 <- data.frame(Est = est, SE = se, SII = est, LCL = est - 1.96 * se, UCL = est + 1.96 * se)
  
  est <- sum(r_c5_full * coef(SII_model))
  se <- sqrt(t(r_c5_full) %*% vcov(SII_model) %*% r_c5_full)
  SII_parms_r_c5 <- data.frame(Est = est, SE = se, SII = est, LCL = est - 1.96 * se, UCL = est + 1.96 * se)
  
  est <- sum(r_c6_full * coef(SII_model))
  se <- sqrt(t(r_c6_full) %*% vcov(SII_model) %*% r_c6_full)
  SII_parms_r_c6 <- data.frame(Est = est, SE = se, SII = est, LCL = est - 1.96 * se, UCL = est + 1.96 * se)
  
  SII_parms <- rbind(SII_parms_r_c1,SII_parms_r_c2,SII_parms_r_c3,SII_parms_r_c4,SII_parms_r_c5,SII_parms_r_c6)
  
  #Get RII and SII parameters
  RII_coefs <- as.numeric(RII_parms$Est)       
  RII_SEs <- as.numeric(RII_parms$SE)            
  RII_values <- as.numeric(RII_parms$RII)      
  RII_CIs <- RII_parms[4:5]  
  RII_trend_coef <- as.numeric(coef(RII_trend_model)["ridit_eqi_ooecd_sex:cycle_cont"])
  RII_trend_SE <- as.numeric(SE(RII_trend_model)["ridit_eqi_ooecd_sex:cycle_cont"])
  RII_trend_pvalue <- summary(RII_trend_model)$coefficients["ridit_eqi_ooecd_sex:cycle_cont", "Pr(>|t|)"]
  RII_trend_sex_coef <- as.numeric(coef(RII_trend_sex_model)["ridit_eqi_ooecd_sex:cycle_cont:sexM"])
  RII_trend_sex_SE <- as.numeric(SE(RII_trend_sex_model)["ridit_eqi_ooecd_sex:cycle_cont:sexM"])
  RII_trend_sex_pvalue <- summary(RII_trend_sex_model)$coefficients["ridit_eqi_ooecd_sex:cycle_cont:sexM", "Pr(>|t|)"]  
  
  SII_coefs <- as.numeric(SII_parms$Est)        
  SII_SEs <- as.numeric(SII_parms$SE)   
  SII_values <- as.numeric(SII_parms$SII)   
  SII_CIs <- SII_parms[4:5]
  SII_trend_coef <- as.numeric(coef(SII_trend_model)["ridit_eqi_ooecd_sex:cycle_cont"])
  SII_trend_SE <- as.numeric(SE(SII_trend_model)["ridit_eqi_ooecd_sex:cycle_cont"])
  SII_trend_pvalue <- summary(SII_trend_model)$coefficients["ridit_eqi_ooecd_sex:cycle_cont", "Pr(>|t|)"]
  SII_trend_sex_coef <- as.numeric(coef(SII_trend_sex_model)["ridit_eqi_ooecd_sex:cycle_cont:sexM"])
  SII_trend_sex_SE <- as.numeric(SE(SII_trend_sex_model)["ridit_eqi_ooecd_sex:cycle_cont:sexM"])
  SII_trend_sex_pvalue <- summary(SII_trend_sex_model)$coefficients["ridit_eqi_ooecd_sex:cycle_cont:sexM", "Pr(>|t|)"]      
  
  RII <- data.frame(
    Age_cat = age_cat_i,
    Sex = sex_i,
    Metric = 'RII',
    Cycle = c(1:6),
    Coef = RII_coefs,
    SE = RII_SEs,
    Value = RII_values,
    Lower_95_CI = as.numeric(RII_CIs$LCL),
    Upper_95_CI = as.numeric(RII_CIs$UCL),
    Trend_coef = RII_trend_coef,
    Trend_SE = RII_trend_SE,
    Trend_pvalue = RII_trend_pvalue,
    Trend_sex_coef = RII_trend_sex_coef,
    Trend_sex_SE = RII_trend_sex_SE,
    Trend_sex_pvalue = RII_trend_sex_pvalue
  )
  
  SII <- data.frame(
    Age_cat = age_cat_i,
    Sex = sex_i,
    Metric = 'SII',
    Cycle = c(1:6),
    Coef = SII_coefs,
    SE = SII_SEs,
    Value = SII_values,
    Lower_95_CI = as.numeric(SII_CIs$LCL),
    Upper_95_CI = as.numeric(SII_CIs$UCL),
    Trend_coef = SII_trend_coef,
    Trend_SE = SII_trend_SE,
    Trend_pvalue = SII_trend_pvalue,
    Trend_sex_coef = SII_trend_sex_coef,
    Trend_sex_SE = SII_trend_sex_SE,
    Trend_sex_pvalue = SII_trend_sex_pvalue
  )
  
  results_list[[i]] <- rbind(RII,SII)
  rm(RII_model,RII_trend_model,RII_trend_sex_model,SII_model,SII_trend_model,SII_trend_sex_model,
     r_c1,r_c2,r_c3,r_c4,r_c5,r_c6,r_c1_full,r_c2_full,r_c3_full,r_c4_full,r_c5_full,r_c6_full,est,se,
     RII_parms_r_c1,RII_parms_r_c2,RII_parms_r_c3,RII_parms_r_c4,RII_parms_r_c5,RII_parms_r_c6,RII_parms,
     SII_parms_r_c1,SII_parms_r_c2,SII_parms_r_c3,SII_parms_r_c4,SII_parms_r_c5,SII_parms_r_c6,SII_parms,
     RII_coefs,RII_SEs,RII_values,RII_CIs,RII_trend_coef,RII_trend_SE,RII_trend_pvalue,RII_trend_sex_coef,RII_trend_sex_SE,RII_trend_sex_pvalue,RII,
     SII_coefs,SII_SEs,SII_values,SII_CIs,SII_trend_coef,SII_trend_SE,SII_trend_pvalue,SII_trend_sex_coef,SII_trend_sex_SE,SII_trend_sex_pvalue,SII)
}

health_inequalities_age_sex_ooecd <- do.call(rbind, results_list)
rownames(health_inequalities_age_sex_ooecd) <- NULL

##############################################################################################################################
##############################################################################################################################
# Sqaure root equivalence scale              
##############################################################################################################################
##############################################################################################################################

##############################################################################################################################
# Stratified by age category (with interactions to parse out cycle-specific effects) - w/o svycontrast
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
  
  #Run RII and SII models
  RII_model <- svyglm(r_poor_cv_health_1 ~ ridit_eqi_sqrt*cycle + age + sex + white + immigrant,
                      design = svy_design,
                      family = quasipoisson(link="log"),
                      subset = age_cat == age_cat_i)
  RII_trend_model <- svyglm(r_poor_cv_health_1 ~ ridit_eqi_sqrt*cycle_cont + age + sex + white + immigrant,
                            design = svy_design,
                            family = quasipoisson(link="log"),
                            subset = age_cat == age_cat_i)
  SII_model <- svyglm(r_poor_cv_health_1 ~ ridit_eqi_sqrt*cycle + age + sex + white + immigrant,
                      design = svy_design,
                      family = gaussian(link="identity"),
                      subset = age_cat == age_cat_i)
  SII_trend_model <- svyglm(r_poor_cv_health_1 ~ ridit_eqi_sqrt*cycle_cont + age + sex + white + immigrant,
                            design = svy_design,
                            family = gaussian(link="identity"),
                            subset = age_cat == age_cat_i)
  
  #Define contrast vectors
  r_c1 <- c("ridit_eqi_sqrt" = 1)
  r_c2 <- c("ridit_eqi_sqrt" = 1, "ridit_eqi_sqrt:cycle2" = 1)
  r_c3 <- c("ridit_eqi_sqrt" = 1, "ridit_eqi_sqrt:cycle3" = 1)
  r_c4 <- c("ridit_eqi_sqrt" = 1, "ridit_eqi_sqrt:cycle4" = 1)
  r_c5 <- c("ridit_eqi_sqrt" = 1, "ridit_eqi_sqrt:cycle5" = 1)
  r_c6 <- c("ridit_eqi_sqrt" = 1, "ridit_eqi_sqrt:cycle6" = 1)
  
  #Fill in zeros for other coefficients and compute contrast estimates, SE, RII/SII, and CIs
  #RII
  r_c1_full <- rep(0, length(coef(RII_model)))
  names(r_c1_full) <- names(coef(RII_model))
  r_c1_full[names(r_c1)] <- r_c1
  est <- sum(r_c1_full * coef(RII_model))
  se <- sqrt(t(r_c1_full) %*% vcov(RII_model) %*% r_c1_full)
  RII_parms_r_c1 <- data.frame(Est = est, SE = se, RII = exp(est), LCL = exp(est - 1.96 * se), UCL = exp(est + 1.96 * se))
  
  r_c2_full <- rep(0, length(coef(RII_model)))
  names(r_c2_full) <- names(coef(RII_model))
  r_c2_full[names(r_c2)] <- r_c2
  est <- sum(r_c2_full * coef(RII_model))
  se <- sqrt(t(r_c2_full) %*% vcov(RII_model) %*% r_c2_full)
  RII_parms_r_c2 <- data.frame(Est = est, SE = se, RII = exp(est), LCL = exp(est - 1.96 * se), UCL = exp(est + 1.96 * se))
  
  r_c3_full <- rep(0, length(coef(RII_model)))
  names(r_c3_full) <- names(coef(RII_model))
  r_c3_full[names(r_c3)] <- r_c3
  est <- sum(r_c3_full * coef(RII_model))
  se <- sqrt(t(r_c3_full) %*% vcov(RII_model) %*% r_c3_full)
  RII_parms_r_c3 <- data.frame(Est = est, SE = se, RII = exp(est), LCL = exp(est - 1.96 * se), UCL = exp(est + 1.96 * se))
  
  r_c4_full <- rep(0, length(coef(RII_model)))
  names(r_c4_full) <- names(coef(RII_model))
  r_c4_full[names(r_c4)] <- r_c4
  est <- sum(r_c4_full * coef(RII_model))
  se <- sqrt(t(r_c4_full) %*% vcov(RII_model) %*% r_c4_full)
  RII_parms_r_c4 <- data.frame(Est = est, SE = se, RII = exp(est), LCL = exp(est - 1.96 * se), UCL = exp(est + 1.96 * se))
  
  r_c5_full <- rep(0, length(coef(RII_model)))
  names(r_c5_full) <- names(coef(RII_model))
  r_c5_full[names(r_c5)] <- r_c5
  est <- sum(r_c5_full * coef(RII_model))
  se <- sqrt(t(r_c5_full) %*% vcov(RII_model) %*% r_c5_full)
  RII_parms_r_c5 <- data.frame(Est = est, SE = se, RII = exp(est), LCL = exp(est - 1.96 * se), UCL = exp(est + 1.96 * se))
  
  r_c6_full <- rep(0, length(coef(RII_model)))
  names(r_c6_full) <- names(coef(RII_model))
  r_c6_full[names(r_c6)] <- r_c6
  est <- sum(r_c6_full * coef(RII_model))
  se <- sqrt(t(r_c6_full) %*% vcov(RII_model) %*% r_c6_full)
  RII_parms_r_c6 <- data.frame(Est = est, SE = se, RII = exp(est), LCL = exp(est - 1.96 * se), UCL = exp(est + 1.96 * se))
  
  RII_parms <- rbind(RII_parms_r_c1,RII_parms_r_c2,RII_parms_r_c3,RII_parms_r_c4,RII_parms_r_c5,RII_parms_r_c6)
  
  #SII
  est <- sum(r_c1_full * coef(SII_model))
  se <- sqrt(t(r_c1_full) %*% vcov(SII_model) %*% r_c1_full)
  SII_parms_r_c1 <- data.frame(Est = est, SE = se, SII = est, LCL = est - 1.96 * se, UCL = est + 1.96 * se)
  
  est <- sum(r_c2_full * coef(SII_model))
  se <- sqrt(t(r_c2_full) %*% vcov(SII_model) %*% r_c2_full)
  SII_parms_r_c2 <- data.frame(Est = est, SE = se, SII = est, LCL = est - 1.96 * se, UCL = est + 1.96 * se)
  
  est <- sum(r_c3_full * coef(SII_model))
  se <- sqrt(t(r_c3_full) %*% vcov(SII_model) %*% r_c3_full)
  SII_parms_r_c3 <- data.frame(Est = est, SE = se, SII = est, LCL = est - 1.96 * se, UCL = est + 1.96 * se)
  
  est <- sum(r_c4_full * coef(SII_model))
  se <- sqrt(t(r_c4_full) %*% vcov(SII_model) %*% r_c4_full)
  SII_parms_r_c4 <- data.frame(Est = est, SE = se, SII = est, LCL = est - 1.96 * se, UCL = est + 1.96 * se)
  
  est <- sum(r_c5_full * coef(SII_model))
  se <- sqrt(t(r_c5_full) %*% vcov(SII_model) %*% r_c5_full)
  SII_parms_r_c5 <- data.frame(Est = est, SE = se, SII = est, LCL = est - 1.96 * se, UCL = est + 1.96 * se)
  
  est <- sum(r_c6_full * coef(SII_model))
  se <- sqrt(t(r_c6_full) %*% vcov(SII_model) %*% r_c6_full)
  SII_parms_r_c6 <- data.frame(Est = est, SE = se, SII = est, LCL = est - 1.96 * se, UCL = est + 1.96 * se)
  
  SII_parms <- rbind(SII_parms_r_c1,SII_parms_r_c2,SII_parms_r_c3,SII_parms_r_c4,SII_parms_r_c5,SII_parms_r_c6)
  
  #Get RII and SII parameters
  RII_coefs <- as.numeric(RII_parms$Est)  
  RII_SEs <- as.numeric(RII_parms$SE)             
  RII_values <- as.numeric(RII_parms$RII)   
  RII_CIs <- RII_parms[4:5]   
  RII_trend_coef <- as.numeric(coef(RII_trend_model)["ridit_eqi_sqrt:cycle_cont"])
  RII_trend_SE <- as.numeric(SE(RII_trend_model)["ridit_eqi_sqrt:cycle_cont"])
  RII_trend_pvalue <- summary(RII_trend_model)$coefficients["ridit_eqi_sqrt:cycle_cont", "Pr(>|t|)"]
  
  SII_coefs <- as.numeric(SII_parms$Est)     
  SII_SEs <- as.numeric(SII_parms$SE) 
  SII_values <- as.numeric(SII_parms$SII)   
  SII_CIs <- SII_parms[4:5] 
  SII_trend_coef <- as.numeric(coef(SII_trend_model)["ridit_eqi_sqrt:cycle_cont"])
  SII_trend_SE <- as.numeric(SE(SII_trend_model)["ridit_eqi_sqrt:cycle_cont"])
  SII_trend_pvalue <- summary(SII_trend_model)$coefficients["ridit_eqi_sqrt:cycle_cont", "Pr(>|t|)"]
  
  RII <- data.frame(
    Age_cat = age_cat_i,
    Metric = 'RII',
    Cycle = c(1:6),
    Coef = RII_coefs,
    SE = RII_SEs,
    Value = RII_values,
    Lower_95_CI = as.numeric(RII_CIs$LCL),
    Upper_95_CI = as.numeric(RII_CIs$UCL),
    Trend_coef = RII_trend_coef,
    Trend_SE = RII_trend_SE,
    Trend_pvalue = RII_trend_pvalue 
  )
  
  SII <- data.frame(
    Age_cat = age_cat_i,
    Metric = 'SII',
    Cycle = c(1:6),
    Coef = SII_coefs,
    SE = SII_SEs,
    Value = SII_values,
    Lower_95_CI = as.numeric(SII_CIs$LCL),
    Upper_95_CI = as.numeric(SII_CIs$UCL),
    Trend_coef = SII_trend_coef,
    Trend_SE = SII_trend_SE,
    Trend_pvalue = SII_trend_pvalue 
  )
  
  results_list[[i]] <- rbind(RII,SII)
  rm(RII_model,RII_trend_model,SII_model,SII_trend_model,
     r_c1,r_c2,r_c3,r_c4,r_c5,r_c6,r_c1_full,r_c2_full,r_c3_full,r_c4_full,r_c5_full,r_c6_full,est,se,
     RII_parms_r_c1,RII_parms_r_c2,RII_parms_r_c3,RII_parms_r_c4,RII_parms_r_c5,RII_parms_r_c6,RII_parms,
     SII_parms_r_c1,SII_parms_r_c2,SII_parms_r_c3,SII_parms_r_c4,SII_parms_r_c5,SII_parms_r_c6,SII_parms,
     RII_coefs,RII_SEs,RII_values,RII_CIs,RII_trend_coef,RII_trend_SE,RII_trend_pvalue,RII,
     SII_coefs,SII_SEs,SII_values,SII_CIs,SII_trend_coef,SII_trend_SE,SII_trend_pvalue,SII)
}

health_inequalities_age_sqrt <- do.call(rbind, results_list)
rownames(health_inequalities_age_sqrt) <- NULL
rm(results_list)

##############################################################################################################################
# Stratified by age category and sex (with interactions to parse out cycle-specific effects) - w/o svycontrast
##############################################################################################################################
#Define possible values for age_cat
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
  
  #Run RII and SII models, stratified by sex
  RII_model <- svyglm(r_poor_cv_health_1 ~ ridit_eqi_sqrt_sex*cycle + age + white + immigrant,
                      design = svy_design,
                      family = quasipoisson(link="log"),
                      subset = age_cat == age_cat_i & sex == sex_i)
  RII_trend_model <- svyglm(r_poor_cv_health_1 ~ ridit_eqi_sqrt_sex*cycle_cont + age + white + immigrant,
                            design = svy_design,
                            family = quasipoisson(link="log"),
                            subset = age_cat == age_cat_i & sex == sex_i)
  RII_trend_sex_model <- svyglm(r_poor_cv_health_1 ~ ridit_eqi_sqrt_sex*cycle_cont*sex + age + white + immigrant,
                                design = svy_design,
                                family = quasipoisson(link="log"),
                                subset = age_cat == age_cat_i)
  SII_model <- svyglm(r_poor_cv_health_1 ~ ridit_eqi_sqrt_sex*cycle + age + white + immigrant,
                      design = svy_design,
                      family = gaussian(link="identity"),
                      subset = age_cat == age_cat_i & sex == sex_i)
  SII_trend_model <- svyglm(r_poor_cv_health_1 ~ ridit_eqi_sqrt_sex*cycle_cont + age + white + immigrant,
                            design = svy_design,
                            family = gaussian(link="identity"),
                            subset = age_cat == age_cat_i & sex == sex_i)
  SII_trend_sex_model <- svyglm(r_poor_cv_health_1 ~ ridit_eqi_sqrt_sex*cycle_cont*sex + age + white + immigrant,
                                design = svy_design,
                                family = gaussian(link="identity"),
                                subset = age_cat == age_cat_i)
  
  #Define contrast vectors
  r_c1 <- c("ridit_eqi_sqrt_sex" = 1)
  r_c2 <- c("ridit_eqi_sqrt_sex" = 1, "ridit_eqi_sqrt_sex:cycle2" = 1)
  r_c3 <- c("ridit_eqi_sqrt_sex" = 1, "ridit_eqi_sqrt_sex:cycle3" = 1)
  r_c4 <- c("ridit_eqi_sqrt_sex" = 1, "ridit_eqi_sqrt_sex:cycle4" = 1)
  r_c5 <- c("ridit_eqi_sqrt_sex" = 1, "ridit_eqi_sqrt_sex:cycle5" = 1)
  r_c6 <- c("ridit_eqi_sqrt_sex" = 1, "ridit_eqi_sqrt_sex:cycle6" = 1)
  
  #Fill in zeros for other coefficients and compute contrast estimates, SE, RII/SII, and CIs
  #RII
  r_c1_full <- rep(0, length(coef(RII_model)))
  names(r_c1_full) <- names(coef(RII_model))
  r_c1_full[names(r_c1)] <- r_c1
  est <- sum(r_c1_full * coef(RII_model))
  se <- sqrt(t(r_c1_full) %*% vcov(RII_model) %*% r_c1_full)
  RII_parms_r_c1 <- data.frame(Est = est, SE = se, RII = exp(est), LCL = exp(est - 1.96 * se), UCL = exp(est + 1.96 * se))
  
  r_c2_full <- rep(0, length(coef(RII_model)))
  names(r_c2_full) <- names(coef(RII_model))
  r_c2_full[names(r_c2)] <- r_c2
  est <- sum(r_c2_full * coef(RII_model))
  se <- sqrt(t(r_c2_full) %*% vcov(RII_model) %*% r_c2_full)
  RII_parms_r_c2 <- data.frame(Est = est, SE = se, RII = exp(est), LCL = exp(est - 1.96 * se), UCL = exp(est + 1.96 * se))
  
  r_c3_full <- rep(0, length(coef(RII_model)))
  names(r_c3_full) <- names(coef(RII_model))
  r_c3_full[names(r_c3)] <- r_c3
  est <- sum(r_c3_full * coef(RII_model))
  se <- sqrt(t(r_c3_full) %*% vcov(RII_model) %*% r_c3_full)
  RII_parms_r_c3 <- data.frame(Est = est, SE = se, RII = exp(est), LCL = exp(est - 1.96 * se), UCL = exp(est + 1.96 * se))
  
  r_c4_full <- rep(0, length(coef(RII_model)))
  names(r_c4_full) <- names(coef(RII_model))
  r_c4_full[names(r_c4)] <- r_c4
  est <- sum(r_c4_full * coef(RII_model))
  se <- sqrt(t(r_c4_full) %*% vcov(RII_model) %*% r_c4_full)
  RII_parms_r_c4 <- data.frame(Est = est, SE = se, RII = exp(est), LCL = exp(est - 1.96 * se), UCL = exp(est + 1.96 * se))
  
  r_c5_full <- rep(0, length(coef(RII_model)))
  names(r_c5_full) <- names(coef(RII_model))
  r_c5_full[names(r_c5)] <- r_c5
  est <- sum(r_c5_full * coef(RII_model))
  se <- sqrt(t(r_c5_full) %*% vcov(RII_model) %*% r_c5_full)
  RII_parms_r_c5 <- data.frame(Est = est, SE = se, RII = exp(est), LCL = exp(est - 1.96 * se), UCL = exp(est + 1.96 * se))
  
  r_c6_full <- rep(0, length(coef(RII_model)))
  names(r_c6_full) <- names(coef(RII_model))
  r_c6_full[names(r_c6)] <- r_c6
  est <- sum(r_c6_full * coef(RII_model))
  se <- sqrt(t(r_c6_full) %*% vcov(RII_model) %*% r_c6_full)
  RII_parms_r_c6 <- data.frame(Est = est, SE = se, RII = exp(est), LCL = exp(est - 1.96 * se), UCL = exp(est + 1.96 * se))
  
  RII_parms <- rbind(RII_parms_r_c1,RII_parms_r_c2,RII_parms_r_c3,RII_parms_r_c4,RII_parms_r_c5,RII_parms_r_c6)
  
  #SII
  est <- sum(r_c1_full * coef(SII_model))
  se <- sqrt(t(r_c1_full) %*% vcov(SII_model) %*% r_c1_full)
  SII_parms_r_c1 <- data.frame(Est = est, SE = se, SII = est, LCL = est - 1.96 * se, UCL = est + 1.96 * se)
  
  est <- sum(r_c2_full * coef(SII_model))
  se <- sqrt(t(r_c2_full) %*% vcov(SII_model) %*% r_c2_full)
  SII_parms_r_c2 <- data.frame(Est = est, SE = se, SII = est, LCL = est - 1.96 * se, UCL = est + 1.96 * se)
  
  est <- sum(r_c3_full * coef(SII_model))
  se <- sqrt(t(r_c3_full) %*% vcov(SII_model) %*% r_c3_full)
  SII_parms_r_c3 <- data.frame(Est = est, SE = se, SII = est, LCL = est - 1.96 * se, UCL = est + 1.96 * se)
  
  est <- sum(r_c4_full * coef(SII_model))
  se <- sqrt(t(r_c4_full) %*% vcov(SII_model) %*% r_c4_full)
  SII_parms_r_c4 <- data.frame(Est = est, SE = se, SII = est, LCL = est - 1.96 * se, UCL = est + 1.96 * se)
  
  est <- sum(r_c5_full * coef(SII_model))
  se <- sqrt(t(r_c5_full) %*% vcov(SII_model) %*% r_c5_full)
  SII_parms_r_c5 <- data.frame(Est = est, SE = se, SII = est, LCL = est - 1.96 * se, UCL = est + 1.96 * se)
  
  est <- sum(r_c6_full * coef(SII_model))
  se <- sqrt(t(r_c6_full) %*% vcov(SII_model) %*% r_c6_full)
  SII_parms_r_c6 <- data.frame(Est = est, SE = se, SII = est, LCL = est - 1.96 * se, UCL = est + 1.96 * se)
  
  SII_parms <- rbind(SII_parms_r_c1,SII_parms_r_c2,SII_parms_r_c3,SII_parms_r_c4,SII_parms_r_c5,SII_parms_r_c6)
  
  #Get RII and SII parameters
  RII_coefs <- as.numeric(RII_parms$Est)       
  RII_SEs <- as.numeric(RII_parms$SE)            
  RII_values <- as.numeric(RII_parms$RII)      
  RII_CIs <- RII_parms[4:5]  
  RII_trend_coef <- as.numeric(coef(RII_trend_model)["ridit_eqi_sqrt_sex:cycle_cont"])
  RII_trend_SE <- as.numeric(SE(RII_trend_model)["ridit_eqi_sqrt_sex:cycle_cont"])
  RII_trend_pvalue <- summary(RII_trend_model)$coefficients["ridit_eqi_sqrt_sex:cycle_cont", "Pr(>|t|)"]
  RII_trend_sex_coef <- as.numeric(coef(RII_trend_sex_model)["ridit_eqi_sqrt_sex:cycle_cont:sexM"])
  RII_trend_sex_SE <- as.numeric(SE(RII_trend_sex_model)["ridit_eqi_sqrt_sex:cycle_cont:sexM"])
  RII_trend_sex_pvalue <- summary(RII_trend_sex_model)$coefficients["ridit_eqi_sqrt_sex:cycle_cont:sexM", "Pr(>|t|)"]  
  
  SII_coefs <- as.numeric(SII_parms$Est)        
  SII_SEs <- as.numeric(SII_parms$SE)   
  SII_values <- as.numeric(SII_parms$SII)   
  SII_CIs <- SII_parms[4:5]
  SII_trend_coef <- as.numeric(coef(SII_trend_model)["ridit_eqi_sqrt_sex:cycle_cont"])
  SII_trend_SE <- as.numeric(SE(SII_trend_model)["ridit_eqi_sqrt_sex:cycle_cont"])
  SII_trend_pvalue <- summary(SII_trend_model)$coefficients["ridit_eqi_sqrt_sex:cycle_cont", "Pr(>|t|)"]
  SII_trend_sex_coef <- as.numeric(coef(SII_trend_sex_model)["ridit_eqi_sqrt_sex:cycle_cont:sexM"])
  SII_trend_sex_SE <- as.numeric(SE(SII_trend_sex_model)["ridit_eqi_sqrt_sex:cycle_cont:sexM"])
  SII_trend_sex_pvalue <- summary(SII_trend_sex_model)$coefficients["ridit_eqi_sqrt_sex:cycle_cont:sexM", "Pr(>|t|)"]      
  
  RII <- data.frame(
    Age_cat = age_cat_i,
    Sex = sex_i,
    Metric = 'RII',
    Cycle = c(1:6),
    Coef = RII_coefs,
    SE = RII_SEs,
    Value = RII_values,
    Lower_95_CI = as.numeric(RII_CIs$LCL),
    Upper_95_CI = as.numeric(RII_CIs$UCL),
    Trend_coef = RII_trend_coef,
    Trend_SE = RII_trend_SE,
    Trend_pvalue = RII_trend_pvalue,
    Trend_sex_coef = RII_trend_sex_coef,
    Trend_sex_SE = RII_trend_sex_SE,
    Trend_sex_pvalue = RII_trend_sex_pvalue
  )
  
  SII <- data.frame(
    Age_cat = age_cat_i,
    Sex = sex_i,
    Metric = 'SII',
    Cycle = c(1:6),
    Coef = SII_coefs,
    SE = SII_SEs,
    Value = SII_values,
    Lower_95_CI = as.numeric(SII_CIs$LCL),
    Upper_95_CI = as.numeric(SII_CIs$UCL),
    Trend_coef = SII_trend_coef,
    Trend_SE = SII_trend_SE,
    Trend_pvalue = SII_trend_pvalue,
    Trend_sex_coef = SII_trend_sex_coef,
    Trend_sex_SE = SII_trend_sex_SE,
    Trend_sex_pvalue = SII_trend_sex_pvalue
  )
  
  results_list[[i]] <- rbind(RII,SII)
  rm(RII_model,RII_trend_model,RII_trend_sex_model,SII_model,SII_trend_model,SII_trend_sex_model,
     r_c1,r_c2,r_c3,r_c4,r_c5,r_c6,r_c1_full,r_c2_full,r_c3_full,r_c4_full,r_c5_full,r_c6_full,est,se,
     RII_parms_r_c1,RII_parms_r_c2,RII_parms_r_c3,RII_parms_r_c4,RII_parms_r_c5,RII_parms_r_c6,RII_parms,
     SII_parms_r_c1,SII_parms_r_c2,SII_parms_r_c3,SII_parms_r_c4,SII_parms_r_c5,SII_parms_r_c6,SII_parms,
     RII_coefs,RII_SEs,RII_values,RII_CIs,RII_trend_coef,RII_trend_SE,RII_trend_pvalue,RII_trend_sex_coef,RII_trend_sex_SE,RII_trend_sex_pvalue,RII,
     SII_coefs,SII_SEs,SII_values,SII_CIs,SII_trend_coef,SII_trend_SE,SII_trend_pvalue,SII_trend_sex_coef,SII_trend_sex_SE,SII_trend_sex_pvalue,SII)
}

health_inequalities_age_sex_sqrt <- do.call(rbind, results_list)
rownames(health_inequalities_age_sex_sqrt) <- NULL
