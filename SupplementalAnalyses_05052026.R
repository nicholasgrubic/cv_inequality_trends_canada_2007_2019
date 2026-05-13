##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
# SUPPLEMENTAL ANALYSES        
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################

#Call packages
.libPaths("Y:/R4.5.1-incomplete-packages_v2") 
library(haven)
library(survey)
library(dplyr)
library(tidyr)
library(tableone)
library(writexl)
library(tibble)
library(emmeans)

#Import analytic data
obj1_analytic_data <- read_sas("Z:/VRDC-PROJ-11125/Grubic_11125/Data/CHMS Linked Source Files/obj1_cohort.sas7bdat")

#Factorize character variables
vars_to_factor <- c("cycle","immigrant","age_cat","sex","race_ethnic_CIHI","white","ei_quartile_moecd","ei_quartile_ooecd","ei_quartile_sqrt","inc_imp_flag","income_imputation")

obj1_analytic_data[vars_to_factor] <- lapply(obj1_analytic_data[vars_to_factor], function(x) factor(x, exclude=NA))
obj1_analytic_data$cycle_cont <- as.numeric(obj1_analytic_data$cycle)
str(obj1_analytic_data)
obj1_analytic_data <- obj1_analytic_data %>% rename(bmi_zscore=`_ZBFA`)

#Define the survey design using replicate weights
svy_design <- svrepdesign(data = obj1_analytic_data, 
                          degf = 68, #Recommended by Statistics Canada for national level analyses of cycles 1-6 
                          #(file:///S:/CHMS%20Cycle%206/CHMS_ECMS_C6_W6_v1/documentation/english/user_guide/CHMS_c1_c6_w2_comb_f1_T15.1_v1_EN.pdf)
                          weights = ~WGT_FULL, 
                          repweights = obj1_analytic_data[, grep("BSW", names(obj1_analytic_data))], 
                          type = "bootstrap") 

##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
# HEALTH INEQUALITIES (INCREMENTAL CONTRASTS BETWEEN ADJACENT CYCLES)    
##############################################################################################################################
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
  RII_model <- svyglm(r_poor_cv_health_1 ~ ridit_eqi_moecd*cycle + age + sex + white + immigrant,
                      design = svy_design,
                      family = quasipoisson(link="log"),
                      subset = age_cat == age_cat_i)
  SII_model <- svyglm(r_poor_cv_health_1 ~ ridit_eqi_moecd*cycle + age + sex + white + immigrant,
                      design = svy_design,
                      family = gaussian(link="identity"),
                      subset = age_cat == age_cat_i)

  RII_emtrends <- emtrends(object=RII_model, specs = ~cycle, var="ridit_eqi_moecd")
  SII_emtrends <- emtrends(object=SII_model, specs = ~cycle, var="ridit_eqi_moecd")
  
  # By default CIs are bonferroni-adjusted, although this adjustment can be removed with adjust=NULL
  # I think this is valid since we are not basing inference on all pairwise comparisons (which is done by default), but only incremental jumps between cycles
  RII_incs <- as.data.frame(pairs(RII_emtrends, reverse=TRUE, adjust=NULL, infer=c(TRUE,TRUE)))
  SII_incs <- as.data.frame(pairs(SII_emtrends, reverse=TRUE, adjust=NULL, infer=c(TRUE,TRUE)))
  
  RII_incs <- RII_incs %>% subset(contrast == "cycle2 - cycle1" | contrast == "cycle3 - cycle2" |  contrast == "cycle4 - cycle3" | contrast == "cycle5 - cycle4" |contrast == "cycle6 - cycle5")
  SII_incs <- SII_incs %>% subset(contrast == "cycle2 - cycle1" | contrast == "cycle3 - cycle2" |  contrast == "cycle4 - cycle3" | contrast == "cycle5 - cycle4" |contrast == "cycle6 - cycle5")
  
  #Get incremental RII and SII parameters
  RII_names <- RII_incs$contrast
  RII_inc_coefs <- as.numeric(RII_incs$estimate)  
  RII_inc_SEs <- as.numeric(RII_incs$SE)             
  RII_inc_CIs <- RII_incs[5:6]   
  RII_inc_critval <- as.numeric(RII_incs$z.ratio) 
  RII_inc_pvalue <- as.numeric(RII_incs$p.value) 

  SII_names <- SII_incs$contrast
  SII_inc_coefs <- as.numeric(SII_incs$estimate)  
  SII_inc_SEs <- as.numeric(SII_incs$SE)             
  SII_inc_CIs <- SII_incs[5:6]   
  SII_inc_critval <- as.numeric(SII_incs$t.ratio) 
  SII_inc_pvalue <- as.numeric(SII_incs$p.value) 
  
  RII <- data.frame(
    Age_cat = age_cat_i,
    Metric = 'RII',
    Inc_Contrast = RII_names,
    Coef = RII_inc_coefs,
    SE = RII_inc_SEs,
    Lower_95_CI = as.numeric(RII_inc_CIs$asymp.LCL),
    Upper_95_CI = as.numeric(RII_inc_CIs$asymp.UCL),
    Crit_val = RII_inc_critval,
    P_value = RII_inc_pvalue
  )
  
  SII <- data.frame(
    Age_cat = age_cat_i,
    Metric = 'SII',
    Inc_Contrast = SII_names,
    Coef = SII_inc_coefs,
    SE = SII_inc_SEs,
    Lower_95_CI = as.numeric(SII_inc_CIs$lower.CL),
    Upper_95_CI = as.numeric(SII_inc_CIs$upper.CL),
    Crit_val = SII_inc_critval,
    P_value = SII_inc_pvalue
  )
  
  results_list[[i]] <- rbind(RII,SII)
  rm(RII_model, SII_model,
     RII_emtrends, SII_emtrends, RII_incs, SII_incs,
     RII_names, RII_inc_coefs, RII_inc_SEs, RII_inc_CIs, RII_inc_critval, RII_inc_pvalue,
     SII_names, SII_inc_coefs, SII_inc_SEs, SII_inc_CIs, SII_inc_critval, SII_inc_pvalue,
     RII, SII)
}

inc_contrasts_age <- do.call(rbind, results_list)
rownames(inc_contrasts_age) <- NULL
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
  RII_model <- svyglm(r_poor_cv_health_1 ~ ridit_eqi_moecd_sex*cycle + age + white + immigrant,
                      design = svy_design,
                      family = quasipoisson(link="log"),
                      subset = age_cat == age_cat_i & sex == sex_i)
  RII_joint_wald_model <- svyglm(r_poor_cv_health_1 ~ ridit_eqi_moecd_sex*cycle*sex + age + white + immigrant,
                                design = svy_design,
                                family = quasipoisson(link="log"),
                                subset = age_cat == age_cat_i)
  SII_model <- svyglm(r_poor_cv_health_1 ~ ridit_eqi_moecd_sex*cycle + age + white + immigrant,
                      design = svy_design,
                      family = gaussian(link="identity"),
                      subset = age_cat == age_cat_i & sex == sex_i)
  SII_joint_wald_model <- svyglm(r_poor_cv_health_1 ~ ridit_eqi_moecd_sex*cycle*sex + age + white + immigrant,
                                design = svy_design,
                                family = gaussian(link="identity"),
                                subset = age_cat == age_cat_i)
  
  RII_emtrends <- emtrends(object=RII_model, specs = ~cycle, var="ridit_eqi_moecd_sex")
  SII_emtrends <- emtrends(object=SII_model, specs = ~cycle, var="ridit_eqi_moecd_sex")
  
  # By default CIs are bonferroni-adjusted, although this adjustment can be removed with adjust=NULL
  # I think this is valid since we are not basing inference on all pairwise comparisons(which is done by default), but only incremental jumps between cycles
  RII_incs <- as.data.frame(pairs(RII_emtrends, reverse=TRUE, adjust=NULL, infer=c(TRUE,TRUE)))
  SII_incs <- as.data.frame(pairs(SII_emtrends, reverse=TRUE, adjust=NULL, infer=c(TRUE,TRUE)))
  
  RII_incs <- RII_incs %>% subset(contrast == "cycle2 - cycle1" | contrast == "cycle3 - cycle2" |  contrast == "cycle4 - cycle3" | contrast == "cycle5 - cycle4" |contrast == "cycle6 - cycle5")
  SII_incs <- SII_incs %>% subset(contrast == "cycle2 - cycle1" | contrast == "cycle3 - cycle2" |  contrast == "cycle4 - cycle3" | contrast == "cycle5 - cycle4" |contrast == "cycle6 - cycle5")
  
  #Get incremental RII and SII parameters
  RII_names <- RII_incs$contrast
  RII_inc_coefs <- as.numeric(RII_incs$estimate)  
  RII_inc_SEs <- as.numeric(RII_incs$SE)             
  RII_inc_CIs <- RII_incs[5:6]   
  RII_inc_critval <- as.numeric(RII_incs$z.ratio) 
  RII_inc_pvalue <- as.numeric(RII_incs$p.value) 
  
  SII_names <- SII_incs$contrast
  SII_inc_coefs <- as.numeric(SII_incs$estimate)  
  SII_inc_SEs <- as.numeric(SII_incs$SE)             
  SII_inc_CIs <- SII_incs[5:6]   
  SII_inc_critval <- as.numeric(SII_incs$t.ratio) 
  SII_inc_pvalue <- as.numeric(SII_incs$p.value) 
  
  #Get joint Wald test p-value for all ridit*cycle*sex interaction terms (using categorical cycle so 5 interaction terms in total)
  jointWald_p_RII <- regTermTest(RII_joint_wald_model, ~ridit_eqi_moecd_sex*cycle*sex, method="Wald")$p
  jointWald_p_SII <- regTermTest(SII_joint_wald_model, ~ridit_eqi_moecd_sex*cycle*sex, method="Wald")$p
  
  RII <- data.frame(
    Age_cat = age_cat_i,
    Sex = sex_i,
    Metric = 'RII',
    Inc_Contrast = RII_names,
    Coef = RII_inc_coefs,
    SE = RII_inc_SEs,
    Lower_95_CI = as.numeric(RII_inc_CIs$asymp.LCL),
    Upper_95_CI = as.numeric(RII_inc_CIs$asymp.UCL),
    Crit_val = RII_inc_critval,
    P_value = RII_inc_pvalue,
    Joint_wald_p_value = jointWald_p_RII
  )
  
  SII <- data.frame(
    Age_cat = age_cat_i,
    Sex = sex_i,
    Metric = 'SII',
    Inc_Contrast = SII_names,
    Coef = SII_inc_coefs,
    SE = SII_inc_SEs,
    Lower_95_CI = as.numeric(SII_inc_CIs$lower.CL),
    Upper_95_CI = as.numeric(SII_inc_CIs$upper.CL),
    Crit_val = SII_inc_critval,
    P_value = SII_inc_pvalue,
    Joint_wald_p_value = jointWald_p_SII
  )
  
  results_list[[i]] <- rbind(RII,SII)
  rm(RII_model, SII_model, RII_joint_wald_model, SII_joint_wald_model,
     RII_emtrends, SII_emtrends, RII_incs, SII_incs,
     RII_names, RII_inc_coefs, RII_inc_SEs, RII_inc_CIs, RII_inc_critval, RII_inc_pvalue,
     SII_names, SII_inc_coefs, SII_inc_SEs, SII_inc_CIs, SII_inc_critval, SII_inc_pvalue,
     RII, SII, jointWald_p_RII, jointWald_p_SII)     
}

inc_contrasts_age_sex <- do.call(rbind, results_list)
rownames(inc_contrasts_age_sex) <- NULL
rm(results_list)