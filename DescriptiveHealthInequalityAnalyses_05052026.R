##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
# MAIN ANALYSES               
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################

#Call packages
.libPaths("//tor-main/Utilities/R/Packages/4.0_top")
library(haven)
library(survey)
library(dplyr)
library(tidyr)
library(tableone)
library(writexl)
library(tibble)

#Import analytic data
obj1_analytic_data <- read_sas("J:/Grubic_11125/Data/CHMS Linked Source Files/obj1_cohort.sas7bdat")

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
# DESCRIPTIVE CHARACTERISTICS        
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################

##############################################################################################################################
# Descriptive characteristics of all respondents (total across survey sample)
##############################################################################################################################
descriptive_chars_total <- svyby(~ age + age_cat + sex + white + race_ethnic_CIHI + immigrant + 
                                   ei_quartile_moecd + eq_moecd_hhi + ei_quartile_ooecd + eq_ooecd_hhi + ei_quartile_sqrt + eq_sqrt_hhi +
                                   r_poor_cv_health_1 + obesity + hbp + dyslipidemia_abnormal + diabetes,
                                   by = ~ pregnant,
                                   design = svy_design,
                                   vartype = c("ci", "cvpct"),
                                   FUN = svymean,
                                   NA.rm = TRUE) #exclude missing

#Sociodemographic characteristics
svy_stats_age <- descriptive_chars_total %>% 
  mutate(age = paste0(round(`age`,1)," (",round(`ci_l.age`,1),"-",round(`ci_u.age`,1),"), CV%= ",round(`cv%.age`,1))) %>%
  select(pregnant,age) %>%
  pivot_longer(cols = -pregnant, names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = pregnant, values_from = value)

svy_stats_age_cat <- descriptive_chars_total %>% 
  mutate(age_cat_6_17_years = paste0(round(`age_cat6-17 years`*100,1)," (",round(`ci_l.age_cat6-17 years`*100,1),"-",round(`ci_u.age_cat6-17 years`*100,1),"), CV%= ",round(`cv%.age_cat6-17 years`,2))) %>%
  mutate(age_cat_18_39_years = paste0(round(`age_cat18-39 years`*100,1)," (",round(`ci_l.age_cat18-39 years`*100,1),"-",round(`ci_u.age_cat18-39 years`*100,1),"), CV%= ",round(`cv%.age_cat18-39 years`,2))) %>%
  mutate(age_cat_40_64_years = paste0(round(`age_cat40-64 years`*100,1)," (",round(`ci_l.age_cat40-64 years`*100,1),"-",round(`ci_u.age_cat40-64 years`*100,1),"), CV%= ",round(`cv%.age_cat40-64 years`,2))) %>%
  select(pregnant,age_cat_6_17_years,age_cat_18_39_years,age_cat_40_64_years) %>%
  pivot_longer(cols = -pregnant, names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = pregnant, values_from = value)

svy_stats_sex <- descriptive_chars_total %>% 
  mutate(sexF = paste0(round(`sexF`*100,1)," (",round(`ci_l.sexF`*100,1),"-",round(`ci_u.sexF`*100,1),"), CV%= ",round(`cv%.sexF`,2))) %>%
  mutate(sexM = paste0(round(`sexM`*100,1)," (",round(`ci_l.sexM`*100,1),"-",round(`ci_u.sexM`*100,1),"), CV%= ",round(`cv%.sexM`,2))) %>%
  select(pregnant,sexF,sexM) %>%
  pivot_longer(cols = -pregnant, names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = pregnant, values_from = value)

svy_stats_white <- descriptive_chars_total %>% 
  mutate(white1 = paste0(round(`white1`*100,1)," (",round(`ci_l.white1`*100,1),"-",round(`ci_u.white1`*100,1),"), CV%= ",round(`cv%.white1`,2))) %>%
  mutate(white0 = paste0(round(`white0`*100,1)," (",round(`ci_l.white0`*100,1),"-",round(`ci_u.white0`*100,1),"), CV%= ",round(`cv%.white0`,2))) %>%
  select(pregnant,white1,white0) %>%
  pivot_longer(cols = -pregnant, names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = pregnant, values_from = value)

svy_stats_race_ethnic <- descriptive_chars_total %>% 
  mutate(race_ethnic_CIHIWhite = paste0(round(`race_ethnic_CIHIWhite`*100,1)," (",round(`ci_l.race_ethnic_CIHIWhite`*100,1),"-",round(`ci_u.race_ethnic_CIHIWhite`*100,1),"), CV%= ",round(`cv%.race_ethnic_CIHIWhite`,2))) %>%
  mutate(race_ethnic_CIHIBlack = paste0(round(`race_ethnic_CIHIBlack`*100,1)," (",round(`ci_l.race_ethnic_CIHIBlack`*100,1),"-",round(`ci_u.race_ethnic_CIHIBlack`*100,1),"), CV%= ",round(`cv%.race_ethnic_CIHIBlack`,2))) %>%
  mutate(race_ethnic_CIHIEastAsian = paste0(round(`race_ethnic_CIHIEast Asian`*100,1)," (",round(`ci_l.race_ethnic_CIHIEast Asian`*100,1),"-",round(`ci_u.race_ethnic_CIHIEast Asian`*100,1),"), CV%= ",round(`cv%.race_ethnic_CIHIEast Asian`,2))) %>%
  mutate(race_ethnic_CIHIIndigenous = paste0(round(`race_ethnic_CIHIIndigenous`*100,1)," (",round(`ci_l.race_ethnic_CIHIIndigenous`*100,1),"-",round(`ci_u.race_ethnic_CIHIIndigenous`*100,1),"), CV%= ",round(`cv%.race_ethnic_CIHIIndigenous`,2))) %>%
  mutate(race_ethnic_CIHILatinAmerican = paste0(round(`race_ethnic_CIHILatin American`*100,1)," (",round(`ci_l.race_ethnic_CIHILatin American`*100,1),"-",round(`ci_u.race_ethnic_CIHILatin American`*100,1),"), CV%= ",round(`cv%.race_ethnic_CIHILatin American`,2))) %>%
  mutate(race_ethnic_CIHIMiddleEastern = paste0(round(`race_ethnic_CIHIMiddle Eastern`*100,1)," (",round(`ci_l.race_ethnic_CIHIMiddle Eastern`*100,1),"-",round(`ci_u.race_ethnic_CIHIMiddle Eastern`*100,1),"), CV%= ",round(`cv%.race_ethnic_CIHIMiddle Eastern`,2))) %>%
  mutate(race_ethnic_CIHISouthAsian = paste0(round(`race_ethnic_CIHISouth Asian`*100,1)," (",round(`ci_l.race_ethnic_CIHISouth Asian`*100,1),"-",round(`ci_u.race_ethnic_CIHISouth Asian`*100,1),"), CV%= ",round(`cv%.race_ethnic_CIHISouth Asian`,2))) %>%
  mutate(race_ethnic_CIHISoutheastAsian = paste0(round(`race_ethnic_CIHISoutheast Asian`*100,1)," (",round(`ci_l.race_ethnic_CIHISoutheast Asian`*100,1),"-",round(`ci_u.race_ethnic_CIHISoutheast Asian`*100,1),"), CV%= ",round(`cv%.race_ethnic_CIHISoutheast Asian`,2))) %>%
  mutate(race_ethnic_CIHIOther = paste0(round(`race_ethnic_CIHIOther`*100,1)," (",round(`ci_l.race_ethnic_CIHIOther`*100,1),"-",round(`ci_u.race_ethnic_CIHIOther`*100,1),"), CV%= ",round(`cv%.race_ethnic_CIHIOther`,2))) %>%
  select(pregnant,race_ethnic_CIHIWhite,race_ethnic_CIHIBlack,race_ethnic_CIHIEastAsian,race_ethnic_CIHIIndigenous,race_ethnic_CIHILatinAmerican,race_ethnic_CIHIMiddleEastern,race_ethnic_CIHISouthAsian,race_ethnic_CIHISoutheastAsian,race_ethnic_CIHIOther) %>%
  pivot_longer(cols = -pregnant, names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = pregnant, values_from = value)

svy_stats_immigrant <- descriptive_chars_total %>% 
  mutate(immigrant1 = paste0(round(`immigrant1`*100,1)," (",round(`ci_l.immigrant1`*100,1),"-",round(`ci_u.immigrant1`*100,1),"), CV%= ",round(`cv%.immigrant1`,2))) %>%
  mutate(immigrant0 = paste0(round(`immigrant0`*100,1)," (",round(`ci_l.immigrant0`*100,1),"-",round(`ci_u.immigrant0`*100,1),"), CV%= ",round(`cv%.immigrant0`,2))) %>%
  select(pregnant,immigrant1,immigrant0) %>%
  pivot_longer(cols = -pregnant, names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = pregnant, values_from = value)

svy_stats_ei_quartile_moecd <- descriptive_chars_total %>% 
  mutate(ei_quartile_moecd1 = paste0(round(`ei_quartile_moecd1`*100,1)," (",round(`ci_l.ei_quartile_moecd1`*100,1),"-",round(`ci_u.ei_quartile_moecd1`*100,1),"), CV%= ",round(`cv%.ei_quartile_moecd1`,2))) %>%
  mutate(ei_quartile_moecd2 = paste0(round(`ei_quartile_moecd2`*100,1)," (",round(`ci_l.ei_quartile_moecd2`*100,1),"-",round(`ci_u.ei_quartile_moecd2`*100,1),"), CV%= ",round(`cv%.ei_quartile_moecd2`,2))) %>%
  mutate(ei_quartile_moecd3 = paste0(round(`ei_quartile_moecd3`*100,1)," (",round(`ci_l.ei_quartile_moecd3`*100,1),"-",round(`ci_u.ei_quartile_moecd3`*100,1),"), CV%= ",round(`cv%.ei_quartile_moecd3`,2))) %>%
  mutate(ei_quartile_moecd4 = paste0(round(`ei_quartile_moecd4`*100,1)," (",round(`ci_l.ei_quartile_moecd4`*100,1),"-",round(`ci_u.ei_quartile_moecd4`*100,1),"), CV%= ",round(`cv%.ei_quartile_moecd4`,2))) %>%
  select(pregnant,ei_quartile_moecd1,ei_quartile_moecd2,ei_quartile_moecd3,ei_quartile_moecd4) %>%
  pivot_longer(cols = -pregnant, names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = pregnant, values_from = value)

svy_stats_eq_moecd_hhi <- descriptive_chars_total %>% 
  mutate(eq_moecd_hhi = paste0(round(`eq_moecd_hhi`,1)," (",round(`ci_l.eq_moecd_hhi`,1),"-",round(`ci_u.eq_moecd_hhi`,1),"), CV%= ",round(`cv%.eq_moecd_hhi`,1))) %>%
  select(pregnant,eq_moecd_hhi) %>%
  pivot_longer(cols = -pregnant, names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = pregnant, values_from = value)

svy_stats_ei_quartile_ooecd <- descriptive_chars_total %>% 
  mutate(ei_quartile_ooecd1 = paste0(round(`ei_quartile_ooecd1`*100,1)," (",round(`ci_l.ei_quartile_ooecd1`*100,1),"-",round(`ci_u.ei_quartile_ooecd1`*100,1),"), CV%= ",round(`cv%.ei_quartile_ooecd1`,2))) %>%
  mutate(ei_quartile_ooecd2 = paste0(round(`ei_quartile_ooecd2`*100,1)," (",round(`ci_l.ei_quartile_ooecd2`*100,1),"-",round(`ci_u.ei_quartile_ooecd2`*100,1),"), CV%= ",round(`cv%.ei_quartile_ooecd2`,2))) %>%
  mutate(ei_quartile_ooecd3 = paste0(round(`ei_quartile_ooecd3`*100,1)," (",round(`ci_l.ei_quartile_ooecd3`*100,1),"-",round(`ci_u.ei_quartile_ooecd3`*100,1),"), CV%= ",round(`cv%.ei_quartile_ooecd3`,2))) %>%
  mutate(ei_quartile_ooecd4 = paste0(round(`ei_quartile_ooecd4`*100,1)," (",round(`ci_l.ei_quartile_ooecd4`*100,1),"-",round(`ci_u.ei_quartile_ooecd4`*100,1),"), CV%= ",round(`cv%.ei_quartile_ooecd4`,2))) %>%
  select(pregnant,ei_quartile_ooecd1,ei_quartile_ooecd2,ei_quartile_ooecd3,ei_quartile_ooecd4) %>%
  pivot_longer(cols = -pregnant, names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = pregnant, values_from = value)

svy_stats_eq_ooecd_hhi <- descriptive_chars_total %>% 
  mutate(eq_ooecd_hhi = paste0(round(`eq_ooecd_hhi`,1)," (",round(`ci_l.eq_ooecd_hhi`,1),"-",round(`ci_u.eq_ooecd_hhi`,1),"), CV%= ",round(`cv%.eq_ooecd_hhi`,1))) %>%
  select(pregnant,eq_ooecd_hhi) %>%
  pivot_longer(cols = -pregnant, names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = pregnant, values_from = value)

svy_stats_ei_quartile_sqrt <- descriptive_chars_total %>% 
  mutate(ei_quartile_sqrt1 = paste0(round(`ei_quartile_sqrt1`*100,1)," (",round(`ci_l.ei_quartile_sqrt1`*100,1),"-",round(`ci_u.ei_quartile_sqrt1`*100,1),"), CV%= ",round(`cv%.ei_quartile_sqrt1`,2))) %>%
  mutate(ei_quartile_sqrt2 = paste0(round(`ei_quartile_sqrt2`*100,1)," (",round(`ci_l.ei_quartile_sqrt2`*100,1),"-",round(`ci_u.ei_quartile_sqrt2`*100,1),"), CV%= ",round(`cv%.ei_quartile_sqrt2`,2))) %>%
  mutate(ei_quartile_sqrt3 = paste0(round(`ei_quartile_sqrt3`*100,1)," (",round(`ci_l.ei_quartile_sqrt3`*100,1),"-",round(`ci_u.ei_quartile_sqrt3`*100,1),"), CV%= ",round(`cv%.ei_quartile_sqrt3`,2))) %>%
  mutate(ei_quartile_sqrt4 = paste0(round(`ei_quartile_sqrt4`*100,1)," (",round(`ci_l.ei_quartile_sqrt4`*100,1),"-",round(`ci_u.ei_quartile_sqrt4`*100,1),"), CV%= ",round(`cv%.ei_quartile_sqrt4`,2))) %>%
  select(pregnant,ei_quartile_sqrt1,ei_quartile_sqrt2,ei_quartile_sqrt3,ei_quartile_sqrt4) %>%
  pivot_longer(cols = -pregnant, names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = pregnant, values_from = value)

svy_stats_eq_sqrt_hhi <- descriptive_chars_total %>% 
  mutate(eq_sqrt_hhi = paste0(round(`eq_sqrt_hhi`,1)," (",round(`ci_l.eq_sqrt_hhi`,1),"-",round(`ci_u.eq_sqrt_hhi`,1),"), CV%= ",round(`cv%.eq_sqrt_hhi`,1))) %>%
  select(pregnant,eq_sqrt_hhi) %>%
  pivot_longer(cols = -pregnant, names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = pregnant, values_from = value)

#Outcomes
svy_stats_r_poor_cv_health <- descriptive_chars_total %>% 
  mutate(r_poor_cv_health_1 = paste0(round(`r_poor_cv_health_1`*100,1)," (",round(`ci_l.r_poor_cv_health_1`*100,1),"-",round(`ci_u.r_poor_cv_health_1`*100,1),"), CV%= ",round(`cv%.r_poor_cv_health_1`,2))) %>%
  select(pregnant,r_poor_cv_health_1) %>%
  pivot_longer(cols = -pregnant, names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = pregnant, values_from = value)

svy_stats_obesity <- descriptive_chars_total %>% 
  mutate(obesity = paste0(round(`obesity`*100,1)," (",round(`ci_l.obesity`*100,1),"-",round(`ci_u.obesity`*100,1),"), CV%= ",round(`cv%.obesity`,2))) %>%
  select(pregnant,obesity) %>%
  pivot_longer(cols = -pregnant, names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = pregnant, values_from = value)

svy_stats_hbp <- descriptive_chars_total %>% 
  mutate(hbp = paste0(round(`hbp`*100,1)," (",round(`ci_l.hbp`*100,1),"-",round(`ci_u.hbp`*100,1),"), CV%= ",round(`cv%.hbp`,2))) %>%
  select(pregnant,hbp) %>%
  pivot_longer(cols = -pregnant, names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = pregnant, values_from = value)

svy_stats_dyslipidemia_abnormal <- descriptive_chars_total %>% 
  mutate(dyslipidemia_abnormal = paste0(round(`dyslipidemia_abnormal`*100,1)," (",round(`ci_l.dyslipidemia_abnormal`*100,1),"-",round(`ci_u.dyslipidemia_abnormal`*100,1),"), CV%= ",round(`cv%.dyslipidemia_abnormal`,2))) %>%
  select(pregnant,dyslipidemia_abnormal) %>%
  pivot_longer(cols = -pregnant, names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = pregnant, values_from = value)

svy_stats_diabetes <- descriptive_chars_total %>% 
  mutate(diabetes = paste0(round(`diabetes`*100,1)," (",round(`ci_l.diabetes`*100,1),"-",round(`ci_u.diabetes`*100,1),"), CV%= ",round(`cv%.diabetes`,2))) %>%
  select(pregnant,diabetes) %>%
  pivot_longer(cols = -pregnant, names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = pregnant, values_from = value)

descriptive_chars_total_table <- rbind(svy_stats_age,svy_stats_age_cat,svy_stats_sex,svy_stats_white,svy_stats_race_ethnic,svy_stats_immigrant,
                                 svy_stats_ei_quartile_moecd,svy_stats_eq_moecd_hhi,svy_stats_ei_quartile_ooecd,svy_stats_eq_ooecd_hhi,svy_stats_ei_quartile_sqrt,svy_stats_eq_sqrt_hhi,
                                 svy_stats_r_poor_cv_health,svy_stats_obesity,svy_stats_hbp,svy_stats_dyslipidemia_abnormal,svy_stats_diabetes)

remove(svy_stats_age,svy_stats_age_cat,svy_stats_sex,svy_stats_white,svy_stats_race_ethnic,svy_stats_immigrant,
       svy_stats_ei_quartile_moecd,svy_stats_eq_moecd_hhi,svy_stats_ei_quartile_ooecd,svy_stats_eq_ooecd_hhi,svy_stats_ei_quartile_sqrt,svy_stats_eq_sqrt_hhi,
       svy_stats_r_poor_cv_health,svy_stats_obesity,svy_stats_hbp,svy_stats_dyslipidemia_abnormal,svy_stats_diabetes)

##############################################################################################################################
# Descriptive characteristics of all respondents by survey cycle
##############################################################################################################################
descriptive_chars <- svyby(~ age + age_cat + sex + white + race_ethnic_CIHI + immigrant + 
                           ei_quartile_moecd + eq_moecd_hhi + ei_quartile_ooecd + eq_ooecd_hhi + ei_quartile_sqrt + eq_sqrt_hhi +
                           r_poor_cv_health_1 + obesity + hbp + dyslipidemia_abnormal + diabetes,
                           by= ~cycle,
                           design = svy_design,
                           vartype = c("ci", "cvpct"),
                           FUN = svymean,
                           NA.rm = TRUE) #exclude missing

#Sociodemographic characteristics
svy_stats_age <- descriptive_chars %>% 
  mutate(age = paste0(round(`age`,1)," (",round(`ci_l.age`,1),"-",round(`ci_u.age`,1),"), CV%= ",round(`cv%.age`,1))) %>%
  select(cycle,age) %>%
  pivot_longer(cols = -cycle, names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = cycle, values_from = value)

svy_stats_age_cat <- descriptive_chars %>% 
  mutate(age_cat_6_17_years = paste0(round(`age_cat6-17 years`*100,1)," (",round(`ci_l.age_cat6-17 years`*100,1),"-",round(`ci_u.age_cat6-17 years`*100,1),"), CV%= ",round(`cv%.age_cat6-17 years`,2))) %>%
  mutate(age_cat_18_39_years = paste0(round(`age_cat18-39 years`*100,1)," (",round(`ci_l.age_cat18-39 years`*100,1),"-",round(`ci_u.age_cat18-39 years`*100,1),"), CV%= ",round(`cv%.age_cat18-39 years`,2))) %>%
  mutate(age_cat_40_64_years = paste0(round(`age_cat40-64 years`*100,1)," (",round(`ci_l.age_cat40-64 years`*100,1),"-",round(`ci_u.age_cat40-64 years`*100,1),"), CV%= ",round(`cv%.age_cat40-64 years`,2))) %>%
  select(cycle,age_cat_6_17_years,age_cat_18_39_years,age_cat_40_64_years) %>%
  pivot_longer(cols = -cycle, names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = cycle, values_from = value)

svy_stats_sex <- descriptive_chars %>% 
  mutate(sexF = paste0(round(`sexF`*100,1)," (",round(`ci_l.sexF`*100,1),"-",round(`ci_u.sexF`*100,1),"), CV%= ",round(`cv%.sexF`,2))) %>%
  mutate(sexM = paste0(round(`sexM`*100,1)," (",round(`ci_l.sexM`*100,1),"-",round(`ci_u.sexM`*100,1),"), CV%= ",round(`cv%.sexM`,2))) %>%
  select(cycle,sexF,sexM) %>%
  pivot_longer(cols = -cycle, names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = cycle, values_from = value)

svy_stats_white <- descriptive_chars %>% 
  mutate(white1 = paste0(round(`white1`*100,1)," (",round(`ci_l.white1`*100,1),"-",round(`ci_u.white1`*100,1),"), CV%= ",round(`cv%.white1`,2))) %>%
  mutate(white0 = paste0(round(`white0`*100,1)," (",round(`ci_l.white0`*100,1),"-",round(`ci_u.white0`*100,1),"), CV%= ",round(`cv%.white0`,2))) %>%
  select(cycle,white1,white0) %>%
  pivot_longer(cols = -cycle, names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = cycle, values_from = value)

svy_stats_race_ethnic <- descriptive_chars %>% 
  mutate(race_ethnic_CIHIWhite = paste0(round(`race_ethnic_CIHIWhite`*100,1)," (",round(`ci_l.race_ethnic_CIHIWhite`*100,1),"-",round(`ci_u.race_ethnic_CIHIWhite`*100,1),"), CV%= ",round(`cv%.race_ethnic_CIHIWhite`,2))) %>%
  mutate(race_ethnic_CIHIBlack = paste0(round(`race_ethnic_CIHIBlack`*100,1)," (",round(`ci_l.race_ethnic_CIHIBlack`*100,1),"-",round(`ci_u.race_ethnic_CIHIBlack`*100,1),"), CV%= ",round(`cv%.race_ethnic_CIHIBlack`,2))) %>%
  mutate(race_ethnic_CIHIEastAsian = paste0(round(`race_ethnic_CIHIEast Asian`*100,1)," (",round(`ci_l.race_ethnic_CIHIEast Asian`*100,1),"-",round(`ci_u.race_ethnic_CIHIEast Asian`*100,1),"), CV%= ",round(`cv%.race_ethnic_CIHIEast Asian`,2))) %>%
  mutate(race_ethnic_CIHIIndigenous = paste0(round(`race_ethnic_CIHIIndigenous`*100,1)," (",round(`ci_l.race_ethnic_CIHIIndigenous`*100,1),"-",round(`ci_u.race_ethnic_CIHIIndigenous`*100,1),"), CV%= ",round(`cv%.race_ethnic_CIHIIndigenous`,2))) %>%
  mutate(race_ethnic_CIHILatinAmerican = paste0(round(`race_ethnic_CIHILatin American`*100,1)," (",round(`ci_l.race_ethnic_CIHILatin American`*100,1),"-",round(`ci_u.race_ethnic_CIHILatin American`*100,1),"), CV%= ",round(`cv%.race_ethnic_CIHILatin American`,2))) %>%
  mutate(race_ethnic_CIHIMiddleEastern = paste0(round(`race_ethnic_CIHIMiddle Eastern`*100,1)," (",round(`ci_l.race_ethnic_CIHIMiddle Eastern`*100,1),"-",round(`ci_u.race_ethnic_CIHIMiddle Eastern`*100,1),"), CV%= ",round(`cv%.race_ethnic_CIHIMiddle Eastern`,2))) %>%
  mutate(race_ethnic_CIHISouthAsian = paste0(round(`race_ethnic_CIHISouth Asian`*100,1)," (",round(`ci_l.race_ethnic_CIHISouth Asian`*100,1),"-",round(`ci_u.race_ethnic_CIHISouth Asian`*100,1),"), CV%= ",round(`cv%.race_ethnic_CIHISouth Asian`,2))) %>%
  mutate(race_ethnic_CIHISoutheastAsian = paste0(round(`race_ethnic_CIHISoutheast Asian`*100,1)," (",round(`ci_l.race_ethnic_CIHISoutheast Asian`*100,1),"-",round(`ci_u.race_ethnic_CIHISoutheast Asian`*100,1),"), CV%= ",round(`cv%.race_ethnic_CIHISoutheast Asian`,2))) %>%
  mutate(race_ethnic_CIHIOther = paste0(round(`race_ethnic_CIHIOther`*100,1)," (",round(`ci_l.race_ethnic_CIHIOther`*100,1),"-",round(`ci_u.race_ethnic_CIHIOther`*100,1),"), CV%= ",round(`cv%.race_ethnic_CIHIOther`,2))) %>%
  select(cycle,race_ethnic_CIHIWhite,race_ethnic_CIHIBlack,race_ethnic_CIHIEastAsian,race_ethnic_CIHIIndigenous,race_ethnic_CIHILatinAmerican,race_ethnic_CIHIMiddleEastern,race_ethnic_CIHISouthAsian,race_ethnic_CIHISoutheastAsian,race_ethnic_CIHIOther) %>%
  pivot_longer(cols = -cycle, names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = cycle, values_from = value)

svy_stats_immigrant <- descriptive_chars %>% 
  mutate(immigrant1 = paste0(round(`immigrant1`*100,1)," (",round(`ci_l.immigrant1`*100,1),"-",round(`ci_u.immigrant1`*100,1),"), CV%= ",round(`cv%.immigrant1`,2))) %>%
  mutate(immigrant0 = paste0(round(`immigrant0`*100,1)," (",round(`ci_l.immigrant0`*100,1),"-",round(`ci_u.immigrant0`*100,1),"), CV%= ",round(`cv%.immigrant0`,2))) %>%
  select(cycle,immigrant1,immigrant0) %>%
  pivot_longer(cols = -cycle, names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = cycle, values_from = value)

svy_stats_ei_quartile_moecd <- descriptive_chars %>% 
  mutate(ei_quartile_moecd1 = paste0(round(`ei_quartile_moecd1`*100,1)," (",round(`ci_l.ei_quartile_moecd1`*100,1),"-",round(`ci_u.ei_quartile_moecd1`*100,1),"), CV%= ",round(`cv%.ei_quartile_moecd1`,2))) %>%
  mutate(ei_quartile_moecd2 = paste0(round(`ei_quartile_moecd2`*100,1)," (",round(`ci_l.ei_quartile_moecd2`*100,1),"-",round(`ci_u.ei_quartile_moecd2`*100,1),"), CV%= ",round(`cv%.ei_quartile_moecd2`,2))) %>%
  mutate(ei_quartile_moecd3 = paste0(round(`ei_quartile_moecd3`*100,1)," (",round(`ci_l.ei_quartile_moecd3`*100,1),"-",round(`ci_u.ei_quartile_moecd3`*100,1),"), CV%= ",round(`cv%.ei_quartile_moecd3`,2))) %>%
  mutate(ei_quartile_moecd4 = paste0(round(`ei_quartile_moecd4`*100,1)," (",round(`ci_l.ei_quartile_moecd4`*100,1),"-",round(`ci_u.ei_quartile_moecd4`*100,1),"), CV%= ",round(`cv%.ei_quartile_moecd4`,2))) %>%
  select(cycle,ei_quartile_moecd1,ei_quartile_moecd2,ei_quartile_moecd3,ei_quartile_moecd4) %>%
  pivot_longer(cols = -cycle, names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = cycle, values_from = value)

svy_stats_eq_moecd_hhi <- descriptive_chars %>% 
  mutate(eq_moecd_hhi = paste0(round(`eq_moecd_hhi`,1)," (",round(`ci_l.eq_moecd_hhi`,1),"-",round(`ci_u.eq_moecd_hhi`,1),"), CV%= ",round(`cv%.eq_moecd_hhi`,1))) %>%
  select(cycle,eq_moecd_hhi) %>%
  pivot_longer(cols = -cycle, names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = cycle, values_from = value)

svy_stats_ei_quartile_ooecd <- descriptive_chars %>% 
  mutate(ei_quartile_ooecd1 = paste0(round(`ei_quartile_ooecd1`*100,1)," (",round(`ci_l.ei_quartile_ooecd1`*100,1),"-",round(`ci_u.ei_quartile_ooecd1`*100,1),"), CV%= ",round(`cv%.ei_quartile_ooecd1`,2))) %>%
  mutate(ei_quartile_ooecd2 = paste0(round(`ei_quartile_ooecd2`*100,1)," (",round(`ci_l.ei_quartile_ooecd2`*100,1),"-",round(`ci_u.ei_quartile_ooecd2`*100,1),"), CV%= ",round(`cv%.ei_quartile_ooecd2`,2))) %>%
  mutate(ei_quartile_ooecd3 = paste0(round(`ei_quartile_ooecd3`*100,1)," (",round(`ci_l.ei_quartile_ooecd3`*100,1),"-",round(`ci_u.ei_quartile_ooecd3`*100,1),"), CV%= ",round(`cv%.ei_quartile_ooecd3`,2))) %>%
  mutate(ei_quartile_ooecd4 = paste0(round(`ei_quartile_ooecd4`*100,1)," (",round(`ci_l.ei_quartile_ooecd4`*100,1),"-",round(`ci_u.ei_quartile_ooecd4`*100,1),"), CV%= ",round(`cv%.ei_quartile_ooecd4`,2))) %>%
  select(cycle,ei_quartile_ooecd1,ei_quartile_ooecd2,ei_quartile_ooecd3,ei_quartile_ooecd4) %>%
  pivot_longer(cols = -cycle, names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = cycle, values_from = value)

svy_stats_eq_ooecd_hhi <- descriptive_chars %>% 
  mutate(eq_ooecd_hhi = paste0(round(`eq_ooecd_hhi`,1)," (",round(`ci_l.eq_ooecd_hhi`,1),"-",round(`ci_u.eq_ooecd_hhi`,1),"), CV%= ",round(`cv%.eq_ooecd_hhi`,1))) %>%
  select(cycle,eq_ooecd_hhi) %>%
  pivot_longer(cols = -cycle, names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = cycle, values_from = value)

svy_stats_ei_quartile_sqrt <- descriptive_chars %>% 
  mutate(ei_quartile_sqrt1 = paste0(round(`ei_quartile_sqrt1`*100,1)," (",round(`ci_l.ei_quartile_sqrt1`*100,1),"-",round(`ci_u.ei_quartile_sqrt1`*100,1),"), CV%= ",round(`cv%.ei_quartile_sqrt1`,2))) %>%
  mutate(ei_quartile_sqrt2 = paste0(round(`ei_quartile_sqrt2`*100,1)," (",round(`ci_l.ei_quartile_sqrt2`*100,1),"-",round(`ci_u.ei_quartile_sqrt2`*100,1),"), CV%= ",round(`cv%.ei_quartile_sqrt2`,2))) %>%
  mutate(ei_quartile_sqrt3 = paste0(round(`ei_quartile_sqrt3`*100,1)," (",round(`ci_l.ei_quartile_sqrt3`*100,1),"-",round(`ci_u.ei_quartile_sqrt3`*100,1),"), CV%= ",round(`cv%.ei_quartile_sqrt3`,2))) %>%
  mutate(ei_quartile_sqrt4 = paste0(round(`ei_quartile_sqrt4`*100,1)," (",round(`ci_l.ei_quartile_sqrt4`*100,1),"-",round(`ci_u.ei_quartile_sqrt4`*100,1),"), CV%= ",round(`cv%.ei_quartile_sqrt4`,2))) %>%
  select(cycle,ei_quartile_sqrt1,ei_quartile_sqrt2,ei_quartile_sqrt3,ei_quartile_sqrt4) %>%
  pivot_longer(cols = -cycle, names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = cycle, values_from = value)

svy_stats_eq_sqrt_hhi <- descriptive_chars %>% 
  mutate(eq_sqrt_hhi = paste0(round(`eq_sqrt_hhi`,1)," (",round(`ci_l.eq_sqrt_hhi`,1),"-",round(`ci_u.eq_sqrt_hhi`,1),"), CV%= ",round(`cv%.eq_sqrt_hhi`,1))) %>%
  select(cycle,eq_sqrt_hhi) %>%
  pivot_longer(cols = -cycle, names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = cycle, values_from = value)

#Outcomes
svy_stats_r_poor_cv_health <- descriptive_chars %>% 
  mutate(r_poor_cv_health_1 = paste0(round(`r_poor_cv_health_1`*100,1)," (",round(`ci_l.r_poor_cv_health_1`*100,1),"-",round(`ci_u.r_poor_cv_health_1`*100,1),"), CV%= ",round(`cv%.r_poor_cv_health_1`,2))) %>%
  select(cycle,r_poor_cv_health_1) %>%
  pivot_longer(cols = -cycle, names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = cycle, values_from = value)

svy_stats_obesity <- descriptive_chars %>% 
  mutate(obesity = paste0(round(`obesity`*100,1)," (",round(`ci_l.obesity`*100,1),"-",round(`ci_u.obesity`*100,1),"), CV%= ",round(`cv%.obesity`,2))) %>%
  select(cycle,obesity) %>%
  pivot_longer(cols = -cycle, names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = cycle, values_from = value)

svy_stats_hbp <- descriptive_chars %>% 
  mutate(hbp = paste0(round(`hbp`*100,1)," (",round(`ci_l.hbp`*100,1),"-",round(`ci_u.hbp`*100,1),"), CV%= ",round(`cv%.hbp`,2))) %>%
  select(cycle,hbp) %>%
  pivot_longer(cols = -cycle, names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = cycle, values_from = value)

svy_stats_dyslipidemia_abnormal <- descriptive_chars %>% 
  mutate(dyslipidemia_abnormal = paste0(round(`dyslipidemia_abnormal`*100,1)," (",round(`ci_l.dyslipidemia_abnormal`*100,1),"-",round(`ci_u.dyslipidemia_abnormal`*100,1),"), CV%= ",round(`cv%.dyslipidemia_abnormal`,2))) %>%
  select(cycle,dyslipidemia_abnormal) %>%
  pivot_longer(cols = -cycle, names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = cycle, values_from = value)

svy_stats_diabetes <- descriptive_chars %>% 
  mutate(diabetes = paste0(round(`diabetes`*100,1)," (",round(`ci_l.diabetes`*100,1),"-",round(`ci_u.diabetes`*100,1),"), CV%= ",round(`cv%.diabetes`,2))) %>%
  select(cycle,diabetes) %>%
  pivot_longer(cols = -cycle, names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = cycle, values_from = value)

descriptive_chars_table <- rbind(svy_stats_age,svy_stats_age_cat,svy_stats_sex,svy_stats_white,svy_stats_race_ethnic,svy_stats_immigrant,
                                 svy_stats_ei_quartile_moecd,svy_stats_eq_moecd_hhi,svy_stats_ei_quartile_ooecd,svy_stats_eq_ooecd_hhi,svy_stats_ei_quartile_sqrt,svy_stats_eq_sqrt_hhi,
                                 svy_stats_r_poor_cv_health,svy_stats_obesity,svy_stats_hbp,svy_stats_dyslipidemia_abnormal,svy_stats_diabetes)

remove(svy_stats_age,svy_stats_age_cat,svy_stats_sex,svy_stats_white,svy_stats_race_ethnic,svy_stats_immigrant,
       svy_stats_ei_quartile_moecd,svy_stats_eq_moecd_hhi,svy_stats_ei_quartile_ooecd,svy_stats_eq_ooecd_hhi,svy_stats_ei_quartile_sqrt,svy_stats_eq_sqrt_hhi,
       svy_stats_r_poor_cv_health,svy_stats_obesity,svy_stats_hbp,svy_stats_dyslipidemia_abnormal,svy_stats_diabetes)

##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
# PREDICTED PROBABILITIES FOR COMPOSITE OUTCOME    
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################

##############################################################################################################################
##############################################################################################################################
# STRATIFIED BY AGE CATEGORY
##############################################################################################################################
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
  pred_prob_model <- svyglm(r_poor_cv_health_1 ~ age + sex + white + immigrant,
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

pred_probs_age <- do.call(rbind, results_list)

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
prev_trend_results_age <- data.frame(
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
  prev_trend_model_age <- svyglm(r_poor_cv_health_1 ~ cycle_cont + age + sex + white + immigrant,
                                    design = svy_design,
                                    family = quasibinomial(link="logit"),
                                    subset = age_cat == age_cat_i)
  
  #Get the p trend
  cycle_coef <- coef(prev_trend_model_age)["cycle_cont"]
  cycle_pval <- summary(prev_trend_model_age)$coefficients["cycle_cont", "Pr(>|t|)"]
  
  #Append the results to the data frame
  prev_trend_results_age <- rbind(prev_trend_results_age, data.frame(
    age_cat = age_cat_i,
    cycle_coef = cycle_coef,
    p_trend = cycle_pval,
    stringsAsFactors = FALSE
  ))
}

prev_trend_results_age <- prev_trend_results_age %>% arrange(age_cat)

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
  pred_prob_model <- svyglm(r_poor_cv_health_1 ~ age + white + immigrant,
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

pred_probs_age_sex <- do.call(rbind, results_list)

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
prev_trend_results_age_sex <- data.frame(
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
  prev_trend_model_age_sex <- svyglm(r_poor_cv_health_1 ~ cycle_cont + age + white + immigrant,
                                     design = svy_design,
                                     family = quasibinomial(link="logit"),
                                     subset = sex == sex_i & age_cat == age_cat_i)
  
  #Get the p trend
  cycle_coef <- coef(prev_trend_model_age_sex)["cycle_cont"]
  cycle_pval <- summary(prev_trend_model_age_sex)$coefficients["cycle_cont", "Pr(>|t|)"]
  
  #Append the results to the data frame
  prev_trend_results_age_sex <- rbind(prev_trend_results_age_sex, data.frame(
    age_cat = age_cat_i,
    sex = sex_i,
    cycle_coef = cycle_coef,
    p_trend = cycle_pval,
    stringsAsFactors = FALSE
  ))
}
prev_trend_results_age_sex <- prev_trend_results_age_sex %>% arrange(age_cat, sex)

remove(results_list,prev_trend_model_age_sex)

##############################################################################################################################
##############################################################################################################################
# STRATIFIED BY AGE CATEGORY AND INCOME QUARTILE
##############################################################################################################################
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
  pred_prob_model <- svyglm(r_poor_cv_health_1 ~ age + sex + white + immigrant,
                            design = svy_design,
                            family = quasibinomial(link="logit"),
                            subset = age_cat == age_cat_i)
  
  #Calculate predicted probabilities
  preds <- as.data.frame(svypredmeans(adjustmodel=pred_prob_model,groupfactor=~interaction(ei_quartile_moecd,cycle))) %>%
    #t-distribution Wald-based logit transformation for CI calculation
    mutate(prob_logit = log(mean/(1 - mean)),
           SE_logit = SE/(mean*(1 - mean)),
           lower_CI_logit = prob_logit - qt(p=0.975, df=68) * SE_logit, 
           upper_CI_logit = prob_logit + qt(p=0.975, df=68) * SE_logit,
           lower_CI = plogis(lower_CI_logit),
           upper_CI = plogis(upper_CI_logit)) %>%
    rownames_to_column(var = "values") %>% mutate(ei_quartile_moecd = sub("\\..*", "", values), cycle = sub(".*\\.", "", values)) %>% rename(prob = mean) %>%
    mutate(age_cat = age_cat_i) %>%
    select(age_cat, ei_quartile_moecd, cycle, prob, SE, lower_CI, upper_CI) %>%
    arrange(ei_quartile_moecd,cycle)
  
  results_list[[i]] <- preds
  rm(preds,pred_prob_model)
  
}

pred_probs_age_ei <- do.call(rbind, results_list)

remove(results_list)

##############################################################################################################################
# Logistic regression for linear trends
# Adjusted for age + sex + white + immigrant
##############################################################################################################################
#Define possible values for age_cat and sex
age_cats <- c("6-17 years", "18-39 years", "40-64 years")
ei_quartiles_moecd <- c("1", "2", "3", "4")
strata_combinations <- expand.grid(age_cat = age_cats, ei_quartile_moecd = ei_quartiles_moecd)
print(strata_combinations)  

#Compute p for trend
prev_trend_results_age_ei <- data.frame(
  age_cat = character(),
  ei_quartile_moecd = character(),
  cycle_coef = numeric(),
  p_trend = numeric(),
  stringsAsFactors = FALSE
)

results_list <- list()
for (i in seq_len(nrow(strata_combinations))) {
  #Extract current values
  age_cat_i <- strata_combinations$age_cat[i]
  ei_quartile_moecd_i <- strata_combinations$ei_quartile_moecd[i]
  
  #Run prevalence trend models
  prev_trend_model_age_ei <- svyglm(r_poor_cv_health_1 ~ cycle_cont + age + sex + white + immigrant,
                                    design = svy_design,
                                    family = quasibinomial(link="logit"),
                                    subset = age_cat == age_cat_i & ei_quartile_moecd == ei_quartile_moecd_i)
  
  #Get the p trend
  cycle_coef <- coef(prev_trend_model_age_ei)["cycle_cont"]
  cycle_pval <- summary(prev_trend_model_age_ei)$coefficients["cycle_cont", "Pr(>|t|)"]
  
  #Append the results to the data frame
  prev_trend_results_age_ei <- rbind(prev_trend_results_age_ei, data.frame(
    age_cat = age_cat_i,
    ei_quartile_moecd = ei_quartile_moecd_i,
    cycle_coef = cycle_coef,
    p_trend = cycle_pval,
    stringsAsFactors = FALSE
  ))
}

prev_trend_results_age_ei <- prev_trend_results_age_ei %>% arrange(age_cat, ei_quartile_moecd)

remove(results_list,prev_trend_model_age_ei)

##############################################################################################################################
##############################################################################################################################
# STRATIFIED BY AGE CATEGORY, SEX, AND INCOME QUARTILE
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
  pred_prob_model <- svyglm(r_poor_cv_health_1 ~ age + white + immigrant,
                            design = svy_design,
                            family = quasibinomial(link="logit"),
                            subset = age_cat == age_cat_i & sex == sex_i)
  
  #Calculate predicted probabilities
  preds <- as.data.frame(svypredmeans(adjustmodel=pred_prob_model,groupfactor=~interaction(ei_quartile_moecd,cycle))) %>%
    #t-distribution Wald-based logit transformation for CI calculation
    mutate(prob_logit = log(mean/(1 - mean)),
           SE_logit = SE/(mean*(1 - mean)),
           lower_CI_logit = prob_logit - qt(p=0.975, df=68) * SE_logit, 
           upper_CI_logit = prob_logit + qt(p=0.975, df=68) * SE_logit,
           lower_CI = plogis(lower_CI_logit),
           upper_CI = plogis(upper_CI_logit)) %>%
    rownames_to_column(var = "values") %>% mutate(ei_quartile_moecd = sub("\\..*", "", values), cycle = sub(".*\\.", "", values)) %>% rename(prob = mean) %>%
    mutate(age_cat = age_cat_i, sex = sex_i) %>%
    select(age_cat, sex, ei_quartile_moecd, cycle, prob, SE, lower_CI, upper_CI) %>%
    arrange(ei_quartile_moecd,sex,cycle)
  
  results_list[[i]] <- preds
  rm(preds,pred_prob_model)
  
}

pred_probs_age_sex_ei <- do.call(rbind, results_list)

remove(results_list)

##############################################################################################################################
# Logistic regression for linear trends
# Adjusted for age + white + immigrant
##############################################################################################################################
#Define possible values for age_cat, sex, and equivalized income quartiles
age_cats <- c("6-17 years", "18-39 years", "40-64 years")
ei_quartiles_moecd <- c("1", "2", "3", "4")
sexes <- c("M", "F")
strata_combinations <- expand.grid(age_cat = age_cats, ei_quartile_moecd = ei_quartiles_moecd, sex = sexes)
print(strata_combinations)  

#Compute p for trend
prev_trend_results_age_sex_ei <- data.frame(
  age_cat = character(),
  sex = character(),
  ei_quartile_moecd = character(),
  cycle_coef = numeric(),
  p_trend = numeric(),
  stringsAsFactors = FALSE
)

results_list <- list()
for (i in seq_len(nrow(strata_combinations))) {
  #Extract current values
  sex_i <- strata_combinations$sex[i]
  age_cat_i <- strata_combinations$age_cat[i]
  ei_quartile_moecd_i <- strata_combinations$ei_quartile_moecd[i]
  
  #Run prevalence trend models
  prev_trend_model_age_sex_ei <- svyglm(r_poor_cv_health_1 ~ cycle_cont + age + white + immigrant,
                                        design = svy_design,
                                        family = quasibinomial(link="logit"),
                                        subset = sex == sex_i & age_cat == age_cat_i & ei_quartile_moecd == ei_quartile_moecd_i)
  
  #Get the p trend
  cycle_coef <- coef(prev_trend_model_age_sex_ei)["cycle_cont"]
  cycle_pval <- summary(prev_trend_model_age_sex_ei)$coefficients["cycle_cont", "Pr(>|t|)"]
  
  #Append the results to the data frame
  prev_trend_results_age_sex_ei <- rbind(prev_trend_results_age_sex_ei, data.frame(
    age_cat = age_cat_i,
    sex = sex_i,
    ei_quartile_moecd = ei_quartile_moecd_i,
    cycle_coef = cycle_coef,
    p_trend = cycle_pval,
    stringsAsFactors = FALSE
  ))
}
prev_trend_results_age_sex_ei <- prev_trend_results_age_sex_ei %>% arrange(age_cat, sex, ei_quartile_moecd)

remove(results_list,prev_trend_model_age_sex_ei)

##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
# HEALTH INEQUALITIES    
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
  RII_trend_model <- svyglm(r_poor_cv_health_1 ~ ridit_eqi_moecd*cycle_cont + age + sex + white + immigrant,
                            design = svy_design,
                            family = quasipoisson(link="log"),
                            subset = age_cat == age_cat_i)
  SII_model <- svyglm(r_poor_cv_health_1 ~ ridit_eqi_moecd*cycle + age + sex + white + immigrant,
                      design = svy_design,
                      family = gaussian(link="identity"),
                      subset = age_cat == age_cat_i)
  SII_trend_model <- svyglm(r_poor_cv_health_1 ~ ridit_eqi_moecd*cycle_cont + age + sex + white + immigrant,
                            design = svy_design,
                            family = gaussian(link="identity"),
                            subset = age_cat == age_cat_i)
  
  #Define contrast vectors
  r_c1 <- c("ridit_eqi_moecd" = 1)
  r_c2 <- c("ridit_eqi_moecd" = 1, "ridit_eqi_moecd:cycle2" = 1)
  r_c3 <- c("ridit_eqi_moecd" = 1, "ridit_eqi_moecd:cycle3" = 1)
  r_c4 <- c("ridit_eqi_moecd" = 1, "ridit_eqi_moecd:cycle4" = 1)
  r_c5 <- c("ridit_eqi_moecd" = 1, "ridit_eqi_moecd:cycle5" = 1)
  r_c6 <- c("ridit_eqi_moecd" = 1, "ridit_eqi_moecd:cycle6" = 1)
  
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
  RII_trend_coef <- as.numeric(coef(RII_trend_model)["ridit_eqi_moecd:cycle_cont"])
  RII_trend_SE <- as.numeric(SE(RII_trend_model)["ridit_eqi_moecd:cycle_cont"])
  RII_trend_pvalue <- summary(RII_trend_model)$coefficients["ridit_eqi_moecd:cycle_cont", "Pr(>|t|)"]
  
  SII_coefs <- as.numeric(SII_parms$Est)     
  SII_SEs <- as.numeric(SII_parms$SE) 
  SII_values <- as.numeric(SII_parms$SII)   
  SII_CIs <- SII_parms[4:5] 
  SII_trend_coef <- as.numeric(coef(SII_trend_model)["ridit_eqi_moecd:cycle_cont"])
  SII_trend_SE <- as.numeric(SE(SII_trend_model)["ridit_eqi_moecd:cycle_cont"])
  SII_trend_pvalue <- summary(SII_trend_model)$coefficients["ridit_eqi_moecd:cycle_cont", "Pr(>|t|)"]
  
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

health_inequalities_age <- do.call(rbind, results_list)
rownames(health_inequalities_age) <- NULL
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
  RII_trend_model <- svyglm(r_poor_cv_health_1 ~ ridit_eqi_moecd_sex*cycle_cont + age + white + immigrant,
                            design = svy_design,
                            family = quasipoisson(link="log"),
                            subset = age_cat == age_cat_i & sex == sex_i)
  RII_trend_sex_model <- svyglm(r_poor_cv_health_1 ~ ridit_eqi_moecd_sex*cycle_cont*sex + age + white + immigrant,
                                design = svy_design,
                                family = quasipoisson(link="log"),
                                subset = age_cat == age_cat_i)
  SII_model <- svyglm(r_poor_cv_health_1 ~ ridit_eqi_moecd_sex*cycle + age + white + immigrant,
                      design = svy_design,
                      family = gaussian(link="identity"),
                      subset = age_cat == age_cat_i & sex == sex_i)
  SII_trend_model <- svyglm(r_poor_cv_health_1 ~ ridit_eqi_moecd_sex*cycle_cont + age + white + immigrant,
                            design = svy_design,
                            family = gaussian(link="identity"),
                            subset = age_cat == age_cat_i & sex == sex_i)
  SII_trend_sex_model <- svyglm(r_poor_cv_health_1 ~ ridit_eqi_moecd_sex*cycle_cont*sex + age + white + immigrant,
                                design = svy_design,
                                family = gaussian(link="identity"),
                                subset = age_cat == age_cat_i)
  
  #Define contrast vectors
  r_c1 <- c("ridit_eqi_moecd_sex" = 1)
  r_c2 <- c("ridit_eqi_moecd_sex" = 1, "ridit_eqi_moecd_sex:cycle2" = 1)
  r_c3 <- c("ridit_eqi_moecd_sex" = 1, "ridit_eqi_moecd_sex:cycle3" = 1)
  r_c4 <- c("ridit_eqi_moecd_sex" = 1, "ridit_eqi_moecd_sex:cycle4" = 1)
  r_c5 <- c("ridit_eqi_moecd_sex" = 1, "ridit_eqi_moecd_sex:cycle5" = 1)
  r_c6 <- c("ridit_eqi_moecd_sex" = 1, "ridit_eqi_moecd_sex:cycle6" = 1)
  
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
  RII_trend_coef <- as.numeric(coef(RII_trend_model)["ridit_eqi_moecd_sex:cycle_cont"])
  RII_trend_SE <- as.numeric(SE(RII_trend_model)["ridit_eqi_moecd_sex:cycle_cont"])
  RII_trend_pvalue <- summary(RII_trend_model)$coefficients["ridit_eqi_moecd_sex:cycle_cont", "Pr(>|t|)"]
  RII_trend_sex_coef <- as.numeric(coef(RII_trend_sex_model)["ridit_eqi_moecd_sex:cycle_cont:sexM"])
  RII_trend_sex_SE <- as.numeric(SE(RII_trend_sex_model)["ridit_eqi_moecd_sex:cycle_cont:sexM"])
  RII_trend_sex_pvalue <- summary(RII_trend_sex_model)$coefficients["ridit_eqi_moecd_sex:cycle_cont:sexM", "Pr(>|t|)"]  
  
  SII_coefs <- as.numeric(SII_parms$Est)        
  SII_SEs <- as.numeric(SII_parms$SE)   
  SII_values <- as.numeric(SII_parms$SII)   
  SII_CIs <- SII_parms[4:5]
  SII_trend_coef <- as.numeric(coef(SII_trend_model)["ridit_eqi_moecd_sex:cycle_cont"])
  SII_trend_SE <- as.numeric(SE(SII_trend_model)["ridit_eqi_moecd_sex:cycle_cont"])
  SII_trend_pvalue <- summary(SII_trend_model)$coefficients["ridit_eqi_moecd_sex:cycle_cont", "Pr(>|t|)"]
  SII_trend_sex_coef <- as.numeric(coef(SII_trend_sex_model)["ridit_eqi_moecd_sex:cycle_cont:sexM"])
  SII_trend_sex_SE <- as.numeric(SE(SII_trend_sex_model)["ridit_eqi_moecd_sex:cycle_cont:sexM"])
  SII_trend_sex_pvalue <- summary(SII_trend_sex_model)$coefficients["ridit_eqi_moecd_sex:cycle_cont:sexM", "Pr(>|t|)"]      
  
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

health_inequalities_age_sex <- do.call(rbind, results_list)
rownames(health_inequalities_age_sex) <- NULL




