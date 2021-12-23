######################################################################
## Title: [Insert full title of paper]
## Short title: [Insert short title of paper (if applicable) here]
## DOI: [Insert DOI of paper here]
## Code author: Steven Kerr steven.kerr@ed.ac.uk
## Description: S-gene negative from 1st Nov 2021 onwards is a good proxy for omicron variant.
##              This code does a logistic regression with s gene positive/negative as the
##              outcome variable
######################################################################

library(tidyverse)
library(mgcv)

setwd('/conf/EAVE/GPanalysis/analyses/B.1.1.529-variant')

# Taken from:
# '/conf/EAVE/GPanalysis/progs/CR/SGene_Omicron/04a_TND_Vacc_Status_17122021.html'
tnd_results <- as.data.frame( matrix( c('uv',	68224,	15757,	0.00,	0.00,	0.00,
                                        'v1_0:3',	1185,	267,	50.11,	42.35,	56.82,
                                        'v1_4+',	13301,	3606,	45.70,	43.01,	48.26,
                                        'v2_0:1',	700,	101,	71.81,	64.90,	77.36,                                        
                                        'v2_2:9',	3782,	270,	86.63,	84.77,	88.26,
                                        'v2_10:14',	12455,	1719,	75.57,	73.97,	77.07,
                                        'v2_15:19',	27565,	5858,	66.90,	65.29,	68.44,
                                        'v2_20:24',	22828,	7030,	54.73,	52.42,	56.92,
                                        'v2_25+',	22594,	6518,	52.13,	49.68,	54.46,
                                        'v3_0',	6240,	1655,	62.36,	59.52,	64.99,
                                        'v3_1',	4046,	373,	89.37,	88.05,	90.54,
                                        'v3_2+',	22902,	1343,	93.17,	92.64,	93.66), 
                          ncol = 6, byrow = TRUE  ) )


# '/conf/EAVE/GPanalysis/progs/CR/SGene_Omicron/04a_TND_Vacc_Status_17122021.html'
tnd_results <- as.data.frame( matrix( c(
                              'uv',	68224,	15757,	0.00,	0.00,	0.00,
                              'v1_0:3_AZ',	12,	3,	59.64,	-52.53,	89.32,
                              'v1_4+_AZ',	1008,	307,	39.22,	29.65,	47.49,
                              'v2_0:1_AZ',	69,	15,	62.51,	31.89,	79.36,
                              'v2_2:9_AZ',	365,	63,	72.30,	63.28,	79.11,
                              'v2_10:14_AZ',	746,	206,	52.88,	44.04,	60.32,
                              'v2_15:19_AZ',	6803,	2316,	45.23,	41.30,	48.89,
                              'v2_20:24_AZ',	15265,	5314,	42.51,	39.25,	45.60,
                              'v2_25+_AZ',	15326,	4972,  41.67,	38.39,	44.78,
                              'v3_0_AZ',	4715,	1365,	53.76,	49.95,	57.29,
                              'v3_1_AZ',	2899,	283,	87.41,	85.62,	88.98,
                              'v3_2+_AZ',	8613,	497,	93.15,	92.35,	93.86,
                              'v1_0:3_Mo',	135,	34,	30.98,	-4.43,	54.39,
                              'v1_4+_Mo',	699,	162,	43.91,	32.39,	53.48,
                              'v2_0:1_Mo',	78,	7,	80.79,	57.53,	91.31,
                              'v2_2:9_Mo',	708,	35,	89.89,	85.71,	92.84,
                              'v2_10:14_Mo',	2750,	273,	82.01,	79.46,	84.24,
                              'v2_15:19_Mo',	2441,	330,	78.32,	75.46,	80.84,
                              'v2_20:24_Mo',	667,	124,	73.12,	67.07,	78.05,
                              'v2_25+_Mo',	77,	9,	80.68,	60.83,	90.47,
                              'v3_0_Mo',	36,	7, 67.90,	25.48,	86.18,
                              'v3_1_Mo',	22,	1,	91.61,	36.45,	98.89,
                              'v3_2+_Mo',	12,	1,	84.09,	-26.14,	97.99,
                              'v1_0:3_PB',	1038,	230,	52.10,	44.10,	58.95,
                              'v1_4+_PB',	11594,	3137,	46.33,	43.53,	48.99,
                              'v2_0:1_PB',	553,	79,	71.82,	63.96,	77.96,
                              'v2_2:9_PB',	2709,	172,	88.04,	85.96,	89.81,
                              'v2_10:14_PB',	8959,	1240, 75.95,	74.17,	77.60,
                              'v2_15:19_PB',	18321,	3212,	71.20,	69.61,	72.71,
                              'v2_20:24_PB',	6896,	1592,	61.07,	58.27,	63.69,
                              'v2_25+_PB',	7191,	1537, 59.93,	57.02,	62.65,
                              'v3_0_PB',	1489,	283,	69.15,	64.54,	73.17,
                              'v3_1_PB',	1125,	89,	89.65,	87.06,	91.72,
                              'v3_2+_PB',	14277,	845,	92.12,	91.43,	92.75), 
                                      ncol = 6, byrow = TRUE  ) )




names(tnd_results) <-c('vs', 'tested',	'positive', 	'OR_TND', 'lcl', 'ucl')

tnd_results <- mutate_at(tnd_results, c(2:6), as.numeric) %>%
  # The OR column is actually vaccine effect expressed as a percentage. Convert it
               mutate_at(c('OR_TND', 'lcl', 'ucl'),  ~1 - ./100) %>%
               mutate( se_ln_OR_TND = (log(ucl)- log(lcl)) /3.92,
                       vs = paste0('vs', vs)) 
               



df <- readRDS('/conf/EAVE/GPanalysis/progs/CR/SGene_Omicron/output/temp/Sgene_all_positive.rds')

# Relevel s_gene so that positive is the reference category.
df <- mutate(df, s_gene = relevel(s_gene, ref= "S_Pos"),
            #vs = relevel(vs, ref = 'uv'))
            vs = relevel(vs, ref = 'v2_2-5'),
            vt = relevel(vt, ref = 'v2_2-5_PB'))

factor_vars <- c('sex', 'simd2020_sc_quintile', 'HB', 'n_risk_gps', 'prev_pos', 'immuno', 'shielding')

df <- mutate_at(df, factor_vars, as.factor)



# Create the n_oth_risk_groups variable
variables_hosp <- c("Q_DIAG_AF", 
                    "Q_DIAG_ASTHMA", 
                    "Q_DIAG_CCF", 
                    "Q_DIAG_CHD", 
                    "Q_DIAG_COPD" , 
                    "Q_DIAG_DIABETES", 
                    "Q_DIAG_EPILEPSY",
                    "Q_DIAG_FRACTURE", 
                    "Q_DIAG_PVD", 
                    "Q_DIAG_RA_SLE", 
                    "Q_DIAG_SEV_MENT_ILL", 
                    "Q_DIAG_STROKE", 
                    "Q_DIAG_VTE", 
                    "Q_DIAG_CKD",
                    "EAVE_Smoke", 
                    "EAVE_BP")


z_vars <- names(df)[grepl("Q_DIAG", names(df))]
z_vars <- z_vars[!(z_vars %in% variables_hosp)]
z_vars <- z_vars[!(z_vars %in% c("Q_DIAG_DIABETES_1"  ,  "Q_DIAG_DIABETES_2" , "Q_DIAG_CKD_LEVEL"))]
z_vars <- c(z_vars, "Q_HOME_CAT","Q_LEARN_CAT")

z <- df %>% 
     dplyr::select_at(c("EAVE_LINKNO", z_vars)) %>%
     mutate_at(z_vars, ~as.numeric(as.character(.))) %>%
     mutate(Q_HOME_CAT = if_else(Q_HOME_CAT >=1, 1, Q_HOME_CAT),
            Q_LEARN_CAT = if_else(Q_LEARN_CAT >=1, 1, Q_LEARN_CAT))%>%
    mutate(N = rowSums(across(all_of(z_vars)))) %>%
     mutate(n_oth_risk_gps = cut(N, breaks=c(-1,0,1,max(N, na.rm=TRUE)), labels=c("0","1","2+")) ) %>% 
     dplyr::select(EAVE_LINKNO, n_oth_risk_gps)

df <- df %>% left_join(z, by="EAVE_LINKNO")






# This function creates a table of ORs and their CIs for the categorical variables
# and OR plots for the spline variables.
# vacc_status_varialbe = 'vt' or 'vs' determines which vaccine status variable to use
# age_lower and age_upper are lower and upper bounds for age
# vacc = 'all', or '1+' or '2+' determines whether all vaccine reicipients, or only
# those who have at least 1/2 doses should be included

fit_model <- function(vacc_status_variable, age_lower, age_upper, dose){
  
  # output subdirectory
  folder = paste0('dose_', dose, '_age_', age_lower, '-', age_upper, '_', vacc_status_variable)
  
  # If output subdirectory doesn't exist, create it
  ifelse(!dir.exists(file.path('./output/', folder)), 
         dir.create(file.path('./output/', folder)), 
         FALSE)
  
  # vs is the name of variable used in analysis. If we want vt instead, rename.
  if (vacc_status_variable == 'vt'){
    df <- dplyr::rename(df, 'vs_old' = 'vs', 'vs' = 'vt' )
  }
  
  subset <- filter(df, s_gene %in% c("S_Pos", "S_Neg"),
                        age_year >= age_lower, 
                        age_year <= age_upper)
  
  if (dose == '1+'){
    subset <- filter(subset, !grepl('uv', vs1 ))
  }
  
  if (dose == '2+'){
    subset <- filter(subset, !grepl('uv', vs2 ))
  }
  
  
  non_spline_predictor_vars <- c('sex',
                                'simd2020_sc_quintile', 
                                'HB',  
                                'n_oth_risk_gps', 
                                'immuno', 
                                'shielding',
                                'prev_pos', 
                                'vs',
                                z_vars)
  
  formula <- as.formula(paste0("s_gene ~ s(age_year) + s(days) + ", 
                               paste(non_spline_predictor_vars ,collapse=" + ")))
  
  
  #logistic Regression Model
  logit <-gam(formula, data=subset, family=binomial)
  
  # Create a table of results
  coef <- summary(logit)$p.coeff
  
  se <- summary(logit)$se[names(coef)]
  
  results <- data.frame( term = names(coef),
                         OR_conditional = round(exp(coef), 2),
                         se_ln_OR_conditional = coef,
                         CI_conditional = paste0( '(', round(exp(coef - 1.96*se),2), '-'
                                            , round(exp(coef + 1.96*se),2), ')'),
                         row.names = 1:length(coef))
  
  ## Calculate unconditional ORs for delta/omicron infection
  # This is approximately equal to 
  # (OR for delta/omicron condition on testing +ve) * (OR for symptomcatic delta infection from TND)
  
  results <- left_join(results, select(tnd_results, vs, OR_TND, lcl, ucl, se_ln_OR_TND), 
                       by = c('term' = 'vs'))
  
  results <- mutate(results, OR_unconditional = OR_conditional * OR_TND,
                    se_ln_OR_unconditional =  sqrt(se_ln_OR_conditional**2 + se_ln_OR_TND**2)) %>%
             mutate(CI_unconditional = ifelse(!is.na(OR_unconditional), 
              paste0( '(', round(exp( log(OR_unconditional) - 1.96*se_ln_OR_unconditional),2), 
               '-', round(exp( log(OR_unconditional) + 1.96*se_ln_OR_unconditional),2), ')'), NA)) %>%
             mutate(OR_unconditional = round(OR_unconditional, 2))
                      
                    
  write.csv(results, paste0('./output/', folder, '/OR_table.csv'))
  
  
  
  
  
  
  # Create a row of test data with all factors set to their baseline level
  test_row <- subset[1, c('age_year', 'days', non_spline_predictor_vars)]
  
  test_row <- mutate_at(test_row, non_spline_predictor_vars,  ~levels(.)[1])

  
  ## Create plot of odds ratio for spline in days
  
  max_days <- max(subset$days)
  
  test_data <- test_row[rep(1, max_days), ]
  
  # Set age to zero and days to the range we wish to predict over
  test_data$age_year <- 0
  test_data$days <- 1:max_days
  
  X <- predict.gam(logit, newdata=test_data, type="lpmatrix")
  
  ln_odds <- X %*% coef(logit)
  
  # Obtain 95% CIs
  se <- sqrt(diag(X %*% vcov(logit, unconditional=TRUE) %*% t(X)))
  lcl <- ln_odds - 1.96*se
  ucl <- ln_odds + 1.96*se
  
  spline_OR <- data.frame( Day = 1:max_days, 
                           OR = exp( ln_odds ),
                           lcl = exp(lcl), 
                           ucl = exp(ucl))
  
  ggplot(data=spline_OR, aes(x=Day, y=OR)) + geom_line() +
        geom_ribbon(aes(ymin=lcl, ymax=ucl), alpha=0.3) + 
        ylab('Odds ratio')
  
  ggsave(paste0('./output/', folder, '/spline_days_OR.png'))
  
  
  
  
  
  ## Create plot of odds ratio for spline in age
  
  max_age <- max(subset$age_year)
  
  test_data <- test_row[rep(1, max_age), ]
  
  # Set age to zero and days to the range we wish to predict over
  test_data$age_year <- 1:max_age
  test_data$days <- 0
  
  X <- predict.gam(logit, newdata=test_data, type="lpmatrix")
  
  ln_odds <- X %*% coef(logit)
  
  # Obtain 95% CIs
  se <- sqrt(diag(X %*% vcov(logit, unconditional=TRUE) %*% t(X)))
  lcl <- ln_odds - 1.96*se
  ucl <- ln_odds + 1.96*se
  
  spline_OR <- data.frame( Age = 1:max_age, 
                           OR = exp( ln_odds ),
                           lcl = exp(lcl), 
                           ucl = exp(ucl))
  
  ggplot(data=spline_OR, aes(x=Age, y=OR)) + geom_line() +
    geom_ribbon(aes(ymin=lcl, ymax=ucl), alpha=0.3) +
    scale_y_continuous(limits = c(0, 0.005))
    ylab('Odds ratio')
  
  
  ggsave(paste0('./output/', folder, '/spline_age_OR.png'))
  
  ## Plot odds ratio for age term
  # This is used when age is not included as a spline
  # plot_OR <- function(model_fit, term){
  #   # plots odds ratios for a single term in a fitted model
  #   
  #   or <- termplot(model_fit, term = term, se = T, plot = F)
  #   
  #   var <- names(or)
  #   
  #   or <- or[[var]]
  #   
  #   or <- mutate(or, ucl = y + 1.96*se,
  #                lcl = y - 1.96*se) %>%
  #     mutate_at(c('y', 'ucl', 'lcl'), exp)
  #   
  #   or <- do.call(data.frame,lapply(or, function(x) replace(x, is.infinite(x),NA)))
  #   
  #   output <- ggplot(data=or, aes(x=x, y=y)) + geom_line() +
  #     geom_ribbon(aes(ymin=lcl, ymax=ucl), linetype=2, alpha=0.3)  + 
  #     ylab("Odds ratio")
  # }
  # 
  # OR_plot <- plot_OR(logit, 'age_year') + xlab('Age')
  # 
  # ggsave('./output/age_OR.png', OR_plot)

}

# Set up a dataframe where each row is a list of arguments to be fed into fit_model
max_age <- max(df$age_year)

arguments <- expand.grid( vacc_status_variable = c('vt', 'vs'),
             age_lower = c(16, 50),
             age_upper = c(50, max_age),
             dose = c('all', '1+', '2+')) %>% 
             filter(age_upper > age_lower) %>%
             mutate_if(is.factor, as.character)


for (row in 1:nrow(arguments)){
  
  fit_model( vacc_status_variable = arguments[row, 'vacc_status_variable'],
            age_lower = arguments[row, 'age_lower'],
            age_upper = arguments[row, 'age_upper'],
            dose = arguments[row, 'dose'])
}



