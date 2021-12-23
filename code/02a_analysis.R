##########################################################
## Title: 1st dose COVID-19 vaccine waning
## Code author(s): Steven Kerr steven.kerr@ed.ac.uk
## Description: 02_descriptive - Descriptive analyses on the
##              baseline cohort
##########################################################


# Libraries
library("finalfit")

source('00_functions.R')

############## 0 Functions #######################

fun.extract <- function(z.fit) {
  #takes a coxph filt using penalised splines and drops off the ps terms
  #make sure no variable begins ps
  z <- summary(z.fit)
  z <- data.frame(z$conf.int)
  z <- z %>% mutate(names = row.names(z)) %>% 
    filter(!(grepl("^ps", names))) %>% 
    dplyr::relocate(names, .before=1) %>% 
    dplyr::select(-exp..coef.)
  names(z) <- c("names","HR","LCL","UCL")
  z
}


plot_HR <- function(model_fit, term){
  # plots hazard ratios for a single term in a fitted model
  
  hr <- termplot(model_fit, term = term, se = T, plot = F)
  
  var <- names(hr)
  
  hr <- hr[[var]]
  
  hr <- mutate(hr, ucl = y + 1.96*se,
               lcl = y - 1.96*se) %>%
    mutate_at(c('y', 'ucl', 'lcl'), exp)
  
  hr <- do.call(data.frame,lapply(hr, function(x) replace(x, is.infinite(x),NA)))
  
  output <- ggplot(data=hr, aes(x=x, y=y)) + geom_line() +
    geom_ribbon(aes(ymin=lcl, ymax=ucl), linetype=2, alpha=0.1, fill = 'steelblue')  + 
    ylab("Hazard Ratio")
  
  if (var == 'ageYear'){
    output <- output + xlab("Age")
  } else if (var == 'days'){
    output <- output + xlab("Days since first specimen collection date")
  }
}

df_pos <- readRDS("output/temp/Sgene_all_positive.rds") 


##### 1 Descriptive tables ####

# Uses function summary_factorlist_wt from 00_functions.R

df_pos$Total = 'Total'

## Full cohort summary tables
rgs <- colnames(df_pos)[startsWith(colnames(df_pos), "Q")]

explanatory <- c("Total",
                 "sex", 
                              "age_gp",  
                 "vs", "n_risk_gps", 'prev_pos', 'flag_covid_symptomatic')

summary_tbl_wt_chrt <- summary_factorlist_wt(df_pos, "Total", explanatory = explanatory, wght="wght_all") 

names(summary_tbl_wt_chrt) <- c('Characteristic', 'Levels', 'Total')
summary_tbl_wt_chrt$Characteristic[duplicated(summary_tbl_wt_chrt$Characteristic)] <- ''
summary_tbl_wt_chrt[1, 'Levels'] <- ''

write.csv(summary_tbl_wt_chrt , paste0("./output/summary_table_weights_cohort_all.csv"), row.names = F)

explanatory <- c("Total", "simd2020_sc_quintile", 
                 "ur6_2016_name", 
                  rgs, "bmi_gp",
                 'EAVE_Smoke',
                 'EAVE_BP')

summary_tbl_wt_chrt <- summary_factorlist_wt(df_pos, "Total", explanatory = explanatory) 
names(summary_tbl_wt_chrt) <- c('Characteristic', 'Levels', 'Total')
summary_tbl_wt_chrt$Characteristic[duplicated(summary_tbl_wt_chrt$Characteristic)] <- ''
summary_tbl_wt_chrt[1, 'Levels'] <- ''

write.csv(summary_tbl_wt_chrt , paste0("./output/summary_table_weights_cohort_eave.csv"), row.names = F)

#by S gene status
explanatory <- c("Total",
                 "sex", 
                                  "age_gp",  
                 "n_risk_gps", "vs", "vt", 'prev_pos')

summary_tbl_wt_chrt <- summary_factorlist_wt(df_pos, "s_gene", explanatory = explanatory, wght="wght_all") 
names(summary_tbl_wt_chrt)[1:2] <- c('Characteristic', 'Levels')
summary_tbl_wt_chrt$Characteristic[duplicated(summary_tbl_wt_chrt$Characteristic)] <- ''
summary_tbl_wt_chrt[1, 'Levels'] <- ''

write.csv(summary_tbl_wt_chrt , paste0("./output/summary_table_weights_cohort_all_sg.csv"), row.names = F)


z_df <- filter(df_pos, s_gene %in% c("S_Pos","S_Neg","Weak_S_Pos"))
z <- summary_factorlist_wt(z_df, "vs", explanatory = "s_gene", wght="wght_all") 
z1 <- z %>% dplyr::select(-characteristic) %>% pivot_longer(cols=names(z)[3:ncol(z)]) %>% pivot_wider(names_from=levels, values_from=value)
colnames(z1)[1] <- "Vaccine Status"
write.csv(z1 , paste0("./output/summary_table_vs_sgene.csv"), row.names = F)

z_df <- filter(df_pos, s_gene %in% c("S_Pos","S_Neg","Weak_S_Pos"))
z <- summary_factorlist_wt(z_df, "vt", explanatory = "s_gene", wght="wght_all") 
z1 <- z %>% dplyr::select(-characteristic) %>% pivot_longer(cols=names(z)[3:ncol(z)]) %>% pivot_wider(names_from=levels, values_from=value)
colnames(z1)[1] <- "Vaccine Status"
write.csv(z1 , paste0("./output/summary_table_vt_sgene.csv"), row.names = F)

z_df <- filter(df_pos, s_gene %in% c("S_Pos","S_Neg","Weak_S_Pos"))
z <- summary_factorlist_wt(z_df, "age_gp", explanatory = "s_gene", wght="wght_all") 
z1 <- z %>% dplyr::select(-characteristic) %>% pivot_longer(cols=names(z)[3:ncol(z)]) %>% pivot_wider(names_from=levels, values_from=value)
colnames(z1)[1] <- "Age Group"
write.csv(z1 , paste0("./output/summary_table_agegp_sgene.csv"), row.names = F)

##################### 2 Descriptive graphs ########################

# Positive tests by day

z <- df_pos %>% select(EAVE_LINKNO, date_ecoss_specimen, lab) %>% 
  group_by(date_ecoss_specimen, lab) %>%
  summarise(N=n()) %>%
  mutate(lab = case_when(lab == "lh" ~ 'Lighthouse',
                               lab == "nhs" ~ 'NHS Labs')) 
z <- z %>% filter(date_ecoss_specimen <= max(z$date_ecoss_specimen) - 2) %>% 
  mutate(days=as.numeric(date_ecoss_specimen - min(z$date_ecoss_specimen)))

z %>%  ggplot(aes(x=date_ecoss_specimen, y=N, colour = lab)) + geom_point() +
  labs(x="Specimen date",y ="Number", colour="Laboratory", title="Positive tests by day") 
ggsave(paste0("./output/pos_tests_by_day_lab.png"), width=14, height=10, unit="cm")



# Emergency covid hospitalisation or covid deaths by day

z <- df_pos %>% filter(hosp_covid_emerg==1 & In_Hosp_At_Test == "no" & lab == 'lh') %>%
  mutate(event_date = admission_date) %>%
  group_by(event_date, s_gene) %>% 
  dplyr::summarise(N=n())  


z_grid <- expand.grid(seq(min(z$event_date), max(z$event_date), by="days"), unique(z$s_gene))
names(z_grid) <- colnames(z)[1:2]  

z_grid <- z_grid %>%
  left_join(z) %>%
  replace_na( list(N = 0) )

z_grid %>%  ggplot(aes(x=event_date, y=N, colour = s_gene)) + geom_point() +
  labs(x="Event date",y ="Number", colour="S Gene", title="Emergency covid hospital admissions by day")
ggsave(paste0("./output/hosp_sgene_day.png"), width=14, height=10, unit="cm")


#sequencing and S gene
z_df <- df_pos %>% filter(date_ecoss_specimen <= a_end_wgs)
z1 <- z_df %>% group_by(s_gene, variant) %>% dplyr::summarise(N=n()) %>% 
  pivot_wider(names_from=variant, values_from=N) %>% 
  mutate(across(where(is.numeric), ~replace_na(.,0)))
write.csv(z1 , paste0("./output/summary_table_sgene_variant.csv"), row.names = F)


# Person years to event by s_gene
z.rv <- "hosp_covid_emerg" 
z.rv.time <- "Time.To.Hosp" 

fmla.plot <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~  s_gene"))
z_df <- df_pos %>% filter(In_Hosp_At_Test == "no" & lab == 'lh')
z_df <- mutate(z_df, Time.To.Hosp = if_else(Time.To.Hosp==0,0.1,Time.To.Hosp))  #add 0.1 to 0 days so survival models work
z.tab <- pyears(fmla.plot, data=z_df, data.frame=TRUE)$data
write.csv(z.tab, paste0("./output/pyears_by_sgene_lh.csv"), row.names = F)

# Person years to event by sg_lab
fmla.plot <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~  sg_lab"))
z_df <- df_pos %>% filter(In_Hosp_At_Test == "no")
z_df <- mutate(z_df, Time.To.Hosp = if_else(Time.To.Hosp==0,0.1,Time.To.Hosp))  #add 0.1 to 0 days so survival models work
z.tab <- pyears(fmla.plot, data=z_df, data.frame=TRUE)$data
write.csv(z.tab, paste0("./output/pyears_by_sglab.csv"), row.names = F)


# cumulative incidence curves
fmla.plot <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~  s_gene"))
z_df <- df_pos %>% filter(In_Hosp_At_Test == "no" & lab == 'lh')
z_df <- mutate(z_df, Time.To.Hosp = if_else(Time.To.Hosp==0,0.1,Time.To.Hosp))  #add 0.1 to 0 days so survival models work
z.survfit <- survfit(fmla.plot, data=z_df)

png("./output/cuminc_sg.png", width=8, height=6, unit="in", res=72)
plot(z.survfit, fun="event", col=1:5, xlab="Days from test to emergency hospital admission",
     ylab="Risk")
legend("topleft",col=1:5, lty=1, legend=levels(z_df$s_gene) )
dev.off()

# cumulative incidence curves
fmla.plot <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~  age_gp"))
z_df <- df_pos %>% filter(In_Hosp_At_Test == "no" & lab == 'lh') %>% filter(s_gene=="S_Pos")
z_df <- mutate(z_df, Time.To.Hosp = if_else(Time.To.Hosp==0,0.1,Time.To.Hosp))  #add 0.1 to 0 days so survival models work
z.survfit <- survfit(fmla.plot, data=z_df)

png("./output/cuminc_spos_age.png", width=8, height=6, unit="in", res=72)
plot(z.survfit, fun="event", col=1:6, xlab="Days from test to emergency hospital admission",
     ylab="Risk")
legend("topleft",col=1:6, lty=1, legend=levels(z_df$age_gp) )
dev.off()

#calculated the expected hospitalisations for s neg based upon s pos
fmla.coxph <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~   days_gp + age_gp + sex +
          simd2020_sc_quintile + n_risk_gps +vt + prev_pos"))
df_pos <- df_pos %>% 
  mutate(days_gp = cut(days, breaks = c(-1,6,13,20,27,34,41,48), labels=paste0("week_",1:7)))
z_df <- df_pos %>% filter(in_eave==1 & lab=="lh" & In_Hosp_At_Test == "no")
z_df <- mutate(z_df, Time.To.Hosp = if_else(Time.To.Hosp==0,0.1,Time.To.Hosp))  #add 0.1 to 0 days so survival models work
z_df$n_risk_gps <- factor(z_df$n_risk_gps)
z <- coxph(fmla.coxph, data=z_df, subset=s_gene=="S_Pos")
summary(z)
saveRDS(z,"./output/coxph_model_s_pos.rds")

z_pred <- predict(z, newdata=z_df, type="expected")
summary(z_pred)
z_df$pred <- z_pred
z_tab_pred <- z_df %>% group_by(s_gene) %>% 
  dplyr::summarise(N=n(), Person_Years = round(sum(Time.To.Hosp)/365.25,1), Covid_Hosp= sum(hosp_covid_emerg),
                   Expected = sum(pred)) %>% 
  as.data.frame() 
z_ci <- epitools::pois.byar(z_tab_pred$Covid_Hosp, z_tab_pred$Expected)
z_tab_pred <- bind_cols(z_tab_pred, z_ci[,c("rate","lower","upper")])
write.csv(z_tab_pred, paste0("./output/pyears_by_sg_exp_hosp_eave.csv"), row.names = F)


z_df_pred <- df_pos %>% filter(lab=="lh" & In_Hosp_At_Test == "no" ) %>% 
  dplyr::select(EAVE_LINKNO, hosp_covid_emerg, Time.To.Hosp, sex, days_gp, age_gp, 
                simd2020_sc_quintile, n_risk_gps, vt,  prev_pos , s_gene)
z_df_pred <- z_df_pred %>% mutate(n_risk_gps = fct_recode(n_risk_gps, "0" = "Unknown"),
                                  simd2020_sc_quintile = fct_explicit_na(simd2020_sc_quintile, na_level = "2"))

z_pred <- predict(z, newdata=z_df_pred, type="expected")
summary(z_pred)
z_df_pred$pred <- z_pred
z_tab_pred <- z_df_pred %>% group_by(s_gene) %>% 
  dplyr::summarise(N=n(), Person_Years = round(sum(Time.To.Hosp)/365.25,1), Covid_Hosp= sum(hosp_covid_emerg),
                   Expected = sum(pred)) %>% 
  as.data.frame() 
z_ci <- epitools::pois.byar(z_tab_pred$Covid_Hosp, z_tab_pred$Expected)
z_tab_pred <- bind_cols(z_tab_pred, z_ci[,c("rate","lower","upper")])
write.csv(z_tab_pred, paste0("./output/pyears_by_sg_exp_hosp_all.csv"), row.names = F)


#Hospitalisation rate by age group
fmla.plot <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~  s_gene + age_gp"))
z_df <- df_pos %>% filter(In_Hosp_At_Test == "no" & lab=="lh")
z_df <- mutate(z_df, Time.To.Hosp = if_else(Time.To.Hosp==0,0.1,Time.To.Hosp))  #add 0.1 to 0 days so survival models work
z.tab <- pyears(fmla.plot, data=z_df, data.frame=TRUE)$data
z.tab <- z.tab %>% mutate(rate_100py=event/(pyears/100))
z_ci <- epitools::pois.byar(z.tab$event, z.tab$pyears/100)
z.tab <- bind_cols(z.tab,z_ci) %>% mutate(lower=if_else(lower<0,0,lower))
write.csv(z.tab, paste0("./output/pyears_by_sg_age.csv"), row.names = F)

z_position <- position_dodge(width=0.2)
z.tab %>% filter(s_gene %in% c("S_Pos","S_Neg")) %>% 
  mutate(s_gene=case_when(s_gene=="S_Pos" ~ "S Positive",
                          s_gene=="S_Neg" ~ "S Negative",
                          TRUE ~ "Unknown")) %>% 
  ggplot(aes(x=age_gp, colour=s_gene)) + geom_point(aes(y=rate_100py), position=z_position) + 
  geom_errorbar(aes(ymin=lower, ymax=upper, width=0.2), position = z_position) +
  labs(x="Age Group", y="Rate Per 100 person years", colour="S Gene")
ggsave(paste0("./output/rate_sg_age.png"), width=14, height=10, unit="cm")


#Sensitivity analysis using individuasl testing positive more than 7 days ago
#calculated the expected hospitalisations for s neg based upon s pos
fmla.coxph <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~   days_gp + age_gp + sex +
                               simd2020_sc_quintile + n_risk_gps +vt + prev_pos"))
df_pos <- df_pos %>% 
  mutate(days_gp = cut(days, breaks = c(-1,6,13,20,27,34,41,48), labels=paste0("week_",1:7)))
z_df <- df_pos %>% filter(in_eave==1 & lab=="lh" & In_Hosp_At_Test == "no") %>% 
  filter(date_ecoss_specimen < max(date_ecoss_specimen)-7)
z_df <- mutate(z_df, Time.To.Hosp = if_else(Time.To.Hosp==0,0.1,Time.To.Hosp))  #add 0.1 to 0 days so survival models work
z_df$n_risk_gps <- factor(z_df$n_risk_gps)
z_df$days_gp <- factor(z_df$days_gp)
z <- coxph(fmla.coxph, data=z_df, subset=s_gene=="S_Pos")
summary(z)

z_pred <- predict(z, newdata=z_df, type="expected")
summary(z_pred)
z_df$pred <- z_pred
z_tab_pred <- z_df %>% group_by(s_gene) %>% 
  dplyr::summarise(N=n(), Person_Years = round(sum(Time.To.Hosp)/365.25,1), Covid_Hosp= sum(hosp_covid_emerg),
                   Expected = sum(pred)) %>% 
  as.data.frame() 
z_ci <- epitools::pois.byar(z_tab_pred$Covid_Hosp, z_tab_pred$Expected)
z_tab_pred <- bind_cols(z_tab_pred, z_ci[,c("rate","lower","upper")])
write.csv(z_tab_pred, paste0("./output/pyears_by_sg_exp_hosp_eave_sens7.csv"), row.names = F)

#Sensitivity analysis using those aged 20 to 59
#calculated the expected hospitalisations for s neg based upon s pos
fmla.coxph <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~   days_gp + age_gp + sex +
                               simd2020_sc_quintile + n_risk_gps +vt + prev_pos"))
df_pos <- df_pos %>% 
  mutate(days_gp = cut(days, breaks = c(-1,6,13,20,27,34,41,48), labels=paste0("week_",1:7)))
z_df <- df_pos %>% filter(in_eave==1 & lab=="lh" & In_Hosp_At_Test == "no") %>% 
  filter(age_year >= 20 & age_year <= 59)
z_df <- mutate(z_df, Time.To.Hosp = if_else(Time.To.Hosp==0,0.1,Time.To.Hosp)) %>%   #add 0.1 to 0 days so survival models work
  mutate(age_gp = cut(age_year, breaks=c(19,24,29,34,39,44,49,54,60), labels=c("20-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59")))
z_df$n_risk_gps <- factor(z_df$n_risk_gps)
z <- coxph(fmla.coxph, data=z_df, subset=s_gene=="S_Pos")
summary(z)

z_pred <- predict(z, newdata=z_df, type="expected")
summary(z_pred)
z_df$pred <- z_pred
z_tab_pred <- z_df %>% group_by(s_gene) %>% 
  dplyr::summarise(N=n(), Person_Years = round(sum(Time.To.Hosp)/365.25,1), Covid_Hosp= sum(hosp_covid_emerg),
                   Expected = sum(pred)) %>% 
  as.data.frame() 
z_ci <- epitools::pois.byar(z_tab_pred$Covid_Hosp, z_tab_pred$Expected)
z_tab_pred <- bind_cols(z_tab_pred, z_ci[,c("rate","lower","upper")])
write.csv(z_tab_pred, paste0("./output/pyears_by_sg_exp_hosp_eave_age2059.csv"), row.names = F)

#Sensitivity analysis using simpler model
#calculated the expected hospitalisations for s neg based upon s pos
fmla.coxph <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~   days_gp + age_gp + sex +
                               n_risk_gps +vs "))
df_pos <- df_pos %>% 
  mutate(days_gp = cut(days, breaks = c(-1,6,13,20,27,34,41,48), labels=paste0("week_",1:7)))
z_df <- df_pos %>% filter(in_eave==1 & lab=="lh" & In_Hosp_At_Test == "no") 
z_df <- mutate(z_df, Time.To.Hosp = if_else(Time.To.Hosp==0,0.1,Time.To.Hosp)) #%>%   #add 0.1 to 0 days so survival models work
z_df$n_risk_gps <- factor(z_df$n_risk_gps)
z <- coxph(fmla.coxph, data=z_df, subset=s_gene=="S_Pos")
summary(z)

z_pred <- predict(z, newdata=z_df, type="expected")
summary(z_pred)
z_df$pred <- z_pred
z_tab_pred <- z_df %>% group_by(s_gene) %>% 
  dplyr::summarise(N=n(), Person_Years = round(sum(Time.To.Hosp)/365.25,1), Covid_Hosp= sum(hosp_covid_emerg),
                   Expected = sum(pred)) %>% 
  as.data.frame() 
z_ci <- epitools::pois.byar(z_tab_pred$Covid_Hosp, z_tab_pred$Expected)
z_tab_pred <- bind_cols(z_tab_pred, z_ci[,c("rate","lower","upper")])
write.csv(z_tab_pred, paste0("./output/pyears_by_sg_exp_hosp_eave_simple.csv"), row.names = F)

#sensitivity stratification of age
#calculated the expected hospitalisations for s neg based upon s pos
fmla.coxph <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~   days_gp + strata(age_gp) + sex +
          simd2020_sc_quintile + n_risk_gps + vs + prev_pos"))
df_pos <- df_pos %>% 
  mutate(days_gp = cut(days, breaks = c(-1,6,13,20,27,34,41,48), labels=paste0("week_",1:7))) 
z_df <- df_pos %>% filter(in_eave==1 & lab=="lh" & In_Hosp_At_Test == "no")
z_df <- mutate(z_df, Time.To.Hosp = if_else(Time.To.Hosp==0,0.1,Time.To.Hosp))  #add 0.1 to 0 days so survival models work
z_df$n_risk_gps <- factor(z_df$n_risk_gps)
z <- coxph(fmla.coxph, data=z_df, subset=s_gene=="S_Pos")
summary(z)

z_pred <- predict(z, newdata=z_df, type="expected")
summary(z_pred)
z_df$pred <- z_pred
z_tab_pred <- z_df %>% group_by(s_gene) %>% 
  dplyr::summarise(N=n(), Person_Years = round(sum(Time.To.Hosp)/365.25,1), Covid_Hosp= sum(hosp_covid_emerg),
                   Expected = sum(pred)) %>% 
  as.data.frame() 
z_ci <- epitools::pois.byar(z_tab_pred$Covid_Hosp, z_tab_pred$Expected)
z_tab_pred <- bind_cols(z_tab_pred, z_ci[,c("rate","lower","upper")])
write.csv(z_tab_pred, paste0("./output/pyears_by_sg_exp_hosp_eave_age_strata.csv"), row.names = F)



