

df_all <- readRDS(paste0("./output/temp/df_all_full_covid_pos.RDS"))

#vaccine uptake by vaccine status and age
z_levs = levels(df_all$vs)
z <- df_all %>% 
  mutate(age_gp = cut(ageYear, breaks = c(-1, 15, 39, 64, 120), labels=c("0-15","16-39","40-64","65+" ))) %>% 
  group_by(age_gp,vs) %>% 
  dplyr::summarise(N=round(sum(eave_weight))) %>%  
  mutate(percent=N/sum(N)*100) %>% 
  mutate(vs=factor(vs,levels=z_levs))

z %>% ggplot(aes(x=vs, y=percent, fill=age_gp)) + geom_col(position="dodge") +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  labs(x="Vaccine Status", y="Percentage", fill="Age\nGroup")
ggsave(paste0("./output/vs_uptake_age.png"), width=14, height=10, unit="cm")


#rates of Spos and Sneg by vaccine status
z_levs = levels(df_all$vs)
z <- df_all %>%  group_by(vs) %>% 
  dplyr::summarise(N=round(sum(eave_weight)), spos=sum(event_s_pos), sneg = sum(event_s_neg)) %>% 
  pivot_longer(cols=spos:sneg) %>% 
  mutate(name=case_when(name=="spos" ~ "S Positive",
                        name=="sneg" ~ "S Negative")) %>% 
  mutate(rate=value/N*100000) %>% 
  mutate(vs=factor(vs,levels=z_levs))

z %>% ggplot(aes(x=vs, y=rate)) + geom_col() + facet_wrap(~name, scales="free_y") +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  labs(x="Vaccine Status", y="Rate per 100,000")
ggsave(paste0("./output/sgene_vs_rate.png"), width=14, height=10, unit="cm")

#rates of Spos and Sneg by age_gp
z_levs = levels(df_all$age_gp)
z <- df_all %>%  group_by(age_gp) %>% 
  dplyr::summarise(N=round(sum(eave_weight)), spos=sum(event_s_pos), sneg = sum(event_s_neg)) %>% 
  pivot_longer(cols=spos:sneg) %>% 
  mutate(name=case_when(name=="spos" ~ "S Positive",
                        name=="sneg" ~ "S Negative")) %>% 
  mutate(rate=value/N*100000) %>% 
  mutate(age_gp=factor(age_gp,levels=z_levs))

z %>% ggplot(aes(x=age_gp, y=rate)) + geom_col() + facet_wrap(~name, scales="free_y") +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  labs(x="Age Group", y="Rate per 100,000")
ggsave(paste0("./output/sgene_age_rate.png"), width=14, height=10, unit="cm")


z <- df_all %>% filter(vs=="uv")  %>% group_by(simd2020_sc_quintile, pos_before_start) %>% 
  dplyr::summarise(N=round(sum(eave_weight)), spos=sum(event_s_pos), sneg = sum(event_s_neg)) %>% 
  pivot_longer(cols=spos:sneg) %>% 
  mutate(rate=value/N*100000) 

z %>% ggplot(aes(x=simd2020_sc_quintile, y=rate, fill=pos_before_start)) + geom_col(position="dodge") + facet_wrap(~name, scales="free_y") +
  theme(axis.text.x=element_text(angle=45,hjust=1))


z <- df_all %>% filter(vs=="uv")  %>% mutate(age_gp =cut(ageYear, breaks=c(-1,11,19,39,59,120), labels=c("0-11","12-19","20-39","40-59","60+"))) %>% 
  group_by(age_gp, pos_before_start) %>% 
  dplyr::summarise(N=round(sum(eave_weight)), spos=sum(event_s_pos), sneg = sum(event_s_neg)) %>% 
  pivot_longer(cols=spos:sneg) %>% 
  mutate(rate=value/N*100000) 

z %>% ggplot(aes(x=age_gp, y=rate, fill=pos_before_start)) + geom_col(position="dodge") + facet_wrap(~name, scales="free_y") +
  theme(axis.text.x=element_text(angle=45,hjust=1))

