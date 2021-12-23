######################################################################
## Title: [Insert full title of paper]
## Short title: [Insert short title of paper (if applicable) here]
## DOI: [Insert DOI of paper here]
## Code author: Steven Kerr steven.kerr@ed.ac.uk
## Description: Descriptives for omicron vaccine effectiveness paper
######################################################################

# 0 Setup ####
#Libraries
library(tidyverse)
library(lubridate)

Location = '/conf/'

setwd(paste0(Location, 'EAVE/GPanalysis/analyses/B.1.1.529-variant'))


# Variant data
wgs <- readRDS(paste0(Location,"EAVE/GPanalysis/data/WGS_latest.rds")) %>% 
  mutate_at(c("Collection_Date","Sequencing_Date","Alignment_Date"), ~ as.Date(. , format="%d/%m/%Y")) %>%
  rename(specimen_date = Collection_Date) %>%
  arrange(EAVE_LINKNO, specimen_date) %>% 
  filter(!duplicated(EAVE_LINKNO))

first_omicron_date <- wgs %>% filter(VariantShorthand == 'B.1.1.529') %>%
                      pull(specimen_date) %>%
                      min()


# All cases starting 1 November 2021
df <- readRDS('/conf/EAVE/GPanalysis/progs/CR/SGene_Omicron/output/temp/Sgene_all_positive.rds')


s_neg <- filter(df, s_gene == 'S_Neg') 

number_s_neg <- s_neg %>% 
                pull(EAVE_LINKNO) %>%
                unique() %>%
                length()


s_neg_count <- s_neg %>%
               group_by(date_ecoss_specimen) %>%
               summarise(n = n())

a_begin <- as.Date("2021-11-01")

a_end <- max(s_neg_count$date_ecoss_specimen)





grid <- expand.grid(seq(a_begin, a_end-1, by="days") )

names(grid) <- 'date_ecoss_specimen'

grid <- grid %>%
  left_join(s_neg_count) %>%
  replace_na( list(n = 0)) 

grid %>%  ggplot(aes(x=date_ecoss_specimen, y=n)) + geom_line(color = 'blue') +
  # Delay of ~2 weeks in sequencing data, so truncate smooth 2 weeks before end
  #geom_smooth(xseq = smooth_start:smooth_end) +
  labs(x="Specimen date",y ="Number", title="S-gene negative tests by day") + 
  theme(legend.title = element_blank()) 

ggsave(paste0("./output/", a_begin, "-", a_end, "_pos_tests_by_day.png"), width=14, height=10, unit="cm")



s_neg_count <- mutate(s_neg_count, ln_cum = log(cumsum(n)),
               days = as.numeric( date_ecoss_specimen - as.Date("2021-11-01")))

model <- lm(ln_cum ~ days, data = s_neg_count)

coefs <- summary(model)$coefficients

doubling_time <- log(2)/coefs['days', 'Estimate']

doubling_time_lcl <- log(2)/( coefs['days', 'Estimate'] - 1.96 * coefs['days', 'Std. Error'] )
doubling_time_ucl <- log(2)/( coefs['days', 'Estimate'] + 1.96 * coefs['days', 'Std. Error'] )







