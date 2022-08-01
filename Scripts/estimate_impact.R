## Multimorbidity interpretation, early withdrawal
library(tidyverse)

ab <- expand.grid(prob0 = seq(0.05, 0.5, 0.1),
                  cc = 0:100,
                  com_mean = c(0.5, 0.75, 1, 1.25, 2)) %>% 
  as_tibble() %>% 
        mutate(
         # com_mean2 = com_mean * com_rr,
         odds0 = prob0/(1-prob0),
         or = 1.11^(cc+1),
         odds = odds0*or,
         probs = odds/(1+odds))

ab <- ab %>% 
  mutate(pmf = dpois(cc, com_mean))

ab_smry <- ab %>% 
  group_by(odds0, com_mean) %>% 
  summarise(sm_pmf = round(sum(pmf),3),
            ew_pmf = weighted.mean(probs, pmf)) %>% 
  ungroup()

## check complete probability distributions
all(ab_smry$sm_pmf ==1)


## spread to wide
ab_smry2 <- ab_smry %>% 
  select(-sm_pmf) %>% 
  mutate(ew_pmf = round(100* ew_pmf) %>% paste0("%"),
         com_ge2 = round(100* (1-ppois(2, com_mean)))%>% paste0("%"),
         withdrawal_base = 100*odds0/(1+odds0),
         cols = paste0(com_mean, " (", com_ge2, ")")) %>% 
  select(withdrawal_base, cols, ew_pmf)

ab_smry_wide <- ab_smry2 %>% 
  # rename("Initial" = com_mean, "Increased to" = com_mean2, "")
  pivot_wider(names_from = cols, values_from = ew_pmf)  
write_csv(ab_smry_wide, "Outputs/change_in_ew_trial.csv")


## individual level
ab <- expand.grid(prob0 = seq(0.05, 0.5, 0.1),
             cc = 0:4) %>% 
  as_tibble()
odds <- function(x) x/(1-x)
prbs <- function(x) x/(1+x)
ab <- ab %>% 
  mutate(prob1 = round(100*prbs(odds(prob0) * 1.1^cc)) %>%  paste0("%"),
         prob0 = round(100*prob0)) %>% 
  spread(cc, prob1)
ab

write_csv(ab, "Outputs/change_in_ew_indiv.csv")
