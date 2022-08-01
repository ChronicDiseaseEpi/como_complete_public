library(tidyverse)
library(brms)
library(tidybayes)

## Read aggregatte-level data
  cfs_csdr <- read_csv("Data/cfs_csdr.csv")
  cfs_yoda <- read_csv("Data/cfs_yoda.csv")
  cfs_comp <- read_csv("Data/numbers_completing.csv")
  linear_models_coefficients_csdr <- read_csv("Data/linear_model_coefficients_csdr.csv")

res <- linear_models_coefficients_csdr 

#### Prepare datasets ----

## join cfs together
cfs_csdr <- cfs_csdr %>%
  select(-statistic, -p.value)

cfs_yoda <- cfs_yoda %>% 
  filter(nct_id != "NCT00168857" & nct_id != "NCT00034840") %>%
  rename(std.error = se) %>% 
  select(-outcome)

cfs <- bind_rows(csdr = cfs_csdr,
                 yoda = cfs_yoda,
                 .id = "repo")

## should be able to use relative paths here
metadata <- read_csv("Supporting/metadata.csv")

## trials need to read in metadata for; condition, treatment 
cfs %>% 
  anti_join(metadata) %>% 
  distinct(nct_id) # None

cfs <- cfs %>% 
  inner_join(metadata) %>%
  mutate(atc5 = str_sub(code, 1, 5))

cfs <- left_join(cfs, cfs_comp, by = "nct_id") %>%
  rename(no_comp = 'Completed', no_term = 'Early termination', prop_term = 'Proportion') %>%
  mutate(no_participants = no_comp + no_term) %>%
  select(-comment) %>%
  filter(nct_id != "NCT00696241" & nct_id != "NCT00696436" & nct_id != "NCT00034840" & nct_id != "NCT00168857", nct_id != "NCT01087762" & 
           nct_id != "NCT01087788" & nct_id != "NCT01032629" & nct_id != "NCT01809327" & nct_id != "NCT01381900" & nct_id != "NCT01989754")
## again, should use relative paths
saveRDS(cfs, "Outputs/cfs.rds")

cfs <- readRDS("Outputs/cfs.rds")

cfs_como = cfs %>% filter(term == "como_cnt" & model_n == 1 & std.error < 5)
cfs_age = cfs %>% filter(term == "age_std" & model_n == 5 & std.error < 5)
cfs_sex = cfs %>% filter((term == "sex" | term == "sexTRUE") & model_n == 5 & std.error < 5)

#### Bayesian hierarchical models ####

# comorbidity 
if(!file.exists("Scratch_data/pooled_model.Rds")){
  mod0_brm_como <- brm(estimate|se(std.error) ~ 1 + (1|nct_id), 
                       data =  compose_data(cfs_como),
                       prior = prior(class = "Intercept", student_t(3, 0, 100)),
                       save_pars = save_pars(all = TRUE), cores = 4,
                       control=list(adapt_delta=0.99))
  mod0_brm_como <- add_criterion(mod0_brm_como, c("loo"), moment_match = TRUE)
  loo1 <- loo(mod0_brm_como)
  reloo1 <- reloo(loo1, mod0_brm_como, chains = 1)
  reloo1
  saveRDS(list(mdl = mod0_brm_como, loo = reloo1), "Scratch_data/pooled_model.Rds")
} 
if(!file.exists("Scratch_data/condition_model.Rds")){
  
mod1_brm_como <- brm(estimate|se(std.error) ~ 1 + (1|nct_id) + (1|sail),
                     data = compose_data(cfs_como),
                     prior = prior(class = "Intercept", student_t(3, 0, 100)),
                     save_pars = save_pars(all = TRUE), cores = 4,
                     control=list(adapt_delta=0.99))
loo1 <- loo(mod1_brm_como, moment_match = TRUE)
reloo1 <- reloo(loo1, mod1_brm_como, chains = 1)
reloo1
saveRDS(list(mdl = mod1_brm_como, loo = reloo1), "Scratch_data/condition_model.Rds")
}

## comorbidity and atc5
if(!file.exists("Scratch_data/drugclass_model.Rds")) {
mod2_brm_como <- brm(estimate|se(std.error) ~ 1 + (1|nct_id) + (1|atc5), 
                     data =  compose_data(cfs_como),
                     prior = prior(class = "Intercept", student_t(3, 0, 100)),                     
                     save_pars = save_pars(all = TRUE), cores = 4,
                     control=list(adapt_delta=0.99))
# mod2_brm_como <- add_criterion(mod2_brm_como, c("loo", "waic"), moment_match = TRUE)
loo1 <- loo(mod2_brm_como, moment_match = TRUE)
reloo1 <- reloo(loo1, mod2_brm_como, chains = 1)
reloo1
saveRDS(list(mdl = mod2_brm_como, loo = reloo1), "Scratch_data/drugclass_model.Rds")
}

if(!file.exists("Scratch_data/drugclasscondition_model.Rds")) {
mod3_brm_como <- brm(estimate|se(std.error) ~ 1 + (1|nct_id) + (1|sail) + (1|atc5), 
                     data =  compose_data(cfs_como),
                     prior = prior(class = "Intercept", student_t(3, 0, 100)),                     
                     save_pars = save_pars(all = TRUE), cores = 4, 
                     control=list(adapt_delta=0.99))
# mod3_brm_como <- add_criterion(mod3_brm_como, c("loo", "waic"), moment_match = TRUE)
loo1 <- loo(mod3_brm_como, moment_match = TRUE)
reloo1 <- reloo(loo1, mod3_brm_como, chains = 1)
reloo1
saveRDS(list(mdl = mod3_brm_como, loo = reloo1), "Scratch_data/drugclasscondition_model.Rds")
}


## re-run with different priors
if(!file.exists("Scratch_data/drugclasscondition_model_altpriors.Rds")) {
  mod3_brm_como <- brm(estimate|se(std.error) ~ 1 + (1|nct_id) + (1|sail) + (1|atc5), 
                       data =  compose_data(cfs_como),
                       prior = c(prior(class = "Intercept", student_t(3, 0, 100)),
                                 prior(class = "sd", student_t(3, 0, 10))),
                       save_pars = save_pars(all = TRUE), cores = 4, 
                       control=list(adapt_delta=0.99))
  # mod3_brm_como <- add_criterion(mod3_brm_como, c("loo", "waic"), moment_match = TRUE)
  loo1 <- loo(mod3_brm_como, moment_match = TRUE)
  reloo1 <- reloo(loo1, mod3_brm_como, chains = 1)
  reloo1
  saveRDS(list(mdl = mod3_brm_como, loo = reloo1), "Scratch_data/drugclasscondition_model_altpriors.Rds")
}


if(!file.exists("Scratch_data/pooled_model_age.Rds")) {
mod0_brm_age <- brm(estimate|se(std.error) ~ 1 + (1|nct_id) , 
                    data =  compose_data(cfs_age), 
                    prior = prior(class = "Intercept", student_t(3, 0, 100)),    
                    save_pars = save_pars(all = TRUE), cores =4) # control=list(adapt_delta=0.99, max_treedepth=13))
# mod0_brm_age <- add_criterion(mod0_brm_age, c("loo", "waic"), moment_match = TRUE)

get_variables(mod0_brm_age)
summary(mod0_brm_age)
saveRDS(mod0_brm_age, "Scratch_data/pooled_model_age.Rds")
}

# age and disease condition
if(!file.exists("Scratch_data/condition_model_age.Rds")) {
mod1_brm_age <- brm(estimate|se(std.error) ~ 1 + (1|nct_id) + (1|sail), 
                    data =  compose_data(cfs_age), 
                    prior = prior(class = "Intercept", student_t(3, 0, 100)),  
                    save_pars = save_pars(all = TRUE), cores =4,
                    control=list(adapt_delta=0.99, max_treedepth=13))
# mod1_brm_age <- add_criterion(mod1_brm_age, c("loo", "waic"), moment_match = TRUE)
saveRDS(mod1_brm_age, "Scratch_data/condition_model_age.Rds")
}


# age and atc5
if(!file.exists("Scratch_data/drugclass_model_age.Rds")){
mod2_brm_age <- brm(estimate|se(std.error) ~ 1 + (1|nct_id) + (1|atc5), 
                    data =  compose_data(cfs_age), 
                    prior = prior(class = "Intercept", student_t(3, 0, 100)),  
                    save_pars = save_pars(all = TRUE), cores =4,
                    control=list(adapt_delta=0.99, max_treedepth=13))
# mod2_brm_age <- add_criterion(mod2_brm_age, c("loo", "waic"), moment_match = TRUE)
saveRDS(mod2_brm_age, "Scratch_data/drugclass_model_age.Rds")
}


if(!file.exists("Scratch_data/drugclass_condition_model_age.Rds")){
mod3_brm_age <- brm(estimate|se(std.error) ~ 1 + (1|nct_id) + (1|sail) + (1|atc5), 
                    data =  compose_data(cfs_age), 
                    prior = prior(class = "Intercept", student_t(3, 0, 100)),  
                    save_pars = save_pars(all = TRUE), cores =4,
                    control=list(adapt_delta=0.99, max_treedepth=13))
# mod3_brm_age <- add_criterion(mod3_brm_age, c("loo", "waic"), moment_match = TRUE)
saveRDS(mod3_brm_age, "Scratch_data/drugclasscondition_model_age.Rds")
}


# sex
if(!file.exists("Scratch_data/pooled_model_sex.Rds")){
mod0_brm_sex <- brm(estimate|se(std.error) ~ 1 + (1|nct_id) , 
                    data =  compose_data(cfs_sex), 
                    prior = prior(class = "Intercept", student_t(3, 0, 100)),  
                    save_pars = save_pars(all = TRUE), cores =4,
                    control=list(adapt_delta=0.99, max_treedepth=13)) 
# mod0_brm_sex <- add_criterion(mod0_brm_sex, c("loo", "waic"), moment_match = TRUE)
saveRDS(mod0_brm_sex, "Scratch_data/pooled_model_sex.Rds")
}

# sex and condition
if(!file.exists("Scratch_data/condition_model_sex.Rds")){
mod1_brm_sex <- brm(estimate|se(std.error) ~ 1 + (1|nct_id) + (1|sail) , 
                    data =  compose_data(cfs_sex), 
                    prior = prior(class = "Intercept", student_t(3, 0, 100)),  
                    save_pars = save_pars(all = TRUE), cores =4,
                    control=list(adapt_delta=0.99, max_treedepth=13)) # no divergent transitions
# mod1_brm_sex <- add_criterion(mod1_brm_sex, c("loo", "waic"), moment_match = TRUE)
saveRDS(mod1_brm_sex, "Scratch_data/condition_model_sex.Rds")
}


# sex and atc5
if(!file.exists("Scratch_data/drugclass_model_sex.Rds")){
  
mod2_brm_sex <- brm(estimate|se(std.error) ~ 1 + (1|nct_id) + (1|atc5), 
                    data =  compose_data(cfs_sex), 
                    prior = prior(class = "Intercept", student_t(3, 0, 100)),  
                    save_pars = save_pars(all = TRUE), cores =4,
                    control=list(adapt_delta=0.99, max_treedepth=13)) # no divergent transitions
# mod2_brm_sex <- add_criterion(mod2_brm_sex, c("loo", "waic"), moment_match = TRUE)
saveRDS(mod2_brm_sex, "Scratch_data/drugclass_model_sex.Rds")
}



# sex, condition and atc5
if(!file.exists("Scratch_data/drugclasscondition_model_sex.Rds")){
mod3_brm_sex <- brm(estimate|se(std.error) ~ 1 + (1|nct_id) + (1|sail) + (1|atc5), 
                    data =  compose_data(cfs_sex), 
                    prior = prior(class = "Intercept", student_t(3, 0, 100)),  
                    save_pars = save_pars(all = TRUE), cores =4,
                    control=list(adapt_delta=0.99, max_treedepth=13))
# mod3_brm_sex <- add_criterion(mod3_brm_sex, c("loo", "waic"), moment_match = TRUE)
saveRDS(mod3_brm_sex, "Scratch_data/drugclasscondition_model_sex.Rds")

}


#### Association between enrolment order and comorbidity count for CSDR trials ####
res %>% filter(term == "(Intercept)") %>% n_distinct() # 63
if(!file.exists("Scratch_data/pooled_model_enrol.Rds")){
mod_como_enrol0 <- brm(estimate|se(std.error) ~ 1 + (1|nct_id) , 
                       data =  compose_data(res %>% filter(term == "como_cnt")), 
                       prior = prior(class = "Intercept", student_t(3, 0, 100)),  
                       save_pars = save_pars(all = TRUE), cores =4) 
saveRDS(mod_como_enrol0, "Scratch_data/pooled_model_enrol.Rds")
}
if(!file.exists("Scratch_data/pooled_model_enrol_age.Rds")){
mod_age_enrol0 <- brm(estimate|se(std.error) ~ 1 + (1|nct_id) , 
                       data =  compose_data(res %>% filter(term == "age_std")), 
                       prior = prior(class = "Intercept", student_t(3, 0, 100)),  
                       save_pars = save_pars(all = TRUE), cores =4) 
saveRDS(mod_age_enrol0, "Scratch_data/pooled_model_enrol_age.Rds")
}
if(!file.exists("Scratch_data/pooled_model_enrol_sex.Rds")){
mod_sex_enrol0 <- brm(estimate|se(std.error) ~ 1 + (1|nct_id) , 
                      data =  compose_data(res %>% filter(term == "sex")), 
                      prior = prior(class = "Intercept", student_t(3, 0, 100)),  
                      save_pars = save_pars(all = TRUE), cores =4) 
saveRDS(mod_sex_enrol0, "Scratch_data/pooled_model_enrol_sex.Rds")
}


#### Association between 1 additional comorbidity count and early withdrawal ####
cfs_como2 = cfs %>% 
  filter(model_n == 4) %>%
  filter(term == "I(como_cnt^2)") %>%
  mutate(como2 = 2*estimate,
         se2 = 2*std.error)
if(!file.exists("Scratch_data/pooled_model_como2.Rds")) {
mod0_como2 <- brm(estimate|se(std.error) ~ 1 + (1|nct_id), 
                       data =  compose_data(cfs_como2),
                       prior = prior(class = "Intercept", student_t(3, 0, 100)),
                       save_pars = save_pars(all = TRUE), cores =4,
                  control=list(adapt_delta=0.99, max_treedepth=13))
# mod0_como2 <- add_criterion(mod0_brm_como, c("loo", "waic"), moment_match = TRUE)
saveRDS(mod0_como2, "Scratch_data/pooled_model_como2.Rds")
}

## by index condition
if(!file.exists("Scratch_data/condition_model_como2.Rds")) {
mod1_como2 <- brm(estimate|se(std.error) ~ 1 + (1|nct_id) + (1|sail), 
                  data =  compose_data(cfs_como2),
                  prior = prior(class = "Intercept", student_t(3, 0, 100)),
                  save_pars = save_pars(all = TRUE), cores =4,
                  control=list(adapt_delta=0.99, max_treedepth=13))
# mod1_como2 <- add_criterion(mod0_brm_como, c("loo", "waic"), moment_match = TRUE)
saveRDS(mod1_como2, "Scratch_data/condition_model_como2.Rds")
}

## by atc drug class
if(!file.exists("Scratch_data/drugclass_model_como2.Rds")) {
mod2_como2 <- brm(estimate|se(std.error) ~ 1 + (1|nct_id) + (1|atc5), 
                  data =  compose_data(cfs_como2),
                  prior = prior(class = "Intercept", student_t(3, 0, 100)),
                  save_pars = save_pars(all = TRUE), cores =4) # control=list(adapt_delta=0.99, max_treedepth=13))
# mod2_como2 <- add_criterion(mod0_brm_como, c("loo", "waic"), moment_match = TRUE)
saveRDS(mod2_como2, "Scratch_data/drugclass_model_como2.Rds")
}


