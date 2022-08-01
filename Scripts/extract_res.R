## extract SD

# calculate bayesian p-values
library(brms)
library(tidyverse)
# library(tidybayes)
# library(stringr)
# library(purrr)

## read in saved models
mdls <- c('condition_model.Rds', 
          'condition_model_age.Rds', 
          'condition_model_sex.Rds', 
          'condition_model_como2.Rds', 
          'drugclass_model.Rds', 
          'drugclass_model_age.Rds', 
          'drugclass_model_sex.Rds', 
          'drugclass_model_como2.Rds', 
          'drugclasscondition_model.Rds', 
          'drugclasscondition_model_altpriors.Rds',
          'drugclasscondition_model_age.Rds', 
          'drugclasscondition_model_sex.Rds', 
          'pooled_model.Rds', 
          'pooled_model_age.Rds', 
          'pooled_model_como2.Rds', 
          'pooled_model_enrol.Rds', 
          'pooled_model_enrol_age.Rds', 
          'pooled_model_enrol_sex.Rds', 
          'pooled_model_sex.Rds')
names(mdls) <- str_sub(mdls, 1, -5)
mdls <- intersect(mdls, list.files("Scratch_data/"))
rv <- map(mdls, ~ readRDS(paste0("Scratch_data/", .x)))

## select model object where have model and loo object inside list
rv <- map(rv, ~ {
  if("mdl" %in% names(.x)) .x$mdl else .x
})
names(rv) <-  str_sub(mdls, 1, -5)

## for effects
getres <- map(rv, ~ {
  a <- summary(.x)
  a$fixed %>% 
    as_tibble() %>% 
    mutate(measure = "effect") %>% 
    select(measure, everything())
})
getres <- bind_rows(getres, .id = "mdl")
getres %>% 
  filter(mdl %>% str_detect("model$")) 

getres2 <- getres %>% 
  mutate(across(Estimate:`u-95% CI`, ~ .x %>%  exp() %>% round(2))) %>% 
  mutate(res = paste0(Estimate, " (", `l-95% CI`, "-", `u-95% CI`, ")"))

getres3 <- getres2 %>% 
  select(mdl, res) %>% 
  separate(mdl, into = c("grps", "model", "var", "other"), fill = "right") %>% 
  mutate(var = if_else(!is.na(other), paste0(var, "_", other), var),
         var = if_else(is.na(var), "como", var)) %>% 
  select(grps, var, res)  %>% 
  filter(!str_detect(var, "enrol")) %>% 
  spread(var, res) %>% 
  select(grps, como, age, sex) %>% 
  arrange(factor(grps, levels = c("pooled", "condition", "drugclass", "drugclasscondition")))
write.csv(getres3, "Outputs/res.csv")

## compare estimates where use wider and standard priors
getres4 <- getres %>% 
  filter(mdl %in% c("drugclasscondition_model_altpriors", "drugclasscondition_model")) %>% 
  mutate(priors = if_else(mdl == "drugclasscondition_model_altpriors", "widersd", "main")) %>% 
  select(priors, Estimate, Est.Error) %>% 
  mutate(across(c(Estimate, Est.Error), ~ round(.x, 3))) %>% 
  mutate(res = paste0(Estimate, " (", Est.Error, ")")) %>% 
  select(priors, res ) %>% 
  spread(priors, res)
write.csv(getres4, "Outputs/compare_res_wider_priors.csv")


