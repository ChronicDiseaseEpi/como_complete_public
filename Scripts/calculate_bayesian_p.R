# calculate bayesian p-values

library(brms)
library(tidyverse)
library(tidybayes)
library(stringr)
library(purrr)

## read in saved models
rv <- map(c("condition_model.Rds", "drugclass_model.Rds", "drugclasscondition_model.Rds", "pooled_model.Rds"), ~ readRDS(paste0("Scratch_data/", .x)))
names(rv) <- c("condition", "drugclass", "drugclasscondition", 
               "pooled")

mdls <- c('condition_model.Rds', 
'condition_model_age.Rds', 
'condition_model_sex.Rds', 
'condition_model_como2.Rds', 
'drugclass_model.Rds', 
'drugclass_model_age.Rds', 
'drugclass_model_sex.Rds', 
'drugclass_model_como2.Rds', 
'drugclasscondition_model.Rds', 
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

rv <- map(mdls, ~ readRDS(paste0("Scratch_data/", .x)))
## select model object where have model and loo object inside list
rv <- map(rv, ~ {
  if("mdl" %in% names(.x)) .x$mdl else .x
})

## for P-values for all models
getp <- map_dbl(rv, ~ {
  a <- posterior_samples(.x, pars = "b_Intercept")
  mean(a$b_Intercept > log(1))
})
getp
write.csv(tibble(model = names(getp), pgtor1 = getp), "Outputs/bayesianP.csv")

xmn <- posterior_predict(rv$pooled$mdl, newdata = data.frame(nct_id = 100, std.error = 100), allow_new_levels = TRUE)
