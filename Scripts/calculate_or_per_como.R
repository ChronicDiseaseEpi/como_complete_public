a <- c(1.11, 1.07, 1.14)
b <- lapply(0:3, function(x) a^x)
library(tidyverse)
d <- map(b, ~ round(.x, digits = 2))
e <- map_chr(d, ~ paste0(.x[1], ", 95% CI ", .x[2], "-", .x[3], ""))
paste(e, collapse = ", ") 


ew <- 0.1/0.9

odds <- map(b, ~ .x *ew)
prbs <- map(odds, ~ 100*.x/(1+.x))
prbs <- map(prbs, round, digits = 1)
p_fin <- map_chr(prbs, ~ paste0(.x[1], "%, 95% CI ", .x[2], "-", .x[3], "%"))
paste(p_fin, collapse = ", ") 


z <- c(est = 0.1197, se = 0.0350, p = 0.0006, ci.lb = 0.0510, cli.ub = 0.1884)
z2 <- exp(z) %>% 
  round(2)
z3 <- z2[c(1, 4, 5)]
z3

