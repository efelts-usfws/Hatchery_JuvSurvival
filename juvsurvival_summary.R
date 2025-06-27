

library(tidyverse)
library(arrow)

# read in nimble parameter estimates

nimble25.parms <- read_rds("outputs/parameter_estimates_2025")

nimble25.chains <- read_rds("outputs/chains_2025")
