#Sarah V. Leavitt
#Boston University
#Pre-chemotherapy TB Severity Analysis

###########################################################################################
# This program performs the mortality Bayesian meta-analysis for the US post-1930s studies
# It takes about 15 to run
###########################################################################################

options(scipen=999)
options(digits = 10)
set.seed(150183)

rm(list = ls())
source("R/utils.R")
reload_source()
source("R/mortality_functions.R")

#Reading in individual mortality data
mortality <- read.csv("data/mortality_data.csv")


#### Set-up ----------------------------------------------------------------------------------------

#Creating interval variables
#interval = 1 implies x1 < t <= x2
#dinterval() isn't working as it should so make all obs have interval = 1
#and set x2 to be 10000 (close enough to infinity) for right censored
mortality <- mortality %>% 
  #Removing people who are censored at the study start
  filter(no_survival == FALSE) %>%
  mutate(time = ifelse(death_tb == 0, NA, interval_l),
         interval = 1,
         x1 = interval_l,
         x2 = ifelse(death_tb == 0, 10000, interval_r))

n.iter <- 31000
n.burnin <- 1000
n.thin <- 30

# n.iter <- 100
# n.burnin <- 10
# n.thin <- 1



#### Running Models --------------------------------------------------------------------------------


#### US Studies: Complete Model ####

#Subsetting and formatting data
mortality_comp <- mortality %>%
  filter(study_id %in% c("1029", "93", "45")) %>%
  mutate(study_sev_num = as.numeric(factor(study_sev)),
         study_id_num = as.numeric(factor(study_id)))

#Running model
output_comp <- run_comp(mortality_comp, 
                        n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin)



#### US Studies: Stratified Model ####

#Subsetting and formatting data
mortality_sev <- mortality %>%
  filter(study_id %in% c("1029", "93", "45"),
         severity != "Unknown") %>%
  mutate(study_sev_num = as.numeric(factor(study_sev)),
         study_id_num = as.numeric(factor(study_id)),
         sev_mod = as.numeric(severity == "Moderate"),
         sev_adv = as.numeric(severity == "Advanced"),
         sev_unk = as.numeric(severity == "Unknown"))

#Running model
output_sev <- run_sev(mortality_sev,
                      n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin)


#### Saving Results------------------------------------------------------------------------------

data_comp <- getData(mortality_comp)
data_sev <- getData(mortality_sev)
res_comp <- output_comp$res
res_sev <- output_sev$res
eval_comp <- output_comp$eval
eval_sev <- output_sev$eval

save(res_comp, res_sev, eval_comp, eval_sev, data_comp, data_sev,
     file = "R/bayesian_severity.RData")


