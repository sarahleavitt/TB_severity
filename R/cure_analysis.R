#Sarah V. Leavitt
#Boston University
#Pre-chemotherapy TB Analysis

##############################################################################
# This program performs the cure Bayesian meta-analysis and creates tables
# to summarize the results
##############################################################################

options(scipen=999)
options(digits = 10)
set.seed(150184)

rm(list = ls())
source("R/utils.R")
reload_source()

#Reading in the study_id correspondence table
studyid <- read.csv("data/study_id.csv")

#Reading in cure data
cureData <- read.csv("data/cure_data.csv")

cureData <- cureData %>%
  mutate(severity = factor(severity, levels = c("Far advanced disease", 
                                                "Moderately advanced disease",
                                                "Minimal disease"),
                           labels = c("Advanced", "Moderate", "Minimal")),
         severityAdv = ifelse(severity != "Advanced", "Min/Mod", severity))



#### Bayesian Logistic Model------------------------------------------------------------------------

#Model
m_cure <- function(){
  
  for (i in 1:nStudy){
    
    #Distribution of cure counts for each severity level
    cAdv[i] ~ dbin(p_0[i], nAdv[i])
    cMod[i] ~ dbin(p_1[i], nMod[i])
    cMin[i] ~ dbin(p_2[i], nMin[i])
    
    #Logit models
    logit(p_0[i]) <- u[i]
    logit(p_1[i]) <- u[i] + bmod
    logit(p_2[i]) <- u[i] + bmin
    
    #Study-level odds for each severity level
    oddsAdv[i] <- exp(u[i])
    oddsMod[i] <- exp(u[i] + bmod)
    oddsMin[i] <- exp(u[i] + bmin)
    
    #Random effects
    u[i] ~ dnorm(alpha, tau)
  }
  
  #Overall odds ratios
  ORmod <- exp(bmod)
  ORmin <- exp(bmin)
  
  #Priors
  bmod ~ dnorm(0, 0.001)
  bmin ~ dnorm(0, 0.001)
  alpha ~ dnorm(0,0.001)
  tau ~ dgamma(1,1)
  
  #Random effect variance
  theta = 1/tau
}

#Parameters to track
par_cure <- c("alpha", "bmod", "bmin", "theta",
              "ORmod", "ORmin", "oddsAdv", "oddsMod", "oddsMin")



#### Three year only analysis, all studies ---------------------------------------------------------

## Running the model
cureData_3 <- cureData %>% filter(study_id != "79_1023")

cureAggregate_3 <- cureData_3 %>%
  group_by(study_id) %>%
  summarize(nMin = sum(severity == "Minimal"),
            nMod = sum(severity == "Moderate"),
            nAdv = sum(severity == "Advanced"),
            cMin = sum(severity == "Minimal" & cure == 1),
            cMod = sum(severity == "Moderate" & cure == 1),
            cAdv = sum(severity == "Advanced" & cure == 1),
            .groups = "drop")

#Data
dt_3 <- list(nStudy = length(unique(cureData_3$study_id)),
             nMin = cureAggregate_3$nMin,
             nMod = cureAggregate_3$nMod,
             nAdv = cureAggregate_3$nAdv,
             cMin = cureAggregate_3$cMin,
             cMod = cureAggregate_3$cMod,
             cAdv = cureAggregate_3$cAdv
)

#Fitting the model
fit_3 <- jags(data = dt_3, model.file = m_cure,
                parameters.to.save = par_cure,
                n.iter = 61000, n.burnin = 1000,
                n.chains = 1, n.thin = 30)

#Extracting data
mcmc_3 <- as.mcmc(fit_3)
eval_3 <- mcmc_3

png("Figures/xyplot_cure_3year.png")
xyplot(eval_3[, c("alpha", "bmod", "bmin", "theta")])
dev.off()
png("Figures/autocorr_cure_3year.png")
autocorr.plot(eval_3[, c("alpha", "bmod", "bmin", "theta")])
dev.off()



#### Running the Analysis -----------------------------------------------------

cureData_3us <- cureData %>% filter(study_id %in% c("1029", "45"))

cureAggregate_3us <- cureData_3us %>%
  group_by(study_id) %>%
  summarize(nMin = sum(severity == "Minimal"),
            nMod = sum(severity == "Moderate"),
            nAdv = sum(severity == "Advanced"),
            cMin = sum(severity == "Minimal" & cure == 1),
            cMod = sum(severity == "Moderate" & cure == 1),
            cAdv = sum(severity == "Advanced" & cure == 1),
            .groups = "drop")

#Data
dt_3us <- list(nStudy = length(unique(cureData_3us$study_id)),
             nMin = cureAggregate_3us$nMin,
             nMod = cureAggregate_3us$nMod,
             nAdv = cureAggregate_3us$nAdv,
             cMin = cureAggregate_3us$cMin,
             cMod = cureAggregate_3us$cMod,
             cAdv = cureAggregate_3us$cAdv
)

#Fitting the model
fit_3us <- jags(data = dt_3us, model.file = m_cure,
              parameters.to.save = par_cure,
              n.iter = 61000, n.burnin = 1000,
              n.chains = 1, n.thin = 30)

#Extracting data
mcmc_3us <- as.mcmc(fit_3us)
eval_3us <- mcmc_3us

png("Figures/xyplot_cure.png")
xyplot(eval_3us[, c("alpha", "bmod", "bmin", "theta")])
dev.off()
png("Figures/autocorr_cure.png")
autocorr.plot(eval_3us[, c("alpha", "bmod", "bmin", "theta")])
dev.off()



#### Tables of results -----------------------------------------------------------------------------

## Three year only, US post-1930s
tabs <- make_cure_tab(eval_3us, cureAggregate_3us)
or <- tabs[[1]]
cureTab <- tabs[[2]]


