#Sarah V. Leavitt
#Boston University
#Pre-chemotherapy TB Severity Analysis

##############################################################################
# This program creates figures and tables of the severity mortality analysis
# with the results
##############################################################################

options(scipen=999)
options(digits = 10)

rm(list = ls())
source("R/utils.R")
reload_source()

#Reading in the study_id correspondence table
studyid <- read.csv("data/study_id.csv")

#Reading in individual mortality data and analysis results
mortality <- read.csv("data/mortality_data.csv")
load('R/bayesian_severity.RData')


#### Formatting results ----------------------------------------------------------------------------

form_comp <- formatBayesian(mortality, res_comp, data_comp, "Combined")
form_sev <- formatBayesian(mortality, res_sev, data_sev, "Severity", fixed = TRUE)


#### Diagnostic Plots ------------------------------------------------------------------------------

png("Figures/xyplot_comp.png")
xyplot(eval_comp)
dev.off()
png("Figures/autocorr_comp.png")
autocorr.plot(eval_comp)
dev.off()

png("Figures/xyplot_sev.png")
xyplot(eval_sev)
dev.off()
png("Figures/autocorr_sev.png")
autocorr.plot(eval_sev)
dev.off()


#### Table of Main Results -------------------------------------------------------------------------

#Combining raw tables from results lists
raw_tab <- bind_rows(form_comp$param, form_sev$param)

#Adding formatted, ordered labels
raw_tab <- raw_tab %>%
  mutate(Severity = ifelse(is.na(severity), "Combined", severity),
         Severity = factor(Severity, level = c("Minimal", "Moderate",
                                               "Advanced", "Combined"))) %>%
  arrange(Severity)

#Extracting the survival probabilities
pred1_tab <- raw_tab %>%
  filter(value == "pred1") %>%
  mutate(`1-Year Survival (95% CI)` = paste0(round(est, 2), " (",
                                             round(cilb, 2), ", ",
                                             round(ciub, 2), ")")) %>%
  select(Severity, `1-Year Survival (95% CI)`)

pred5_tab <- raw_tab %>%
  filter(value == "pred5") %>%
  mutate(`5-Year Survival (95% CI)` = paste0(round(est, 2), " (",
                                             round(cilb, 2), ", ",
                                             round(ciub, 2), ")")) %>%
  select(Severity, `5-Year Survival (95% CI)`)

pred10_tab <- raw_tab %>%
  filter(value == "pred10") %>%
  mutate(`10-Year Survival (95% CI)` = paste0(round(est, 2), " (",
                                              round(cilb, 2), ", ",
                                              round(ciub, 2), ")")) %>%
  select(Severity, `10-Year Survival (95% CI)`)

med_tab <- raw_tab %>%
  filter(value == "median") %>%
  mutate(`Median Survival Time (95% CI)` = paste0(round(est, 2), " (",
                                              round(cilb, 2), ", ",
                                              round(ciub, 2), ")")) %>%
  select(Severity, `Median Survival Time (95% CI)`)


#Extracting the distribution parameters
sdlog <- raw_tab %>%
  filter(value == "sdlog") %>%
  select(label, sdlog = est)

dist_tab <- raw_tab %>%
  filter(value == "meanlog") %>%
  full_join(sdlog, by = "label") %>%
  mutate(`Survival Distribution` = paste0("lognormal(", round(est, 2), ", ",
                                          round(sdlog, 2), ")")) %>%
  select(Severity, `Survival Distribution`)


#Finding number of papers, cohorts, individual for each analysis
data_comp2 <- as.data.frame(t(data_comp[[2]]))
data_comp2$Severity <- "Combined"

counts_comp <- data_comp2

#Counts stratified by severity for US post-1930s subset
mortality_sub <- mortality %>% filter(severity != "Unknown",
                                      study_id %in% c("1029", "93", "45"))

counts <- mortality_sub %>%
  group_by(severity) %>%
  summarize(nStudies = length(unique(study_id)),
            nCohorts = length(unique(cohort_id)),
            nIndividuals = n(),
            .groups = "drop") %>%
  mutate(Severity = paste0(severity, ""))


#Combining the tables
final_tab <- dist_tab %>%
  full_join(pred1_tab, by = c("Severity")) %>%
  full_join(pred5_tab, by = c("Severity")) %>%
  full_join(pred10_tab, by = c("Severity")) %>%
  full_join(med_tab, by = c("Severity")) %>%
  full_join(counts, by = c("Severity"))


#Variance of frailty terms
theta <- raw_tab %>%
  filter(value == "theta") %>%
  mutate(theta = round(est, 2)) %>%
  select(label, theta)

