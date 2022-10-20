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


#### Formatting results --------------------------------------------------------

form_comp <- formatBayesian(mortality, res_comp, data_comp, "Combined")
form_sev <- formatBayesian(mortality, res_sev, data_sev, "Severity",
                           fixed = TRUE)


#### Diagnostic Plots ----------------------------------------------------------

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


#### Table of Main Results -----------------------------------------------------

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

#Counts stratified by severity for baseline sample size
mortality_sub1 <- mortality %>% filter(severity != "Unknown",
                                      study_id %in% c("1029", "93", "45"))

#Counts stratified by severity for modeling sample size
mortality_sub <- mortality_sub1 %>% filter(no_survival == FALSE)

counts <- mortality_sub %>%
  group_by(severity) %>%
  summarize(nStudies = length(unique(study_id)),
            nCohorts = length(unique(cohort_id)),
            nIndividuals = n(),
            .groups = "drop") %>%
  mutate(Severity = paste0(severity, "")) %>%
  bind_rows(data_comp2) %>%
  select(-severity, -nSeverity)

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


#### Summary Survival Curves ---------------------------------------------------

ggplot(form_comp$surv_dens) +
  geom_line(aes(x = x, y = surv),
            color = "black", size = 1, linetype = "solid") +
  geom_smooth(aes(x = x, y = surv_est, ymin = cilb, ymax = ciub),
              stat = "identity", linetype = 0, alpha = 0.25, na.rm = TRUE) +
  scale_y_continuous(name = "Survival, 1 - F(t)", limits = c(0, 1)) +
  scale_x_continuous(name = "Years", limits = c(0, 30)) +
  theme_bw()

ggsave("Figures/summary_curves_comp.png", width = 5, height = 4.5)

ggplot(form_sev$surv_dens) +
  geom_line(aes(x = x, y = surv, color = severity),
            size = 1, linetype = "solid") +
  geom_smooth(aes(x = x, y = surv_est, ymin = cilb, ymax = ciub,
                  fill = severity),
              stat = "identity", linetype = 0, alpha = 0.15, na.rm = TRUE) +
  scale_y_continuous(name = "Survival, 1 - F(t)", limits = c(0, 1)) +
  scale_x_continuous(name = "Years", limits = c(0, 30)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_color_manual("", values = c("Minimal" = "seagreen",
                                    "Moderately advanced" = "goldenrod1",
                                    "Far advanced" = "firebrick2")) +
  scale_fill_manual("",
                    values = c("Minimal" = "seagreen",
                               "Moderately advanced" = "goldenrod1",
                               "Far advanced" = "firebrick2"))

ggsave("Figures/summary_curves_sev.png", width = 5, height = 5)



#### Individual Survival Curves ------------------------------------------------

#Survival curves for complete model
p1 <- ggplot(form_comp$ind_surv) +
  geom_line(aes(x = x, y = surv, group = study_sev, color = severity),
            size = 0.7, alpha = 0.3) +
  geom_line(data = form_comp$surv_dens, aes(x = x, y = surv),
            color = "black", size = 1, linetype = "longdash") +
  scale_y_continuous(name = "Survival, 1 - F(t)", limits = c(0, 1)) +
  scale_x_continuous(name = "Years", limits = c(0, 30)) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_color_manual("Disease Severity",
                     values = c("Minimal" = "seagreen", 
                                "Moderately advanced" = "goldenrod1",
                                "Far advanced" = "firebrick2"))

#Survival curves for stratified model
p2 <- ggplot(form_sev$ind_surv) +
  geom_line(aes(x = x, y = surv, group = study_sev, color = severity),
            size = 0.7, alpha = 0.3) +
  geom_line(data = form_sev$surv_dens,
            aes(x = x, y = surv, color = severity),
            linetype = "longdash", size = 1) +
  scale_y_continuous(name = "Survival, 1 - F(t)", limits = c(0, 1)) +
  scale_x_continuous(name = "Years", limits = c(0, 30)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_color_manual("", drop = FALSE,
                     values = c("Minimal" = "seagreen",
                                "Moderately advanced" = "goldenrod1",
                                "Far advanced" = "firebrick2"))

p_comb <- arrangeGrob(p1, p2, nrow = 2)
ggsave("Figures/individual_curves.png", p_comb, width = 5.5, height = 8)



##### Forest plots -------------------------------------------------------------

#TB survival for full model
ggplot(form_comp$pred_comb %>% filter(value != "median"),
       aes(x = est, y = first_author, xmin = cilb, xmax = ciub,
           shape = shape)) +
  facet_grid(severity ~ pred_label, scales = "free", space = "free") +
  geom_point() +
  geom_point(data = form_comp$pred_comb %>% filter(shape == "Overall" & value != "median"),
             color = 'black', shape = 18, size = 3) +
  geom_errorbar(width = 0.3) +
  scale_x_continuous(name = "1-Year Survival Probability", limits = c(0, 1),
                     breaks = c(0, 0.5, 1)) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        legend.position = "bottom") +
  scale_shape_discrete(guide = "none")

ggsave("Figures/forest_full.png", width = 7, height = 5)

#TB survival for stratified model
ggplot(form_sev$pred_comb %>% filter(value != "median"),
       aes(x = est, y = first_author, xmin = cilb, xmax = ciub, shape = shape)) +
  geom_point() +
  geom_point(data = form_sev$pred_comb %>% filter(shape == "Overall" & value != "median"),
             color = 'black', shape = 18, size = 3) +
  geom_errorbar(width = 0.3) +
  scale_x_continuous(name = "1-Year Survival Probability", limits = c(0, 1),
                     breaks = c(0, 0.5, 1)) +
  facet_grid(severity ~ pred_label, scales = "free", space = "free") +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        legend.position = "bottom") +
  scale_shape_discrete(guide = "none")

ggsave("Figures/forest_stratified.png", width = 7, height = 5)
