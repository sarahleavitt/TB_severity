#Laura F White
#Boston University
#Pre-chemotherapy TB Analysis

##############################################################################
# This program plots to overall cure data that is not stratified by severity
##############################################################################
setwd("G:/My Drive/Work/TB Projects/Duration Research/pre_chemo_tb")

options(scipen=999)
options(digits = 10)

rm(list = ls())
source("R/utils.R")
reload_source()

#Reading in the study_id correspondence table
studyid <- read.csv("data/study_id.csv")

#Reading in all data
dataList <- read_excel_allsheets("data/pre_chemo_data.xlsx")
dataList$`Data dictionary` <- NULL

# get data we want from US with severity
cureList <- list(dataList$`1029_1055`, dataList$`1029_1056`, dataList$`48_1000_1029`,
                 dataList$`45`)

# make into a data frame
cureData1 <- do.call(rbind.data.frame, cureList)
cureData1$study_id <- as.character((cureData1$study_id))

# add in author and date of study
cureData2 <- as.data.frame(inner_join(cureData1,studyid,by="study_id"))

# calculate the cure rate
cureData2$cureRate <- cureData2$c2/cureData2$n

# basic plot of cure rate by year, author and severity; 
ggplot(cureData2[cureData2$severity!="None",],aes(x=interval_r,y=cureRate,group=cohort_id,color=severity))+
  geom_point(aes(shape=first_author))+geom_line()+xlim(0,10)+
  labs(x="Year since diagnosis",y="Probability of Natural Recovery",color="Severity")

  