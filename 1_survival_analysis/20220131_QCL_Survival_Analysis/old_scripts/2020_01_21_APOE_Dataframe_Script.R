###--------------------------------------------------------------------------###
### Project: ApoE in ALS                                                       -
### Purpose: SBA data extraction                                               -
### Authors: FS Sanders, QC Lin                                                -
### Last Edit: 2022-02-10, QC Lin                                              -
### Contact: q.c.lin@students.uu.nl or sebastian.lewandowski@ki.se             -
###--------------------------------------------------------------------------###
### Comment:                                                                   -
### This script generates new dataframes from the objects I received.          -
### It only uses Utrecht cohort data for APOE analyses                         -
### Date last modified: 2020/01/21 # F.S.Sanders                               -
### Most of the code to generate specific columns are adapted from original    -
### script made by Julia and Anna.                                             -
### ~ F.S. Sanders
###--------------------------------------------------------------------------###

### Dependencies and general settings ------------------------------------------

# Load libraries
library(tidyverse)
library(SBA)
library(survival)
library(survminer)
library(testthat)

# Set working directory and seed for reproducibility
setwd("~/Desktop/20220131_QCL_Survival_Analysis/")
set.seed(747)

# Create output directory for pdfs
ifelse(!dir.exists(file.path(getwd(), "generated_data_R")),
       dir.create(file.path(getwd(), "generated_data_R")),
       FALSE)

### Data load-in ---------------------------------------------------------------

# Load data
load(file="R_data/ALS_8 Utrecht data file pure ALS and over 5 years after diagnosis removed.Rdata")
als.ma.utrecht <- als.utrecht
# Remove original datasource
rm(als.utrecht)

### Isolation of bead info -----------------------------------------------------

bead_info_utrecht <- als.ma.utrecht@bead[,1:3]

### Utrecht data (ApoE filtered) -----------------------------------------------

# Adding an ApoE column
apoe <- rep(NA,nrow(als.ma.utrecht@sample))
apoe[which(als.ma.utrecht@sample$Diagnosis=="ALS" & als.ma.utrecht@sample$common_name=="Apo-e2e2")] <- "ALS Apo-e2e2"
apoe[which(als.ma.utrecht@sample$Diagnosis=="ALS" & als.ma.utrecht@sample$common_name=="Apo-e3e3")] <- "ALS Apo-e3e3"
apoe[which(als.ma.utrecht@sample$Diagnosis=="ALS" & als.ma.utrecht@sample$common_name=="Apo-e3e4")] <- "ALS Apo-e3e4"
apoe[which(als.ma.utrecht@sample$Diagnosis=="ALS" & als.ma.utrecht@sample$common_name=="Apo-e4e4")] <- "ALS Apo-e4e4"
apoe[which(als.ma.utrecht@sample$class=="Control" & als.ma.utrecht@sample$common_name=="Apo-e2e2")] <- "Control Apo-e2e2"
apoe[which(als.ma.utrecht@sample$class=="Control" & als.ma.utrecht@sample$common_name=="Apo-e3e3")] <- "Control Apo-e3e3"
apoe[which(als.ma.utrecht@sample$class=="Control" & als.ma.utrecht@sample$common_name=="Apo-e3e4")] <- "Control Apo-e3e4"
apoe[which(als.ma.utrecht@sample$class=="Control" & als.ma.utrecht@sample$common_name=="Apo-e4e4")] <- "Control Apo-e4e4"
als.ma.utrecht@sample$apoe <- apoe

# Survival calculations
# Calculating survival months for alive patients, last check-up in April 2018
year <- substr(as.character(als.ma.utrecht@sample$Date_of_onset), 1, 4)
year.onset <- as.numeric(year)
month <- substr(as.character(als.ma.utrecht@sample$Date_of_onset), 6, 7)
all.months <- c("01","02","03","04","05","06","07","08","09","10","11","12")
month.onset <- match(month,all.months)
year.check <- rep(2018,length(year.onset))
month.check <- rep(4,length(month.onset))
survival.months <- (year.check-year.onset)*12+month.check-month.onset # Calculate months, disregard days
death.year <- substr(as.character(als.ma.utrecht@sample$Date.of.Death), 1, 4)
death.year <- as.numeric(death.year)
death.month <- substr(as.character(als.ma.utrecht@sample$Date.of.Death), 6, 7)
death.month <- match(death.month,all.months)
death.survival <- (death.year-year.onset)*12+death.month-month.onset # Calculate months, disregard days
survival.months[which(is.na(death.survival)==FALSE)] <- death.survival[which(is.na(death.survival)==FALSE)] # Replace the known deaths
als.ma.utrecht@sample$survival_months <- survival.months
died <- rep(FALSE,length(year))
died[which(is.na(death.survival)==FALSE)] <- TRUE
als.ma.utrecht@sample$died <- died
als.ma.utrecht@sample$survobj <- with(als.ma.utrecht@sample, Surv(survival_months, died))

# Site on onset at Bulbar (Boolean)
onset <- rep(FALSE,length(als.ma.utrecht@sample$Site_of_Onset))
onset[which(als.ma.utrecht@sample$Site_of_Onset=="Bulbar")] <- TRUE
als.ma.utrecht@sample$onset_bulbar <- onset

# Get protein MFI values from als.ma.utrecht@X
raw.apoev <- als.ma.utrecht@X[,1:193]
raw.apoev <- as.data.frame(raw.apoev)
colnames(raw.apoev) <- paste("RawB.",colnames(raw.apoev),sep = "")

# Creating new APOE dataframe
# Subsetting columns of interest, adding protein values, include cohort signature.
# Initially only interested in apoe, hence the name of the new dataframe.
als.ma.apoe <- als.ma.utrecht@sample[,c("sample_name","class", "gender", "age_at_onset",
                                        "sampling_age","Date_of_onset", "sampling_date",
                                        "survival_months", "died","survobj","D50", "dx", "onset_bulbar", "apoe", "Date.of.Death")]

#Adding gap (sampling delay information), "immortal timeline"
#This gap is more precise, might use it in the future. For Utrecht the differences do not seem to be huge.
# als.ma.apoe <- mutate(als.ma.apoe, gap = as.Date(als.ma.apoe$sampling_date) - as.Date(als.ma.apoe$Date_of_onset))
# als.ma.apoe$gap <- as.numeric(als.ma.apoe$gap/30.5)
# als.ma.apoe$gap <- round(als.ma.apoe$gap/12)
# als.ma.apoe$gap[als.ma.apoe$gap == -1] <- 0

# Bind together
als.ma.apoe <- cbind(als.ma.apoe, raw.apoev)
als.ma.apoe$cohort <- "Utrecht"

### Data transfer QC -----------------------------------------------------------

# Following function to test if sample, MFI values and survival months are still correctly linked
# Note: Rownames ID are transferred from original dataframe, 
# BUT the sample_id column in original dataframe contains an empty additional attribute.
# Consequently test is without checking attributes so -> check attributes is FALSE.
# The function gives an error if they are not equal
test_that("New DF patients and values are similar to original DF",{
  x <- cbind(as.character(als.ma.utrecht@sample$sample_name), als.ma.utrecht@sample$survival_months ,als.ma.utrecht@X[,1])
  y <- cbind(as.character(als.ma.apoe$sample_name), als.ma.apoe$survival_months ,als.ma.apoe$RawB.1)
  expect_equal(x, y, check.attributes = FALSE)
})

### Restricted data filtering (Case) -------------------------------------------

# Add Raw identifier to bead_info dataframe
bead_info_utrecht$Raw_name <- colnames(als.ma.apoe)[16:208]

# Adding identifiable names for beads to als.ma.apoe
colnames(als.ma.apoe)[16:208] <- c(as.character(bead_info_utrecht$Bead_name))

# Only cases dataframe
apoeALS <- subset(als.ma.apoe, als.ma.apoe$class == "Case")
# Gap calculation
apoeALS <- mutate(apoeALS, gap = apoeALS$sampling_age - apoeALS$age_at_onset)
apoeALS$gap[apoeALS$gap == -1] <- 0

# Export dataframe as rds
saveRDS(bead_info_utrecht, file = "generated_data_R/bead_info_utrecht.rds")
saveRDS(als.ma.apoe, file = "generated_data_R/als.ma.apoe.rds")
saveRDS(apoeALS, file = "generated_data_R/apoeALS.rds")

### Restricted data filtering --------------------------------------------------
# Filtering on genotype info, only e3e3 and e3e4
genotype_df <- subset(apoeALS, apoeALS$apoe == "ALS Apo-e3e3" | apoeALS$apoe == "ALS Apo-e3e4" )

# Even smaller dataframe with info needed
apoecox <- genotype_df %>%
  transmute(survobj,
            APOE90 = RawB.90,
            APOE182 = RawB.182,
            survival_months,
            died,
            Genotype = apoe,
            Gender = gender,
            Onset = factor(onset_bulbar, labels = c("Thoracic/Spinal", "Bulbar")),
            Sampling_Age = sampling_age,
            Sampling_delay = gap,
            cohort)

# Quartiled based on filtered data
apoecox <- mutate(apoecox, QB.APOE90 = as.factor(ntile(apoecox$APOE90, 4)))
apoecox <- mutate(apoecox, QB.APOE182 = as.factor(ntile(apoecox$APOE182, 4)))

# Only thorasic
thor <- subset(apoecox, apoecox$Onset == "Thoracic/Spinal")

# Making cutoffs, in data filtered on Thorasic Onset
single.cut <- surv_cutpoint(thor, time = "survival_months", event = "died", variables = c("APOE90"),minprop = 0.1)
max90 <- single.cut[["APOE90"]]
max90
cut90 <- max90$estimate
single.cut <- surv_cutpoint(thor, time = "survival_months", event = "died", variables = c("APOE182"),minprop = 0.1)
max182 <- single.cut[["APOE182"]]
max182
cut182 <- max182$estimate

# Stratifying on earlier establish cutoffs.
thor <- thor %>% mutate(APOE90_Level = case_when(
  thor$APOE90 <= cut90 ~ "Lower",
  thor$APOE90 > cut90 ~ "Upper"))
thor <- thor %>% mutate(APOE182_Level = case_when(
  thor$APOE182 <= cut182 ~ "Lower",
  thor$APOE182 > cut182 ~ "Upper"))

# Only e3e3 genotype in thoracic dataset
thore3 <- subset(thor, thor$Genotype == "ALS Apo-e3e3")

# Quartiled based on filtered data
thore3 <- mutate(thore3, QB.APOE90 = as.factor(ntile(thore3$APOE90, 4)))
thore3 <- mutate(thore3, QB.APOE182 = as.factor(ntile(thore3$APOE182, 4)))
thor <- mutate(thor, QB.APOE90 = as.factor(ntile(thor$APOE90, 4)))
thor <- mutate(thor, QB.APOE182 = as.factor(ntile(thor$APOE182, 4)))
