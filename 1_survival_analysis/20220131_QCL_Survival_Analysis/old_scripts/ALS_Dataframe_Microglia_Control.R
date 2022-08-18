###### 
# This script generates new dataframes from the objects I received. 
# Making these allowed me to get to know accustom myself with the data and create a DF I know. 
# It makes a new dataframe for Utrecht, Ulm, Leuven and combines them into a big dataframe for combined analyses. 
# Tests are included. 
# Date modified: 2019/10/31 # F.S.Sanders
# Date modified: 2019/11/08 # FSSanders

# Microglia beads:
# B.183_AIF1_HPA049234
# B.34_CD11b ITGAM_550282
# B.29_CD68 ED1_MCA341GA
# B.186_AIF1_HPA062949

# Control beads
# B.1_Anti-IgG_control
# B.4_Empty_control
# B.3_Rabbit IgG_control

###### 

###### 
# Load workspace and files
######
setwd("~/Documents/ALS-research/WorkSpace/")

load(file="VASC2 ALS Ulm Leuven and Utrecht data median scaled.Rdata")

als.ma.utrecht <- als.utrecht.scale
als.ma.ulm <- als.ulm
als.ma.leuven <- als.leuven.scale

rm(als.utrecht.scale)
rm(als.ulm)
rm(als.leuven.scale)

# Compared unscaled dataframe with scaled one using the following:
# all.equal(scaled, old)
# Results only show difference in MFI mean relative differences for the proteins, which is as expected.

###### I run all packages I generally use in other scripts working with 
library(SBA)
library(dplyr)
library(gplots)
library(beeswarm)
library(ggbeeswarm)
library(lubridate)
library(survival)
library(survminer)
library(coxme)
# library(GGally)
library(ggpubr)
#library(timeROC)
#library(cutpointr)
library(ggplot2)
library(gridExtra)
library(testthat)

########################################
# Utrecht Dataframe
# Most of the code to generate specific columns are adapted from original script made by Julia and Anna

# Adding an apoe column
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

# Site on onset at Bulbar true or false
onset <- rep(FALSE,length(als.ma.utrecht@sample$Site_of_Onset))
onset[which(als.ma.utrecht@sample$Site_of_Onset=="Bulbar")] <- TRUE
als.ma.utrecht@sample$onset_bulbar <- onset



#Get protein MFI values from als.ma.utrecht@X
#Spp1, Apoe, Col6a1, Apoe#2, 4 kinds of NEFL/H
raw.apoev <- als.ma.utrecht@X[,c("1","3","4","29","34","53","90","123","182","88","94","97","183","186","193")]
raw.apoev <- as.data.frame(raw.apoev)
raw.apoev <- raw.apoev %>%
  rename("RawB.1" ="1","RawB.3" ="3","RawB.4" ="4","RawB.29" ="29","RawB.34" ="34","RawB.53" = "53",
         "RawB.88" = "88","RawB.90" = "90", "RawB.94" = "94","RawB.97" = "97","RawB.123" = "123",
         "RawB.182" = "182","RawB.183" = "183","RawB.186" = "186","RawB.193" = "193")

#Creating new APOE dataframe
#Subsetting columns of interest, adding protein values, inculde cohort signature.
#Initially only interested in apoe, hence the name of the new dataframe.
als.ma.apoe <- als.ma.utrecht@sample[,c("sample_name","class", "gender", "age_at_onset",
                                        "sampling_age","Date_of_diagnosis", "sampling_date",
                                        "survival_months", "died","survobj","D50", "dx", "onset_bulbar", "apoe")]
als.ma.apoe <- cbind(als.ma.apoe, raw.apoev)
als.ma.apoe$cohort <- "Utrecht"

#Only cases dataframe
apoeALS <- subset(als.ma.apoe, als.ma.apoe$class == "Case")

#Adding gap (sampling delay information), "immortal timeline"
apoeALS <- mutate(apoeALS, gap = apoeALS$sampling_age - apoeALS$age_at_onset)
apoeALS$gap[apoeALS$gap == -1] <- 0
#The following can be added to generate 4 equal quartiles of MFI values
# apoeALS <- mutate(apoeALS, QB.53 = ntile(apoeALS$RawB.53, 4))
# apoeALS <- mutate(apoeALS, QB.90 = ntile(apoeALS$RawB.90, 4))
# apoeALS <- mutate(apoeALS, QB.123 = ntile(apoeALS$RawB.123, 4))
# apoeALS <- mutate(apoeALS, QB.182 = ntile(apoeALS$RawB.182, 4))

#Create dataframes only containing specific genotypes
# als.ma.apoeGT <- subset(als.ma.apoe, als.ma.apoe$apoe == "ALS Apo-e3e3" | als.ma.apoe$apoe == "ALS Apo-e3e4" )
# alsE3 <- subset(als.ma.apoe, als.ma.apoe$apoe == "ALS Apo-e3e3")
# alsE4 <- subset(als.ma.apoe, als.ma.apoe$apoe == "ALS Apo-e3e4")

# Following function to test if sample, MFI values and survival months are still correctly linked
# Note: Rownames ID are transferred from original dataframe, BUT the sample_id column in original dataframe contains an empty additional attribute.
# Consequently test is without checking attributes so -> check attributes is FALSE.
# The function gives an error if they are not equal
test_that("New DF patients and values are similar to original DF",{
  x <- cbind(as.character(als.ma.utrecht@sample$sample_name), als.ma.utrecht@sample$survival_months ,als.ma.utrecht@X[,53])
  y <- cbind(as.character(als.ma.apoe$sample_name), als.ma.apoe$survival_months ,als.ma.apoe$RawB.53)
  expect_equal(x, y, check.attributes = FALSE)
})

# used this to compare attributes 
# x <- cbind(as.character(als.ma.ulm@sample$sample_name), als.ma.ulm@sample$survival_months ,als.ma.ulm@X[,53])
# y <- cbind(as.character(als.ulm.apoe$sample_name), als.ulm.apoe$survival_months ,als.ulm.apoe$RawB.53)
# str(x)
# str(y)

########################################
#ULM dataframe

apoe <- rep(NA,nrow(als.ma.ulm@sample))
als.ma.ulm@sample$apoe <- apoe

# Adding survival months for alive patients, last check-up in April 2018
als.ma.ulm@sample$died <- als.ma.ulm@sample$died_update_2018
survival_months <- als.ma.ulm@sample$survival_months
died <- als.ma.ulm@sample$died
als.ma.ulm@sample$survobj <- with(als.ma.ulm@sample, Surv(survival_months, died))

# Site on onset at Bulbar true or false
onset <- als.ma.ulm@sample$init_symp_bulbar_update_2018
als.ma.ulm@sample$onset_bulbar <- onset

#G#Get protein MFI values from als.ma.ulm@X
raw.apoev <- als.ma.ulm@X[,c("1","3","4","29","34","53","90","123","182","88","94","97","183","186","193")]
raw.apoev <- as.data.frame(raw.apoev)
raw.apoev <- raw.apoev %>%
  rename("RawB.1" ="1","RawB.3" ="3","RawB.4" ="4","RawB.29" ="29","RawB.34" ="34","RawB.53" = "53","RawB.88" = "88","RawB.90" = "90",
         "RawB.94" = "94","RawB.97" = "97",
         "RawB.123" = "123", "RawB.182" = "182","RawB.183" = "183","RawB.186" = "186","RawB.193" = "193")
#Creating new APOE dataframe
#Subsetting columns of interest, adding protein values, inculde cohort signature.
#Initially only interested in apoe, hence the name of the new dataframe.
als.ulm.apoe <- als.ma.ulm@sample[,c("sample_name","class", "gender", "age_at_onset","sampling_age",
                                     "survival_months", "died","survobj" , "D50", "dx", "onset_bulbar", "apoe")]
als.ulm.apoe <- cbind(als.ulm.apoe, raw.apoev)
als.ulm.apoe$cohort <- "Ulm"

#Select patients only, into new DF
ulmALS <- subset(als.ulm.apoe, als.ulm.apoe$class == "Case")

#Remove NA survival months, which is present in Ulm data
ulmALS <- ulmALS[!is.na(ulmALS$survival_months),]
#20 rows are removed.

# ulmALS <- mutate(ulmALS, QB.53 = ntile(ulmALS$RawB.53, 4))
# ulmALS <- mutate(ulmALS, QB.90 = ntile(ulmALS$RawB.90, 4))
# ulmALS <- mutate(ulmALS, QB.123 = ntile(ulmALS$RawB.123, 4))
# ulmALS <- mutate(ulmALS, QB.182 = ntile(ulmALS$RawB.182, 4))

# Add gap variable to account for the delay between onset and sampling.
ulmALS <- mutate(ulmALS, gap = ulmALS$sampling_age - ulmALS$age_at_onset)
ulmALS$gap[ulmALS$gap == -1] <- 0

# Same test as before
test_that("New DF patients and values are similar to original DF",{
  x <- cbind(as.character(als.ma.ulm@sample$sample_name), als.ma.ulm@sample$survival_months ,als.ma.ulm@X[,53])
  y <- cbind(as.character(als.ulm.apoe$sample_name), als.ulm.apoe$survival_months ,als.ulm.apoe$RawB.53)
  expect_equal(x, y, check.attributes = FALSE)
})

########################################
#Leuven dataframe

apoe <- rep(NA,nrow(als.ma.leuven@sample))
als.ma.leuven@sample$apoe <- apoe

# Adding survival months for alive patients, last check-up in April 2018
als.ma.leuven@sample$died <- als.ma.leuven@sample$died
survival_months <- als.ma.leuven@sample$survival_months
died <- als.ma.leuven@sample$died
als.ma.leuven@sample$survobj <- with(als.ma.leuven@sample, Surv(survival_months, died))

# Site on onset at Bulbar true or false
onset <- rep(FALSE,nrow(als.ma.leuven@sample))
onset[which(als.ma.leuven@sample$Site_of_Onset_for_analysis=="Bulbar")] <- TRUE
als.ma.leuven@sample$onset_bulbar <- onset

#Selecting MFI values from als.ma.leuven@X
raw.apoev <- als.ma.leuven@X[ ,c("B.1_Anti-IgG_control","B.3_Rabbit IgG_control", "B.4_Empty_control","B.29_CD68 ED1_MCA341GA",
                                 "B.34_CD11b ITGAM_550282", "B.53_SPP1_HPA005562","B.88_NEFH_HPA061615", "B.90_APOE_HPA065539",
                                 "B.94_NEFL_Novus NBP237525","B.97_NEFH_AMAb91025",
                                 "B.123_COL6A1_HPA019142","B.182_APOE_HPA068768","B.183_AIF1_HPA049234","B.186_AIF1_HPA062949","B.193_NEFL_Abnova 4747001")]
raw.apoev <- as.data.frame(raw.apoev)
raw.apoev <- raw.apoev %>%
  rename("RawB.1" ="B.1_Anti-IgG_control","RawB.3" ="B.3_Rabbit IgG_control","RawB.4" ="B.4_Empty_control",
         "RawB.29" ="B.29_CD68 ED1_MCA341GA","RawB.34" = "B.34_CD11b ITGAM_550282",
         "RawB.53" = "B.53_SPP1_HPA005562","RawB.88" = "B.88_NEFH_HPA061615", "RawB.90" = "B.90_APOE_HPA065539",
         "RawB.94" = "B.94_NEFL_Novus NBP237525","RawB.97" = "B.97_NEFH_AMAb91025", "RawB.123" = "B.123_COL6A1_HPA019142","RawB.182" = "B.182_APOE_HPA068768",
         "RawB.183" ="B.183_AIF1_HPA049234","RawB.186" = "B.186_AIF1_HPA062949","RawB.193" = "B.193_NEFL_Abnova 4747001")

#Creating new APOE dataframe
als.leuven.apoe <- als.ma.leuven@sample[,c("sample_name","class3", "gender", "Age_of_onset","Age_at_sampling",
                                           "survival_months", "died","survobj" , "D50", "dx", "onset_bulbar", "apoe")]

# als.leuven.apoe <- cbind(als.leuven.apoe, quartiles)
als.leuven.apoe <- cbind(als.leuven.apoe, raw.apoev)
als.leuven.apoe <- als.leuven.apoe %>%
  rename("age_at_onset" = "Age_of_onset", "class" = "class3","sampling_age" = "Age_at_sampling")
als.leuven.apoe$cohort <- "Leuven"

#Select cases for patient only DF
leuvenALS <- subset(als.leuven.apoe, als.leuven.apoe$class == "ALS")

#Remove NA survival months, which is present in leuven data
leuvenALS <- leuvenALS[!is.na(leuvenALS$survival_months),]

#Gap
leuvenALS <- mutate(leuvenALS, gap = leuvenALS$sampling_age - leuvenALS$age_at_onset)
leuvenALS$gap[leuvenALS$gap == -1] <- 0
# Same test as before
test_that("New DF patients and values are similar to original DF",{
  x <- cbind(as.character(als.ma.leuven@sample$sample_name), als.ma.leuven@sample$survival_months ,als.ma.leuven@X[,53])
  y <- cbind(as.character(als.leuven.apoe$sample_name), als.leuven.apoe$survival_months ,als.leuven.apoe$RawB.53)
  expect_equal(x, y, check.attributes = FALSE)
})

# Quartiles can be added running the following
# leuvenALS <- mutate(leuvenALS, QB.53 = ntile(leuvenALS$RawB.53, 4))
# leuvenALS <- mutate(leuvenALS, QB.90 = ntile(leuvenALS$RawB.90, 4))
# leuvenALS <- mutate(leuvenALS, QB.123 = ntile(leuvenALS$RawB.123, 4))
# leuvenALS <- mutate(leuvenALS, QB.182 = ntile(leuvenALS$RawB.182, 4))

########################################
# For further analyses we want to combine datasets into one big DF 
# Remove columns not present in Ulm and Leuven
apoeALS <- select(apoeALS,-c(Date_of_diagnosis, sampling_date))
als.ma.apoe <- select(als.ma.apoe, -c(Date_of_diagnosis, sampling_date))

# Combine the three cohorts
### NOTE Survobj does not transfer after rbind. Need to recompute the column with Surv function. 
combinedALSinc <- rbind(apoeALS,ulmALS,leuvenALS)
combinedALSinc$survobj <- with(combinedALSinc, Surv(survival_months, died))
cohorts <- rbind(als.ma.apoe,als.ulm.apoe,als.leuven.apoe)

# Only Utrecht and Ulm
# combinedALS <- rbind(apoeALS,ulmALS)

#Clean global environment
rm(als.ma.utrecht)
rm(als.ma.ulm)
rm(als.ma.leuven)
cat("Script done:\nUtrecht, Ulm, Leuven and cohorts combined loaded in Global Environment.")
