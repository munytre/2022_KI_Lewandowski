#### Correlation plots APOE ####
# This script generates all spearman correlation plots and boxplots
# Last edited: 2020-01-21
# F.S.Sanders

library(tidyverse)
library(forestmodel)
library(timeROC)

setwd("~/Documents/ALS-research/APOE_2020_01_20/")

# See if there is a correlation with age
# ggplot(data= genotype_df, aes(age_at_onset, RawB.90),
#        color = apoe) +
#   geom_point(aes(color = apoe))+ 
#   geom_smooth(aes(color = apoe),method = "lm")+
#   stat_cor(aes(color =  apoe),method = "spearman", label.x = 75)+
#   theme_bw()
# 
# ggplot(data= genotype_df, aes(age_at_onset, RawB.182),
#        color = apoe) +
#   geom_point(aes(color = apoe))+ 
#   geom_smooth(aes(color = apoe),method = "lm")+
#   stat_cor(aes(color =  apoe),method = "spearman", label.x = 75)+
#   theme_bw()

#No correlation found 

# Gender differences
# ggboxplot(apoeALS, x = "gender", y = "RawB.182", xlab = "Case vs Control",
#                ylab = "MFI APOE (B.182)",title = "Gender differences",
#                add = "jitter" )+
#   stat_compare_means(label.y = 10)
# 
# ggboxplot(apoeALS, x = "gender", y = "RawB.90", xlab = "Case vs Control",
#                ylab = "MFI APOE (B.90)",title = "Gender differences",
#                add = "jitter" )+
#   stat_compare_means(label.y = 10)

# No differences in between male and female patients. 

# Genotype onset differences
# Had to recreate to include missing genotype data for overview

APOEgen <- apoeALS %>% 
  transmute(survobj,
            APOE90 = RawB.90,
            APOE182 = RawB.182,
            NEFL = RawB.193,
            Gender = gender,
            Genotype = apoe,
            Onset = factor(onset_bulbar, labels = c("Thoracic/Spinal", "Bulbar")),
            Sampling_Age = sampling_age,
            Sampling_delay = gap,
            cohort)

APOEgen <- APOEgen %>% mutate(Genotype = replace_na(Genotype, "Missing"))

genoplot <- ggboxplot(APOEgen, x = "Onset", y = "APOE90", 
              xlab = "Thoracic/Spinal vs Bulbar Onset",
                   ylab = "MFI APOE (B.90)",title = "Onset differences",)+
              geom_jitter(aes(color = Genotype))+
              stat_compare_means(label.y = 9.5)

pdf("APOE_Genotype_Onset_Comparison.pdf", width = 8, height = 6)
print(genoplot)
dev.off()

# Noteworthy: 5 patients with e4e4 are all in the lower quartile of MFI values.

# Survival correlations

a <- ggplot(data= genotype_df, aes(survival_months, RawB.182),
       color = apoe) +
  geom_point( )+ 
  geom_smooth( method = "lm")+
  stat_cor(method = "spearman", label.x = 80, size  = 3)+
  labs(y= "MFI (B.182)", x="Survival Time in Months")+
  theme_bw()


b <- ggplot(data= genotype_df, aes(survival_months, RawB.90),
       color = apoe) +
  geom_point( )+ 
  geom_smooth( method = "lm")+
  stat_cor(method = "spearman", label.x = 80, size  = 3)+
  labs(y= "MFI (B.90)", x="Survival Time in Months")+
  theme_bw()

c <- ggplot(data= genotype_df, aes(survival_months, RawB.182),
       color = apoe) +
  geom_point(aes(color = apoe))+ 
  geom_smooth(aes(color = apoe),method = "lm")+
  stat_cor(aes(color =  apoe),method = "spearman", label.x = 80, size  = 3)+
  labs(y= "MFI (B.182)", x="Survival Time in Months")+
  theme_bw()

d <- ggplot(data= genotype_df, aes(survival_months, RawB.90),
       color = apoe) +
  geom_point(aes(color = apoe))+ 
  geom_smooth(aes(color = apoe),method = "lm")+
  stat_cor(aes(color =  apoe),method = "spearman", label.x = 80, size  = 3)+
  labs(y= "MFI (B.90)", x="Survival Time in Months")+
  theme_bw()

# Run lay, plotlist and gridExtra if you want to manually arrange plots.

# lay <- rbind(c(1,1,2,2),
#              c(3,3,4,4))
# plotlist <<- list(b,a,d,c)
pdf(file = "APOE_Survival_Protein_Correlations.pdf", width = 10,height = 8)
# gridExtra::grid.arrange(grobs = plotlist, layout_matrix = lay)
ggarrange(b,a,d,c, labels = c("A", "B", "C","D"), ncol = 2, nrow = 2, legend="top")
dev.off()

#Thoracic and Bulbar whole cohort
# Bulbar <- subset(apoeALS, apoeALS$onset_bulbar == TRUE)
# Thoracic <- subset(apoeALS, apoeALS$onset_bulbar == FALSE)
# 
# a <- ggplot(data= Thoracic, aes(survival_months, RawB.90)) +
#   geom_point()+
#   geom_smooth(method = "lm")+
#   stat_cor(method = "spearman", label.x = 80, size = 3)+
#   labs(y= "MFI (B.90)", x="Survival Time in Months")+
#   theme_bw()+
#   ggtitle("Thoracic only")
# 
# b <- ggplot(data= Bulbar, aes(survival_months, RawB.90)) +
#   geom_point()+
#   geom_smooth(method = "lm")+
#   stat_cor(method = "spearman", label.x = 60, size = 3)+
#   labs(y= "MFI (B.90)", x="Survival Time in Months")+
#   theme_bw() +
#   ggtitle("Bulbar only")
# a
# b


#Thoracic, Bulbar after genotype selection
Bulbar <- subset(genotype_df, genotype_df$onset_bulbar == TRUE)
Thoracic <- subset(genotype_df, genotype_df$onset_bulbar == FALSE)

a <- ggplot(data= Thoracic, aes(survival_months, RawB.90)) +
  geom_point()+ 
  geom_smooth(method = "lm")+
  stat_cor(method = "spearman", label.x = 80, size = 3)+
  labs(y= "MFI (B.90)", x="Survival Time in Months")+
  theme_bw()+
  ggtitle("Thoracic only")

b <- ggplot(data= Bulbar, aes(survival_months, RawB.90)) +
  geom_point()+ 
  geom_smooth(method = "lm")+
  stat_cor(method = "spearman", label.x = 60, size = 3)+
  labs(y= "MFI (B.90)", x="Survival Time in Months")+
  theme_bw() +
  ggtitle("Bulbar only")


c <- ggplot(data= Thoracic, aes(survival_months, RawB.90),
       color = apoe) +
  geom_point(aes(color = apoe))+ 
  geom_smooth(aes(color = apoe),method = "lm")+
  stat_cor(aes(color =  apoe),method = "spearman", label.x = 80, size = 3)+
  labs(y= "MFI (B.90)", x="Survival Time in Months")+
  theme_bw()+
  ggtitle("Thoracic only")

d <- ggplot(data= Bulbar, aes(survival_months, RawB.90),
       color = apoe) +
  geom_point(aes(color = apoe))+ 
  geom_smooth(aes(color = apoe),method = "lm")+
  stat_cor(aes(color =  apoe),method = "spearman", label.x = 60, size = 3)+
  labs(y= "MFI (B.90)", x="Survival Time in Months")+
  theme_bw() +
  ggtitle("Bulbar only")

pdf(file = "APOE_BulbThor_Correlations.pdf", width = 10,height = 8)
ggarrange(b,a,d,c, labels = c("A", "B", "C","D"), ncol = 2, nrow = 2, legend="top")
dev.off()

# Genotype plots
# Need tot include controls
GT_E3 <- subset(als.ma.apoe, als.ma.apoe$apoe == "ALS Apo-e3e3" | als.ma.apoe$apoe == "Control Apo-e3e3")
GT_E4 <- subset(als.ma.apoe, als.ma.apoe$apoe == "ALS Apo-e3e4" | als.ma.apoe$apoe == "Control Apo-e3e4")

alsE3 <- subset(als.ma.apoe, als.ma.apoe$apoe == "ALS Apo-e3e3")
alsE4 <- subset(als.ma.apoe, als.ma.apoe$apoe == "ALS Apo-e3e4")

ggboxplot(data= GT_E3, x = "class", y = "RawB.182", xlab = "Case vs Control",
          ylab = "MFI APOE (B.182)",title = "APO-e3e3 Control vs Case",
          add = "jitter", fill = "class" ,palette = c("#00AFBB", "#E7B800")) +
  stat_compare_means()  + theme_bw()

apoegt <- subset(als.ma.apoe, (!is.na(als.ma.apoe["apoe"]))) #this does nothing, Column is filled without "NA" instead if NA
apoefilter <- apoegt %>% 
  filter(!apoe %in% c("ALS Apo-e2e2", "ALS Apo-e4e4", "Control Apo-e2e2", "Control Apo-e4e4"))
compare <-list(c("Control Apo-e3e3","ALS Apo-e3e3" ),c("Control Apo-e3e4", "ALS Apo-e3e4"))

a <- ggboxplot(apoefilter, x = "apoe", y = "RawB.182", color  = "apoe",xlab = "Case vs Control",
          ylab = "MFI APOE (B.182)",title = "APOE Control vs Case",
          add = "jitter", order = c("ALS Apo-e3e3", "Control Apo-e3e3", "ALS Apo-e3e4","Control Apo-e3e4" ))+
  stat_compare_means(comparisons =  compare)+
  stat_compare_means(label.y = 10)

b <- ggboxplot(apoefilter, x = "apoe", y = "RawB.90", color  = "apoe",xlab = "Case vs Control",
          ylab = "MFI APOE (B.90)",title = "APOE Control vs Case",
          add = "jitter", order = c("ALS Apo-e3e3", "Control Apo-e3e3", "ALS Apo-e3e4","Control Apo-e3e4" ))+
  stat_compare_means(comparisons =  compare)+
  stat_compare_means(label.y = 10)
b
pdf(file = "Genotype_Comparison_APOE.pdf", width = 14, height = 5)
lay <- rbind(c(1,2))
plotlist <- list(b,a)
gridExtra::grid.arrange(grobs = plotlist, layout_matrix = lay)
dev.off()

#### Final cox model ####

# Another dataframe with all patients, disregarding genotypes
# Dataframe with info needed, not stratifying
allcox <- apoeALS %>%
  transmute(survobj,
            APOE90 = RawB.90,
            APOE182 = RawB.182,
            survival_months,
            Genotype = apoe,
            Gender = gender,
            Onset = factor(onset_bulbar, labels = c("Thorasic/Spinal", "Bulbar")),
            Sampling_Age = sampling_age,
            Sampling_delay = gap,
            cohort)

# Stratifying on earlier establish cutoffs.
allcox <- allcox %>% mutate(APOE90_Level = case_when(
  allcox$APOE90 <= cut90 ~ "Lower",
  allcox$APOE90 > cut90 ~ "Upper"))
allcox <- allcox %>% mutate(APOE182_Level = case_when(
  allcox$APOE182 <= cut182 ~ "Lower",
  allcox$APOE182 > cut182 ~ "Upper"))

#Quartiled
allcox <- mutate(allcox, QB.APOE90 = as.factor(ntile(allcox$APOE90, 4)))
allcox <- mutate(allcox, QB.APOE182 = as.factor(ntile(allcox$APOE182, 4)))

# Models and plotting
cox1 <- coxph(survobj ~ QB.APOE90 +
                Gender + Sampling_Age + strata(Sampling_delay), data = thor)
cox1
pdf(file = "ThorassicCox90.pdf", height = 6, width = 10)
forestmodel::forest_model(cox1)
dev.off()

cox2 <- coxph(survobj ~ QB.APOE90 +
                Gender + Sampling_Age + strata(Sampling_delay), data = thore3)
cox2
pdf(file = "ThorassicE3E3Cox90.pdf", height = 6, width = 10)
forestmodel::forest_model(cox2)
dev.off()

cox3 <- coxph(survobj ~ QB.APOE182 +
                Gender + Sampling_Age + strata(Sampling_delay), data = thor)
cox3
pdf(file = "ThorassicCox182.pdf", height = 6, width = 10)
forestmodel::forest_model(cox3)
dev.off()

cox4 <- coxph(survobj ~ QB.APOE182 +
                Gender + Sampling_Age + strata(Sampling_delay), data = thore3)
cox4
pdf(file = "ThorassicE3E3Cox182.pdf", height = 6, width = 10)
forestmodel::forest_model(cox4)
dev.off()

cox5 <- coxph(survobj ~ QB.APOE90 +
                Gender + Sampling_Age + strata(Sampling_delay), data = allcox)
cox5
pdf(file = "APOE90Cox.pdf", height = 6, width = 10)
forestmodel::forest_model(cox5)
dev.off()

cox6 <- coxph(survobj ~ QB.APOE182 +
                Gender + Sampling_Age + strata(Sampling_delay), data = allcox)
cox6
pdf(file = "APOE182Cox.pdf", height = 10, width = 10)
forestmodel::forest_model(cox6)
dev.off()

### For genotype hazard ###

gtcox <- coxph(survobj ~ QB.APOE90 + Genotype + Onset + Gender + Sampling_Age + strata(Sampling_delay), data = apoecox)
pdf(file = "APOEGenotypehazardCox.pdf")
forestmodel::forest_model(gtcox, format_options = list(text_size = 3))
forestmodel::forest_model(gtcox)
dev.off()

#### PLOT ROCs Combined #####

ROC.all <-  timeROC(T=apoeALS$survival_months,
                    delta=apoeALS$died,
                    marker=apoeALS$RawB.182,
                    cause=TRUE,weighting="marginal",
                    times=quantile(apoeALS$survival_months))
ROC.all
ROC.Thorassic <- timeROC(T=thor$survival_months,
                         delta=thor$died,
                         marker=thor$APOE182,
                         cause=TRUE,weighting="marginal",
                         times=quantile(apoeALS$survival_months))
ROC.Thorassic
# ROC.Thorassic90 <- timeROC(T=thor$survival_months,
#                          delta=thor$died,
#                          marker=thor$APOE90,
#                          cause=TRUE,weighting="marginal",
#                          times=quantile(apoeALS$survival_months))
# ROC.Thorassic90
ROC.ThorassicE3 <- timeROC(T=thore3$survival_months,
                           delta=thore3$died,
                           marker=thore3$APOE182,
                           cause=TRUE,weighting="marginal",
                           times=quantile(apoeALS$survival_months))
ROC.ThorassicE3

pdf("ROC_APOE_182.pdf", width = 8, height = 8)
plot(ROC.all, time = 32.5 ,title = FALSE)
plot(ROC.Thorassic, time = 32.5, add = TRUE, col = "blue", lty = 5,title = FALSE)
plot(ROC.ThorassicE3, time = 32.5, add = TRUE, col = "darkgreen", lty = 10,title = FALSE)
legend(0.6,0.45, c('All,       AUC: 59.34','Thor/Spin, AUC: 70.21','Thor/e3e3, AUC: 74.19'),lty=c(1,1),
       lwd=c(2,2),col=c('red','blue','darkgreen'))
title("C/D Time Dependent ROC's")
dev.off()
# For comparing objects using the package
# plotAUCcurveDiff(ROC.Thorassic,ROC.Thorassic90,conf.int=TRUE,conf.band=TRUE,ylim=c(-0.2,0.5))

