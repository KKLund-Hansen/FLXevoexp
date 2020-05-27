################################################################################################
########################### STANDARD FITNESS ASSAY GENERATION 15 & 18 ##########################
################################################################################################

#Set up environment
library(car)
library(Hmisc)
library(lme4)
library(lmerTest)

#Read in csv file with data
SFmid.data <- read.table(file = "StandardFitnessG1518.csv", h = T, sep = ",")
SFmidM.data <- read.table(file = "StandardFitnessG1518Mean.csv", h = T, sep = ",")

##########################################  STATISTIC  #########################################


##############################
####THE MODEL FOR ALL DATA####
##############################


### FEMALE TYPE ###

### T-TEST ###
#First test if the types are significant different within FLX & CFM
FLXsfm <- subset(SFmid.data, regime == "aFLX" & sex == "f")
CFMsfm <- subset(SFmid.data, regime == "bCFM" & sex == "f")

#Do a t-test for FLX
t.test(FLXsfm$standard ~ FLXsfm$type)
#FLX 2X and 1X are not significant different: P = 0.129

#Do a t-test for CFM
t.test(CFMsfm$standard ~ CFMsfm$type)
#CFM 2X and 1X are not significant different: P = 0.6116


##2X and 1X can be tested together


### LINEAT MODEL ###

#Linear mixed model testing if regime has a significant effect on the standard fitness
model.SFm <- lmer(standard ~ regime * sex + (1|regime:rep_population), data = SFmid.data)
#ANOVA
anova(model.SFm, test = "F")
#Regime is not significant, P = 0.8138
#Sex is not, P = 0.8199
#And the interaction is not, P = 0.1779
ranova(model.SFm)
# Nested factor is not significant, P = 0.451


#Homogeneity of variances
leveneTest(SFmid.data$standard ~ SFmid.data$regime * SFmid.data$sex)
#Not significant, assumption of homogeneity of variances is met

#Test the model by looking at the residuls
resid.SFm <- residuals(model.SFm)
hist(resid.SFm)
shapiro.test(resid.SFm)
qqnorm(resid.SFm)
qqline(resid.SFm)
#The residuls are more or less normally distrubuted, meaning that the model is valid 



###################################
####THE MODEL FOR THE MEAN DATA####
###################################


### FEMALE TYPE ###

### T-TEST ###
#First test if the types are significant different within FLX & CFM
FLXsfm <- subset(SFmidM.data, regime == "aFLX" & sex == "f")
CFMsfm <- subset(SFmidM.data, regime == "bCFM" & sex == "f")

#Do a t-test for FLX
t.test(FLXsfm$standard ~ FLXsfm$type)
#FLX 2X and 1X are not significant different: P = 0.1113

#Do a t-test for CFM
t.test(CFMsfm$standard ~ CFMsfm$type)
#CFM 2X and 1X are not significant different: P = 0.5263


##2X and 1X can be tested together


### LINEAT MODEL ###

#Linear mixed model testing if regime has a significant effect on the standard fitness
model.SFmM <- lm(standard ~ regime * sex, data = SFmidM.data)
#ANOVA
anova(model.SFmM, test = "F")
#Regime is not significant, P = 0.4814
#Sex is not, P = 0.9979
#And the interaction is not, P = 0.1121


#Homogeneity of variances
leveneTest(SFmidM.data$standard ~ SFmidM.data$regime * SFmidM.data$sex)
#Not significant, assumption of homogeneity of variances is met


#Test the model by looking at the residuls
resid.SFmM <- residuals(model.SFmM)
hist(resid.SFmM)
shapiro.test(resid.SFmM)
qqnorm(resid.SFmM)
qqline(resid.SFmM)
#The residuls are more or less normally distrubuted, meaning that the model is valid 


#########################################  PLOT DATA  #########################################


#Colours
RP1 <- "#ff766f"
RP2 <- "#ffe34d"
RP3 <- "#6da8ff"
RP4 <- "#a5f584"


### PLOT ALL DATA ###

#Add fitted values
SFmid.data$fitted <- fitted(model.SFm)

#Caluclate mean
SFmid.M <- aggregate(fitted ~ regime + sex, mean, data = SFmid.data)
#SE
SFmid.SE <- (summary(model.SFm)$coefficients[,2])


## The plot

pdf("StandardFitnessG1518.pdf", 7, 7)
par(mar = c(3, 5, 2, 2))
#First make the plot
plot(NULL, xlim = c(0.2, 3), xlab = "", xaxt = "n", ylim = c(-2.5, 2.5), ylab = "Standard Fitness", yaxt = "n",
     cex.lab = 1.8, las = 1, main = "Generation 15-18", cex.main = 1.8)
#Add points
points(jitter(rep(0.4, 80), -0.07, 0.07), SFmid.data$standard[SFmid.data$regimesex == "aFLXf"],
       pch = 16, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(0.8, 80), -0.07, 0.07), SFmid.data$standard[SFmid.data$regimesex == "bCFMf"],
       pch = 15, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(1.2, 40), -0.07, 0.07), SFmid.data$standard[SFmid.data$regimesex == "cCwtf"],
       pch = 17, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(2, 40), -0.07, 0.07), SFmid.data$standard[SFmid.data$regimesex == "dFLXm"],
       pch = 16, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(2.4, 40), -0.07, 0.07), SFmid.data$standard[SFmid.data$regimesex == "eCFMm"],
       pch = 15, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(2.8, 40), -0.07, 0.07), SFmid.data$standard[SFmid.data$regimesex == "fCwtm"],
       pch = 17, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
#Add errorbar
xSFm <- c(0.4, 0.8, 1.2, 2, 2.4, 2.8)
errbar(xSFm, SFmid.M$fitted, SFmid.M$fitted + SFmid.SE, SFmid.M$fitted - SFmid.SE,
       pch = c(16, 15, 17), cex = 2.5, lwd = 2.5, add = T)
#Add axis
axis(1, at = c(0.8, 2.4), cex.axis = 1.6, labels = c("Female", "Male"))
axis(2, at = seq(-2.5, 2.5,  by = 0.5), cex.axis = 1.3, las = 1)
#Add legend
legend("bottomleft", c("FLX", "CFM", "Cwt"), pch = c(16, 15, 17), pt.cex = 2)
dev.off()



### PLOT MEAN ###

#Add fitted values
SFmidM.data$fitted <- fitted(model.SFmM)

#Caluclate mean
SFmidM.M <- aggregate(fitted ~ regime + sex, mean, data = SFmidM.data)
#SE
BFmidM.SE <- (summary(model.SFmM)$coefficients[,2])


## The plot

pdf("StandardFitnessG1518mean.pdf", 7, 7)
par(mar = c(3, 5, 2, 2))
#First make the plot
plot(NULL, xlim = c(0.2, 3), xlab = "", xaxt = "n", ylim = c(-2.5, 2.5), ylab = "Standard Fitness", yaxt = "n",
     cex.lab = 1.8, las = 1, main = "Generation 15-18", cex.main = 1.8)
#Add points
points(jitter(rep(0.4, 8), -0.07, 0.07), SFmidM.data$standard[SFmidM.data$regimesex == "aFLXf"],
       pch = 16, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(0.8, 8), -0.07, 0.07), SFmidM.data$standard[SFmidM.data$regimesex == "bCFMf"],
       pch = 15, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(1.2, 4), -0.07, 0.07), SFmidM.data$standard[SFmidM.data$regimesex == "cCwtf"],
       pch = 17, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(2, 4), -0.07, 0.07), SFmidM.data$standard[SFmidM.data$regimesex == "dFLXm"],
       pch = 16, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(2.4, 4), -0.07, 0.07), SFmidM.data$standard[SFmidM.data$regimesex == "eCFMm"],
       pch = 15, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(2.8, 4), -0.07, 0.07), SFmidM.data$standard[SFmidM.data$regimesex == "fCwtm"],
       pch = 17, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
#Add errorbar
xSFmM <- c(0.4, 0.8, 1.2, 2, 2.4, 2.8)
errbar(xSFmM, SFmidM.M$fitted, SFmidM.M$fitted + SFmidM.SE, SFmidM.M$fitted - SFmidM.SE, 
       pch = c(16, 15, 17), cex = 2.5, lwd = 2.5, add = T)
#Add axis
axis(1, at = c(0.8, 2.4), cex.axis = 1.6, labels = c("Female", "Male"))
axis(2, at = seq(-2.5, 2.5,  by = 0.5), cex.axis = 1.3, las = 1)
#Add legend
legend("bottomleft", c("FLX", "CFM", "Cwt"), pch = c(16, 15, 17), pt.cex = 2)
dev.off()

