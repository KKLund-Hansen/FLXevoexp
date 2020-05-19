################################################################################################
############################### STANDARD FITNESS GENERATION 39-41 ##############################
################################################################################################

#Set working directory
setwd("~/Google Drev/Work/PhD/Drosophila/FLX/1.FLX/R/FitnessAssay")

#Set up environment
library(car)
library(Hmisc)
library(lme4)
library(lmerTest)

#Read in csv file with data
SF.data <- read.table(file = "StandardFitnessG3941.csv", h = T, sep = ",")
SFM.data <- read.table(file = "StandardFitnessG3941Mean.csv", h = T, sep = ",")


##########################################  STATISTIC  #########################################


##############################
####THE MODEL FOR ALL DATA####
##############################


### FEMALE TYPE ###

### T-TEST ###
#First test if the types are significant different within FLX & CFM
FLXsf <- subset(SF.data, regime == "aFLX" & sex == "f")
CFMsf <- subset(SF.data, regime == "bCFM" & sex == "f")

#Do a t-test for FLX
t.test(FLXsf$standard ~ FLXsf$type)
#FLX 2X and 1X are not significant different: P = 0.2268

#Do a t-test for CFM
t.test(CFMsf$standard ~ CFMsf$type)
#CFM 2X and 1X are not significant different: P = 0.02979


##2X and 1X can not be tested together


### LINEAT MODEL ###

#Linear mixed model testing if regime has a significant effect on standard fitness
model.SF <- lmer(standard ~ regimetype + block + (1|regimetype:rep_population), data = SF.data)
#ANOVA
anova(model.SF, test = "F")
#RegimeType is not significant, P = 0.17895
#Block is, P = 0.02253
ranova(model.SF)
#Nested factor is not significant, P = 0.09538


#Homogeneity of variances
leveneTest(SF.data$standard, SF.data$regimetype)
#Not significant, assumption of homogeneity of variances is met

#Test the model by looking at the residuls
resid.SF <- residuals(model.SF)
hist(resid.SF)
shapiro.test(resid.SF)
qqnorm(resid.SF)
qqline(resid.SF)
#The residuls are more or less normally distrubuted, meaning that the model is valid 



###################################
####THE MODEL FOR THE MEAN DATA####
###################################


### FEMALE TYPE ###

### T-TEST ###
#First test if the types are significant different within FLX & CFM
FLXsfM <- subset(SFM.data, regime == "aFLX" & sex == "f")
CFMsfM <- subset(SFM.data, regime == "bCFM" & sex == "f")

#Do a t-test for FLX
t.test(FLXsfM$standard ~ FLXsfM$type)
#FLX 2X and 1X are not significant different: P = 0.4006

#Do a t-test for CFM
t.test(CFMsfM$standard ~ CFMsfM$type)
#CFM 2X and 1X are not significant different: P = 0.2012


##2X and 1X can be tested together


### LINEAT MODEL ###

#Linear mixed model testing if regime has a significant effect on standard fitness
model.SFM <- lm(standard ~ regime * sex + block, data = SFM.data)
#ANOVA
anova(model.SFM, test = "F")
#regime is not significant, P = 0.1447
#Sex is not, P = 0.6068
#And the interaction is not, P = 0.6770
#Block is not, P = 0.1224

#Homogeneity of variances
leveneTest(SFM.data$standard ~ SFM.data$regime * SFM.data$sex)
#Not significant, assumption of homogeneity of variances is met

#Test the model by looking at the residuls
resid.SFM <- residuals(model.SFM)
hist(resid.SFM)
shapiro.test(resid.SFM)
qqnorm(resid.SFM)
qqline(resid.SFM)
#The residuls are more or less normally distrubuted, meaning that the model is valid 


#########################################  PLOT DATA  #########################################


#Colours
RP1 <- "#ff766f"
RP2 <- "#ffe34d"
RP3 <- "#6da8ff"
RP4 <- "#a5f584"


### PLOT ALL DATA ###

#Add fitted values
SF.data$fitted <- fitted(model.SF)

#Caluclate mean
SF.M <- tapply(SF.data$fitted, SF.data$regimetype, mean)
#SE
SF.SE <- (summary(model.SF)$coefficients[1:8,2])


## The plot

pdf("StandardFitnessG3941.pdf", 7, 7)
par(mar = c(3, 5, 2, 2))
#First make the plot
plot(NULL, xlim = c(0.3, 7.2), xlab = "", xaxt = "n", ylim = c(-3, 4), ylab = "Standard Fitness",
     cex.lab = 1.8, cex.axis = 1.3, las = 1, main = "Generation 39-41", cex.main = 1.8)
#Add points
points(jitter(rep(0.5, 48), -0.07, 0.07), SF.data$standard[SF.data$regimetype == "aFLX2X"],
       pch = 16, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(1, 48), -0.07, 0.07), SFM.data$standard[SF.data$regimetype == "bFLX1X"],
       pch = 1, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(2, 48), -0.07, 0.07), SF.data$standard[SF.data$regimetype == "cCFM2X"],
       pch = 15, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(2.5, 48), -0.07, 0.07), SF.data$standard[SF.data$regimetype == "dCFM1X"],
       pch = 0, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(3.5, 48), -0.07, 0.07), SF.data$standard[SF.data$regimetype == "eCwtf"],
       pch = 17, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(5, 48), -0.07, 0.07), SF.data$standard[SF.data$regimetype == "fFLXm"],
       pch = 16, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(6, 48), -0.07, 0.07), SF.data$standard[SF.data$regimetype == "gCFMm"],
       pch = 15, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(7, 48), -0.07, 0.07), SF.data$standard[SF.data$regimetype == "hCwtm"],
       pch = 17, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
#Add errorbar
xSF <- c(0.5, 1, 2, 2.5, 3.5, 5, 6, 7)
errbar(xSF, SF.M, SF.M + SF.SE, SF.M - SF.SE, 
       pch = c(16, 1, 15, 0, 17, 16, 15, 17), cex = 2.5, lwd = 2.5, add = T)
#Add axis
axis(1, at = c(2, 6), cex.axis = 1.6, labels = c("Female", "Male"))
#Add legend
legend("topleft", c("FLX", "CFM", "Cwt"), pch = c(16, 15, 17), pt.cex = 2)
dev.off()



### PLOT MEAN ###

#Add fitted values
SFM.data$fitted <- fitted(model.SFM)

#Caluclate mean
SFM.M <- aggregate(fitted ~ regime + sex, mean, data = SFM.data)
#SE
SFM.SE <- (summary(model.SFM)$coefficients[c(1:4, 6:7),2])


## The plot

pdf("StandardFitnessG3941mean.pdf", 7, 7)
par(mar = c(3, 5, 2, 2))
#First make the plot
plot(NULL, xlim = c(0.2, 3), xlab = "", xaxt = "n", ylim = c(-2.5, 2.5), ylab = "Standard Fitness", yaxt = "n",
     cex.lab = 1.8, las = 1, main = "Generation 39-41", cex.main = 1.8)
#Add points
points(jitter(rep(0.4, 8), -0.07, 0.07), SFM.data$standard[SFM.data$regimesex == "aFLXf"],
       pch = 16, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(0.8, 8), -0.07, 0.07), SFM.data$standard[SFM.data$regimesex == "bCFMf"],
       pch = 15, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(1.2, 4), -0.07, 0.07), SFM.data$standard[SFM.data$regimesex == "cCwtf"],
       pch = 17, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(2, 8), -0.07, 0.07), SFM.data$standard[SFM.data$regimesex == "dFLXm"],
       pch = 16, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(2.4, 8), -0.07, 0.07), SFM.data$standard[SFM.data$regimesex == "eCFMm"],
       pch = 15, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(2.8, 8), -0.07, 0.07), SFM.data$standard[SFM.data$regimesex == "fCwtm"],
       pch = 17, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
#Add errorbar
xSFM <- c(0.4, 0.8, 1.2, 2, 2.4, 2.8)
errbar(xSFM, SFM.M$fitted, SFM.M$fitted + SFM.SE, SFM.M$fitted - SFM.SE, 
       pch = c(16, 15, 17), cex = 2.5, lwd = 2.5, add = T)
#Add axis
axis(1, at = c(0.8, 2.4), cex.axis = 1.6, labels = c("Female", "Male"))
axis(2, at = seq(-2.5, 2.5,  by = 0.5), cex.axis = 1.3, las = 1)
#Add legend
legend("bottomleft", c("FLX", "CFM", "Cwt"), pch = c(16, 15, 17), pt.cex = 2)
dev.off()

