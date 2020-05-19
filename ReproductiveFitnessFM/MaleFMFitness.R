################################################################################################
################## FM BALANCER MALE REPRODUCTIVE FITNESS ASSAY GENERATION 126 ##################
################################################################################################

#Set working directory
setwd("~/Google Drev/Work/PhD/Drosophila/FLX/1.FLX/R/FMFitnessAssay")

#Set up environment
library(car)
library(Hmisc)
library(lme4)
library(lmerTest)
library(multcomp)

#Read in csv file with data
MFMF.data <- read.table(file = "MaleFMFitness.csv", h = T, sep = ",")
MFMFm.data <- read.table(file = "MaleFMFitnessMean.csv", h = T, sep = ",")


##########################################  STATISTIC  #########################################


##############################
####THE MODEL FOR ALL DATA####
##############################

### CALCULATE RELATIVE FITNESS ###
#Calculate the proportion of red eyed offspring
MFMF.data$prop_red <- MFMF.data$red / MFMF.data$total
#Now divid each proportion red by 5 to get the number for each male
MFMF.data$prop_red_male <- MFMF.data$prop_red / 5
#Find maxium
MFMF.max <- max(MFMF.data$prop_red_male)
#Calculate relative fitness by dividing each proportion by the maximum
MFMF.data$relative_fit <- MFMF.data$prop_red_male / MFMF.max


### LINEAT MODEL ###

#Linear mixed model testing if regime type has a significant effect on the relative fitness
model.MFMF <- lmer(relative_fit ~ regimetype + block + (1|regimetype:rep_population), data = MFMF.data)
#ANOVA
anova(model.MFMF, test = "F")
#Regime/Type is significant, P = 3.979e-09
#Block is, P = 0.0388
ranova(model.MFMF)
#And the nested factor is, P = 0.001203


#Homogeneity of variances
leveneTest(MFMF.data$relative_fit, MFMF.data$regimetype)
#Significant, assumption of homogeneity of variances is not met

#Test the model by looking at the residuls
resid.MFMF <- residuals(model.MFMF)
hist(resid.MFMF)
shapiro.test(resid.MFMF)
qqnorm(resid.MFMF)
qqline(resid.MFMF)
#The residuls are more or less normally distrubuted, meaning that the model is valid 


#I do a post hoc test to see which regimes are significantly different
RESULT.MFMF <- glht(model.MFMF, mcp(regimetype = "Tukey"))
#Look at the results
summary(RESULT.MFMF)
#FM is significant different from none-FM


### PERCENTAGE CHANGE ###

#First calculate mean

FMmean <- tapply(MFMF.data$relative_fit, MFMF.data$regimetype, mean)

#FLX
((FMmean[1]-FMmean[2])/FMmean[1])*100
#64.33545

#CFM
((FMmean[3]-FMmean[4])/FMmean[3])*100
#80.78323

#Cwt
((FMmean[5]-FMmean[6])/FMmean[5])*100
#84.38201



###################################
####THE MODEL FOR THE MEAN DATA####
###################################

### CALCULATE RELATIVE FITNESS ###
#Calculate the proportion of red eyed offspring
MFMFm.data$prop_red <- MFMFm.data$red / MFMFm.data$total
#Now divid each proportion red by 5 to get the number for each male
MFMFm.data$prop_red_male <- MFMFm.data$prop_red / 5
#Find maxium
MFMFm.max <- max(MFMFm.data$prop_red_male)
#Calculate relative fitness by dividing each proportion by the maximum
MFMFm.data$relative_fit <- MFMFm.data$prop_red_male / MFMFm.max


### LINEAT MODEL ###

#Linear model testing if regime type has a significant effect on the relative fitness
model.MFMFm <- lm(relative_fit ~ regimetype + block, data = MFMFm.data)
#ANOVA
anova(model.MFMFm, test = "F")
#Regime/type is significant, P < 2.2e-16
#Block is not, P = 0.1551


#Homogeneity of variances
leveneTest(MFMFm.data$relative_fit, MFMFm.data$regimetype)
#Significant, assumption of homogeneity of variances is not met

#Test the model by looking at the residuls
resid.MFMFm <- residuals(model.MFMFm)
hist(resid.MFMFm)
shapiro.test(resid.MFMFm)
qqnorm(resid.MFMFm)
qqline(resid.MFMFm)
#The residuls are more or less normally distrubuted, meaning that the model is valid 

#I do a post hoc test to see which regimes are significantly different
RESULT.MFMFm <- glht(model.MFMFm, mcp(regimetype = "Tukey"))
#Look at the results
summary(RESULT.MFMFm)
#FM is significant different from none-FM


### PERCENTAGE CHANGE ###

#First calculate mean

FMmmean <- tapply(MFMFm.data$relative_fit, MFMFm.data$regimetype, mean)

#FLX
((FMmmean[1]-FMmmean[2])/FMmmean[1])*100
#64.8554

#CFM
((FMmmean[3]-FMmmean[4])/FMmmean[3])*100
#80.39575

#Cwt
((FMmmean[5]-FMmmean[6])/FMmmean[5])*100
#84.51666


#########################################  PLOT DATA  #########################################


#Colours
RP1 <- "#ff766f"
RP2 <- "#ffe34d"
RP3 <- "#6da8ff"
RP4 <- "#a5f584"


### PLOT ALL DATA ###

#Add fitted values
MFMF.data$fitted <- fitted(model.MFMF)

#Caluclate mean
MFMF.M <- tapply(MFMF.data$fitted, MFMF.data$regimetype, mean)
#SE
MFMF.SE <- (summary(model.MFMF)$coefficients[1:6, 2])


## The plot

pdf("MFMF.pdf", 7, 7)
par(mar = c(3, 5, 2, 2))
#First make the plot
plot(NULL, xlim = c(0, 4.5), xlab = "", xaxt = "n", ylim = c(-0.1, 1.1), ylab = "Relative fitness", yaxt = "n",
     cex.lab = 1.8, cex.axis = 1.3, las = 1, main = "Males", cex.main = 1.8)
#Add points
points(jitter(rep(0.5, 40), -0.07, 0.07), MFMF.data$relative_fit[MFMF.data$regimetype == "aFLXXY"],
       pch = 16, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(1, 40), -0.07, 0.07), MFMF.data$relative_fit[MFMF.data$regimetype == "bFLXFMY"],
       pch = 1, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(2, 40), -0.07, 0.07), MFMF.data$relative_fit[MFMF.data$regimetype == "cCFMXY"],
       pch = 15, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(2.5, 40), -0.07, 0.07), MFMF.data$relative_fit[MFMF.data$regimetype == "dCFMFMY"],
       pch = 0, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(3.5, 40), -0.07, 0.07), MFMF.data$relative_fit[MFMF.data$regimetype == "eCwtXY"],
       pch = 17, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(4, 40), -0.07, 0.07), MFMF.data$relative_fit[MFMF.data$regimetype == "fCwtFMY"],
       pch = 2, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
#Add errorbar
xMFMF <- c(0.5, 1, 2, 2.5, 3.5, 4)
errbar(xMFMF, MFMF.M, MFMF.M + MFMF.SE, MFMF.M - MFMF.SE,
       pch = c(16, 1, 15, 0, 17, 2), cex = 2.5, lwd = 2.5, add = T)
#Add axis
axis(1, at = c(0.75, 2.25, 3.75), cex.axis = 1.5, labels = c("FLX", "Control FM", "Control wt"))
axis(2, at = seq(-0.1, 1.1,  by = 0.2), cex.axis = 1.3, las = 1)
#Add legend
legend("bottomleft", c("Closed symbols: X/X", "Open symbols: FM/X"), adj = c(0.1, 0.5), cex = 1.1)
#Add significance
mtext("a", side = 3, line = -2, at = 0.5, cex = 1.5)
mtext("b", side = 3, line = -2, at = 1, cex = 1.5)
mtext("a", side = 3, line = -2, at = 2, cex = 1.5)
mtext("b", side = 3, line = -2, at = 2.5, cex = 1.5)
mtext("a", side = 3, line = -2, at = 3.5, cex = 1.5)
mtext("b", side = 3, line = -2, at = 4, cex = 1.5)
dev.off()



### PLOT MEAN ###

#Add fitted values
MFMFm.data$fitted <- fitted(model.MFMFm)

#Caluclate mean
MFMFm.M <- tapply(MFMFm.data$fitted, MFMFm.data$regimetype, mean)
#SE
MFMFm.SE <- (summary(model.MFMFm)$coefficients[1:6, 2])


## The plot

pdf("MFMFm.pdf", 7, 7)
par(mar = c(3, 5, 2, 2))
#First make the plot
plot(NULL, xlim = c(0, 4.5), xlab = "", xaxt = "n", ylim = c(-0.1, 1.1), ylab = "Relative fitness", yaxt = "n",
     cex.lab = 1.8, cex.axis = 1.3, las = 1, main = "Males", cex.main = 1.8)
#Add points
points(jitter(rep(0.5, 8), -0.07, 0.07), MFMFm.data$relative_fit[MFMFm.data$regimetype == "aFLXXY"],
       pch = 16, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(1, 8), -0.07, 0.07), MFMFm.data$relative_fit[MFMFm.data$regimetype == "bFLXFMY"],
       pch = 1, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(2, 8), -0.07, 0.07), MFMFm.data$relative_fit[MFMFm.data$regimetype == "cCFMXY"],
       pch = 15, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(2.5, 8), -0.07, 0.07), MFMFm.data$relative_fit[MFMFm.data$regimetype == "dCFMFMY"],
       pch = 0, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(3.5, 8), -0.07, 0.07), MFMFm.data$relative_fit[MFMFm.data$regimetype == "eCwtXY"],
       pch = 17, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(4, 8), -0.07, 0.07), MFMFm.data$relative_fit[MFMFm.data$regimetype == "fCwtFMY"],
       pch = 2, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
#Add errorbar
xMFMFm <- c(0.5, 1, 2, 2.5, 3.5, 4)
errbar(xMFMFm, MFMFm.M, MFMFm.M + MFMFm.SE, MFMFm.M - MFMFm.SE,
       pch = c(16, 1, 15, 0, 17, 2), cex = 2.5, lwd = 2.5, add = T)
#Add axis
axis(1, at = c(0.75, 2.25, 3.75), cex.axis = 1.5, labels = c("FLX", "Control FM", "Control wt"))
axis(2, at = seq(-0.1, 1.1,  by = 0.2), cex.axis = 1.3, las = 1)
#Add legend
legend("bottomleft", c("Closed symbols: X/X", "Open symbols: FM/X"), adj = c(0.1, 0.5), cex = 1.1)
#Add significance
mtext("a", side = 3, line = -2, at = 0.5, cex = 1.5)
mtext("b", side = 3, line = -2, at = 1, cex = 1.5)
mtext("a", side = 3, line = -2, at = 2, cex = 1.5)
mtext("b", side = 3, line = -2, at = 2.5, cex = 1.5)
mtext("a", side = 3, line = -2, at = 3.5, cex = 1.5)
mtext("b", side = 3, line = -2, at = 4, cex = 1.5)
dev.off()

