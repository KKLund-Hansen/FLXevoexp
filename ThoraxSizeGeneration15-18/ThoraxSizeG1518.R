################################################################################################
############################### THORAX LENGTH GENERATION 15 & 18 ###############################
################################################################################################

#Set working directory
setwd("~/Google Drev/Work/PhD/Drosophila/FLX/1.FLX/R/MidwayAssays/ThoraxSize")

#Set up environment
library(car)
library(emmeans)
library(Hmisc)
library(lme4)
library(lmerTest)
library(multcomp)

#Read in csv file with data
TSmid.data <- read.table(file = "ThoraxSizeG1518.csv", h = T, sep = ",")
TSmidM.data <- read.table(file = "ThoraxSizeG1518Mean.csv", h = T, sep = ",")


##########################################  STATISTIC  #########################################


##############################
####THE MODEL FOR ALL DATA####
##############################


### FEMALE TYPE ###

### T-TEST ###
#First test if the types are significant different within FLX & CFM
FLXtsm <- subset(TSmid.data, regime == "aFLX" & sex == "f")
CFMtsm <- subset(TSmid.data, regime == "bCFM" & sex == "f")

#Do a t-test for FLX
t.test(FLXtsm$size_mm ~ FLXtsm$type)
#FLX 2X and 1X are significant different: P = 0.003386

#Do a t-test for CFM
t.test(CFMbtm$size_mm ~ CFMbtm$type)
#CFM 2X and 1X are not significant different: P = 0.7067


##2X and 1X can not be tested together


### LINEAT MODEL ###

#Linear mixed model testing if regime has a significant effect on size
model.TSmid <- lmer(size_mm ~ regimetype + (1|regimetype:rep_population), data = TSmid.data)
#ANOVA
anova(model.TSmid, test = "F")
#RegimeType is significant, P < 2.2e-16
ranova(model.TSmid)
#Nested factor is significant, P < 2.2e-16


#Homogeneity of variances
leveneTest(TSmid.data$size_mm, TSmid.data$regimetype)
#Significant, assumption of homogeneity of variances is not met

#Test the model by looking at the residuls
resid.TSmid <- residuals(model.TSmid)
hist(resid.TSmid)
shapiro.test(resid.TSmid)
qqnorm(resid.TSmid)
qqline(resid.TSmid)
#The residuls are more or less normally distrubuted, meaning that the model is valid 


#Do a Tukey test
RESULT.TSmid <- glht(model.TSmid, linfct = mcp(regimetype = "Tukey"))
#Look at the results
summary(RESULT.TSmid)
#All females and males are significant different from each other
#FLX males are significant different from the other two
#FLX 2X and 1X are significant different from the other female types, but not each other



###################################
####THE MODEL FOR THE MEAN DATA####
###################################


### FEMALE TYPE ###

### T-TEST ###
#First test if the types are significant different within FLX & CFM
FLXtsmM <- subset(TSmidM.data, regime == "aFLX" & sex == "f")
CFMtsmM <- subset(TSmidM.data, regime == "bCFM" & sex == "f")

#Do a t-test for FLX
t.test(FLXtsmM$size_mm ~ FLXtsmM$type)
#FLX 2X and 1X are not significant different: P = 0.3036

#Do a t-test for CFM
t.test(CFMtsmM$size_mm ~ CFMtsmM$type)
#CFM 2X and 1X are not significant different: P = 0.8002


##2X and 1X can be tested together


### LINEAT MODEL ###

#Linear mixed model testing if regime has a significant effect on size
model.TSmidM <- lm(size_mm ~ regime * sex, data = TSmidM.data)
#ANOVA
anova(model.TSmidM, test = "F")
#Regime is significant, P = 3.223e-14
#Sex is significant, P < 2.2e-16
#And the interaction is not, P = 0.1118


#Homogeneity of variances
leveneTest(TSmidM.data$size_mm ~ TSmidM.data$regime * TSmidM.data$sex)
#Not significant, assumption of homogeneity of variances is met


#Test the model by looking at the residuls
resid.TSmidM <- residuals(model.TSmidM)
hist(resid.TSmidM)
shapiro.test(resid.TSmidM)
qqnorm(resid.TSmidM)
qqline(resid.TSmidM)
#The residuls are more or less normally distrubuted, meaning that the model is valid 


#Use emmeans to test if any of treatments are significant different from each other
emmeans(model.TSmidM, pairwise ~ regime * sex)
#All females and males are significant different from each other
#FLX males are significant different from Cwt
#FLX females are significant different from the other females



#########################################  PLOT DATA  #########################################


#Colours
RP1 <- "#ff766f"
RP2 <- "#ffe34d"
RP3 <- "#6da8ff"
RP4 <- "#a5f584"


### PLOT ALL DATA ###

#Add fitted values
TSmid.data$fitted <- fitted(model.TSmid)

#Caluclate mean
TSmid.M <- tapply(TSmid.data$fitted, TSmid.data$regimetype, mean)
#SE
TSmid.SE <- (summary(model.TSmid)$coefficients[,2])


## The plot

pdf("ThoraxSizeG1518.pdf", 7, 7)
par(mar = c(3, 5, 2, 2))
#First make the plot
plot(NULL, xlim = c(0.3, 7.2), xlab = "", xaxt = "n", ylim = c(0.75, 1.1), ylab = "Thorax size/mm", 
     cex.lab = 1.8, cex.axis = 1.3, las = 1, main = "Generation 15-18", cex.main = 1.8)
#Add points
points(jitter(rep(0.5, 134), -0.07, 0.07), TSmid.data$size_mm[TSmid.data$regimetype == "aFLX2X"],
       pch = 16, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(1, 127), -0.07, 0.07), TSmid.data$size_mm[TSmid.data$regimetype == "bFLX1X"],
       pch = 1, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(2, 132), -0.07, 0.07), TSmid.data$size_mm[TSmid.data$regimetype == "cCFM2X"],
       pch = 15, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(2.5, 125), -0.07, 0.07), TSmid.data$size_mm[TSmid.data$regimetype == "dCFM1X"],
       pch = 0, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(3.5, 133), -0.07, 0.07), TSmid.data$size_mm[TSmid.data$regimetype == "eCwtf"],
       pch = 17, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(5, 151), -0.07, 0.07), TSmid.data$size_mm[TSmid.data$regimetype == "fFLXm"],
       pch = 16, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(6, 145), -0.07, 0.07), TSmid.data$size_mm[TSmid.data$regimetype == "gCFMm"],
       pch = 15, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(7, 150), -0.07, 0.07), TSmid.data$size_mm[TSmid.data$regimetype == "hCwtm"],
       pch = 17, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
#Add errorbar
xTSm <- c(0.5, 1, 2, 2.5, 3.5, 5, 6, 7)
errbar(xTSm, TSmid.M, TSmid.M + TSmid.SE, TSmid.M - TSmid.SE, 
       cex = 3, lwd = 3, pch = c(16, 1, 15, 0, 17, 16, 15, 17), add = T)
#Add axis
axis(1, at = c(2, 6), cex.axis = 1.6, labels = c("Female", "Male"))
#Mark significant between the treatments with letters
mtext("a", side = 3, line = -2, at = 0.5, cex = 1.5)
mtext("a", side = 3, line = -2, at = 1, cex = 1.5)
mtext("b", side = 3, line = -2, at = 2, cex = 1.5)
mtext("b", side = 3, line = -2, at = 2.5, cex = 1.5)
mtext("b", side = 3, line = -2, at = 3.5, cex = 1.5)
mtext("c", side = 3, line = -2, at = 5, cex = 1.5)
mtext("d", side = 3, line = -2, at = 6, cex = 1.5)
mtext("d", side = 3, line = -2, at = 7, cex = 1.5)
#Add legend
legend("bottomleft", c("FLX", "CFM", "Cwt"), pch = c(16, 15, 17), pt.cex = 2)
dev.off()



### PLOT MEAN ###

#Add fitted values
TSmidM.data$fitted <- fitted(model.TSmidM)

#Caluclate mean
TSmidM.M <- aggregate(fitted ~ regime + sex, mean, data = TSmidM.data)
#SE
TSmidM.SE <- (summary(model.TSmidM)$coefficients[,2])


## The plot

pdf("ThoraxSizeG1518mean.pdf", 7, 7)
par(mar = c(3, 5, 2, 2))
#First make the plot
plot(NULL, xlim = c(0.2, 3), xlab = "", xaxt = "n", ylim = c(0.8, 1.05), ylab = "Thorax size/mm", 
     cex.lab = 1.8, cex.axis = 1.3, las = 1, main = "Generation 15-18", cex.main = 1.8)
#Add points
points(jitter(rep(0.4, 8), -0.07, 0.07), TSmidM.data$size_mm[TSmidM.data$regimesex == "aFLXf"],
       pch = 16, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(0.8, 8), -0.07, 0.07), TSmidM.data$size_mm[TSmidM.data$regimesex == "bCFMf"],
       pch = 15, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(1.2, 4), -0.07, 0.07), TSmidM.data$size_mm[TSmidM.data$regimesex == "cCwtf"],
       pch = 17, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(2, 4), -0.07, 0.07), TSmidM.data$size_mm[TSmidM.data$regimesex == "dFLXm"],
       pch = 16, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(2.4, 4), -0.07, 0.07), TSmidM.data$size_mm[TSmidM.data$regimesex == "eCFMm"],
       pch = 15, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(2.8, 4), -0.07, 0.07), TSmidM.data$size_mm[TSmidM.data$regimesex == "fCwtm"],
       pch = 17, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
#Add errorbar
xTSmM <- c(0.4, 0.8, 1.2, 2, 2.4, 2.8)
errbar(xTSmM, TSmidM.M$fitted, TSmidM.M$fitted + TSmidM.SE, TSmidM.M$fitted - TSmidM.SE, 
       pch = c(16, 15, 17), cex = 2.5, lwd = 2.5, add = T)
#Add axis
axis(1, at = c(0.8, 2.4), cex.axis = 1.6, labels = c("Female", "Male"))
#Add legend
legend("bottomleft", c("FLX", "CFM", "Cwt"), pch = c(16, 15, 17), pt.cex = 2)
#Mark significant between the treatments with letters
mtext("a", side = 3, line = -2, at = 0.4, cex = 1.5)
mtext("b", side = 3, line = -2, at = 0.8, cex = 1.5)
mtext("b", side = 3, line = -2, at = 1.2, cex = 1.5)
mtext("c", side = 3, line = -2, at = 2, cex = 1.5)
mtext("cd", side = 3, line = -2, at = 2.4, cex = 1.5)
mtext("d", side = 3, line = -2, at = 2.8, cex = 1.5)
dev.off()

