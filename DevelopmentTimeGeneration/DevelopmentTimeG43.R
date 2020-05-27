################################################################################################
################################# DEVELOPMENT TIME GENERATION 43 ###############################
################################################################################################

#Set up environment
library(car)
library(emmeans)
library(Hmisc)
library(lme4)
library(lmerTest)

#Read in csv file with data
DT.data <- read.table(file = "DevelopmentTimeG43.csv", h = T, sep = ",")
DTM.data <- read.table(file = "DevelopmentTimeG43mean.csv", h = T, sep = ",")


##########################################  STATISTIC  #########################################


##############################
####THE MODEL FOR ALL DATA####
##############################


### FEMALE TYPE ###

### T-TEST ###
#First test if the types are significant different within FLX & CFM
FLXdt <- subset(DT.data, regime == "aFLX" & sex == "f")
CFMdt <- subset(DT.data, regime == "bCFM" & sex == "f")

#Do a t-test for FLX
t.test(FLXdt$dev_time ~ FLXdt$type)
#FLX 2X and 1X are not significant different: P = 0.6636

#Do a t-test for CFM
t.test(CFMdt$dev_time ~ CFMdt$type)
#CFM 2X and 1X are not significant different: P = 0.9617

##2X and 1X can be tested together


### LINEAT MODEL ###

#Linear mixed model testing if regime has a significant effect on development time
model.dt <- lmer(dev_time ~ regime * sex + (1|regime:rep_population), data = DT.data)
#ANOVA
anova(model.dt, test = "F")
#Regime is not significant, P = 0.06846
#Sex is not, P = 0.06844
#And the interaction is not, P = 0.99116 


#Homogeneity of variances
leveneTest(DT.data$dev_time ~ DT.data$regime * DT.data$sex)
#Not significant, assumption of homogeneity of variances is met


#Test the model by looking at the residuls
resid.dt <- residuals(model.dt)
hist(resid.dt)
shapiro.test(resid.dt)
qqnorm(resid.dt)
qqline(resid.dt)
#The residuls are more or less normally distrubuted, meaning that the model is valid 



###################################
####THE MODEL FOR THE MEAN DATA####
###################################


### FEMALE TYPE ###

### T-TEST ###
#First test if the types are significant different within FLX & CFM
FLXdtM <- subset(DTM.data, regime == "aFLX" & sex == "f")
CFMdtM <- subset(DTM.data, regime == "bCFM" & sex == "f")

#Do a t-test for FLX
t.test(FLXdtM$dev_time ~ FLXdtM$type)
#FLX 2X and 1X are not significant different: P = 0.7891

#Do a t-test for CFM
t.test(CFMdtM$dev_time ~ CFMdtM$type)
#CFM 2X and 1X are not significant different: P = 0.9689

##2X and 1X can be tested together


### LINEAT MODEL ###

#Linear mixed model testing if regime has a significant effect on development time
model.dtM <- lm(dev_time ~ regime * sex, data = DTM.data)
#ANOVA
anova(model.dtM, test = "F")
#Regime is significant, P = 0.0006584
#But sex is not, P = 0.2789451
#And the interaction is not, P = 0.9968467

#Homogeneity of variances
leveneTest(DTM.data$dev_time ~ DTM.data$regime * DTM.data$sex)
#Significant, assumption of homogeneity of variances is not met


#Test the model by looking at the residuls
resid.dtM <- residuals(model.dtM)
hist(resid.dtM)
shapiro.test(resid.dtM)
qqnorm(resid.dtM)
qqline(resid.dtM)
#The residuls are more or less normally distrubuted, meaning that the model is valid 


#Use emmeans to test if any of regimes are significant different from each other
emmeans(model.dtM, pairwise ~ regime * sex)
#FLX and CFM females are significant from Cwt males


#########################################  PLOT DATA  #########################################


#Colours
RP1 <- "#ff766f"
RP2 <- "#ffe34d"
RP3 <- "#6da8ff"
RP4 <- "#a5f584"


### PLOT ALL DATA ###

#Add fitted values
DT.data$fitted <- fitted(model.dt)

#Caluclate mean
DT.M <- aggregate(fitted ~ regime + sex, mean, data = DT.data)
#SE
DT.SE <- (summary(model.dt)$coefficients[,2])


## The plot
pdf("DevelopmentTimeG43.pdf", 7, 7)
par(mar = c(3, 5, 2, 2))
#First make the plot
plot(NULL, xlim = c(0.2, 3), xlab = "", xaxt = "n", ylim = c(200, 280), ylab = "Development time/h",
     cex.lab = 1.8, cex.axis = 1.3, las = 1)
#Add points
points(jitter(rep(0.4, 40), -0.07, 0.07), DT.data$dev_time[DT.data$regimesex == "aFLXf"],
       pch = 16, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(0.8, 40), -0.07, 0.07), DT.data$dev_time[DT.data$regimesex == "bCFMf"],
       pch = 15, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(1.2, 20), -0.07, 0.07), DT.data$dev_time[DT.data$regimesex == "cCwtf"],
       pch = 17, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(2, 20), -0.07, 0.07), DT.data$dev_time[DT.data$regimesex == "dFLXm"],
       pch = 16, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(2.4, 20), -0.07, 0.07), DT.data$dev_time[DT.data$regimesex == "eCFMm"],
       pch = 15, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(2.8, 20), -0.07, 0.07), DT.data$dev_time[DT.data$regimesex == "fCwtm"],
       pch = 17, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
#Add errorbars
xDT <- c(0.4, 0.8, 1.2, 2, 2.4, 2.8)
errbar(xDT, DT.M$fitted, DT.M$fitted + DT.SE, DT.M$fitted - DT.SE,
       pch = c(16, 15, 17), cex = 2.5, lwd = 2.5, add = T)
#Add axis
axis(1, at = c(0.8, 2.4), cex.axis = 1.6, labels = c("Female", "Male"))
#Add legend
legend("bottomright", c("FLX", "CFM", "Cwt"), pch = c(16, 15, 17), pt.cex = 2)
dev.off()



### PLOT MEAN ###

#Add fitted values
DTM.data$fitted <- fitted(model.dtM)

#Caluclate mean
DTM.M <- aggregate(fitted ~ regime + sex, mean, data = DTM.data)
#SE
DTM.SE <- (summary(model.dtM)$coefficients[,2])


## The plot

pdf("DevelopmentTimeG43mean.pdf", 7, 7)
par(mar = c(3, 5, 2, 2))
#First make the plot
plot(NULL, xlim = c(0.2, 3), xlab = "", xaxt = "n", ylim = c(215, 255), ylab = "Development time/h", yaxt = "n",
     cex.lab = 1.8, las = 1)
#Add points
points(jitter(rep(0.4, 8), -0.07, 0.07), DTM.data$dev_time[DTM.data$regimesex == "aFLXf"],
       pch = 16, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(0.8, 8), -0.07, 0.07), DTM.data$dev_time[DTM.data$regimesex == "bCFMf"],
       pch = 15, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(1.2, 4), -0.07, 0.07), DTM.data$dev_time[DTM.data$regimesex == "cCwtf"],
       pch = 17, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(2, 4), -0.07, 0.07), DTM.data$dev_time[DTM.data$regimesex == "dFLXm"],
       pch = 16, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(2.4, 4), -0.07, 0.07), DTM.data$dev_time[DTM.data$regimesex == "eCFMm"],
       pch = 15, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(2.8, 4), -0.07, 0.07), DTM.data$dev_time[DTM.data$regimesex == "fCwtm"],
       pch = 17, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
#Add errorbar
xDTM <- c(0.4, 0.8, 1.2, 2, 2.4, 2.8)
errbar(xDTM, DTM.M$fitted, DTM.M$fitted + DTM.SE, DTM.M$fitted - DTM.SE, 
       pch = c(16, 15, 17), cex = 2.5, lwd = 2.5, add = T)
#Add axis
axis(1, at = c(0.8, 2.4), cex.axis = 1.6, labels = c("Female", "Male"))
axis(2, at = seq(215, 255,  by = 5), cex.axis = 1.3, las = 1)
#Add legend
legend("bottomright", c("FLX", "CFM", "Cwt"), pch = c(16, 15, 17), pt.cex = 2)
#Mark significant between the treatments with letters
mtext("b", side = 3, line = -2, at = 0.4, cex = 1.5)
mtext("b", side = 3, line = -2, at = 0.8, cex = 1.5)
mtext("ab", side = 3, line = -2, at = 1.2, cex = 1.5)
mtext("ab", side = 3, line = -2, at = 2, cex = 1.5)
mtext("ab", side = 3, line = -2, at = 2.4, cex = 1.5)
mtext("a", side = 3, line = -2, at = 2.8, cex = 1.5)
dev.off()

