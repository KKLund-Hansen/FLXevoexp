################################################################################################
############################## LOCOMOTION ACTIVITY GENERATION 123 ##############################
################################################################################################

#Set up environment
library(blmeco)
library(car)
library(emmeans)
library(Hmisc)
library(lme4)

#Read in csv file with data
LocoG123.data <- read.table(file = "LocoAssayG123.csv", h = T, sep = ",")
LocoG123m.data <- read.table(file = "LocoAssayG123mean.csv", h = T, sep = ",")


##########################################  STATISTIC  #########################################


##############################
####THE MODEL FOR ALL DATA####
##############################


### FEMALE TYPE ###

### T-TEST ###
#First test if the types are significant different within FLX & CFM
FLXlocog123 <- subset(LocoG123.data, regime == "aFLX" & sex == "f")
CFMlocog123 <- subset(LocoG123.data, regime == "bCFM" & sex == "f")

#Do a t-test for FLX
t.test((cbind(FLXlocog123$active, FLXlocog123$inactive)) ~ FLXlocog123$type)
#FLX 2X and 1X are not significant different: P = 1

#Do a t-test for CFM
t.test((cbind(CFMlocog123$active, CFMlocog123$inactive)) ~ CFMlocog123$type)
#CFM 2X and 1X are not significant different: P = 1

##2X and 1X can be tested together


### LINEAT MODEL ###

#General linear mixed model testing if regime has a significant effect on locomotion
model.locog123 <- glmer((cbind(active, inactive)) ~ regime * sex + (1|regime:rep_population), 
                          data = LocoG123.data, family = "binomial")

#Test for overdispersion
dispersion_glmer(model.locog123)
#If between 0.75 and 1.4, there may not be an overdispersion problem
#1.267359, could work, but I'll to add an extra observation
LocoG123.data$obs <- 1:nrow(LocoG123.data)

#General linear mixed model testing if regime has a significant effect on locomotion
model.locog123 <- glmer((cbind(active, inactive)) ~ regime * sex + (1|regime:rep_population) + (1|obs),
                         data = LocoG123.data, family = "binomial")
#ANOVA
Anova(model.locog123, test = "Chisq")
#Regime is not significant, P = 0.533194
#Sex is significant, P = 0.001639 
#Interaction is not significany, P = 0.199746


#Homogeneity of variances
leveneTest(LocoG123.data$locoactive ~ LocoG123.data$regime * LocoG123.data$sex)
#Significant, assumption of homogeneity of variances is not met


#Test the model by looking at the residuls
resid.locog123 <- residuals(model.locog123)
hist(resid.locog123)
shapiro.test(resid.locog123)
qqnorm(resid.locog123)
qqline(resid.locog123)
#The residuls are more or less normally distrubuted, meaning that the model is valid 


#Use emmeans to test if any of regimes are significant different from each other
emmeans(model.locog123, pairwise ~ regime * sex)
#FLX females are significant different from FLX and Cwt males



###################################
####THE MODEL FOR THE MEAN DATA####
###################################


### FEMALE TYPE ###

### T-TEST ###
#First test if the types are significant different within FLX & CFM
FLXlocog123m <- subset(LocoG123m.data, regime == "aFLX" & sex == "f")
CFMlocog123m <- subset(LocoG123m.data, regime == "bCFM" & sex == "f")

#Do a t-test for FLX
t.test((cbind(FLXlocog123m$active, FLXlocog123m$inactive)) ~ FLXlocog123m$type)
#FLX 2X and 1X are not significant different: P = 1

#Do a t-test for CFM
t.test((cbind(CFMlocog123m$active, CFMlocog123m$inactive)) ~ CFMlocog123m$type)
#CFM 2X and 1X are not significant different: P = 1

##2X and 1X can be tested together


### LINEAT MODEL ###

#General linear mixed model testing if regime has a significant effect on locomotion
model.locog123m <- glm((cbind(active, inactive)) ~ regime * sex, data = LocoG123m.data, family = "binomial")

#First look at the residuals
summary(model.locog123m)
#Residual deviance: 21.186 on 26 degrees of freedom - okay 

#ANOVA
Anova(model.locog123m, test = "Wald")
#Treatment is not significant, P = 0.3573341
#Sex is significant, P = 0.0001102
#Interaction is significany, P = 0.0911764


#Homogeneity of variances
leveneTest(LocoG123m.data$locoactive ~ LocoG123m.data$regime * LocoG123m.data$sex)
#Not significant, assumption of homogeneity of variances is met


#Test the model by looking at the residuls
resid.locog123m <- residuals(model.locog123m)
hist(resid.locog123m)
shapiro.test(resid.locog123m)
qqnorm(resid.locog123m)
qqline(resid.locog123m)
#The residuls are more or less normally distrubuted, meaning that the model is valid 


#Use emmeans to test if any of regimes are significant different from each other
emmeans(model.locog123m, pairwise ~ regime * sex)
#FLX/CFM females are significant different from FLX and Cwt males



#########################################  PLOT DATA  #########################################


#Colours
RP1 <- "#ff766f"
RP2 <- "#ffe34d"
RP3 <- "#6da8ff"
RP4 <- "#a5f584"
### PLOT ALL DATA ###


#Add fitted values
LocoG123.data$fitted <- fitted(model.locog123)

#Caluclate mean
LocoG123.M <- aggregate(fitted ~ regime + sex, mean, data = LocoG123.data)
#SE
LocoG123.SE <- (summary(model.locog123)$coefficients[,2])


## The plot

pdf("LocoG123.pdf", 7, 7)
par(mar = c(3, 5, 2, 2))
#First make the plot
plot(NULL, xlim = c(0.2, 3), xlab = "", xaxt = "n", ylim = c(-0.1, 1.1), 
     ylab = "Locomotion Activity",  yaxt = "n", cex.lab = 1.8, las = 1)
#Add points
points(jitter(rep(0.4, 40), -0.07, 0.07), LocoG123.data$locoactive[LocoG123.data$regimesex == "aFLXf"],
       pch = 16, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(0.8, 40), -0.07, 0.07), LocoG123.data$locoactive[LocoG123.data$regimesex == "bCFMf"],
       pch = 15, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(1.2, 20), -0.07, 0.07), LocoG123.data$locoactive[LocoG123.data$regimesex == "cCwtf"],
       pch = 17, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(2, 20), -0.07, 0.07), LocoG123.data$locoactive[LocoG123.data$regimesex == "dFLXm"],
       pch = 16, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(2.4, 20), -0.07, 0.07), LocoG123.data$locoactive[LocoG123.data$regimesex == "eCFMm"],
       pch = 15, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(2.8, 20), -0.07, 0.07), LocoG123.data$locoactive[LocoG123.data$regimesex == "fCwtm"],
       pch = 17, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
#Add errorbar
xLocoG123 <- c(0.4, 0.8, 1.2, 2, 2.4, 2.8)
errbar(xLocoG123, LocoG123.M$fitted, LocoG123.M$fitted + LocoG123.SE, LocoG123.M$fitted - LocoG123.SE, 
       cex = 3, lwd = 3, pch = c(16, 15, 17), add = T)
#Add axis
axis(1, at = c(0.8, 2.4), cex.axis = 1.5, labels = c("Female", "Male"))
axis(2, at = seq(-0.1, 1.1,  by = 0.2), cex.axis = 1.3, las = 1)
#Add legend
legend("bottomleft", c("FLX", "CFM", "Cwt"), pch = c(16, 15, 17), pt.cex = 2)
#Mark significant between the treatments with letters
mtext("b", side = 3, line = -1.8, at = 0.4, cex = 1.5)
mtext("ab", side = 3, line = -1.8, at = 0.8, cex = 1.5)
mtext("ab", side = 3, line = -1.8, at = 1.2, cex = 1.5)
mtext("a", side = 3, line = -1.8, at = 2, cex = 1.5)
mtext("ab", side = 3, line = -1.8, at = 2.4, cex = 1.5)
mtext("a", side = 3, line = -1.8, at = 2.8, cex = 1.5)
dev.off()



### PLOT MEAN ###

#Add fitted values
LocoG123m.data$fitted <- fitted(model.locog123m)

#Caluclate mean
LocoG123m.M <- aggregate(fitted ~ regime + sex, mean, data = LocoG123m.data)
#SE
LocoG123m.SE <- (summary(model.locog123m)$coefficients[,2])


## The plot

pdf("LocoG123mean.pdf", 7, 7)
par(mar = c(3, 5, 2, 2))
#First make the plot
plot(NULL, xlim = c(0.2, 3), xlab = "", xaxt = "n", ylim = c(0, 0.7), ylab = "Locomotion Activity",
     cex.lab = 1.8, cex.axis = 1.3, las = 1)
#Add points
points(jitter(rep(0.4, 8), -0.07, 0.07), LocoG123m.data$locoactive[LocoG123m.data$regimesex == "aFLXf"],
       pch = 16, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(0.8, 8), -0.07, 0.07), LocoG123m.data$locoactive[LocoG123m.data$regimesex == "bCFMf"],
       pch = 15, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(1.2, 4), -0.07, 0.07), LocoG123m.data$locoactive[LocoG123m.data$regimesex == "cCwtf"],
       pch = 17, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(2, 4), -0.07, 0.07), LocoG123m.data$locoactive[LocoG123m.data$regimesex == "dFLXm"],
       pch = 16, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(2.4, 4), -0.07, 0.07), LocoG123m.data$locoactive[LocoG123m.data$regimesex == "eCFMm"],
       pch = 15, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(2.8, 4), -0.07, 0.07), LocoG123m.data$locoactive[LocoG123m.data$regimesex == "fCwtm"],
       pch = 17, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
#Add errorbar
xLocoG123m <- c(0.4, 0.8, 1.2, 2, 2.4, 2.8)
errbar(xLocoG123m, LocoG123m.M$fitted, LocoG123m.M$fitted + LocoG123m.SE, LocoG123m.M$fitted - LocoG123m.SE, 
       pch = c(16, 15, 17), cex = 2.5, lwd = 2.5, add = T)
#Add axis
axis(1, at = c(0.8, 2.4), cex.axis = 1.6, labels = c("Female", "Male"))
#Add legend
legend("bottomleft", c("FLX", "CFM", "Cwt"), pch = c(16, 15, 17), pt.cex = 2)
#Mark significant between the treatments with letters
mtext("b", side = 3, line = -1.8, at = 0.4, cex = 1.5)
mtext("b", side = 3, line = -1.8, at = 0.8, cex = 1.5)
mtext("ab", side = 3, line = -1.8, at = 1.2, cex = 1.5)
mtext("a", side = 3, line = -1.8, at = 2, cex = 1.5)
mtext("ab", side = 3, line = -1.8, at = 2.4, cex = 1.5)
mtext("a", side = 3, line = -1.8, at = 2.8, cex = 1.5)
dev.off()

