################################################################################################
################################# THORAX LENGTH GENERATION 72 ##################################
################################################################################################

#Set up environment
library(car)
library(emmeans)
library(Hmisc)
library(lme4)
library(lmerTest)
library(multcomp)

#Read in csv file with data
BT.data <- read.table(file = "BothThorax.csv", h = T, sep = ",")
BTM.data <- read.table(file = "BothThoraxMean.csv", h = T, sep = ",")


##########################################  STATISTIC  #########################################


##############################
####THE MODEL FOR ALL DATA####
##############################


### FEMALE TYPE ###

### T-TEST ###
#First test if the types are significant different within FLX & CFM
FLXbt <- subset(BT.data, treatment == "aFLX" & sex == "f")
CFMbt <- subset(BT.data, treatment == "bCFM" & sex == "f")

#Do a t-test for FLX
t.test(FLXbt$size_mm ~ FLXbt$type)
#FLX 2X and 1X are not significant different: P = 0.0005261

#Do a t-test for CFM
t.test(CFMbt$size_mm ~ CFMbt$type)
#CFM 2X and 1X are not significant different: P = 1


##2X and 1X can not be tested together


### LINEAT MODEL ###

#Linear mixed model testing if treatment has a significant effect on size
model.BT <- lmer(size_mm ~ treattype + (1|treattype:rep_population), data = BT.data)
#ANOVA
anova(model.BT, test = "F")
#Treattrype is significant, P = 4.433e-15
ranova(model.BT)
#Nested factor is significant, P = 3.58e-15


#Homogeneity of variances
leveneTest(BT.data$size_mm, BT.data$treattype)
#Significant, assumption of homogeneity of variances is not met

#Test the model by looking at the residuls
resid.BT <- residuals(model.BT)
hist(resid.BT)
shapiro.test(resid.BT)
qqnorm(resid.BT)
qqline(resid.BT)
#The residuls are more or less normally distrubuted, meaning that the model is valid 


#Do a Tukey test
RESULT.BT <- glht(model.BT, linfct = mcp(treattype = "Tukey"))
#Look at the results
summary(RESULT.BT)
#All females and males are significant different from each other
#FLX 2X are significant different from Cwt



###################################
####THE MODEL FOR THE MEAN DATA####
###################################


### FEMALE TYPE ###

### T-TEST ###
#First test if the types are significant different within FLX & CFM
FLXbtM <- subset(BTM.data, treatment == "aFLX" & sex == "f")
CFMbtM <- subset(BTM.data, treatment == "bCFM" & sex == "f")

#Do a t-test for FLX
t.test(FLXbtM$size_mm ~ FLXbtM$type)
#FLX 2X and 1X are not significant different: P = 0.15

#Do a t-test for CFM
t.test(CFMbtM$size_mm ~ CFMbtM$type)
#CFM 2X and 1X are not significant different: P = 1


##2X and 1X can be tested together


### LINEAT MODEL ###

#Linear mixed model testing if treatment has a significant effect on size
model.BTM <- lm(size_mm ~ treatment * sex, data = BTM.data)
#ANOVA
anova(model.BTM, test = "F")
#Treatment is significant, P = 3.801e-06
#Sex is significant, P < 2.2e-16
#And the interaction is not, P = 0.6576


#Homogeneity of variances
leveneTest(size_mm ~ treatment * sex, data = BTM.data)
#Not significant, assumption of homogeneity of variances is not met

#Test the model by looking at the residuls
resid.BTM <- residuals(model.BTM)
hist(resid.BTM)
shapiro.test(resid.BTM)
qqnorm(resid.BTM)
qqline(resid.BTM)
#The residuls are more or less normally distrubuted, meaning that the model is valid 


#Use emmeans to test if any of treatments are significant different from each other
emmeans(model.BTM, pairwise ~ treatment * sex)
#All females and males are significant different from each other



#########################################  PLOT DATA  #########################################

### PLOT ALL DATA ###

#Create new dataframe for plotting
FLXBTplot <- BT.data
#Remove NAs
FLXBTplot <- FLXBTplot[which(is.na(FLXBTplot$size_mm) == F),]


#Add fitted values
FLXBTplot$fitted <- fitted(model.BT)

#Caluclate mean
BT.M <- tapply(FLXBTplot$fitted, FLXBTplot$treattype, mean)
#SE
BT.SE <- (summary(model.BT)$coefficients[,2])


## The plot

pdf("BT.pdf", 7, 7)
par(mar = c(3, 5, 2, 2))
#First make the plot
plot(NULL, xlim = c(0, 5), xlab = "", xaxt = "n", ylim = c(0.7, 1.1), ylab = "Thorax size/mm",
     cex.axis = 1.2, cex.lab = 1.5, las = 1)
#Add points
stripchart(BT.data$size_mm ~ BT.data$treattype, add = T, vertical = T, 
           method = "jitter", pch = c(16, 1, 15, 0, 17, 16, 15, 17), cex = 1.4, col = "gray80",
           at = c(0.5, 1, 1.5, 2, 2.5, 3.5, 4, 4.5))
#Add errorbar
xBT <- c(0.5, 1, 1.5, 2, 2.5, 3.5, 4, 4.5)
errbar(xBT, BT.M, BT.M + BT.SE, BT.M - BT.SE, 
       cex = 3, lwd = 3, pch = c(16, 1, 15, 0, 17, 16, 15, 17), add = T)
#Add axis
axis(1, at = c(1.5, 4), cex.axis = 1.5, labels = c("Female", "Male"))
#Mark significant between the treatments with letters
mtext("a", side = 3, line = -2, at = 0.5, cex = 1.5)
mtext("ab", side = 3, line = -2, at = 1, cex = 1.5)
mtext("ab", side = 3, line = -2, at = 1.5, cex = 1.5)
mtext("ab", side = 3, line = -2, at = 2, cex = 1.5)
mtext("b", side = 3, line = -2, at = 2.5, cex = 1.5)
mtext("c", side = 3, line = -2, at = 3.5, cex = 1.5)
mtext("c", side = 3, line = -2, at = 4, cex = 1.5)
mtext("c", side = 3, line = -2, at = 4.5, cex = 1.5)
dev.off()



### PLOT MEAN ###

#Add fitted values
BTM.data$fitted <- fitted(model.BTM)

#Caluclate mean
BTM.M <- aggregate(fitted ~ treatment + sex, mean, data = BTM.data)
#SE
BTM.SE <- (summary(model.BTM)$coefficients[,2])


## The plot

pdf("BTM.pdf", 7, 7)
par(mar = c(3, 5, 2, 2))
#First make the plot
plot(NULL, xlim = c(0.2, 3), xlab = "", xaxt = "n", ylim = c(0.75, 1), ylab = "Thorax size/mm",
     cex.axis = 1.2, cex.lab = 1.5, las = 1)
#Add points
stripchart(BTM.data$size_mm ~ BTM.data$treatsex, add = T, vertical = T, 
           method = "jitter", pch = c(16, 15, 17), cex = 1.4, col = "gray80",
           at = c(0.4, 0.8, 1.2, 2, 2.4, 2.8))
#Add errorbar
xBTM <- c(0.4, 0.8, 1.2, 2, 2.4, 2.8)
errbar(xBTM, BTM.M$fitted, BTM.M$fitted + BTM.SE, BTM.M$fitted - BTM.SE, 
       cex = 3, lwd = 3, pch = c(16, 15, 17), add = T)
#Add axis
axis(1, at = c(0.8, 2.4), cex.axis = 1.5, labels = c("Female", "Male"))
#Add legend
legend("bottomleft", c("FLX", "CFM", "Cwt"), pch = c(16, 15, 17), pt.cex = 2)
#Mark significant between the treatments with letters
mtext("a", side = 3, line = -2, at = 0.4, cex = 1.5)
mtext("a", side = 3, line = -2, at = 0.8, cex = 1.5)
mtext("a", side = 3, line = -2, at = 1.2, cex = 1.5)
mtext("b", side = 3, line = -2, at = 2, cex = 1.5)
mtext("b", side = 3, line = -2, at = 2.4, cex = 1.5)
mtext("b", side = 3, line = -2, at = 2.8, cex = 1.5)
dev.off()

