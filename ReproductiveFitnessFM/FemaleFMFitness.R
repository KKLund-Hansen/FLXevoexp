################################################################################################
################# FM BALANCER FEMALE REPRODUCTIVE FITNESS ASSAY GENERATION 127 #################
################################################################################################

#Set working directory
setwd("~/Google Drev/Work/PhD/Drosophila/FLX/1.FLX/R/FMFitnessAssay")

#Set up environment
library(car)
library(Hmisc)
library(lme4)
library(lmerTest)

#Read in csv file with data
FFMF.data <- read.table(file = "FemaleFMFitness.csv", h = T, sep = ",")
FFMFm.data <- read.table(file = "FemaleFMFitnessMean.csv", h = T, sep = ",")

##########################################  STATISTIC  #########################################


##############################
####THE MODEL FOR ALL DATA####
##############################

### CALCULATE RELATIVE FECUNDITY ###
#Find the maximum
FFMF.max <- max(FFMF.data$offspring, na.rm = T)
#calculate relative fecudity
FFMF.data$relative_fec <- FFMF.data$offspring / FFMF.max


### LINEAT MODEL ###

#Linear mixed model testing if regime type has a significant effect on the relative fitness
model.FFMF <- lmer(relative_fec ~ regimetype + block + (1|regimetype:rep_population), data = FFMF.data)
#ANOVA
anova(model.FFMF, test = "F")
#Regime type is not significant, P = 0.66453
#Block is, P = 0.00142
ranova(model.FFMF)
#But the nested factor is not, P = 0.1646


#Homogeneity of variances
leveneTest(FFMF.data$relative_fec, FFMF.data$regimetype)
#Not significant, assumption of homogeneity of variances is met

#Test the model by looking at the residuls
resid.FFMF <- residuals(model.FFMF)
hist(resid.FFMF)
shapiro.test(resid.FFMF)
qqnorm(resid.FFMF)
qqline(resid.FFMF)
#The residuls are more or less normally distrubuted, meaning that the model is valid 



###################################
####THE MODEL FOR THE MEAN DATA####
###################################

### CALCULATE RELATIVE FECUNDITY ###
#Find the maximum
FFMFm.max <- max(FFMFm.data$offspring)
#calculate relative fecudity
FFMFm.data$relative_fec <- FFMFm.data$offspring / FFMFm.max


### LINEAT MODEL ###

#Linear model testing if regime type has a significant effect on the relative fitness
model.FFMFm <- lm(relative_fec ~ regimetype + block, data = FFMFm.data)
#ANOVA
anova(model.FFMFm, test = "F")
#Treatment type is not significant, P = 0.75177
#Block is, P = 0.03931


#Homogeneity of variances
leveneTest(FFMFm.data$relative_fec, FFMFm.data$regimetype)
#Not significant, assumption of homogeneity of variances is met

#Test the model by looking at the residuls
resid.FFMFm <- residuals(model.FFMFm)
hist(resid.FFMFm)
shapiro.test(resid.FFMFm)
qqnorm(resid.FFMFm)
qqline(resid.FFMFm)
#The residuls are more or less normally distrubuted, meaning that the model is valid 



#########################################  PLOT DATA  #########################################


#Colours
RP1 <- "#ff766f"
RP2 <- "#ffe34d"
RP3 <- "#6da8ff"
RP4 <- "#a5f584"


### PLOT ALL DATA ###

#Create new dataframe for plotting
FFMFplot <- FFMF.data[c(1:2, 6)]
#Remove NAs
FFMFplot <- FFMFplot[which(is.na(FFMF.data$relative_fec) == F),]

#Add fitted values
FFMFplot$fitted <- fitted(model.FFMF)

#Caluclate mean
FFMF.M <- tapply(FFMFplot$fitted, FFMFplot$regimetype, mean)
#SE
FFMF.SE <- (summary(model.FFMF)$coefficients[1:6, 2])


## The plot

pdf("FFMF.pdf", 7, 7)
par(mar = c(3, 5, 2, 2))
#First make the plot
plot(NULL, xlim = c(0, 4.5), xlab = "", xaxt = "n", ylim = c(0.4, 1), ylab = "Relative fitness",
     cex.lab = 1.8, cex.axis = 1.3, las = 1, main = "Females", cex.main = 1.8)
#Add points
points(jitter(rep(0.5, 40), -0.07, 0.07), FFMFplot$relative_fec[FFMFplot$regimetype == "aFLX2X"],
       pch = 16, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(1, 39), -0.07, 0.07), FFMFplot$relative_fec[FFMFplot$regimetype == "bFLXFMX"],
       pch = 1, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(2, 40), -0.07, 0.07), FFMFplot$relative_fec[FFMFplot$regimetype == "cCFM2X"],
       pch = 15, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(2.5, 40), -0.07, 0.07), FFMFplot$relative_fec[FFMFplot$regimetype == "dCFMFMX"],
       pch = 0, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(3.5, 40), -0.07, 0.07), FFMFplot$relative_fec[FFMFplot$regimetype == "eCwt"],
       pch = 17, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(4, 40), -0.07, 0.07), FFMFplot$relative_fec[FFMFplot$regimetype == "fCwtFMX"],
       pch = 2, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
#Add errorbar
xFFMF <- c(0.5, 1, 2, 2.5, 3.5, 4)
errbar(xFFMF, FFMF.M, FFMF.M + FFMF.SE, FFMF.M - FFMF.SE,
       pch = c(16, 1, 15, 0, 17, 2), cex = 2.5, lwd = 2.5, add = T)
#Add axis
axis(1, at = c(0.75, 2.25, 3.75), cex.axis = 1.5, labels = c("FLX", "Control FM", "Control wt"))
#Add legend
legend("bottomright", c("Closed symbols: X/X", "Open symbols: FM/X"), adj = c(0.1, 0.5), cex = 1.1)
dev.off()



### PLOT MEAN ###

#Add fitted values
FFMFm.data$fitted <- fitted(model.FFMFm)

#Caluclate mean
FFMFm.M <- tapply(FFMFm.data$fitted, FFMFm.data$regimetype, mean)
#SE
FFMFm.SE <- (summary(model.FFMFm)$coefficients[1:6, 2])


## The plot

pdf("FFMFm.pdf", 7, 7)
par(mar = c(3, 5, 2, 2))
#First make the plot
plot(NULL, xlim = c(0, 4.5), xlab = "", xaxt = "n", ylim = c(0.6, 1), ylab = "Relative fitness",
     cex.lab = 1.8, cex.axis = 1.3, las = 1, main = "Females", cex.main = 1.8)
#Add points
points(jitter(rep(0.5, 8), -0.07, 0.07), FFMFm.data$relative_fec[FFMFm.data$regimetype == "aFLX2X"],
       pch = 16, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(1, 8), -0.07, 0.07), FFMFm.data$relative_fec[FFMFm.data$regimetype == "bFLXFMX"],
       pch = 1, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(2, 8), -0.07, 0.07), FFMFm.data$relative_fec[FFMFm.data$regimetype == "cCFM2X"],
       pch = 15, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(2.5, 8), -0.07, 0.07), FFMFm.data$relative_fec[FFMFm.data$regimetype == "dCFMFMX"],
       pch = 0, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(3.5, 8), -0.07, 0.07), FFMFm.data$relative_fec[FFMFm.data$regimetype == "eCwt"],
       pch = 17, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
points(jitter(rep(4, 8), -0.07, 0.07), FFMFm.data$relative_fec[FFMFm.data$regimetype == "fCwtFMX"],
       pch = 2, cex = 1.4, col = c(RP1, RP2, RP3, RP4))
#Add errorbar
xFFMFm <- c(0.5, 1, 2, 2.5, 3.5, 4)
errbar(xFFMFm, FFMFm.M, FFMFm.M + FFMFm.SE, FFMFm.M - FFMFm.SE,
       pch = c(16, 1, 15, 0, 17, 2), cex = 2.5, lwd = 2.5, add = T)
#Add axis
axis(1, at = c(0.75, 2.25, 3.75), cex.axis = 1.5, labels = c("FLX", "Control FM", "Control wt"))
#Add legend
legend("bottomleft", c("Closed symbols: X/X", "Open symbols: FM/X"), adj = c(0.1, 0.5), cex = 1.1)
dev.off()

