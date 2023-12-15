#SICKO Statistical Testing

#install.packages("purrr")
#install.packages("coin")
#install.packages("dplyr")
#install.packages("tidyr")
#install.packages("ggpubr")
#install.packages("tidyverse")
#install.packages("vioplot")
#install.packages("sm")
#install.packages("zoo")

#Load Libraries
library(lattice)
library(ggplot2)
library(survival)
library(dplyr)
library(tidyr)
library(purrr)
library(tidyverse)
library(broom)
library(vioplot)
library(coin)
# library(ade4)
# library(RcppEigen)
# library(lme4)
# library(car)
# library(FactoMineR)
# library(mixOmics)
# library(RVAideMemoire)
# library(rstatix)
# library(ggpubr)
# library(stats)
# library(rlang)


# ****************** READ INSTRUCTIONS BELOW FIRST***************************************************

# Press "Ctrl" + "F" to begin changing IDs and Labels for experiment
  # Search "Control Sample ID" and replace all with ID for control conditions title used in experimental folder saving (ex: N2)
  # Search "Control Label" and replace all with label desired for graphical outputs (ex: WT)
  # Search "Test Sample ID" and replace all with ID for test condition used titles used in experimental folder saving (ex: KU25 [strain name])
  # Search "Test Label" and replace all with label desired for graphical outputs (ex: pmk-1 [gene KO])
# Change working directory (line 51) to "outputs" folder created in the image processing and compiling
# Set "Data.set"(line 53) to csv denoted "_compiled_analyzed.csv" in "outputs" folder to feed data from SICKO processing
# Fill in numbers of worms alive and died during the washing process of both the Control and Test conditions (lines: 56, 57, 61, 62)

# ***************************************************************************************************

# set working directory
setwd("********_*****************_****************_********************_***************_outputs")

Data.set <- read.csv("******************************_compiled_analyzed.csv")

#Input numbers to account for worms died in process (Numbers in SICKO Coefficient Values spread sheet)
Control.Living.Population <- 464 #How many living worms survived wash plates (selection population)
Control.Died.Transfer <- 19 #How many worms died during total washes (missing from selection population)
Control.Total.Pop <- Control.Living.Population + Control.Died.Transfer #Effective population (includes population lost before observation)
Control.Died.Proportion <- Control.Died.Transfer/Control.Total.Pop
Control.Observed.Proportion <- 1 - Control.Died.Proportion
Test.Living.Population <- 387
Test.Died.Transfer <- 60
Test.Total.Pop <- Test.Living.Population + Test.Died.Transfer
Test.Died.Proportion <- Test.Died.Transfer/Test.Total.Pop
Test.Observed.Proportion <- 1 - Test.Died.Proportion

#Make labels columns for future figure making in proper notation
Data.set$Labels <- NA

#Fill labels column based off conditions
for (i in 1:length(Data.set$Condition)){
  if(Data.set$Condition[i] == "Control Sample ID"){#________Enter Data______Control experiment
    Data.set$Labels[i] <- "Control Label"#________Enter Data______Control label
  } else if(Data.set$Condition[i] != "Control Sample ID"){#________Enter Data______Control experiment
    Data.set$Labels[i] <- "Test Label"#________Enter Data______Test
  }
}

#Statistical Output for SICKO Tests Run (Can find line in script by searching SICKO.Stat)
SICKO.Statistics <- data.frame(matrix(NA,
                               ncol = 8,
                               nrow = 35))
colnames(SICKO.Statistics) <- c("Comparison", "Statistical Test","p-value", "Statistic Type",  "Statistic", "Parameter Type", "Parameter", "Associated Figure")
#Add means to table
SICKO.Statistics$Mean_1 <- NA
SICKO.Statistics$Mean_2 <- NA
SICKO.Statistics$SD <- NA
SICKO.Statistics$Median_1 <- NA
SICKO.Statistics$Median_2 <- NA


#Make Day of Death column from last day of observation +1
Data.set$Day.of.Death <- NA

#Fill column
for (i in 1:length(Data.set$Last.Day.of.Observation)){
  if(Data.set$Is.Dead[i] == 1){
    Data.set$Day.of.Death[i] <- Data.set$Last.Day.of.Observation[i] + 1
  } else {Data.set$Day.of.Death[i] <- NaN
  }
}

Conditions <- unique(Data.set$Condition)
Conditions <- factor(Conditions,levels = c("Control Sample ID","Test Sample ID")) #________Enter Data______Control then Test experiement

#Create and factor Infection status column and conditions___________________________________________________________________________________________________________________
Is.Infected <- Data.set#Need to fix to new column with infection status
Is.Infected$Infection.Status <- NA#Add new column to data frame
for (i in 1:length(Is.Infected$Infection.Status)){ #Make data from "first day of non-zero data" into infection status
if(is.nan(Is.Infected$First.Day.of.nonzero.data[i]) ){
  Is.Infected$Infection.Status[i] <- "Not Infected"
  } else {
    Is.Infected$Infection.Status[i] <- "Infected"
  }
}

Infection.State <- unique(Is.Infected$Infection.Status) #Define infection conditions
Infection.State <- factor(Is.Infected$Infection.Status, levels = c("Not Infected", "Infected"))

Is.Infected$Condition <- factor(Is.Infected$Condition, levels = c("Control Sample ID","Test Sample ID")) #Turn conditions into factor#________Enter Data______Control then Test experiment
Is.Infected$Infection.Status <- factor(Is.Infected$Infection.Status, levels = c("Not Infected", "Infected"))
Is.Infected$Labels <- factor(Is.Infected$Labels, levels = c("Control Label","Test Label"))

#Create Infection grouping by Condition and Infection state and create factoring
Is.Infected$Infection.by.Group <- NA
for (i in 1:length(Is.Infected$Infection.Status)) {
  Is.Infected$Infection.by.Group[i] <- paste(Is.Infected$Labels[i], Is.Infected$Infection.Status[i])
}

#Factor groups and infection for final graphs to make comparisons between groups and infection state
Is.Infected$Infection.by.Group <- factor(Is.Infected$Infection.by.Group, levels = c("Control Label Not Infected", "Control Label Infected", "Test Label Not Infected", "Test Label Infected"))#________Input Data______

#Create individual data frames for Control condition and test condition (will be dependent on experiment)
Control.DF <- Is.Infected[Is.Infected$Labels == "Control Label",] #________Input Data______
Test.DF <- Is.Infected[Is.Infected$Labels != "Control Label",] #________Input Data______

#Chi Square Test if Infection is Related to Alive or Dead Status in Control Animals_________________________________________________________________________________________

#Number of infected animals vs Total animals
Control.Not.Infected.DF <- Control.DF[Control.DF$Infection.Status=="Not Infected",]
Control.Not.Infected.Num <- nrow(Control.Not.Infected.DF)
Control.Not.Infected.Alive.Num <- length(which(Control.Not.Infected.DF$Is.Dead==0))
Control.Not.Infected.Dead.Num <- length(which(Control.Not.Infected.DF$Is.Dead==1))

Control.Infected.DF <- Control.DF[Control.DF$Infection.Status=="Infected",]
Control.Infected.Num <- nrow(Control.Infected.DF)
Control.Infected.Alive.Num <-length(which(Control.Infected.DF$Is.Dead==0))
Control.Infected.Dead.Num <-length(which(Control.Infected.DF$Is.Dead==1))


#Make Table for Chi Square Test
DF.Infected.Dead.Control <- data.frame(Alive = NA, Dead = NA)
DF.Infected.Dead.Control[1,1] <- Control.Not.Infected.Alive.Num
DF.Infected.Dead.Control[1,2] <- Control.Not.Infected.Dead.Num
DF.Infected.Dead.Control[2,1] <- Control.Infected.Alive.Num
DF.Infected.Dead.Control[2,2] <- Control.Infected.Dead.Num
rownames(DF.Infected.Dead.Control) <- c("Uninfected","Infected")

DF.Infected.Dead.Control

#Chi square Infection v Death test in Control
Chisq.Control.Infection.V.Death <-  chisq.test(DF.Infected.Dead.Control, correct = FALSE)

#SICKO.Stat 1
SICKO.Statistics[1,1] <- "Control Condition Infection v Death"
SICKO.Statistics[1,2] <- Chisq.Control.Infection.V.Death$method
SICKO.Statistics[1,3] <- Chisq.Control.Infection.V.Death$p.value
SICKO.Statistics[1,4] <- "t"
SICKO.Statistics[1,5] <- Chisq.Control.Infection.V.Death$statistic
SICKO.Statistics[1,6] <- "df"
SICKO.Statistics[1,7] <- Chisq.Control.Infection.V.Death$parameter
SICKO.Statistics[1,8] <- "Heat Map"
 
#Chi Square Test if Infection is Related to Alive or Dead Status in Test Animals_________________________________Still in work: Really shows survival past median because Test are far past________________________________________________________
  
#Number of infected animals vs Total animals
Test.Not.Infected.DF <- Test.DF[Test.DF$Infection.Status=="Not Infected",]
Test.Not.Infected.Num <- nrow(Test.Not.Infected.DF)
Test.Not.Infected.Alive.Num <- length(which(Test.Not.Infected.DF$Is.Dead==0))
Test.Not.Infected.Dead.Num <- length(which(Test.Not.Infected.DF$Is.Dead==1))
  
Test.Infected.DF <- Test.DF[Test.DF$Infection.Status=="Infected",]
Test.Infected.Num <- nrow(Test.Infected.DF)
Test.Infected.Alive.Num <-length(which(Test.Infected.DF$Is.Dead==0))
Test.Infected.Dead.Num <-length(which(Test.Infected.DF$Is.Dead==1))
  
  
#Make Table for Chi Square Test
DF.Infected.Dead.Test <- data.frame(Alive = NA, Dead = NA)
DF.Infected.Dead.Test[1,1] <- Test.Not.Infected.Alive.Num
DF.Infected.Dead.Test[1,2] <- Test.Not.Infected.Dead.Num
DF.Infected.Dead.Test[2,1] <- Test.Infected.Alive.Num
DF.Infected.Dead.Test[2,2] <- Test.Infected.Dead.Num
rownames(DF.Infected.Dead.Test) <- c("Uninfected","Infected")

DF.Infected.Dead.Test
  
#Chi square Infection v Death test
Chisq.Test.Infection.V.Death <- chisq.test(DF.Infected.Dead.Test, correct = FALSE)

#SICKO.Stat 2
SICKO.Statistics[2,1] <- "Test Condition Infection v Death"
SICKO.Statistics[2,2] <- Chisq.Test.Infection.V.Death$method
SICKO.Statistics[2,3] <- Chisq.Test.Infection.V.Death$p.value
SICKO.Statistics[2,4] <- "t"
SICKO.Statistics[2,5] <- Chisq.Test.Infection.V.Death$statistic
SICKO.Statistics[2,6] <- "df"
SICKO.Statistics[2,7] <- Chisq.Test.Infection.V.Death$parameter
SICKO.Statistics[2,8] <- "Heat Map"


#Obtain Values to Chi test Infection vs Condition_______________________________________________________________________________________________________
Control.Total.Num <- nrow(Control.DF)
Test.Total.Num <- nrow(Test.DF)

Test.Not.Infected.DF <- Test.DF[Test.DF$First.Day.of.nonzero.data==0,]
Test.Not.Infected.Num <- nrow(Test.Not.Infected.DF)
Test.Infected.Num <- Test.Total.Num - Test.Not.Infected.Num

#DF for Chi Test of Infection by Condition
DF.Infected.Condition <- data.frame(Not.Infected = NA, Infected = NA)
DF.Infected.Condition[1,1] <- Control.Not.Infected.Num
DF.Infected.Condition[1,2] <- Control.Infected.Num
DF.Infected.Condition[2,1] <- Test.Not.Infected.Num
DF.Infected.Condition[2,2] <- Test.Infected.Num

rownames(DF.Infected.Condition) <- c("Control Label","Test Label") #________Input Data______

DF.Infected.Condition

#Chi square test Infection v Condition
chisq.Test.Stat <- chisq.test(DF.Infected.Condition, correct = FALSE) 

#SICKO.Stat 3
SICKO.Statistics[3,1] <- "Condition v Infection"
SICKO.Statistics[3,2] <- chisq.Test.Stat$method
SICKO.Statistics[3,3] <- chisq.Test.Stat$p.value
SICKO.Statistics[3,4] <- "t"
SICKO.Statistics[3,5] <- chisq.Test.Stat$statistic
SICKO.Statistics[3,6] <- "df"
SICKO.Statistics[3,7] <- chisq.Test.Stat$parameter
SICKO.Statistics[3,8] <- "Heat Map"

#Kaplan Meyer Infection vs Dead Control and Test (Only Dead Animals)_____________________________Not the best meteric due to survivor bias________________________________________________________________________________

Is.Dead.Control <-Control.DF[Control.DF$Is.Dead == 1,]

survfit.Life.Infection.Control <- survfit(Surv(Last.Day.of.Observation)~Infection.Status, data=Is.Dead.Control)

plot(survfit.Life.Infection.Control, xlab = "Time(days)",
     ylab = "Fraction Surviving",
     col = c("black","red","blue","forestgreen"),
     conf.int = F,
     mark.time = T,
     main = "Control Label Survival Post Challenge") #_________________________________Input Data_________________________

legend("bottomleft",
       lty = 1,
       col = c("black", "red", "blue", "forestgreen"),
       legend = levels(Is.Dead.Control$Labels))


survdiff(Surv(Last.Day.of.Observation)~Infection.Status, data = Is.Dead.Control)


#Test Condition
Is.Dead.Test <-Test.DF[Test.DF$Is.Dead == 1,]

survfit.Life.Infection.Test <- survfit(Surv(Last.Day.of.Observation)~Infection.Status, data=Is.Dead.Test)

plot(survfit.Life.Infection.Test, xlab = "Time(days)",
     ylab = "Fraction Surviving",
     col = c("black","red","blue","forestgreen"),
     conf.int = F,
     mark.time = T,
     main = "Test Label -/- Survival Post Challenge") #_________________________________Input Data_________________________

legend("bottomleft",
       lty = 1,
       col = c("black", "red", "blue", "forestgreen"),
       legend = levels(Infection.State))


survdiff(Surv(Last.Day.of.Observation)~First.Day.of.nonzero.data, data = Is.Dead.Test)

#Perform Kaplan meyer and log rank test on just conditions #Come back #Need to calculate proportion lost before sampling, could start at total proportion lost then kaplan does proportion of those that esisted

survfit.Life.Conditions <- survfit(Surv(Last.Day.of.Observation, Is.Dead == 1)~Labels, data=Is.Infected)

plot(survfit.Life.Conditions)

#Kaplan Meyer Infection vs Alive Control and Test (All Animals Alive and Dead)_____________________________________________________________________________________________________________

#Adujust for population lost prior to observation

#Adujust for population lost prior to observation
Control.Sample.Size <- nrow(Control.DF) #Population observed in SICKO
Control.Pop.Adjust.Lifespan <- Control.Sample.Size/Control.Observed.Proportion #cross multiply to include proportion lost into observered proportion
Control.Round.Pop <- round(Control.Pop.Adjust.Lifespan) # Round number to make rows to fill
Control.Missing.Rows.Adjustment <- Control.Round.Pop - Control.Sample.Size #Row amount that must be added to account for those that die in SICKO process
Control.Num.Died.Observed <- Control.Infected.Dead.Num + Control.Not.Infected.Dead.Num # Number of animals that died in observation so estimate infect v not infected death prior to observation
Control.Prop.Died.Infected <- Control.Infected.Dead.Num/Control.Num.Died.Observed #Proportion of those that were observed that died and were infected
Control.Rows.Died.Missing.Infected <- Control.Missing.Rows.Adjustment*Control.Prop.Died.Infected #apply proportion of death to infection to population lost before observation
Control.Rows.Died.Missing.Infected.Rounded <- round(Control.Rows.Died.Missing.Infected) #round number to create proper row numbers
Control.Rows.Died.Missing.Not.Infected.Rounded <- Control.Missing.Rows.Adjustment - Control.Rows.Died.Missing.Infected.Rounded
Test.Sample.Size <- nrow(Test.DF) 
Test.Pop.Adjust.Lifespan <- Test.Sample.Size/Test.Observed.Proportion 
Test.Round.Pop <- round(Test.Pop.Adjust.Lifespan) 
Test.Missing.Rows.Adjustment <- Test.Round.Pop - Test.Sample.Size 
Test.Num.Died.Observed <- Test.Infected.Dead.Num + Test.Not.Infected.Dead.Num 
Test.Prop.Died.Infected <- Test.Infected.Dead.Num/Test.Num.Died.Observed
Test.Rows.Died.Missing.Infected <- Test.Missing.Rows.Adjustment*Test.Prop.Died.Infected
Test.Rows.Died.Missing.Infected.Rounded <- round(Test.Rows.Died.Missing.Infected)
Test.Rows.Died.Missing.Not.Infected.Rounded <- Test.Missing.Rows.Adjustment - Test.Rows.Died.Missing.Infected.Rounded

#Create dataframe to append observed data for infected animals that died in the process
Control.Lifespan.DF.Supp.Infected <- data.frame(matrix(NA,
                                                       nrow = Control.Rows.Died.Missing.Infected.Rounded,
                                                       ncol = ncol(Control.DF),
))
colnames(Control.Lifespan.DF.Supp.Infected) = colnames(Control.DF[c(1:37)])

Control.Lifespan.DF.Supp.Not.Infected <- data.frame(matrix(NA,
                                                           nrow = Control.Rows.Died.Missing.Not.Infected.Rounded,
                                                           ncol = ncol(Control.DF),
))
colnames(Control.Lifespan.DF.Supp.Not.Infected) = colnames(Control.DF[c(1:37)])

#Now with Test Condition
Test.Lifespan.DF.Supp.Infected <- data.frame(matrix(NA,
                                                    nrow = Test.Rows.Died.Missing.Infected.Rounded,
                                                    ncol = ncol(Test.DF),
))
colnames(Test.Lifespan.DF.Supp.Infected) = colnames(Test.DF[c(1:37)])

Test.Lifespan.DF.Supp.Not.Infected <- data.frame(matrix(NA,
                                                        nrow = Test.Rows.Died.Missing.Not.Infected.Rounded,
                                                        ncol = ncol(Test.DF),
))
colnames(Test.Lifespan.DF.Supp.Not.Infected) = colnames(Test.DF[c(1:37)])

#Fill ammending df with data needed for lifespan analysis______________Enter Data_______
if (Control.Rows.Died.Missing.Not.Infected.Rounded > 0) {
Control.Lifespan.DF.Supp.Not.Infected$Is.Dead <- 1
Control.Lifespan.DF.Supp.Not.Infected$Is.Last.Day.Censored <- 0
Control.Lifespan.DF.Supp.Not.Infected$Last.Day.of.Observation <- 0
Control.Lifespan.DF.Supp.Not.Infected$Condition <- "Control Sample ID"
Control.Lifespan.DF.Supp.Not.Infected$Labels <- "Control Label"
Control.Lifespan.DF.Supp.Not.Infected$Infection.Status <- "Not Infected"
Control.Lifespan.DF.Supp.Not.Infected$Infection.by.Group <- "Control Label Not Infected"
}

if (Control.Rows.Died.Missing.Infected.Rounded > 0) {
Control.Lifespan.DF.Supp.Infected$Is.Dead <- 1
Control.Lifespan.DF.Supp.Infected$Is.Last.Day.Censored <- 0
Control.Lifespan.DF.Supp.Infected$Last.Day.of.Observation <- 0
Control.Lifespan.DF.Supp.Infected$Condition <- "Control Sample ID"
Control.Lifespan.DF.Supp.Infected$Labels <- "Control Label"
Control.Lifespan.DF.Supp.Infected$Infection.Status <- "Infected"
Control.Lifespan.DF.Supp.Infected$Infection.by.Group <- "Control Label Infected"
}

#Fill Test amending df with data needed for lifespan analysis
if (Test.Rows.Died.Missing.Not.Infected.Rounded > 0) {
Test.Lifespan.DF.Supp.Not.Infected$Is.Dead <- 1
Test.Lifespan.DF.Supp.Not.Infected$Is.Last.Day.Censored <- 0
Test.Lifespan.DF.Supp.Not.Infected$Last.Day.of.Observation <- 0
Test.Lifespan.DF.Supp.Not.Infected$Condition <- "Test Sample ID"
Test.Lifespan.DF.Supp.Not.Infected$Labels <- "Test Label"
Test.Lifespan.DF.Supp.Not.Infected$Infection.Status <- "Not Infected"
Test.Lifespan.DF.Supp.Not.Infected$Infection.by.Group <- "Test Label Not Infected"
}

if (Test.Rows.Died.Missing.Infected.Rounded > 0) {
Test.Lifespan.DF.Supp.Infected$Is.Dead <- 1
Test.Lifespan.DF.Supp.Infected$Is.Last.Day.Censored <- 0
Test.Lifespan.DF.Supp.Infected$Last.Day.of.Observation <- 0
Test.Lifespan.DF.Supp.Infected$Condition <- "Test Sample ID"
Test.Lifespan.DF.Supp.Infected$Labels <- "Test Label"
Test.Lifespan.DF.Supp.Infected$Infection.Status <- "Infected"
Test.Lifespan.DF.Supp.Infected$Infection.by.Group <- "Test Label Infected"
}

Control.DF.Lifespan <- rbind(Control.DF,Control.Lifespan.DF.Supp.Not.Infected,Control.Lifespan.DF.Supp.Infected) #Amend observed data to account for those lost in infection and wash process
Test.DF.Lifespan <- rbind(Test.DF,Test.Lifespan.DF.Supp.Not.Infected,Test.Lifespan.DF.Supp.Infected)


#Export Figures
pdf(file = "SICKO_Figures.pdf", width = 5, height = 5)

par(
  bty = "l", #box around the graph
  mar = c(5,5,5,3), #Margins around the graph(standard is (5,4,4,2))
  tck = -0.03, #tick marks on axis (standard is ~ -0.05)
  lwd = 3, #line with
  las = 1
)

#Control Kaplan Meyer Infected vs Non-Infected (Figure 1)

survfit.Life.Infection.Control <- survfit(Surv(Last.Day.of.Observation, Is.Dead == 1 )~Infection.Status, data=Control.DF.Lifespan)

plot(survfit.Life.Infection.Control, 
     xlab = "Time Post Challenge (days)",
     ylab = "Fraction Surviving",
     col = c("black","red","blue","forestgreen"),
     lwd = 3,
     conf.int = F,
     mark.time = T,
     cex.lab = 1.2,
     cex.main = 1.2,
     xaxt='n',
     main = "Control Label Infection and Survival Post Challenge") #_________________________________Input Data_________________________

axis(side = 1, at=seq(0,9, by = 1))

legend("bottomleft",
       lty = 1,
       col = c("black", "red", "blue", "forestgreen"),
       legend = levels(Infection.State))

Kaplan.Control.Stat <- survdiff(Surv(Last.Day.of.Observation,Is.Dead == 1 )~Infection.Status, data = Control.DF.Lifespan)

Control.Lifespan.Summary <- summary(survfit.Life.Infection.Control) #pull our medial lifespan with this then use ###$table to access

#SICKO.Stat 4
SICKO.Statistics[4,1] <- "Control Infected V Non-Infected"
SICKO.Statistics[4,2] <- "Log-Rank Test"
SICKO.Statistics[4,3] <- Kaplan.Control.Stat$pvalue
SICKO.Statistics[4,4] <- "Chi-squared"
SICKO.Statistics[4,5] <- Kaplan.Control.Stat$chisq
SICKO.Statistics[4,6] <- "df"
SICKO.Statistics[4,7] <- 1
SICKO.Statistics[4,8] <- "Pg 1"
SICKO.Statistics[4,9] <- Control.Lifespan.Summary$table[1,5]
SICKO.Statistics[4,10] <- Control.Lifespan.Summary$table[2,5]
SICKO.Statistics[4,12] <- Control.Lifespan.Summary$table[1,7]
SICKO.Statistics[4,13] <- Control.Lifespan.Summary$table[2,7]

#Test Survfit
survfit.Life.Infection.Test <- survfit(Surv(Last.Day.of.Observation,Is.Dead == 1)~Infection.Status, data=Test.DF.Lifespan)

#Test Kaplan Meyer Infected vs Non-Infected (Figure 2)
plot(survfit.Life.Infection.Test, 
     xlab = "Time Post Challenge (days)",
     ylab = "Fraction Surviving",
     col = c("black","red","blue","forestgreen"),
     lwd = 3,
     conf.int = F,
     mark.time = T,
     cex.lab = 1.2,
     cex.main = 1.2,
     font.main = 2,
     xaxt='n',
     main = substitute(paste(bolditalic('Test Label '), bold(' Infection and Survival Post Challenge')))) #_________________________________Input Data_________________________

axis(side = 1, at=seq(0,9, by = 1))

legend("bottomleft",
       lty = 1,
       col = c("black", "red", "blue", "forestgreen"),
       legend = levels(Infection.State))


Kaplan.Test.Stat <- survdiff(Surv(Last.Day.of.Observation, Is.Dead == 1)~Infection.Status, data = Test.DF.Lifespan)

Test.Lifespan.Summary <- summary(survfit.Life.Infection.Test)

#SICKO.Stat 5
SICKO.Statistics[5,1] <- "Test Infected V Non-Infected"
SICKO.Statistics[5,2] <- "Log-Rank Test"
SICKO.Statistics[5,3] <- Kaplan.Test.Stat$pvalue
SICKO.Statistics[5,4] <- "Chi-squared"
SICKO.Statistics[5,5] <- Kaplan.Test.Stat$chisq
SICKO.Statistics[5,6] <- "df"
SICKO.Statistics[5,7] <- 1
SICKO.Statistics[5,8] <- "Pg 2"
SICKO.Statistics[5,9] <- Test.Lifespan.Summary$table[1,5]
SICKO.Statistics[5,10] <- Test.Lifespan.Summary$table[2,5]
SICKO.Statistics[5,12] <- Test.Lifespan.Summary$table[1,7]
SICKO.Statistics[5,13] <- Test.Lifespan.Summary$table[2,7]


#Now compare infected v non-infected in both conditions control and test
All.DF.Lifespan <- rbind(Control.DF.Lifespan,Test.DF.Lifespan)

survfit.Life.Infection.Compiled <- survfit(Surv(Last.Day.of.Observation, Is.Dead == 1)~Infection.by.Group, data=All.DF.Lifespan)

#Plot Control vs Test in both Infected and non-infected (Figure 3)
plot(survfit.Life.Infection.Compiled, 
     xlab = "Time Post Challenge (days)",
     ylab = "Fraction Surviving",
     col = c("black","red","blue","forestgreen"),
     lwd = 3, #line thickness
     conf.int = F,
     mark.time = T,
     cex.lab = 1.2, #axis label size
     cex.main = 1.2, #Title size
     xaxt='n',
     main = substitute(paste(bolditalic('Test Label'),bold(' v Control Label Infection and Survival Post Challenge')))) #_________________________________Input Data_________________________

axis(side = 1, at=seq(0,9, by = 1))

legend("bottomleft",
       cex = 0.9,
       lty = 1,
       col = c("black", "red", "blue", "forestgreen"),
       legend = levels(All.DF.Lifespan$Infection.by.Group))


Infection.Conditions.Stats <- survdiff(Surv(Last.Day.of.Observation,Is.Dead == 1 )~Infection.by.Group, data = All.DF.Lifespan)

#SICKO.Stat 6
SICKO.Statistics[6,1] <- "Control v Test; Infected V Non-Infected"
SICKO.Statistics[6,2] <- "Log-Rank Test"
SICKO.Statistics[6,3] <- Infection.Conditions.Stats$pvalue
SICKO.Statistics[6,4] <- "Chi-squared"
SICKO.Statistics[6,5] <- Infection.Conditions.Stats$chisq
SICKO.Statistics[6,6] <- "df"
SICKO.Statistics[6,7] <- 1
SICKO.Statistics[6,8] <- "Pg 3"

#Create Dataframe for survival comparison of infection/non between text and control

Infected.Worms.Lifespan <- All.DF.Lifespan[All.DF.Lifespan$Infection.Status=="Infected",]
Not.Infected.Worms.Lifespan <- All.DF.Lifespan[All.DF.Lifespan$Infection.Status=="Not Infected",]

#Graph and compare Infected animals in both conditions (Figure 4)
survfit.Life.All.Infected <- survfit(Surv(Last.Day.of.Observation, Is.Dead == 1)~Infection.by.Group, data=Infected.Worms.Lifespan)

plot(survfit.Life.All.Infected, 
     xlab = "Time Post Challenge (days)",
     ylab = "Fraction Surviving",
     col = c("black","red","blue","forestgreen"),
     lwd = 3, #line thickness
     conf.int = F,
     mark.time = T,
     cex.lab = 1.2, #axis label size
     cex.main = 1.2, #Title size
     xaxt='n',
     main = substitute(paste(bolditalic('Test Label'),bold(' v Control Label Infected Survival Post Challenge')))) #_________________________________Input Data_________________________

axis(side = 1, at=seq(0,9, by = 1))

legend("bottomleft",
       cex = 0.9,
       lty = 1,
       col = c("black", "red", "blue", "forestgreen"),
       legend = levels(Infected.Worms.Lifespan$Labels))


Infected.Animals.Conditions.Stats <- survdiff(Surv(Last.Day.of.Observation,Is.Dead == 1 )~Infection.by.Group, data = Infected.Worms.Lifespan)

Infected.Lifespan.Summary <- summary(survfit.Life.All.Infected)

#SICKO.Stat 7
SICKO.Statistics[7,1] <- "Infected Control v Test"
SICKO.Statistics[7,2] <- "Log-Rank Test"
SICKO.Statistics[7,3] <- Infected.Animals.Conditions.Stats$pvalue
SICKO.Statistics[7,4] <- "Chi-squared"
SICKO.Statistics[7,5] <- Infected.Animals.Conditions.Stats$chisq
SICKO.Statistics[7,6] <- "df"
SICKO.Statistics[7,7] <- 1
SICKO.Statistics[7,8] <- "Pg 4"
SICKO.Statistics[7,9] <- Infected.Lifespan.Summary$table[1,5]
SICKO.Statistics[7,10] <- Infected.Lifespan.Summary$table[2,5]
SICKO.Statistics[7,12] <- Infected.Lifespan.Summary$table[1,7]
SICKO.Statistics[7,13] <- Infected.Lifespan.Summary$table[2,7]

#Compare non-infected animals in both conditions (Figure 5)
survfit.Life.All.Not.Infected <- survfit(Surv(Last.Day.of.Observation, Is.Dead == 1)~Infection.by.Group, data=Not.Infected.Worms.Lifespan)

plot(survfit.Life.All.Not.Infected, 
     xlab = "Time Post Challenge (days)",
     ylab = "Fraction Surviving",
     col = c("black","red","blue","forestgreen"),
     lwd = 3, #line thickness
     conf.int = F,
     mark.time = T,
     cex.lab = 1.2, #axis label size
     cex.main = 1.2, #Title size
     xaxt='n',
     main = substitute(paste(bolditalic('Test Label'),bold(' v Control Label Not Infected Survival Post Challenge')))) #_________________________________Input Data_________________________

axis(side = 1, at=seq(0,9, by = 1))

legend("bottomleft",
       cex = 0.9,
       lty = 1,
       col = c("black", "red", "blue", "forestgreen"),
       legend = levels(Not.Infected.Worms.Lifespan$Labels))


Not.Infected.Animals.Conditions.Stats <- survdiff(Surv(Last.Day.of.Observation,Is.Dead == 1 )~Infection.by.Group, data = Not.Infected.Worms.Lifespan)

Not.Infected.Lifespan.Summary <- summary(survfit.Life.All.Not.Infected)

#SICKO.Stat 8
SICKO.Statistics[8,1] <- "Non-Infected Control v Test"
SICKO.Statistics[8,2] <- "Log-Rank Test"
SICKO.Statistics[8,3] <- Not.Infected.Animals.Conditions.Stats$pvalue
SICKO.Statistics[8,4] <- "Chi-squared"
SICKO.Statistics[8,5] <- Not.Infected.Animals.Conditions.Stats$chisq
SICKO.Statistics[8,6] <- "df"
SICKO.Statistics[8,7] <- 1
SICKO.Statistics[8,8] <- "Pg 5"
SICKO.Statistics[8,9] <- Not.Infected.Lifespan.Summary$table[1,5]
SICKO.Statistics[8,10] <- Not.Infected.Lifespan.Summary$table[2,5]
SICKO.Statistics[8,12] <- Not.Infected.Lifespan.Summary$table[1,7]
SICKO.Statistics[8,13] <- Not.Infected.Lifespan.Summary$table[2,7]

#Kaplan Meyer Tracking Onset of Infection Between Conditions (All Animals Alive and Dead)_____________________________________________________________________________________________________________

#Create Column to capture if worm infected and alive to properly create kaplan-meyer
Is.Infected$Infection.Start <- NA #create column in dataset

#Fill column with proper data
for (i in 1:length(Is.Infected$Infection.Status)){
  if(Is.Infected$Infection.Status[i] == "Infected") {
    Is.Infected$Infection.Start[i] <- Is.Infected$First.Day.of.nonzero.data[i]
  } else if(Is.Infected$Infection.Status[i] == "Not Infected" 
            && Is.Infected$Is.Dead[i] == 0) {
    Is.Infected$Infection.Start[i] <- 9
  } else if(Is.Infected$Infection.Status[i] == "Not Infected"
            && Is.Infected$Is.Dead[i] == 1) {
    Is.Infected$Infection.Start[i] <- Is.Infected$Last.Day.of.Observation[i] #### Come Back when not tired_______This doesn't look right
  }
}

survfit.Infection.Onset <- survfit(Surv(Infection.Start, Infection.Status == "Infected")~Condition, data=Is.Infected)

#Plot Latent Infections In Control Label and Test Conditions (Figure 6)
plot(survfit.Infection.Onset, xlab = "Time(days)",
     ylab = "Fraction Not Infected",
     col = c("black","red","blue","forestgreen"),
     lwd = 3,
     conf.int = F,
     mark.time = T,
     cex.lab = 1.2,
     cex.main = 1.2,
     xaxt = 'n',
     main = "Onset of Infection")

axis(side = 1, at=seq(0,9, by = 1))

legend("bottomleft",
       lty = 1,
       col = c("black", "red", "blue", "forestgreen"),
       legend = levels(Is.Infected$Labels))


Infection.Onset <- survdiff(Surv(Infection.Start, Infection.Status == "Infected")~Condition, data=Is.Infected)

Infect.Onset.Summary <- summary(survfit.Infection.Onset)

#SICKO.Stat 9
SICKO.Statistics[9,1] <- "Infection Onset Control v Test"
SICKO.Statistics[9,2] <- "Log-Rank Test"
SICKO.Statistics[9,3] <- Infection.Onset$pvalue
SICKO.Statistics[9,4] <- "Chi-squared"
SICKO.Statistics[9,5] <- Infection.Onset$chisq
SICKO.Statistics[9,6] <- "df"
SICKO.Statistics[9,7] <- 1
SICKO.Statistics[9,8] <- "Pg 6"
SICKO.Statistics[9,9] <- Infect.Onset.Summary$table[1,5]
SICKO.Statistics[9,10] <- Infect.Onset.Summary$table[2,5]
SICKO.Statistics[9,12] <- Infect.Onset.Summary$table[1,7]
SICKO.Statistics[9,13] <- Infect.Onset.Summary$table[2,7]

#Analyze regression of individual animals and compare regression scores across conditions________________________________________________________________________________________
#Fix Is.Infected dataframe to get rid of animals with no usable data for calculating Infection progression

#Find and remove worthless worms with no usable data from animal (will create issues for regression)
Is.Infected$Worthless.Worm <- NA

for (i in 1:nrow(Is.Infected)) {
  Worm.c <- Is.Infected[i,16:24]
  a <- Worm.c == -1
  a[is.na(a)] =FALSE
  b <- is.nan(unlist(Worm.c, use.names = FALSE))
  if(sum(a+b) == 9) {
    Is.Infected$Worthless.Worm[i] <- 1
  } else {
    Is.Infected$Worthless.Worm[i] <- 0
  }
}

#Remove Worthless worms
Infection.Progression.DF <- Is.Infected[!(Is.Infected$Worthless.Worm == 1),]

#Create smaller DF with only infection area values for LM
Infection.Progression.Area <-  Infection.Progression.DF[c('Area_data1', 'Area_data2', 'Area_data3', 'Area_data4', 'Area_data5', 'Area_data6', 'Area_data7', 'Area_data8', 'Area_data9')]

#Get rid of -1 values from death or censorship so not calculated in regression
Infection.Progression.Area[Infection.Progression.Area == -1] <- NaN

#Calculate linear groControl Labelh of infection in every individual animal 
i <- 1
for(i in 1:dim(Infection.Progression.Area)[1]) {
  area.c <- as.numeric(Infection.Progression.Area[i,])
  days.c <- 1:9

  lm.c <- lm(area.c~days.c)
  sum.c <- summary(lm.c)

  # print(sum.c$coef[2])  
  Infection.Progression.DF$LM.Slope.Area[i] <- sum.c$coef[2]
}


#Repeat same process with integrated intensity data

#Create smaller DF with only infection intensity values for LM
Infection.Progression.Intensity <-  Infection.Progression.DF[c('Intensity_data1', 'Intensity_data2', 'Intensity_data3', 'Intensity_data4', 'Intensity_data5', 'Intensity_data6', 'Intensity_data7', 'Intensity_data8', 'Intensity_data9')]

#Get rid of -1 values so not calculated in regression
Infection.Progression.Intensity[Infection.Progression.Intensity == -1] <- NaN


i <- 1
for(i in 1:dim(Infection.Progression.Intensity)[1]) {
  Intensity.c <- as.numeric(Infection.Progression.Intensity[i,])
  days.c <- 1:9
  
  lm.c <- lm(Intensity.c~days.c)
  sum.c <- summary(lm.c)
  
  # print(sum.c$coef[2])  
  Infection.Progression.DF$LM.Slope.Intensity[i] <- sum.c$coef[2]
}

#Single out dead animals
Infection.Progression.DF.Dead <- Infection.Progression.DF[Infection.Progression.DF$Is.Dead == 1,]
Infection.Progression.DF.Dead.Control <- Infection.Progression.DF.Dead[Infection.Progression.DF.Dead$Labels == "Control Label",] #_________________________________Input Data_________________________
Infection.Progression.DF.Dead.Test <- Infection.Progression.DF.Dead[Infection.Progression.DF.Dead$Labels != "Control Label",] #_________________________________Input Data_________________________

#Progression of Infection (Sloped of LM) vs Survival in Conditions_____________________________________________________________________________________________________________________________________

N.Graph.Y.Max <- max(Infection.Progression.DF.Dead$LM.Slope.Area, na.rm = TRUE)
N.Graph.Y.Min <- min(Infection.Progression.DF.Dead$LM.Slope.Area, na.rm = TRUE)

#Perform LM for infection AREA progression vs death for Control (Figure 7)
Progression.V.Death.Area.Control <- lm(formula = LM.Slope.Area ~ Last.Day.of.Observation, data = Infection.Progression.DF.Dead.Control)
LM.Control.Area <- summary(Progression.V.Death.Area.Control)
Summary.LM.Control.Area <- summary(LM.Control.Area)
plot(formula = LM.Slope.Area/1000 ~ Day.of.Death, data = Infection.Progression.DF.Dead.Control,
     main = "Control Label Infection Area Progression v Death",
     cex.main = 1.2,
     ylab = "Infection Area Slope (1K pixels/day)",
     xlab = "Survival Post Challenge (days)",
     ylim = c(N.Graph.Y.Min/1000,N.Graph.Y.Max/1000), #######___Dependent on Data_____
     cex = 2,
     lwd = 2,
     cex.lab = 1.2)
abline(lm(formula = LM.Slope.Area/1000 ~ Last.Day.of.Observation, data = Infection.Progression.DF.Dead.Control),
       col = "red",
       lty = "dashed",
       lwd = 4)

#SICKO.Stat 10
SICKO.Statistics[10,1] <- "Control Infection Area Progression v Death"
SICKO.Statistics[10,2] <- "Linear Regression"
SICKO.Statistics[10,3] <- LM.Control.Area$coefficients[2,4]
SICKO.Statistics[10,4] <- "t"
SICKO.Statistics[10,5] <- LM.Control.Area$coefficients[2,3]
SICKO.Statistics[10,6] <- "df"
SICKO.Statistics[10,7] <- LM.Control.Area$df[2]
SICKO.Statistics[10,8] <- "Pg 7"

#Perform LM for infection AREA progression vs death in Test (Figure 8)
Progression.V.Death.Area.Test <- lm(formula = LM.Slope.Area ~ Last.Day.of.Observation, data = Infection.Progression.DF.Dead.Test)
LM.Test.Area <- summary(Progression.V.Death.Area.Test)
plot(formula = LM.Slope.Area/1000 ~ Last.Day.of.Observation, data = Infection.Progression.DF.Dead.Test,
     main = substitute(paste(bolditalic('Test Label'),bold('Area Progression v Death'))), #______Input Data____
     cex.main = 1.2,
     ylab = "Infection Area Slope (1K pixels/day)",
     xlab = "Life Remaining (days)",
     ylim = c(N.Graph.Y.Min/1000,N.Graph.Y.Max/1000), #######___Dependent on Data_____
     cex = 2,
     lwd = 2,
     cex.lab = 1.2)
abline(lm(formula = LM.Slope.Area/1000 ~ Last.Day.of.Observation, data = Infection.Progression.DF.Dead.Test),
       col = "red",
       lty = "dashed",
       lwd = 4)

#SICKO.Stat 11
SICKO.Statistics[11,1] <- "Test Infection Area Progression v Death"
SICKO.Statistics[11,2] <- "Linear Regression"
SICKO.Statistics[11,3] <- LM.Test.Area$coefficients[2,4]
SICKO.Statistics[11,4] <- "t"
SICKO.Statistics[11,5] <- LM.Test.Area$coefficients[2,3]
SICKO.Statistics[11,6] <- "df"
SICKO.Statistics[11,7] <- LM.Test.Area$df[2]
SICKO.Statistics[11,8] <- "Pg 8"

#Perform LM for infection progression vs death INTENSITY in Control

P.Graph.Y.Max <- max(Infection.Progression.DF.Dead$LM.Slope.Intensity, na.rm = TRUE)
P.Graph.Y.Min <- min(Infection.Progression.DF.Dead$LM.Slope.Intensity, na.rm = TRUE)

#Control Label Infection Intensity Progression (figure 9)
Progression.V.Death.Intensity.Control <- lm(formula = LM.Slope.Intensity ~ Last.Day.of.Observation, data = Infection.Progression.DF.Dead.Control)
LM.Control.Intensity <- summary(Progression.V.Death.Intensity.Control)
plot(formula = LM.Slope.Intensity/10000000 ~ Last.Day.of.Observation, data = Infection.Progression.DF.Dead.Control,
     main = "Control Label Intensity Progression v Death",
     cex.main = 1.2,
     ylab = "Infection Intensity Slope (a.u./day)",
     xlab = "Life Remaining (days)",
     ylim = c(P.Graph.Y.Min/10000000,P.Graph.Y.Max/10000000), #######___Dependent on Data_____
     cex = 2,
     lwd = 2,
     cex.lab = 1.2)
abline(lm(formula = LM.Slope.Intensity/10000000 ~ Last.Day.of.Observation, data = Infection.Progression.DF.Dead.Control),
       col = "red",
       lty = "dashed",
       lwd = 4)

#SICKO.Stat 12
SICKO.Statistics[12,1] <- "Control Infection Intensity Progression v Death"
SICKO.Statistics[12,2] <- "Linear Regression"
SICKO.Statistics[12,3] <- LM.Control.Intensity$coefficients[2,4]
SICKO.Statistics[12,4] <- "t"
SICKO.Statistics[12,5] <- LM.Control.Intensity$coefficients[2,3]
SICKO.Statistics[12,6] <- "df"
SICKO.Statistics[12,7] <- LM.Control.Intensity$df[2]
SICKO.Statistics[12,8] <- "Pg 9"

#Perform LM for infection Intensity progression vs in Test (figure 10)
Progression.V.Death.Intensity.Test <- lm(formula = LM.Slope.Intensity ~ Last.Day.of.Observation, data = Infection.Progression.DF.Dead.Test)
LM.Test.Intensity <- summary(Progression.V.Death.Intensity.Test)
plot(formula = LM.Slope.Intensity/10000000 ~ Last.Day.of.Observation, data = Infection.Progression.DF.Dead.Test,
     main = substitute(paste(bolditalic('Test Label'),bold('Intensity Progression v Death'))), #______Input Data____
     cex.main = 1.2,
     ylab = "Infection Intensity Slope (a.u./day)",
     xlab = "Life Remaining (days)",
#     ylim = c(P.Graph.Y.Min/10000000,P.Graph.Y.Max/10000000), #######___Dependent on Data_____ ##Come back_____ comment out only for bacteria strain comparison
     cex = 2,
     lwd = 2,
     cex.lab = 1.2)
abline(lm(formula = LM.Slope.Intensity/10000000 ~ Last.Day.of.Observation, data = Infection.Progression.DF.Dead.Test),
       col = "red",
       lty = "dashed",
       lwd = 4)

#SICKO.Stat 13
SICKO.Statistics[13,1] <- "Test Infection Intensity Progression v Death"
SICKO.Statistics[13,2] <- "Linear Regression"
SICKO.Statistics[13,3] <- LM.Test.Intensity$coefficients[2,4]
SICKO.Statistics[13,4] <- "t"
SICKO.Statistics[13,5] <- LM.Test.Intensity$coefficients[2,3]
SICKO.Statistics[13,6] <- "df"
SICKO.Statistics[13,7] <- LM.Test.Intensity$df[2]
SICKO.Statistics[13,8] <- "Pg 10"

#Test if regression data for Conditions is significantly different_____________________________________________________________________________________________________________________________________
Infection.Progression.DF.Infected <- Infection.Progression.DF[!(Infection.Progression.DF$LM.Slope.Area == 0.0000),]

T_Test_Area <- t.test(LM.Slope.Area ~ Labels, data = Infection.Progression.DF.Infected)

T_Test_Intensity <- t.test(LM.Slope.Intensity ~ Labels, data = Infection.Progression.DF.Infected)

# Plot Slope of Infection Area Progression by populations (Figure 11)
vioplot(LM.Slope.Area/100 ~ Labels, data = Infection.Progression.DF.Infected,
        main = "Infection Area Progression", 
        cex.main = 1.2,
        ylab = "Infection Area Slope (1K pixels/day)",
        xlab = "Conditions",
        cex = 2,
        lwd = 2,
        col = c("grey","red"),
        cex.lab = 1.2)

#SICKO.Stat 14
SICKO.Statistics[14,1] <- "Infection Area Progression Control v Test"
SICKO.Statistics[14,2] <- "t-test"
SICKO.Statistics[14,3] <- T_Test_Area$p.value
SICKO.Statistics[14,4] <- "t"
SICKO.Statistics[14,5] <- T_Test_Area$statistic
SICKO.Statistics[14,6] <- "df"
SICKO.Statistics[14,7] <- T_Test_Area$parameter
SICKO.Statistics[14,8] <- "Pg 11"
SICKO.Statistics[14,9] <- T_Test_Area$estimate[1]
SICKO.Statistics[14,10] <- T_Test_Area$estimate[2]
SICKO.Statistics[14,11] <- T_Test_Area$stderr

# Plot Slope of Infection Intensity Progression by populations (Figure 12)
vioplot(LM.Slope.Intensity/10000000 ~ Labels, data = Infection.Progression.DF.Infected,
        main = "Infection Intensity Progression", 
        cex.main = 1.2,
        ylab = "Infection Intensity Slope (a.u./day)",
        xlab = "Conditions",
        cex = 2,
        lwd = 2,
        col = c("grey","red"),
        cex.lab = 1.2)

#SICKO.Stat 15
SICKO.Statistics[15,1] <- "Infection Intensity Progression Control v Test"
SICKO.Statistics[15,2] <- "t-test"
SICKO.Statistics[15,3] <- T_Test_Intensity$p.value
SICKO.Statistics[15,4] <- "t"
SICKO.Statistics[15,5] <- T_Test_Intensity$statistic
SICKO.Statistics[15,6] <- "df"
SICKO.Statistics[15,7] <- T_Test_Intensity$parameter
SICKO.Statistics[15,8] <- "Pg 12"
SICKO.Statistics[15,9] <- T_Test_Intensity$estimate[1]
SICKO.Statistics[15,10] <- T_Test_Intensity$estimate[2]
SICKO.Statistics[15,11] <- T_Test_Intensity$stderr

 # #Perform LM for infection progression vs death__________________________________________________________________________________For All animals all conditions**********
 # Progression.V.Death.Intensity <- lm(formula = LM.Slope.Intensity ~ Last.Day.of.Observation, fill = Infection.Progression.DF$Labels, data = Infection.Progression.DF)
# summary(Progression.V.Death.Intensity)
# plot(formula = LM.Slope.Intensity ~ Last.Day.of.Observation, data = Infection.Progression.DF)
# abline(lm(formula = LM.Slope.Intensity ~ Last.Day.of.Observation, data = Infection.Progression.DF))



# #Infection quantity last day of ovservation vs death___________________________________________________________________________Interesting distribution but infection intensity vs days remaining graph (later) may be better metric
# Control.DF.Infected.Dead <- Control.Infected.DF[Control.Infected.DF$Is.Dead == 1,]
# Final.Intensity.V.Death <- lm(formula = Intensity.at.Last.Day.of.Observation ~ Last.Day.of.Observation, data = Control.DF.Infected.Dead)
# plot(formula = Intensity.at.Last.Day.of.Observation ~ Last.Day.of.Observation, data = Control.DF.Infected.Dead)
# summary(Final.Intensity.V.Death)

#Time of Initial Infection correlation with death_____________________________________________________________________________________________________________________________________

Q.Graph.Y.Max <- max(Is.Infected$First.Day.of.nonzero.data, na.rm =TRUE)
Q.Graph.Y.Min <- min(Is.Infected$First.Day.of.nonzero.data, na.rm =TRUE)
Q.Graph.X.Max <- max(Is.Infected$Last.Day.of.Observation, na.rm = TRUE)
Q.Graph.X.Min <- min(Is.Infected$Last.Day.of.Observation, na.rm = TRUE)

#Does initial infection correlate to death? Control  (Figure 13)______Many warning messages
Control.DF.Infected.Dead <- Control.Infected.DF[Control.Infected.DF$Is.Dead == 1,]
Control.First.Infection.V.Death <- lm(formula = First.Day.of.nonzero.data ~ Last.Day.of.Observation, data = Control.DF.Infected.Dead)
plot(y=jitter(Control.DF.Infected.Dead$First.Day.of.nonzero.data, 0.5),  x=jitter(Control.DF.Infected.Dead$Day.of.Death, 0.5), data = Control.DF.Infected.Dead,
     main = "Control Label Initial Infection v Survival", #_____Input Data_____
     cex.main = 1.2,
     ylab = "First Day Detected Infection",
     xlab = "Day of Death",
     # ylim = c(0,Q.Graph.Y.Max),
     # xlim = c(Q.Graph.X.Min, Q.Graph.X.Max),
     cex = 2,
     lwd = 2,
     xaxt = 'n',
     cex.lab = 1.2)
axis(side = 1, at=seq(0,9, by = 1))

abline(lm(formula = First.Day.of.nonzero.data ~ Day.of.Death, data = Control.DF.Infected.Dead),
       col = "red",
       lty = "dashed",
       lwd = 4)
Control.LM.Infect.Death.Summary <- summary(Control.First.Infection.V.Death)

#SICKO.Stat 16
SICKO.Statistics[16,1] <- "Control Infection Initiation v Death"
SICKO.Statistics[16,2] <- "Linear Regression"
SICKO.Statistics[16,3] <- Control.LM.Infect.Death.Summary$coefficients[2,4]
SICKO.Statistics[16,4] <- "t"
SICKO.Statistics[16,5] <- Control.LM.Infect.Death.Summary$coefficients[2,3]
SICKO.Statistics[16,6] <- "df"
SICKO.Statistics[16,7] <- Control.LM.Infect.Death.Summary$df[2]
SICKO.Statistics[16,8] <- "Pg 13"

#With test Condition (Figure 14
Test.DF.Infected.Dead <- Test.Infected.DF[Test.Infected.DF$Is.Dead == 1,]
Test.First.Infection.V.Death <- lm(formula = First.Day.of.nonzero.data ~ Last.Day.of.Observation, data = Test.DF.Infected.Dead)
plot(y=jitter(Test.DF.Infected.Dead$First.Day.of.nonzero.data, 0.5),  x=jitter(Test.DF.Infected.Dead$Day.of.Death, 0.5), data = Test.DF.Infected.Dead,
     main = substitute(paste(bolditalic('Test Label'),bold('Initial Infection v Survival'))), 
     cex.main = 1.2,
     ylab = "First Day Detected Infection",
     xlab = "Day of Death",
     ylim = c(0,Q.Graph.Y.Max),
     xlim = c(Q.Graph.X.Min, Q.Graph.X.Max),
     cex = 2,
     lwd = 2,
     xaxt = 'n',
     cex.lab = 1.2)
axis(side = 1, at=seq(0,9, by = 1))

abline(lm(formula = First.Day.of.nonzero.data ~ Day.of.Death, data = Test.DF.Infected.Dead),
       col = "red",
       lty = "dashed",
       lwd = 4)
Test.LM.Infect.Death.Summary <- summary(Test.First.Infection.V.Death)

#SICKO.Stat 17
SICKO.Statistics[17,1] <- "Test Infection Initiation v Death"
SICKO.Statistics[17,2] <- "Linear Regression"
SICKO.Statistics[17,3] <- Test.LM.Infect.Death.Summary$coefficients[2,4]
SICKO.Statistics[17,4] <- "t"
SICKO.Statistics[17,5] <- Test.LM.Infect.Death.Summary$coefficients[2,3]
SICKO.Statistics[17,6] <- "df"
SICKO.Statistics[17,7] <- Test.LM.Infect.Death.Summary$df[2]
SICKO.Statistics[17,8] <- "Pg 14"

#Look at infection detection to death times compared between the two conditions
Is.Infected.Animal <- Is.Infected[Is.Infected$Infection.Status == "Infected",]
Is.Infected.Dead <- Is.Infected.Animal[Is.Infected.Animal$Is.Dead == 1,]

Is.Infected.Dead$Infection.To.Death <- NA
for (i in 1:length(Is.Infected.Dead$Infection.To.Death)){
  x <- Is.Infected.Dead$Last.Day.of.Observation[i]
  y <- Is.Infected.Dead$First.Day.of.nonzero.data[i]
  Is.Infected.Dead$Infection.To.Death[i] <-  x-y+1
}

#Violin plot of length of infection to death (Figure 15)

#Make DF containing infection to death along with infection to last day observation
Is.Infected.Animal$Infection.To.Final.Observation <- NA
for (i in 1:length(Is.Infected.Animal$Infection.To.Final.Observation)){
  x <- Is.Infected.Animal$Last.Day.of.Observation[i]
  y <- Is.Infected.Animal$First.Day.of.nonzero.data[i]
  Is.Infected.Animal$Infection.To.Final.Observation[i] <- x-y+1
}

#Sort out if animals fled
Is.Infected.Animal$Did.Flea <- NA
Is.Infected.Animal.Present <- Is.Infected.Animal[(Is.Infected.Animal$Is.Dead == 1 & Is.Infected.Animal$Last.Day.of.Observation < 9)
                                                 |(Is.Infected.Animal$Last.Day.of.Observation == 9 & Is.Infected.Animal$Is.Dead == 0),]


T_Test_Infection_To_Death <- t.test(Is.Infected.Animal.Present$Infection.To.Final.Observation ~ Labels, data = Is.Infected.Animal.Present)
vioplot(Is.Infected.Animal.Present$Infection.To.Final.Observation ~ Labels, data = Is.Infected.Animal.Present,
        main = "Infection Detection to Death", 
        cex.main = 1.2,
        ylab = "Time Infection to Death (days)",
        xlab = "Conditions",
        cex = 2,
        lwd = 2,
        col = c("grey","red"),
        cex.lab = 1.2)

#SICKO.Stat 18
SICKO.Statistics[18,1] <- "Infection Detection to Death Control v Test"
SICKO.Statistics[18,2] <- "t-test"
SICKO.Statistics[18,3] <- T_Test_Infection_To_Death$p.value
SICKO.Statistics[18,4] <- "t"
SICKO.Statistics[18,5] <- T_Test_Infection_To_Death$statistic
SICKO.Statistics[18,6] <- "df"
SICKO.Statistics[18,7] <- T_Test_Infection_To_Death$parameter
SICKO.Statistics[18,8] <- "Pg 15"
SICKO.Statistics[18,9] <- T_Test_Infection_To_Death$estimate[1]
SICKO.Statistics[18,10] <- T_Test_Infection_To_Death$estimate[2]
SICKO.Statistics[18,11] <- T_Test_Infection_To_Death$stderr

#Test if Cumulative Infection is related to earlier death___________________________________Problematic of course worms that live longer will have more accumulated intensity (maybe divide by days?)____________________________________________________________
Dead.Worms <- Data.set[Data.set$Is.Dead == 1,]
Dead.Worms.Control <- Dead.Worms[Dead.Worms$Labels == "Control Label",]
Dead.Worms.Test <- Dead.Worms[Dead.Worms$Labels != "Control Label",]

R.Graph.Y.Max <- max(Dead.Worms$Intensity.Integrated.Across.Time, na.rm = TRUE)
R.Graph.Y.Min <- min(Dead.Worms$Intensity.Integrated.Across.Time, na.rm = TRUE)

#Plot Cumulative infection Intensity of Control Condition against survival (Figure 16)
Control.Total.Intensity.V.Death <- lm(formula = Intensity.Integrated.Across.Time ~ Last.Day.of.Observation, data = Dead.Worms.Control)
plot(y = Dead.Worms.Control$Intensity.Integrated.Across.Time/10000000, x = Dead.Worms.Control$Last.Day.of.Observation + 1,
     main = "Control Label Cumulative Infection Intensity v Survival", #_____Input Data_____
     cex.main = 1.2,
     ylab = "Cumulative Infection Intensity (a.u.)",
     xlab = "Day of Death",
     ylim = c(R.Graph.Y.Min/10000000, R.Graph.Y.Max/10000000),
     cex = 2,
     lwd = 2,
     cex.lab = 1.2)
Control.LM.Cumulative.Intensity.V.Death <- summary(Control.Total.Intensity.V.Death)

#SICKO.Stat 19
SICKO.Statistics[19,1] <- "Control Cumulative Infection Intensity v Death"
SICKO.Statistics[19,2] <- "Linear Regression"
SICKO.Statistics[19,3] <- Control.LM.Cumulative.Intensity.V.Death$coefficients[2,4]
SICKO.Statistics[19,4] <- "t"
SICKO.Statistics[19,5] <- Control.LM.Cumulative.Intensity.V.Death$coefficients[2,3]
SICKO.Statistics[19,6] <- "df"
SICKO.Statistics[19,7] <- Control.LM.Cumulative.Intensity.V.Death$df[2]
SICKO.Statistics[19,8] <- "Pg 16"


#Plot Cumulative infection Intensity of Test Condition against survival (Figure 17)
Test.Total.Intensity.V.Death <- lm(formula = Intensity.Integrated.Across.Time ~ Last.Day.of.Observation, data = Dead.Worms.Test)
plot(y = Dead.Worms.Test$Intensity.Integrated.Across.Time/10000000, x = Dead.Worms.Test$Last.Day.of.Observation + 1,
     main = substitute(paste(bolditalic('Test Label'),bold('Cumulative Infection Intensity v Survival'))), #_____Input Data_____
     cex.main = 1.2,
     ylab = "Cumulative Infection Intensity (a.u.)",
     xlab = "Day of Death",
     # ylim = c(R.Graph.Y.Min/10000000, R.Graph.Y.Max/10000000),   ####Come back only change for bacteria comparison
     cex = 2,
     lwd = 2,
     cex.lab = 1.2)
Test.LM.Cumulative.Intensity.V.Death <- summary(Test.Total.Intensity.V.Death)

#SICKO.Stat 20
SICKO.Statistics[20,1] <- "Test Cumulative Infection Intensity v Death"
SICKO.Statistics[20,2] <- "Linear Regression"
SICKO.Statistics[20,3] <- Test.LM.Cumulative.Intensity.V.Death$coefficients[2,4]
SICKO.Statistics[20,4] <- "t"
SICKO.Statistics[20,5] <- Test.LM.Cumulative.Intensity.V.Death$coefficients[2,3]
SICKO.Statistics[20,6] <- "df"
SICKO.Statistics[20,7] <- Test.LM.Cumulative.Intensity.V.Death$df[2]
SICKO.Statistics[20,8] <- "Pg 17"

S.Graph.Y.Max <- max(Dead.Worms$Area.Integrated.Across.Time, na.rm = TRUE)
S.Graph.Y.Min <- min(Dead.Worms$Area.Integrated.Across.Time, na.rm = TRUE)

#Plot Cumulative infection Area of Control Condition against survival (Figure 18)
Control.Total.Area.V.Death <- lm(formula = Area.Integrated.Across.Time ~ Last.Day.of.Observation, data = Dead.Worms.Control)
plot(y = Dead.Worms.Control$Area.Integrated.Across.Time/1000, x = Dead.Worms.Control$Last.Day.of.Observation + 1,
     main = "Control Label Cumulative Infection Area v Survival", #_____Input Data_____
     cex.main = 1.2,
     ylab = "Cumulative Infection Area (a.u.)",
     xlab = "Day of Death",
     ylim = c(S.Graph.Y.Min/1000,S.Graph.Y.Max/1000),
     cex = 2,
     lwd = 2,
     cex.lab = 1.2)
Control.LM.Cumulative.Area.V.Death <- summary(Control.Total.Area.V.Death)

#SICKO.Stat 21
SICKO.Statistics[21,1] <- "Control Cumulative Infection Area v Death"
SICKO.Statistics[21,2] <- "Linear Regression"
SICKO.Statistics[21,3] <- Control.LM.Cumulative.Area.V.Death$coefficients[2,4]
SICKO.Statistics[21,4] <- "t"
SICKO.Statistics[21,5] <- Control.LM.Cumulative.Area.V.Death$coefficients[2,3]
SICKO.Statistics[21,6] <- "df"
SICKO.Statistics[21,7] <- Control.LM.Cumulative.Area.V.Death$df[2]
SICKO.Statistics[21,8] <- "Pg 18"

#Plot Cumulative infection Area of Test Condition against survival (Figure 19)
Test.Total.Area.V.Death <- lm(formula = Area.Integrated.Across.Time ~ Last.Day.of.Observation, data = Dead.Worms.Test)
plot(y = Dead.Worms.Test$Area.Integrated.Across.Time/1000, x = Dead.Worms.Test$Last.Day.of.Observation + 1,
     main = substitute(paste(bolditalic('Test Label'),bold(' Cumulative Infection Area v Survival'))), #_____Input Data_____
     cex.main = 1.2,
     ylab = "Cumulative Infection Area (a.u.)",
     xlab = "Day of Death",
     ylim = c(S.Graph.Y.Min/1000,S.Graph.Y.Max/1000),
     cex = 2,
     lwd = 2,
     cex.lab = 1.2)
Test.LM.Cumulative.Area.V.Death <- summary(Test.Total.Area.V.Death)

#SICKO.Stat 22
SICKO.Statistics[22,1] <- "Test Cumulative Infection Area v Death"
SICKO.Statistics[22,2] <- "Linear Regression"
SICKO.Statistics[22,3] <- Test.LM.Cumulative.Area.V.Death$coefficients[2,4]
SICKO.Statistics[22,4] <- "t"
SICKO.Statistics[22,5] <- Test.LM.Cumulative.Area.V.Death$coefficients[2,3]
SICKO.Statistics[22,6] <- "df"
SICKO.Statistics[22,7] <- Test.LM.Cumulative.Area.V.Death$df[2]
SICKO.Statistics[22,8] <- "Pg 19"

#Early vs Late Infection and SURVIVAL____________________________________________________________________________________________

#Specifically look at animals who lived longer than the first 2 days as animal infection may be in plateau

Is.Infected.Dead$Infection.Period <- NA
for (i in 1:length(Is.Infected.Dead$Infection.Period)){
  if (Is.Infected.Dead$First.Day.of.nonzero.data[i] > 2) {
    Is.Infected.Dead$Infection.Period[i] <- "Late"
  } else if(Is.Infected.Dead$First.Day.of.nonzero.data[i] < 3){
    Is.Infected.Dead$Infection.Period[i] <- "Early"
  }
}

Control.Is.Infected.Dead <- Is.Infected.Dead[Is.Infected.Dead$Labels == "Control Label",]#________Insert Data__________
Test.Is.Infected.Dead <- Is.Infected.Dead[Is.Infected.Dead$Labels != 'Control Label',]#________Insert Data__________

#Plot Late and early infection vs length infection to death in Control (Figure 20)
vioplot(Infection.To.Death ~ Infection.Period, data = Control.Is.Infected.Dead,
        main = "Control Label Early v Late Infection", #________Insert Data__________
        cex.main = 1.2,
        ylab = "Time Infection to Death (Days)",
        xlab = "Infection Detection",
        cex = 2,
        lwd = 2,
        col = c("grey","red"),
        cex.lab = 1.2)

#Plot Late and early infection vs length infection to death in Test (Figure 21)
vioplot(Infection.To.Death ~ Infection.Period, data = Test.Is.Infected.Dead,
        main = substitute(paste(bolditalic('Test Label'),bold(' Early v Late Infection'))), #________Insert Data__________
        cex.main = 1.2,
        ylab = "Time Infection to Death (Days)",
        xlab = "Infection Detection",
        cex = 2,
        lwd = 2,
        col = c("grey","red"),
        cex.lab = 1.2)

Control.Is.Infected.Dead$Infection.To.Death

Condition.Variable.Control <- table(Control.Is.Infected.Dead$Infection.Period)

if(!any(Condition.Variable.Control < 3) & length(unique(Control.Is.Infected.Dead$Infection.Period)) >1){
  Control_T_Test_Early_V_Late_Infection <- t.test(Infection.To.Death ~ Infection.Period, data = Control.Is.Infected.Dead) #Come back not normlly distributed wilcox
  
  #SICKO.Stat 23
  SICKO.Statistics[23,1] <- "Control Early v Late Infection Initiation to Death"
  SICKO.Statistics[23,2] <- "t-test"
  SICKO.Statistics[23,3] <- Control_T_Test_Early_V_Late_Infection$p.value
  SICKO.Statistics[23,4] <- "t"
  SICKO.Statistics[23,5] <- Control_T_Test_Early_V_Late_Infection$statistic
  SICKO.Statistics[23,6] <- "df"
  SICKO.Statistics[23,7] <- Control_T_Test_Early_V_Late_Infection$parameter
  SICKO.Statistics[23,8] <- "Pg 20"
  SICKO.Statistics[23,9] <- Control_T_Test_Early_V_Late_Infection$estimate[1]
  SICKO.Statistics[23,10] <- Control_T_Test_Early_V_Late_Infection$estimate[2]
  SICKO.Statistics[23,11] <- Control_T_Test_Early_V_Late_Infection$stderr

} else {
  #SICKO.Stat 23
  SICKO.Statistics[23,1] <- "Control Early v Late Infection Initiation to Death"
  SICKO.Statistics[23,2] <- "t-test"
  SICKO.Statistics[23,3] <- NA
  SICKO.Statistics[23,4] <- "t"
  SICKO.Statistics[23,5] <- NA
  SICKO.Statistics[23,6] <- "df"
  SICKO.Statistics[23,7] <- NA
  SICKO.Statistics[23,8] <- "Pg 20"
  SICKO.Statistics[23,9] <- NA
  SICKO.Statistics[23,10] <- NA
  SICKO.Statistics[23,11] <- NA
}

Condition.Variable.Test <- table(Test.Is.Infected.Dead$Infection.Period)


if(!any(Condition.Variable.Test < 3) & length(unique(Test.Is.Infected.Dead$Infection.Period)) >1) {
  Test_T_Test_Early_V_Late_Infection <- t.test(Infection.To.Death ~ Infection.Period, data = Test.Is.Infected.Dead)
  
  #SICKO.Stat 24
  SICKO.Statistics[24,1] <- "Test Early v Late Infection Initiation to Death"
  SICKO.Statistics[24,2] <- "t-test"
  SICKO.Statistics[24,3] <- Test_T_Test_Early_V_Late_Infection$p.value
  SICKO.Statistics[24,4] <- "t"
  SICKO.Statistics[24,5] <- Test_T_Test_Early_V_Late_Infection$statistic
  SICKO.Statistics[24,6] <- "df"
  SICKO.Statistics[24,7] <- Test_T_Test_Early_V_Late_Infection$parameter
  SICKO.Statistics[24,8] <- "Pg 21"
  SICKO.Statistics[24,9] <- Test_T_Test_Early_V_Late_Infection$estimate[1]
  SICKO.Statistics[24,10] <- Test_T_Test_Early_V_Late_Infection$estimate[2]
  SICKO.Statistics[24,11] <- Test_T_Test_Early_V_Late_Infection$stderr
} else {

  #SICKO.Stat 24
  SICKO.Statistics[24,1] <- "Test Early v Late Infection Initiation to Death"
  SICKO.Statistics[24,2] <- "t-test"
  SICKO.Statistics[24,3] <- NA
  SICKO.Statistics[24,4] <- "t"
  SICKO.Statistics[24,5] <- NA
  SICKO.Statistics[24,6] <- "df"
  SICKO.Statistics[24,7] <- NA
  SICKO.Statistics[24,8] <- "Pg 21"
  SICKO.Statistics[24,9] <- NA
  SICKO.Statistics[24,10] <- NA
  SICKO.Statistics[24,11] <- NA
}


#Regression/progression of Infection Compared to length from infection to death_____________________________________________________________________________________________________________________________
Infection.Progression.DF.Dead <- Infection.Progression.DF[Infection.Progression.DF$Is.Dead == 1,]
Infection.Progression.DF.Dead.Infected <- Infection.Progression.DF.Dead[Infection.Progression.DF.Dead$Infection.Status == "Infected",]
Infection.Progression.DF.Dead.Infected$Infection.To.Death <- NA
for (i in 1:length(Infection.Progression.DF.Dead.Infected$Infection.To.Death)){
  x <- Infection.Progression.DF.Dead.Infected$Day.of.Death[i]
  y <- Infection.Progression.DF.Dead.Infected$First.Day.of.nonzero.data[i]
  Infection.Progression.DF.Dead.Infected$Infection.To.Death[i] <-  x-y
}

LM.Control.Infection.Resiliance <- Infection.Progression.DF.Dead.Infected[Infection.Progression.DF.Dead.Infected$Labels == "Control Label",]#_____Insert Data______
LM.Test.Infection.Resiliance <- Infection.Progression.DF.Dead.Infected[Infection.Progression.DF.Dead.Infected$Labels != "Control Label",]#_____Insert Data______

#Linear Regression of Infection AREA progression v length infection to death

T.Graph.Y.Max <- max(Infection.Progression.DF.Dead.Infected$LM.Slope.Area, na.rm = TRUE)
T.Graph.Y.Min <- min(Infection.Progression.DF.Dead.Infected$LM.Slope.Area, na.rm = TRUE)

#Linear Regression of Infection Progression Area versus time infection to death In Control Condition(Figure 22)
LM.Control.Infection.Area.Progression.Resiliance <- lm(formula = LM.Slope.Area ~ Infection.To.Death, data = LM.Control.Infection.Resiliance )
Summary.LM.Control.Infection.Area.Progression.Resiliance <- summary(LM.Control.Infection.Area.Progression.Resiliance)
plot(formula = LM.Slope.Area/1000 ~ Infection.To.Death, data = LM.Control.Infection.Resiliance,
     main = "Control Label Infection Progression Resiliance", #_____Input Data_____
     cex.main = 1.2,
     ylab = "Infection Area Slope (1K pixels/day)",
     xlab = "Infection Initiation to Death (days)",
     ylim = c(T.Graph.Y.Min/1000, T.Graph.Y.Max/1000),
     cex = 2,
     lwd = 2,
     cex.lab = 1.2)
abline(lm(formula = LM.Slope.Area/1000 ~ Infection.To.Death, data = LM.Control.Infection.Resiliance),
       col = "red",
       lty = "dashed",
       lwd = 4)

#SICKO.Stat 25
SICKO.Statistics[25,1] <- "Control Infection Area Progression v Infection initiation to Death"
SICKO.Statistics[25,2] <- "Linear Regression"
SICKO.Statistics[25,3] <- Summary.LM.Control.Infection.Area.Progression.Resiliance$coefficients[2,4]
SICKO.Statistics[25,4] <- "t"
SICKO.Statistics[25,5] <- Summary.LM.Control.Infection.Area.Progression.Resiliance$coefficients[2,3]
SICKO.Statistics[25,6] <- "df"
SICKO.Statistics[25,7] <- Summary.LM.Control.Infection.Area.Progression.Resiliance$df[2]
SICKO.Statistics[25,8] <- "Pg 22"

#Linear Regression of Infection Progression Area versus time infection to death in Test Condition (Figure 23)
LM.Test.Infection.Area.Progression.Resiliance <- lm(formula = LM.Slope.Area ~ Infection.To.Death, data = LM.Test.Infection.Resiliance )
Summary.LM.Test.Infection.Area.Progression.Resiliance <- summary(LM.Test.Infection.Area.Progression.Resiliance)
plot(formula = LM.Slope.Area/1000 ~ Infection.To.Death, data = LM.Test.Infection.Resiliance,
     main = substitute(paste(bolditalic('Test Label'),bold(' Infection Progression Resiliance'))), #_____Input Data_____
     cex.main = 1.2,
     ylab = "Infection Area Slope (1K pixels/day)",
     xlab = "Infection Initiation to Death (days)",
     ylim = c(T.Graph.Y.Min/1000, T.Graph.Y.Max/1000),
     cex = 2,
     lwd = 2,
     cex.lab = 1.2)
abline(lm(formula = LM.Slope.Area/1000 ~ Infection.To.Death, data = LM.Test.Infection.Resiliance),
       col = "red",
       lty = "dashed",
       lwd = 4)

#SICKO.Stat 26
SICKO.Statistics[26,1] <- "Test Infection Area Progression v Infection initiation to Death"
SICKO.Statistics[26,2] <- "Linear Regression"
SICKO.Statistics[26,3] <- Summary.LM.Test.Infection.Area.Progression.Resiliance$coefficients[2,4]
SICKO.Statistics[26,4] <- "t"
SICKO.Statistics[26,5] <- Summary.LM.Test.Infection.Area.Progression.Resiliance$coefficients[2,3]
SICKO.Statistics[26,6] <- "df"
SICKO.Statistics[26,7] <- Summary.LM.Test.Infection.Area.Progression.Resiliance$df[2]
SICKO.Statistics[26,8] <- "Pg 23"

#Linear Regression of Infection Intensity progression v length infection to death

#Obtain graph limits
U.Graph.Y.Max <- max(Infection.Progression.DF.Dead.Infected$LM.Slope.Intensity, na.rm = TRUE)
U.Graph.Y.Min <- min(Infection.Progression.DF.Dead.Infected$LM.Slope.Intensity, na.rm = TRUE)

#Linear Regression of Infection Progression Intensity versus time infection to death in Control Condition (Figure 24)
LM.Control.Infection.Intensity.Progression.Resiliance <- lm(formula = LM.Slope.Intensity ~ Infection.To.Death, data = LM.Control.Infection.Resiliance )
Summary.LM.Control.Infection.Intensity.Progression.Resiliance <- summary(LM.Control.Infection.Intensity.Progression.Resiliance)
plot(formula = LM.Slope.Intensity/10000000 ~ Infection.To.Death, data = LM.Control.Infection.Resiliance,
     main = "Control Label Infection Progression Resiliance", #_____Input Data_____
     cex.main = 1.2,
     ylab = "Infection Intensity Slope (a.u./day)",
     xlab = "Infection Initiation to Death (days)",
     ylim = c(U.Graph.Y.Min/10000000,U.Graph.Y.Max/10000000),
     cex = 2,
     lwd = 2,
     cex.lab = 1.2)
abline(lm(formula = LM.Slope.Intensity/10000000 ~ Infection.To.Death, data = LM.Control.Infection.Resiliance),
       col = "red",
       lty = "dashed",
       lwd = 4)

#SICKO.Stat 27
SICKO.Statistics[27,1] <- "Control Infection Intensity Progression v Infection initiation to Death"
SICKO.Statistics[27,2] <- "Linear Regression"
SICKO.Statistics[27,3] <- Summary.LM.Control.Infection.Intensity.Progression.Resiliance$coefficients[2,4]
SICKO.Statistics[27,4] <- "t"
SICKO.Statistics[27,5] <- Summary.LM.Control.Infection.Intensity.Progression.Resiliance$coefficients[2,3]
SICKO.Statistics[27,6] <- "df"
SICKO.Statistics[27,7] <- Summary.LM.Control.Infection.Intensity.Progression.Resiliance$df[2]
SICKO.Statistics[27,8] <- "Pg 24"

#Linear Regression of Infection Progression Intensity versus time infection to death in Test Condition (Figure 25)
LM.Test.Infection.Intensity.Progression.Resiliance <- lm(formula = LM.Slope.Intensity ~ Infection.To.Death, data = LM.Test.Infection.Resiliance )
Summary.LM.Test.Infection.Intensity.Progression.Resiliance <- summary(LM.Test.Infection.Intensity.Progression.Resiliance)
plot(formula = LM.Slope.Intensity/10000000 ~ Infection.To.Death, data = LM.Test.Infection.Resiliance,
     main = substitute(paste(bolditalic('Test Label'),bold(' Infection Progression Resiliance'))), #_____Input Data_____
     cex.main = 1.2,
     ylab = "Infection Intensity Slope (a.u./day)",
     xlab = "Infection Initiation to Death (days)",
     # ylim = c(U.Graph.Y.Min/10000000,U.Graph.Y.Max/10000000),     ####Come back
     cex = 2,
     lwd = 2,
     cex.lab = 1.2)
abline(lm(formula = LM.Slope.Intensity/10000000 ~ Infection.To.Death, data = LM.Test.Infection.Resiliance),
       col = "red",
       lty = "dashed",
       lwd = 4)

#SICKO.Stat 28
SICKO.Statistics[28,1] <- "Test Infection Intensity Progression v Infection initiation to Death"
SICKO.Statistics[28,2] <- "Linear Regression"
SICKO.Statistics[28,3] <- Summary.LM.Test.Infection.Intensity.Progression.Resiliance$coefficients[2,4]
SICKO.Statistics[28,4] <- "t"
SICKO.Statistics[28,5] <- Summary.LM.Test.Infection.Intensity.Progression.Resiliance$coefficients[2,3]
SICKO.Statistics[28,6] <- "df"
SICKO.Statistics[28,7] <- Summary.LM.Test.Infection.Intensity.Progression.Resiliance$df[2]
SICKO.Statistics[28,8] <- "Pg 25"

#Regression values by condition___________________________________________________________________________________________________________

#Make DF for conditions and infection state
Infection.Progression.DF.Control <- Infection.Progression.DF[Infection.Progression.DF$Condition == "Control Label",]#_____Input Data_____
Infection.Progression.DF.Control.Infected <- Infection.Progression.DF.Control[Infection.Progression.DF.Control$Infection.Status == "Infected",]
Infection.Progression.DF.Test <- Infection.Progression.DF[Infection.Progression.DF$Condition != "Control Label",]#_____Input Data_____
Infection.Progression.DF.Test.Infected <- Infection.Progression.DF.Test[Infection.Progression.DF.Test$Infection.Status =="Infected",]
Combo.Infected <- rbind(Infection.Progression.DF.Control.Infected,Infection.Progression.DF.Test.Infected)

Avg.Progression.Area.Control <- median(Infection.Progression.DF.Control.Infected$LM.Slope.Area, na.rm=TRUE)
Avg.Progression.Area.Test <- median(Infection.Progression.DF.Test.Infected$LM.Slope.Area, na.rm=TRUE)

wilcox.test(LM.Slope.Area ~ Condition, data=Combo.Infected)

Avg.Progression.Intensity.Control <- median(Infection.Progression.DF.Control.Infected$LM.Slope.Intensity, na.rm=TRUE)
Avg.Progression.Intensity.Test <- median(Infection.Progression.DF.Test.Infected$LM.Slope.Intensity, na.rm=TRUE)


#Make Comparable Data to classic run by cross-section design_____________________________________________________________________May use this in different script as should not need these calculations besides validation.... But all the previous variables here are used
Infected.D1 <- Is.Infected[Is.Infected$First.Day.of.nonzero.data == "1",]
#Create Columns that can give quantifiable intensity
Infected.D1$Infect.Area.Norm <- NA
Infected.D1$Infect.Intensity.Norm <- NA



#Create DF for Control condition
Infected.D1.Control <- Infected.D1[Infected.D1$Labels == "Control Label",] #________Input Data______

#Create dataframe for each Control repeat
Infected.D1.Control.R1 <- Infected.D1.Control[Infected.D1.Control$Biological.Replicate == "Repeat 1",]
Infected.D1.Control.R2 <- Infected.D1.Control[Infected.D1.Control$Biological.Replicate == "Repeat 2",]
Infected.D1.Control.R3 <- Infected.D1.Control[Infected.D1.Control$Biological.Replicate == "Repeat 3",]

#Calculate number of infected animals seen on Day 1 in each repeat of control condition
Infected.D1.Control.R1.Num <- nrow(Infected.D1.Control.R1)
Infected.D1.Control.R2.Num <- nrow(Infected.D1.Control.R2)
Infected.D1.Control.R3.Num <- nrow(Infected.D1.Control.R3)

#Calculate total number of animals observed during each biological repeat within control condition
R1.Control.DF <- Control.DF[Control.DF$Biological.Replicate == "Repeat 1",]
R2.Control.DF <- Control.DF[Control.DF$Biological.Replicate == "Repeat 2",]
R3.Control.DF <- Control.DF[Control.DF$Biological.Replicate == "Repeat 3",]

#Calculate proportion of the population seen infected on Day 1 of control condition
Control.Prop.Infected.D1.R1 <- Infected.D1.Control.R1.Num/nrow(R1.Control.DF)
Control.Prop.Infected.D1.R2 <- Infected.D1.Control.R2.Num/nrow(R2.Control.DF)
Control.Prop.Infected.D1.R3 <- Infected.D1.Control.R3.Num/nrow(R3.Control.DF)

#Create list of infected proportion values for Control group
D1.Control.Infected.Prop.Repeats <- c(Control.Prop.Infected.D1.R1, Control.Prop.Infected.D1.R2, Control.Prop.Infected.D1.R3)

#Repeat all calculations with test condition

Infected.D1 <- Is.Infected[Is.Infected$First.Day.of.nonzero.data == "1",]
Infected.D1.Test <- Infected.D1[Infected.D1$Labels != "Control Label",] #________Input Data______


#Create dataframe for each Test repeat
Infected.D1.Test.R1 <- Infected.D1.Test[Infected.D1.Test$Biological.Replicate == "Repeat 1",]
Infected.D1.Test.R2 <- Infected.D1.Test[Infected.D1.Test$Biological.Replicate == "Repeat 2",]
Infected.D1.Test.R3 <- Infected.D1.Test[Infected.D1.Test$Biological.Replicate == "Repeat 3",]

#Calculate number of infected animals seen on Day 1 in each repeat of Test condition
Infected.D1.Test.R1.Num <- nrow(Infected.D1.Test.R1)
Infected.D1.Test.R2.Num <- nrow(Infected.D1.Test.R2)
Infected.D1.Test.R3.Num <- nrow(Infected.D1.Test.R3)

#Calculate total number of animals observed during each biological repeat within Test condition
R1.Test.DF <- Test.DF[Test.DF$Biological.Replicate == "Repeat 1",]
R2.Test.DF <- Test.DF[Test.DF$Biological.Replicate == "Repeat 2",]
R3.Test.DF <- Test.DF[Test.DF$Biological.Replicate == "Repeat 3",]

#Calculate proportion of the population seen infected on Day 1 of Test condition
Test.Prop.Infected.D1.R1 <- Infected.D1.Test.R1.Num/nrow(R1.Test.DF)
Test.Prop.Infected.D1.R2 <- Infected.D1.Test.R2.Num/nrow(R2.Test.DF)
Test.Prop.Infected.D1.R3 <- Infected.D1.Test.R3.Num/nrow(R3.Test.DF)

#Create list of infected proportion values for Test group
D1.Test.Infected.Prop.Repeats <- c(Test.Prop.Infected.D1.R1, Test.Prop.Infected.D1.R2, Test.Prop.Infected.D1.R3)

#Plot Infected Proportion by group (Figure 26)
boxplot(D1.Control.Infected.Prop.Repeats, D1.Test.Infected.Prop.Repeats,
        main = "Infected Population Day 1", #________Input Data______
        xlab = "Condition",
        ylab = "Fraction Infected",
        col = c("grey","red"),
        cex.lab = 1.2,
        cex.main = 1.2,
        names = c('Control Label',"Test Label")) #_____Input Data____

q <- ifelse(D1.Control.Infected.Prop.Repeats >= 0,1,D1.Control.Infected.Prop.Repeats)
m <- ifelse(is.na(q),0, q)

g <- ifelse(D1.Test.Infected.Prop.Repeats >= 0, 1, D1.Test.Infected.Prop.Repeats)
h <- ifelse(is.na(g),0, g)

if(sum(m) >= 3&
   sum(h) >= 3){
  Results.D1.Infected.Proportions <- t.test(D1.Control.Infected.Prop.Repeats,D1.Test.Infected.Prop.Repeats)
  
  #SICKO.Stat 29
  SICKO.Statistics[29,1] <- "Day 1 Proportion Populations Infected"
  SICKO.Statistics[29,2] <- "t-test"
  SICKO.Statistics[29,3] <- Results.D1.Infected.Proportions$p.value
  SICKO.Statistics[29,4] <- "t"
  SICKO.Statistics[29,5] <- Results.D1.Infected.Proportions$statistic
  SICKO.Statistics[29,6] <- "df"
  SICKO.Statistics[29,7] <- Results.D1.Infected.Proportions$parameter
  SICKO.Statistics[29,8] <- "Pg 26"
  SICKO.Statistics[29,9] <- Results.D1.Infected.Proportions$estimate[1]
  SICKO.Statistics[29,10] <- Results.D1.Infected.Proportions$estimate[2]
  SICKO.Statistics[29,11] <- Results.D1.Infected.Proportions$stderr

} else {
  #SICKO.Stat 29
  SICKO.Statistics[29,1] <- "Day 1 Proportion Populations Infected"
  SICKO.Statistics[29,2] <- "t-test"
  SICKO.Statistics[29,3] <- NA
  SICKO.Statistics[29,4] <- "t"
  SICKO.Statistics[29,5] <- NA
  SICKO.Statistics[29,6] <- "df"
  SICKO.Statistics[29,7] <- NA
  SICKO.Statistics[29,8] <- "Pg 26"
  SICKO.Statistics[29,9] <- NA
  SICKO.Statistics[29,10] <- NA
  SICKO.Statistics[29,11] <- NA
}


#Boxplot D1 infection Intensity by group  (Figure 27)                    #***Need George script for ~prism like graph
boxplot(Intensity_data1/10000000 ~ Labels, data = Infected.D1,
        main = "Infection Intensity Day 1",
        xlab = "Condition",
        ylab = "Infection Intensity (a.u.)",
        cex.lab = 1.2,
        cex.main = 1.2,
        col = c("grey", "red"))  


if(sum(m) >= 3&
   sum(h) >= 3) {
 Results.D1.Infection.Intensity <- t.test(Intensity_data1 ~ Labels, data = Infected.D1)

 #SICKO.Stat 30
 SICKO.Statistics[30,1] <- "Day 1 Average Infection Intensity by Group"
 SICKO.Statistics[30,2] <- "t-test"
 SICKO.Statistics[30,3] <- Results.D1.Infection.Intensity$p.value
 SICKO.Statistics[30,4] <- "t"
 SICKO.Statistics[30,5] <- Results.D1.Infection.Intensity$statistic
 SICKO.Statistics[30,6] <- "df"
 SICKO.Statistics[30,7] <- Results.D1.Infection.Intensity$parameter
 SICKO.Statistics[30,8] <- "Pg 27"
 SICKO.Statistics[30,9] <- Results.D1.Infection.Intensity$estimate[1]
 SICKO.Statistics[30,10] <- Results.D1.Infection.Intensity$estimate[2]
 SICKO.Statistics[30,11] <- Results.D1.Infection.Intensity$stderr
} else {
  #SICKO.Stat 30
  SICKO.Statistics[30,1] <- "Day 1 Average Infection Intensity by Group"
  SICKO.Statistics[30,2] <- "t-test"
  SICKO.Statistics[30,3] <- NA
  SICKO.Statistics[30,4] <- "t"
  SICKO.Statistics[30,5] <- NA
  SICKO.Statistics[30,6] <- "df"
  SICKO.Statistics[30,7] <- NA
  SICKO.Statistics[30,8] <- "Pg 27"
  SICKO.Statistics[30,9] <- NA
  SICKO.Statistics[30,10] <- NA
  SICKO.Statistics[30,11] <- NA
}

#Boxplot D1 Infection Area by Group (Figure 28)
boxplot(Area_data1/1000 ~ Labels, data = Infected.D1,
        main = "Infection Area Day 1",
        xlab = "Condition",
        ylab = "Infection Area (1K pixels)",
        cex.lab = 1.2,
        cex.main = 1.2,
        col = c("grey", "red"))

if(sum(m) >= 3&
   sum(h) >= 3){
 Results.D1.Infection.Area <- t.test(Area_data1 ~ Labels, data = Infected.D1)

 #SICKO.Stat 31
 SICKO.Statistics[31,1] <- "Day 1 Average Infection Area by Group"
 SICKO.Statistics[31,2] <- "t-test"
 SICKO.Statistics[31,3] <- Results.D1.Infection.Area$p.value
 SICKO.Statistics[31,4] <- "t"
 SICKO.Statistics[31,5] <- Results.D1.Infection.Area$statistic
 SICKO.Statistics[31,6] <- "df"
 SICKO.Statistics[31,7] <- Results.D1.Infection.Area$parameter
 SICKO.Statistics[31,8] <- "Pg 28"
 SICKO.Statistics[31,9] <- Results.D1.Infection.Area$estimate[1]
 SICKO.Statistics[31,10] <- Results.D1.Infection.Area$estimate[2]
 SICKO.Statistics[31,11] <- Results.D1.Infection.Area$stderr
} else {
  #SICKO.Stat 31
  SICKO.Statistics[31,1] <- "Day 1 Average Infection Area by Group"
  SICKO.Statistics[31,2] <- "t-test"
  SICKO.Statistics[31,3] <- NA
  SICKO.Statistics[31,4] <- "t"
  SICKO.Statistics[31,5] <- NA
  SICKO.Statistics[31,6] <- "df"
  SICKO.Statistics[31,7] <- NA
  SICKO.Statistics[31,8] <- "Pg 28"
  SICKO.Statistics[31,9] <- NA
  SICKO.Statistics[31,10] <- NA
  SICKO.Statistics[31,11] <- NA
}

#Graph intensity by how many days remaining in animal lifespan________Death Clock_____________________________________________________________________________________________________________

#Account for Control Animals that stayed
Control.Infected.DF.Stayed <- Control.Infected.DF[Control.Infected.DF$Is.Last.Day.Censored != 1,]

#Single out the Control intensity data
Control.Intensity.Data <- Control.Infected.DF.Stayed[c('Intensity_data1', 'Intensity_data2', 'Intensity_data3', 'Intensity_data4', 'Intensity_data5', 'Intensity_data6', 'Intensity_data7', 'Intensity_data8', 'Intensity_data9')]

#************Need to deal with NaN in DF
is.nan.data.frame <- function(x)
do.call(cbind, lapply(x, is.nan))
Control.Intensity.Data[is.nan(Control.Intensity.Data)] <- 0


#Create Control DF for Intensity measurements and length until death
Intensity.Clock.Control <- data.frame(matrix(NA, nrow = nrow(Control.DF)*9, ncol = 2))
names(Intensity.Clock.Control) <- c("Intensity_Value", "Remaining_Days")


#Fill Control Condition dataframe with forloop
j=0
for (i in 1:length(Control.Intensity.Data$Intensity_data1)){
  for (k in 1:ncol(Control.Intensity.Data)){
    j = j+1
    if (Control.Intensity.Data[i,k] > 0){
      Intensity.Clock.Control$Intensity_Value[j] <- Control.Intensity.Data[i,k]
      Intensity.Clock.Control$Remaining_Days[j] <- Control.Infected.DF.Stayed$Last.Day.of.Observation[i] - k
    }
  }
}

#Repeat for Test condition
Test.Infected.DF.Stayed <- Test.Infected.DF[Test.Infected.DF$Is.Last.Day.Censored != 1,]

#Single out the intensity data for Test Condition
Test.Intensity.Data <- Test.Infected.DF.Stayed[c('Intensity_data1', 'Intensity_data2', 'Intensity_data3', 'Intensity_data4', 'Intensity_data5', 'Intensity_data6', 'Intensity_data7', 'Intensity_data8', 'Intensity_data9')]

#************Need to deal with NaN in DF for Test Condtion
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
Test.Intensity.Data[is.nan(Test.Intensity.Data)] <- 0


#Create DF for Intensity measurements and length until death for Test Condition
Intensity.Clock.Test <- data.frame(matrix(NA, nrow = nrow(Test.DF)*9, ncol = 2))
names(Intensity.Clock.Test) <- c("Intensity_Value", "Remaining_Days")


#Fill Test Condition dataframe with forloop
j=0
for (i in 1:length(Test.Intensity.Data$Intensity_data1)){
  for (k in 1:ncol(Test.Intensity.Data)){
    j = j+1
    if (Test.Intensity.Data[i,k] > 0){
      Intensity.Clock.Test$Intensity_Value[j] <- Test.Intensity.Data[i,k]
      Intensity.Clock.Test$Remaining_Days[j] <- Test.Infected.DF.Stayed$Last.Day.of.Observation[i] - k
    }
  }
}

#Combine dataframes to attain maximum limits for the graphs
Max.Intensity.Value.Finder <- rbind(Intensity.Clock.Control, Intensity.Clock.Test)

#Obtain graph limits
V.Graph.Y.Max <- max(Max.Intensity.Value.Finder$Intensity_Value, na.rm = TRUE)
V.Graph.Y.Min <- min(Max.Intensity.Value.Finder$Intensity_Value, na.rm = TRUE)
VW.Graph.X.Max <- max(Max.Intensity.Value.Finder$Remaining_Days, na.rm = TRUE)

#Regression of Infection intensity vs due death in Control Condition (Figure 29)
Control.Intensity.Clock.LM <- lm(formula = Intensity_Value ~ Remaining_Days, data = Intensity.Clock.Control)
Death.Clock.Control.Intensity <- summary(Control.Intensity.Clock.LM)
plot(formula = Intensity_Value/10000000 ~ Remaining_Days, data = Intensity.Clock.Control,
     main = "Control Label Infection Intensity v Life Remaining",
     xlab = "Life Remaining(days)",
     ylab = "Infection Intensity (a.u.)",
     ylim = c(V.Graph.Y.Min/10000000,V.Graph.Y.Max/10000000),
     xlim = c(0,VW.Graph.X.Max),
     lwd = 2,
     cex = 2,
     cex.lab = 1.2,
     cex.main = 1.2)
abline(lm(formula = Intensity_Value/10000000 ~ Remaining_Days, data = Intensity.Clock.Control),
       col = "red",
       lty = "dashed",
       lwd = 4 )

#SICKO.Stat 32
SICKO.Statistics[32,1] <- "Control Infection Intensity v Time Until Death"
SICKO.Statistics[32,2] <- "Linear Regression"
SICKO.Statistics[32,3] <- Death.Clock.Control.Intensity$coefficients[2,4]
SICKO.Statistics[32,4] <- "t"
SICKO.Statistics[32,5] <- Death.Clock.Control.Intensity$coefficients[2,3]
SICKO.Statistics[32,6] <- "df"
SICKO.Statistics[32,7] <- Death.Clock.Control.Intensity$df[2]
SICKO.Statistics[32,8] <- "Pg 29"

#Regression of Infection Intensity vs due death In Test Condition (Figure 30)
Test.Intensity.Clock.LM <- lm(formula = Intensity_Value ~ Remaining_Days, data = Intensity.Clock.Test)
Death.Clock.Test.Intensity <- summary(Test.Intensity.Clock.LM)
plot(formula = Intensity_Value/10000000 ~ Remaining_Days, data = Intensity.Clock.Test,
     main = substitute(paste(bolditalic('Test Label'),bold(' Infection Intensity v Life Remaining'))),
     xlab = "Life Remaining(days)",
     ylab = "Infection Intensity (a.u.)",
#     ylim = c(V.Graph.Y.Min/10000000,V.Graph.Y.Max/10000000),    ### Come back
     xlim = c(0,VW.Graph.X.Max),
     lwd = 2,
     cex = 2,
     cex.lab = 1.2,
     cex.main = 1.2)
abline(lm(formula = Intensity_Value/10000000 ~ Remaining_Days, data = Intensity.Clock.Test),
       col = "red",
       lty = "dashed",
       lwd = 4)

#SICKO.Stat 33
SICKO.Statistics[33,1] <- "Test Infection Intensity v Time Until Death"
SICKO.Statistics[33,2] <- "Linear Regression"
SICKO.Statistics[33,3] <- Death.Clock.Test.Intensity$coefficients[2,4]
SICKO.Statistics[33,4] <- "t"
SICKO.Statistics[33,5] <- Death.Clock.Test.Intensity$coefficients[2,3]
SICKO.Statistics[33,6] <- "df"
SICKO.Statistics[33,7] <- Death.Clock.Test.Intensity$df[2]
SICKO.Statistics[33,8] <- "Pg 30"

#Repeat Death Clock for infection area

Control.Area.Data <- Control.Infected.DF.Stayed[c('Area_data1', 'Area_data2', 'Area_data3', 'Area_data4', 'Area_data5', 'Area_data6', 'Area_data7', 'Area_data8', 'Area_data9')]

#************Need to deal with NaN in DF Control
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
Control.Area.Data[is.nan(Control.Area.Data)] <- 0


#Create DF for Area measurements and length until death in Control
Area.Clock.Control <- data.frame(matrix(NA, nrow = nrow(Control.DF)*9, ncol = 2))
names(Area.Clock.Control) <- c("Area_Value", "Remaining_Days")


#Fill dataframe with forloop for Control
j=0
for (i in 1:length(Control.Area.Data$Area_data1)){
  for (k in 1:ncol(Control.Area.Data)){
    j = j+1
    if (Control.Area.Data[i,k] > 0){
      Area.Clock.Control$Area_Value[j] <- Control.Area.Data[i,k]
      Area.Clock.Control$Remaining_Days[j] <- Control.Infected.DF.Stayed$Last.Day.of.Observation[i] - k
    }
  }
}

#Repeat for Test condition

Test.Infected.DF.Stayed <- Test.Infected.DF[Test.Infected.DF$Is.Last.Day.Censored != 1,]

#Single out the Area data in Test
Test.Area.Data <- Test.Infected.DF.Stayed[c('Area_data1', 'Area_data2', 'Area_data3', 'Area_data4', 'Area_data5', 'Area_data6', 'Area_data7', 'Area_data8', 'Area_data9')]

#************Need to deal with NaN in DF in Test
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
Test.Area.Data[is.nan(Test.Area.Data)] <- 0


#Create DF for Area measurements and length until death for Test
Area.Clock.Test <- data.frame(matrix(NA, nrow = nrow(Test.DF)*9, ncol = 2))
names(Area.Clock.Test) <- c("Area_Value", "Remaining_Days")


#Fill dataframe with forloop for Test
j=0
for (i in 1:length(Test.Area.Data$Area_data1)){
  for (k in 1:ncol(Test.Area.Data)){
    j = j+1
    if (Test.Area.Data[i,k] > 0){
      Area.Clock.Test$Area_Value[j] <- Test.Area.Data[i,k]
      Area.Clock.Test$Remaining_Days[j] <- Test.Infected.DF.Stayed$Last.Day.of.Observation[i] - k
    }
  }
}

#Combine dataframes to attain maximum limits for the graphs
Max.Area.Value.Finder <- rbind(Area.Clock.Control, Area.Clock.Test)

#Find max limits for graph
W.Graph.Y.Max <- max(Max.Area.Value.Finder$Area_Value, na.rm = TRUE)
W.Graph.Y.Min <- min(Max.Area.Value.Finder$Area_Value, na.rm = TRUE)

#Calculate statistics using linear regression and graph Infection Area by life Remaining in Control (Figure 31)
Control.Area.Clock.LM <- lm(formula = Area_Value ~ Remaining_Days, data = Area.Clock.Control)
Death.Clock.Control.Area <- summary(Control.Area.Clock.LM)
plot(formula = Area_Value/1000 ~ Remaining_Days, data = Area.Clock.Control,
     main = "Control Label Infection Area v Life Remaining",
     xlab = "Life Remaining (days)",
     ylab = "Infection Area (1K pixels)",
     ylim = c(0,W.Graph.Y.Max/1000),
     xlim = c(0,VW.Graph.X.Max),
     lwd = 2,
     cex = 2,
     cex.lab = 1.2,
     cex.main = 1.2)
abline(lm(formula = Area_Value/1000 ~ Remaining_Days, data = Area.Clock.Control),
       col = "red",
       lty = "dashed",
       lwd = 4)

#SICKO.Stat 34
SICKO.Statistics[34,1] <- "Control Infection Area v Time Until Death"
SICKO.Statistics[34,2] <- "Linear Regression"
SICKO.Statistics[34,3] <- Death.Clock.Control.Area$coefficients[2,4]
SICKO.Statistics[34,4] <- "t"
SICKO.Statistics[34,5] <- Death.Clock.Control.Area$coefficients[2,3]
SICKO.Statistics[34,6] <- "df"
SICKO.Statistics[34,7] <- Death.Clock.Control.Area$df[2]
SICKO.Statistics[34,8] <- "Pg 31"

#Now perform regression on Area vs due death in Test Condition (Figure 32)
Test.Area.Clock.LM <- lm(formula = Area_Value ~ Remaining_Days, data = Area.Clock.Test)
Death.Clock.Test.Area <- summary(Test.Area.Clock.LM)
plot(formula = Area_Value/1000 ~ Remaining_Days, data = Area.Clock.Test,
     main = substitute(paste(bolditalic('Test Label'),bold(' Infection Area v Life Remaining'))),
     xlab = "Life Remaining (days)",
     ylab = "Infection Area (1K pixels)",
     # ylim = c(0,W.Graph.Y.Max/1000),
     xlim = c(0,VW.Graph.X.Max),
     lwd = 2,
     cex = 2,
     cex.lab = 1.2,
     cex.main = 1.2)
abline(lm(formula = Area_Value/1000 ~ Remaining_Days, data = Area.Clock.Test),
       col = "red",
       lty = "dashed",
       lwd = 4)

#SICKO.Stat 35
SICKO.Statistics[35,1] <- "Test Infection Area v Time Until Death"
SICKO.Statistics[35,2] <- "Linear Regression"
SICKO.Statistics[35,3] <- Death.Clock.Test.Area$coefficients[2,4]
SICKO.Statistics[35,4] <- "t"
SICKO.Statistics[35,5] <- Death.Clock.Test.Area$coefficients[2,3]
SICKO.Statistics[35,6] <- "df"
SICKO.Statistics[35,7] <-Death.Clock.Test.Area$df[2]
SICKO.Statistics[35,8] <- "Pg 32"

#Finish Graphing Material
dev.off()

write.csv(SICKO.Statistics, "SICKO_Statistics.csv", quote = FALSE, row.names = FALSE)
#still need to make table of statistical outputs


#Log Transform intensity data for better statistiscal testing
Intensity.Clock.Control$Intensity_Value <- log(Intensity.Clock.Control$Intensity_Value)
qqnorm(y=(log)(Intensity.Clock.Control$Intensity_Value),x=Intensity.Clock.Control$Remaining_Days, ylim = c(0,max((log)(Intensity.Clock.Control$Intensity_Value), na.rm = TRUE)))
Intensity.Clock.LM.Temp <- lm(formula = (log)(Intensity_Value) ~ Remaining_Days, data = Intensity.Clock.Control)
summary(Intensity.Clock.LM.Temp)
plot(formula = (log)(Intensity_Value) ~ Remaining_Days, data = Intensity.Clock.Control)
abline(lm(formula = (log)(Intensity_Value) ~ Remaining_Days, data = Intensity.Clock.Control))







