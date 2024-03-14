setwd("C:/Users/iberd/OneDrive/Desktop")
FullData = read.csv("FullData.csv")
library(tidyr)
library(dplyr)
library(ggplot2)

#DATAWRANGLING
#BASAL OCR
exp_conc_ocr <- subset(FullData, select = c('Dex_Exposure','Dex_Concentration', 'Basal_OCR.Av') )


#Control values at each measure of time
control_ocr <- subset(exp_conc_ocr, Dex_Concentration == "Control")
control_ocr$Dex_Exposure[control_ocr$Dex_Exposure == "15min"] <- .25
control_ocr$Dex_Exposure[control_ocr$Dex_Exposure == "30min"] <- .5
control_ocr$Dex_Exposure[control_ocr$Dex_Exposure == "1h"] <- 1
control_ocr$Dex_Exposure[control_ocr$Dex_Exposure == "2h"] <- 2
control_ocr$Dex_Exposure[control_ocr$Dex_Exposure == "4h"] <- 4
control_ocr$Dex_Exposure[control_ocr$Dex_Exposure == "8h"] <- 8


#Concentration values at each measure of time
concentration_ocr <- subset(exp_conc_ocr, Dex_Concentration != "Control")
concentration_ocr$Dex_Exposure[concentration_ocr$Dex_Exposure == "15min"] <- .25
concentration_ocr$Dex_Exposure[concentration_ocr$Dex_Exposure == "30min"] <- .5
concentration_ocr$Dex_Exposure[concentration_ocr$Dex_Exposure == "1h"] <- 1
concentration_ocr$Dex_Exposure[concentration_ocr$Dex_Exposure == "2h"] <- 2
concentration_ocr$Dex_Exposure[concentration_ocr$Dex_Exposure == "4h"] <- 4
concentration_ocr$Dex_Exposure[concentration_ocr$Dex_Exposure == "8h"] <- 8

#BASAL ECAR
exp_conc_ecar <- subset(FullData, select = c('Dex_Exposure','Dex_Concentration', 'Basal_ECAR.Av') )

control_ecar <- subset(exp_conc_ecar, Dex_Concentration == "Control")
concentration_ecar <- subset(exp_conc_ecar, Dex_Concentration != "Control")


#_____________________________________________________________________________________________
#Step 1: Test Statistic - AUC (median & plots)

#sort through data concentrations (10nM, 25nM, 50nM, 100nM, 500nM, 1uM, 10uM) &
#exposures (15min, 30min, 1h, 2h, 4h, 8h)



#calculate the median for control for each exposure
medians0.combos <- control_ocr %>%
  group_by(Dex_Exposure) %>%
  summarize(median_OCR = median(Basal_OCR.Av))


#calculate median for each concentration & exposure
medians1.combos <- concentration_ocr %>%
  group_by(Dex_Exposure, Dex_Concentration) %>%
  summarize(median_OCR = median(Basal_OCR.Av))


#PLOT - CONTROL
exposures <- c(.25, .5, 1, 2, 4, 8) 
concentrations <- c("10nm", "25nm", "50nm", "100nm", "500nm", "1um", "10um")

controlOCRvals <- c(medians0.combos$median_OCR)

conc10nm <- c(medians1.combos$median_OCR[medians1.combos$Dex_Concentration == "10nM"])
conc25nm <- c(medians1.combos$median_OCR[medians1.combos$Dex_Concentration == "25nM"])
conc50nm <- c(medians1.combos$median_OCR[medians1.combos$Dex_Concentration == "50nM"])
conc100nm <- c(medians1.combos$median_OCR[medians1.combos$Dex_Concentration == "100nM"])
conc500nm <- c(medians1.combos$median_OCR[medians1.combos$Dex_Concentration == "500nM"])
conc1um <- c(medians1.combos$median_OCR[medians1.combos$Dex_Concentration == "1uM"])
conc10um <- c(medians1.combos$median_OCR[medians1.combos$Dex_Concentration == "10uM"])


#control

plot(exposures, controlOCRvals, type = "o", ylim = range(controlOCRvals),
     xlim = range(exposures), main = "Control treatment plot",
     xlab = "Exposure", ylab = "OCR values")

#10nM
plot(exposures, conc10nm, type = "o", ylim = range(conc10nm),
     xlim = range(exposures), main = "10nM treatment plot",
     xlab = "Exposure", ylab = "OCR values")

plot(exposures, controlOCRvals, type = "o", ylim = range(c(controlOCRvals, conc100nm)),
     xlim = range(exposures), main = "Control vs 100nM treatment plot",
     xlab = "Exposure", ylab = "OCR values")

# Plot  100nM treatment on the same graph
lines(exposures, conc100nm, type = "o", col = "red")

polygon(c(exposures, rev(exposures)), c(controlOCRvals, rev(conc100nm)), col = "gray", border = NA)

#25nM
plot(exposures, conc25nm, type = "o", ylim = range(conc25nm),
     xlim = range(exposures), main = "25nM treatment plot",
     xlab = "Exposure", ylab = "OCR values")
#50nM
plot(exposures, conc50nm, type = "o", ylim = range(conc50nm),
     xlim = range(exposures), main = "50nM treatment plot",
     xlab = "Exposure", ylab = "OCR values")
#100nM
plot(exposures, conc100nm, type = "o", ylim = range(conc100nm),
     xlim = range(exposures), main = "100nM treatment plot",
     xlab = "Exposure", ylab = "OCR values")
#500nM
plot(exposures, conc500nm, type = "o", ylim = range(conc500nm),
     xlim = range(exposures), main = "500nM treatment plot",
     xlab = "Exposure", ylab = "OCR values")
#1uM
plot(exposures, conc1um, type = "o", ylim = range(conc1um),
     xlim = range(exposures), main = "10nM treatment plot",
     xlab = "Exposure", ylab = "OCR values")
#10uM
plot(exposures, conc10um, type = "o", ylim = range(conc10um),
     xlim = range(exposures), main = "10nM treatment plot",
     xlab = "Exposure", ylab = "OCR values")


#_________________________________________________________________________________


#difference in lines

#establish functions for each concentration 
line_control <- approxfun(exposures, controlOCRvals)
line_10nm <- approxfun(exposures, conc10nm)
line_25nm <- approxfun(exposures, conc25nm)
line_50nm <- approxfun(exposures, conc50nm)
line_100nm <- approxfun(exposures, conc100nm)
line_500nm <- approxfun(exposures, conc500nm)
line_1um <- approxfun(exposures, conc1um)
line_10um <- approxfun(exposures, conc10um)

x_values <- c(.25, .5, 1, 2, 4, 8)


#use in loop later
#difference in 10nM
difference10 <- function(x){ line_control(x) - line_10nm(x)}
area10 <- integrate(difference10, lower = min(x_values), upper = max(x_values))$value

#difference in 25nM
difference25 <- function(x){ line_control(x) - line_25nm(x)}
area25 <- integrate(difference10, lower = min(x_values), upper = max(x_values))$value

#difference in 50nM
difference50 <- function(x){ line_control(x) - line_50nm(x)}
area50 <- integrate(difference10, lower = min(x_values), upper = max(x_values))$value

#difference in 100nM
difference100 <- function(x){ line_control(x) - line_100nm(x)}
area100 <- integrate(difference10, lower = min(x_values), upper = max(x_values))$value

#difference in 500nM
difference500 <- function(x){ line_control(x) - line_500nm(x)}
area500 <- integrate(difference10, lower = min(x_values), upper = max(x_values))$value

#difference in 1uM
difference1um <- function(x){ line_control(x) - line_1um(x)}
area1um <- integrate(difference10, lower = min(x_values), upper = max(x_values))$value

#difference in 10uM
difference10um <- function(x){ line_control(x) - line_10um(x)}
area10um <- integrate(difference10, lower = min(x_values), upper = max(x_values))$value



conc_lines1 <- list(line_control, line_10nm, line_25nm, line_50nm, 
                   line_100nm, line_500nm, line_1um, line_10um)

areas1 <- numeric(length(conc_lines1) - 1)



# Loop through each concentration line except the control line
for (i in 1:length(conc_lines1)) {
 
  
  difference_values1 <- conc_lines1[[1]](x_values) - conc_lines1[[i]](x_values)
  

  difference_fun1 <- function(x) conc_lines1[[1]](x) - conc_lines1[[i]](x)
  
 
  areas1[i - 1] <- integrate(difference_fun1, lower = min(x_values), upper = max(x_values))$value
}


areas1
#[1] -25.69073 -27.60759 -35.13256 -36.69737 -39.99983 -47.28749
#[7] -45.44549


result_df <- data.frame(Difference_10nm = areas1[1],
                        Difference_25nm = areas1[2],
                        Difference_50nm = areas1[3],
                        Difference_100nm = areas1[4],
                        Difference_500nm = areas1[5],
                        Difference_1um = areas1[6],
                        Difference_10um = areas1[7])

View(result_df)
#_______________________________________________________________________________________

#Calculating Fisher with AUC

View(BinaryData)
BinaryData <- FullData
BinaryData$Concentration <- BinaryData$Dex_Concentration 
BinaryData$Dex_Concentration<-ifelse(BinaryData$Dex_Concentration=="Control",0,1)
BinaryData$Dex_Exposure[BinaryData$Dex_Exposure == "15min"] <- .25
BinaryData$Dex_Exposure[BinaryData$Dex_Exposure == "30min"] <- .5
BinaryData$Dex_Exposure[BinaryData$Dex_Exposure == "1h"] <- 1
BinaryData$Dex_Exposure[BinaryData$Dex_Exposure == "2h"] <- 2
BinaryData$Dex_Exposure[BinaryData$Dex_Exposure == "4h"] <- 4
BinaryData$Dex_Exposure[BinaryData$Dex_Exposure == "8h"] <- 8




#dataset for random  matrix to compare 10nM  vs control
num_units10nm <- sum(FullData$Dex_Concentration == "10nM") + sum(FullData$Dex_Concentration == "Control") #3680
num_treatments <- 2
num_W_rep <- choose(num_units10nm, num_treatments)
num_loop <- 10000
randommatrix <- matrix(nrow = num_units10nm, ncol = num_loop)
for (i in 1:num_loop) {
  randommatrix[, i] <- sample(BinaryData$Dex_Concentration[BinaryData$Concentration == "10nM"| 
                                                                   BinaryData$Concentration == "Control"])
}
randommatrix <- unique(randommatrix)
tstats10 <- rep(NA, length = length(randommatrix[1,]))



#dataset for random  matrix to compare 25nM  vs control
num_units25nm <- sum(FullData$Dex_Concentration == "25nM") + sum(FullData$Dex_Concentration == "Control") #3680
num_W_rep25 <- choose(num_units25nm, num_treatments)
randommatrix25 <- matrix(nrow = num_units25nm, ncol = num_loop)
for (i in 1:num_loop) {
  randommatrix25[, i] <- sample(BinaryData$Dex_Concentration[BinaryData$Concentration == "25nM"| 
                                                             BinaryData$Concentration == "Control"])
}
randommatrix25 <- unique(randommatrix25)
tstats25 <- rep(NA, length = length(randommatrix25[1,]))



#dataset for random  matrix to compare 50nM  vs control
num_units50nm <- sum(FullData$Dex_Concentration == "500nM") + sum(FullData$Dex_Concentration == "Control") #3680
num_W_rep50 <- choose(num_units50nm, num_treatments)
randommatrix50 <- matrix(nrow = num_units50nm, ncol = num_loop)
for (i in 1:num_loop) {
  randommatrix50[, i] <- sample(BinaryData$Dex_Concentration[BinaryData$Concentration == "50nM"| 
                                                             BinaryData$Concentration == "Control"])
}
randommatrix50 <- unique(randommatrix50)
tstats50 <- rep(NA, length = length(randommatrix50[1,]))



#dataset for random  matrix to compare 100nM  vs control
num_units100nm <- sum(FullData$Dex_Concentration == "100nM") + sum(FullData$Dex_Concentration == "Control") #3680
num_W_rep100 <- choose(num_units100nm, num_treatments)
randommatrix100 <- matrix(nrow = num_units100nm, ncol = num_loop)
for (i in 1:num_loop) {
  randommatrix100[, i] <- sample(BinaryData$Dex_Concentration[BinaryData$Concentration == "100nM"| 
                                                             BinaryData$Concentration == "Control"])
}
randommatrix100 <- unique(randommatrix100)
tstats100 <- rep(NA, length = length(randommatrix100[1,]))


#dataset for random  matrix to compare 500nM  vs control
num_units500nm <- sum(FullData$Dex_Concentration == "500nM") + sum(FullData$Dex_Concentration == "Control") #3680
num_W_rep500 <- choose(num_units500nm, num_treatments)
randommatrix500 <- matrix(nrow = num_units500nm, ncol = num_loop)
for (i in 1:num_loop) {
  randommatrix500[, i] <- sample(BinaryData$Dex_Concentration[BinaryData$Concentration == "500nM"| 
                                                             BinaryData$Concentration == "Control"])
}
randommatrix500 <- unique(randommatrix500)
tstats500 <- rep(NA, length = length(randommatrix500[1,]))



#dataset for random  matrix to compare 1uM  vs control
num_units1um <- sum(FullData$Dex_Concentration == "1uM") + sum(FullData$Dex_Concentration == "Control") #3680
num_W_rep1um <- choose(num_units1um, num_treatments)
randommatrix1um <- matrix(nrow = num_units1um, ncol = num_loop)
for (i in 1:num_loop) {
  randommatrix1um[, i] <- sample(BinaryData$Dex_Concentration[BinaryData$Concentration == "1uM"| 
                                                             BinaryData$Concentration == "Control"])
}
randommatrix1um <- unique(randommatrix1um)
tstats1um <- rep(NA, length = length(randommatrix1um[1,]))


#dataset for random  matrix to compare 10uM  vs control
num_units10um <- sum(FullData$Dex_Concentration == "10uM") + sum(FullData$Dex_Concentration == "Control") #3680
num_W_rep10um <- choose(num_units10um, num_treatments)
randommatrix10um <- matrix(nrow = num_units10um, ncol = num_loop)
for (i in 1:num_loop) {
  randommatrix10um[, i] <- sample(BinaryData$Dex_Concentration[BinaryData$Concentration == "10uM"| 
                                                                BinaryData$Concentration == "Control"])
}
randommatrix10um <- unique(randommatrix10um)
tstats10um <- rep(NA, length = length(randommatrix10um[1,]))




#______________________________________________________________________________________
View(BinaryTemp10)



# Create list to store test statistics for 10nM
for(i in 1:num_loop) { 
  BinaryTemp10 <- BinaryData[BinaryData$Concentration == "10nM" | BinaryData$Concentration == "Control", ]
  BinaryTemp10$Dex_Concentration <- randommatrix[, i]
  control0 <- BinaryTemp10[BinaryTemp10$Dex_Concentration == 0,]
  concentration10 <- BinaryTemp10[BinaryTemp10$Dex_Concentration == 1,]
  
  
  # Calculate median for each coherence within control
  coh25_control <- median(control0$Basal_OCR.Av[control0$Dex_Exposure == .25])
  coh50_control <- median(control0$Basal_OCR.Av[control0$Dex_Exposure == .50])
  coh1_control <- median(control0$Basal_OCR.Av[control0$Dex_Exposure == 1])
  coh2_control <- median(control0$Basal_OCR.Av[control0$Dex_Exposure == 2])
  coh4_control <- median(control0$Basal_OCR.Av[control0$Dex_Exposure == 4])
  coh8_control <- median(control0$Basal_OCR.Av[control0$Dex_Exposure == 8])
  
  
  # Calculate median for each coherence within 10nM condition
  coh25_conc <- median(concentration10$Basal_OCR.Av[concentration10$Dex_Exposure == .25])
  coh50_conc <- median(concentration10$Basal_OCR.Av[concentration10$Dex_Exposure == .50])
  coh1_conc <- median(concentration10$Basal_OCR.Av[concentration10$Dex_Exposure == 1])
  coh2_conc <- median(concentration10$Basal_OCR.Av[concentration10$Dex_Exposure == 2])
  coh4_conc <- median(concentration10$Basal_OCR.Av[concentration10$Dex_Exposure == 4])
  coh8_conc <- median(concentration10$Basal_OCR.Av[concentration10$Dex_Exposure == 8])
  
  line0points <- c(coh25_control, coh50_control, coh1_control, coh2_control,  coh4_control,  coh8_control)
  line1points <- c(coh25_conc, coh50_conc, coh1_conc, coh2_conc,  coh4_conc,  coh8_conc)
  
  line0 <- approxfun(as.numeric(x_values), line0points)
  line1 <- approxfun(as.numeric(x_values), line1points)
  
  difference <- function(x) line0(x) - line1(x)
  
  areas <- integrate(difference, lower = min(x_values), upper = max(x_values))$value
  
  tstats10[i]<-abs(areas)
  
  
}
TestStat10 = -25.60793
pval10 <- (sum(tstats10 <= abs(TestStat10)) / length(tstats10))

#25
for(i in 1:num_loop) { 
  BinaryTemp25 <- subset(BinaryData, Concentration == "25nM" | Concentration == "Control")
  BinaryTemp25$Dex_Concentration <- randommatrix25[, i]
  control0 <- BinaryTemp25[BinaryTemp25$Dex_Concentration == 0,]
  concentration25 <- BinaryTemp25[BinaryTemp25$Dex_Concentration == 1,]
  
  
  # Calculate median for each coherence within control
  coh25_control <- median(control0$Basal_OCR.Av[control0$Dex_Exposure == .25])
  coh50_control <- median(control0$Basal_OCR.Av[control0$Dex_Exposure == .50])
  coh1_control <- median(control0$Basal_OCR.Av[control0$Dex_Exposure == 1])
  coh2_control <- median(control0$Basal_OCR.Av[control0$Dex_Exposure == 2])
  coh4_control <- median(control0$Basal_OCR.Av[control0$Dex_Exposure == 4])
  coh8_control <- median(control0$Basal_OCR.Av[control0$Dex_Exposure == 8])
  
  
  # Calculate median for each coherence within 10nM condition
  coh25_conc25 <- median(concentration25$Basal_OCR.Av[concentration25$Dex_Exposure == .25])
  coh50_conc25 <- median(concentration25$Basal_OCR.Av[concentration25$Dex_Exposure == .50])
  coh1_conc25 <- median(concentration25$Basal_OCR.Av[concentration25$Dex_Exposure == 1])
  coh2_conc25 <- median(concentration25$Basal_OCR.Av[concentration25$Dex_Exposure == 2])
  coh4_conc25 <- median(concentration25$Basal_OCR.Av[concentration25$Dex_Exposure == 4])
  coh8_conc25 <- median(concentration25$Basal_OCR.Av[concentration25$Dex_Exposure == 8])
  
  line0points <- c(coh25_control, coh50_control, coh1_control, coh2_control,  coh4_control,  coh8_control)
  line1points25 <- c(coh25_conc25, coh50_conc25, coh1_conc25, coh2_conc25,  coh4_conc25,  coh8_conc25)
  
  line0 <- approxfun(as.numeric(x_values), line0points)
  line1.25 <- approxfun(as.numeric(x_values), line1points25)
  
  difference25 <- function(x) line0(x) - line1.25(x)
  
  areas25 <- integrate(difference25, lower = min(x_values), upper = max(x_values))$value
  
  tstats25[i]<-abs(areas25)
  
  
}
TestStat25 = -27.60759
pval25 <- (sum(tstats25 <= abs(TestStat25)) / length(tstats25))

#50nM
for(i in 1:num_loop) { 
  BinaryTemp50 <- subset(BinaryData, Concentration == "50nM" | Concentration == "Control")
  BinaryTemp50$Dex_Concentration <- randommatrix50[, i]
  control0 <- BinaryTemp50[BinaryTemp50$Dex_Concentration == 0,]
  concentration50 <- BinaryTemp50[BinaryTemp50$Dex_Concentration == 1,]
  
  
  # Calculate median for each coherence within control
  coh25_control <- median(control0$Basal_OCR.Av[control0$Dex_Exposure == .25])
  coh50_control <- median(control0$Basal_OCR.Av[control0$Dex_Exposure == .50])
  coh1_control <- median(control0$Basal_OCR.Av[control0$Dex_Exposure == 1])
  coh2_control <- median(control0$Basal_OCR.Av[control0$Dex_Exposure == 2])
  coh4_control <- median(control0$Basal_OCR.Av[control0$Dex_Exposure == 4])
  coh8_control <- median(control0$Basal_OCR.Av[control0$Dex_Exposure == 8])
  
  
  # Calculate median for each coherence within 10nM condition
  coh25_conc50 <- median(concentration50$Basal_OCR.Av[concentration50$Dex_Exposure == .25])
  coh50_conc50 <- median(concentration50$Basal_OCR.Av[concentration50$Dex_Exposure == .50])
  coh1_conc50 <- median(concentration50$Basal_OCR.Av[concentration50$Dex_Exposure == 1])
  coh2_conc50 <- median(concentration50$Basal_OCR.Av[concentration50$Dex_Exposure == 2])
  coh4_conc50 <- median(concentration50$Basal_OCR.Av[concentration50$Dex_Exposure == 4])
  coh8_conc50 <- median(concentration50$Basal_OCR.Av[concentration50$Dex_Exposure == 8])
  
  line0points <- c(coh25_control, coh50_control, coh1_control, coh2_control,  coh4_control,  coh8_control)
  line1points50 <- c(coh25_conc50, coh50_conc50, coh1_conc50, coh2_conc50,  coh4_conc50,  coh8_conc50)
  
  line0 <- approxfun(as.numeric(x_values), line0points)
  line1.50 <- approxfun(as.numeric(x_values), line1points50)
  
  difference50 <- function(x) line0(x) - line1.50(x)
  
  areas50 <- integrate(difference50, lower = min(x_values), upper = max(x_values))$value
  
  tstats50[i]<-abs(areas50)
  
  
}
TestStat50 = -35.13256
pval50 <- (sum(tstats50 <= abs(TestStat50)) / length(tstats50))

#100nM
for(i in 1:num_loop) { 
  BinaryTemp100 <- subset(BinaryData, Concentration == "100nM" | Concentration == "Control")
  BinaryTemp100$Dex_Concentration <- randommatrix100[, i]
  control0 <- BinaryTemp100[BinaryTemp100$Dex_Concentration == 0,]
  concentration100 <- BinaryTemp100[BinaryTemp100$Dex_Concentration == 1,]
  
  
  # Calculate median for each coherence within control
  coh25_control <- median(control0$Basal_OCR.Av[control0$Dex_Exposure == .25])
  coh50_control <- median(control0$Basal_OCR.Av[control0$Dex_Exposure == .50])
  coh1_control <- median(control0$Basal_OCR.Av[control0$Dex_Exposure == 1])
  coh2_control <- median(control0$Basal_OCR.Av[control0$Dex_Exposure == 2])
  coh4_control <- median(control0$Basal_OCR.Av[control0$Dex_Exposure == 4])
  coh8_control <- median(control0$Basal_OCR.Av[control0$Dex_Exposure == 8])
  
  
  # Calculate median for each coherence within 10nM condition
  coh25_conc100 <- median(concentration100$Basal_OCR.Av[concentration100$Dex_Exposure == .25])
  coh50_conc100 <- median(concentration100$Basal_OCR.Av[concentration100$Dex_Exposure == .50])
  coh1_conc100 <- median(concentration100$Basal_OCR.Av[concentration100$Dex_Exposure == 1])
  coh2_conc100 <- median(concentration100$Basal_OCR.Av[concentration100$Dex_Exposure == 2])
  coh4_conc100 <- median(concentration100$Basal_OCR.Av[concentration100$Dex_Exposure == 4])
  coh8_conc100 <- median(concentration100$Basal_OCR.Av[concentration100$Dex_Exposure == 8])
  
  line0points <- c(coh25_control, coh50_control, coh1_control, coh2_control,  coh4_control,  coh8_control)
  line1points100 <- c(coh25_conc100, coh50_conc100, coh1_conc100, coh2_conc100,  coh4_conc100,  coh8_conc100)
  
  line0 <- approxfun(as.numeric(x_values), line0points)
  line1.100 <- approxfun(as.numeric(x_values), line1points100)
  
  difference100 <- function(x) line0(x) - line1.100(x)
  
  areas100 <- integrate(difference100, lower = min(x_values), upper = max(x_values))$value
  
  tstats100[i]<-abs(areas100)
  
  
}
TestStat100 = 25.60793
pval100 <- (sum(tstats100 >= (TestStat100)) / length(tstats100))


# Create a histogram of the test statistics
ggplot(data.frame(Test_Stats = tstats100), aes(x = Test_Stats)) +
  geom_histogram(bins = 20, fill = "skyblue", color = "black") +
  xlim(-30, 30) +
  labs(x = "Test Statistic", y = "Frequency", title = "Histogram of Test Statistics") +
  geom_vline(xintercept = TestStat100, linetype = "dashed", color = "red") +
  geom_text(x = TestStat100 + 1, y = max(hist(tstats100)$counts) * 0.9,
            label = paste("Observed Test Statistic =", round(TestStat100, 1)), color = "red") +
  geom_text(x = -29, y = max(hist(tstats100)$counts) * 0.8,
            label = paste("p-value:", round(pval100, 4)), color = "red") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())



#500uM
for(i in 1:num_loop) { 
  BinaryTemp500 <- subset(BinaryData, Concentration == "500nM" | Concentration == "Control")
  BinaryTemp500$Dex_Concentration <- randommatrix500[, i]
  control0 <- BinaryTemp500[BinaryTemp10$Dex_Concentration == 0,]
  concentration500 <- BinaryTemp500[BinaryTemp500$Dex_Concentration == 1,]
  
  
  # Calculate median for each exposure within control
  coh25_control <- median(control0$Basal_OCR.Av[control0$Dex_Exposure == .25])
  coh50_control <- median(control0$Basal_OCR.Av[control0$Dex_Exposure == .50])
  coh1_control <- median(control0$Basal_OCR.Av[control0$Dex_Exposure == 1])
  coh2_control <- median(control0$Basal_OCR.Av[control0$Dex_Exposure == 2])
  coh4_control <- median(control0$Basal_OCR.Av[control0$Dex_Exposure == 4])
  coh8_control <- median(control0$Basal_OCR.Av[control0$Dex_Exposure == 8])
  
  
  # Calculate median for each exposure within 10nM condition
  coh25_conc500 <- median(concentration500$Basal_OCR.Av[concentration500$Dex_Exposure == .25])
  coh50_conc500 <- median(concentration500$Basal_OCR.Av[concentration500$Dex_Exposure == .50])
  coh1_conc500 <- median(concentration500$Basal_OCR.Av[concentration500$Dex_Exposure == 1])
  coh2_conc500 <- median(concentration500$Basal_OCR.Av[concentration500$Dex_Exposure == 2])
  coh4_conc500 <- median(concentration500$Basal_OCR.Av[concentration500$Dex_Exposure == 4])
  coh8_conc500 <- median(concentration500$Basal_OCR.Av[concentration500$Dex_Exposure == 8])
  
  line0points <- c(coh25_control, coh50_control, coh1_control, coh2_control,  coh4_control,  coh8_control)
  line1points500 <- c(coh25_conc500, coh50_conc500, coh1_conc500, coh2_conc500,  coh4_conc500,  coh8_conc500)
  
  line0 <- approxfun(as.numeric(x_values), line0points)
  line1.500 <- approxfun(as.numeric(x_values), line1points500)
  
  difference500 <- function(x) line0(x) - line1.500(x)
  
  areas500 <- integrate(difference500, lower = min(x_values), upper = max(x_values))$value
  
  tstats500[i]<-(areas500)
  
  
}
TestStat500 = -39.99983
pval500 <- (sum(tstats500 >= (TestStat500)) / length(tstats500))

#1uM
for(i in 1:num_loop) { 
  BinaryTemp1um <- subset(BinaryData, Concentration == "1uM" | Concentration == "Control")
  BinaryTemp1um$Dex_Concentration <- randommatrix1um[, i]
  control0 <- BinaryTemp1um[BinaryTemp1um$Dex_Concentration == 0,]
  concentration1um <- BinaryTemp1um[BinaryTemp1um$Dex_Concentration == 1,]

  
  # Calculate median for each exposure within control
  coh25_control <- median(control0$Basal_OCR.Av[control0$Dex_Exposure == .25])
  coh50_control <- median(control0$Basal_OCR.Av[control0$Dex_Exposure == .50])
  coh1_control <- median(control0$Basal_OCR.Av[control0$Dex_Exposure == 1])
  coh2_control <- median(control0$Basal_OCR.Av[control0$Dex_Exposure == 2])
  coh4_control <- median(control0$Basal_OCR.Av[control0$Dex_Exposure == 4])
  coh8_control <- median(control0$Basal_OCR.Av[control0$Dex_Exposure == 8])
  
  
  # Calculate median for each exposure within 10nM condition
  coh25_conc1um <- median(concentration1um$Basal_OCR.Av[concentration1um$Dex_Exposure == .25])
  coh50_conc1um <- median(concentration1um$Basal_OCR.Av[concentration1um$Dex_Exposure == .50])
  coh1_conc1um <- median(concentration1um$Basal_OCR.Av[concentration1um$Dex_Exposure == 1])
  coh2_conc1um <- median(concentration1um$Basal_OCR.Av[concentration1um$Dex_Exposure == 2])
  coh4_conc1um <- median(concentration1um$Basal_OCR.Av[concentration1um$Dex_Exposure == 4])
  coh8_conc1um <- median(concentration1um$Basal_OCR.Av[concentration1um$Dex_Exposure == 8])
  
  line0points <- c(coh25_control, coh50_control, coh1_control, coh2_control,  coh4_control,  coh8_control)
  line1points1um <- c(coh25_conc1um, coh50_conc1um, coh1_conc1um, coh2_conc1um,  coh4_conc1um,  coh8_conc1um)
  
  line0 <- approxfun(as.numeric(x_values), line0points)
  line1.1um <- approxfun(as.numeric(x_values), line1points1um)
  
  difference1um <- function(x) line0(x) - line1.1um(x)
  
  areas1um <- integrate(difference1um, lower = min(x_values), upper = max(x_values))$value
  
  tstats1um[i]<-(areas1um)
  
  
}
TestStat1um = -47.28749
pval1um <- (sum(tstats1um >= (TestStat1um)) / length(tstats1um))


#10uM
for(i in 1:num_loop) { 
  BinaryTemp10um <- subset(BinaryData, Concentration == "10uM" | Concentration == "Control")
  BinaryTemp10um$Dex_Concentration <- randommatrix10um[, i]
  control0 <- BinaryTemp10um[BinaryTemp10um$Dex_Concentration == 0,]
  concentration10um <- BinaryTemp10um[BinaryTemp10um$Dex_Concentration == 1,]
  
  
  # Calculate median for each exposure within control
  coh25_control <- median(control0$Basal_OCR.Av[control0$Dex_Exposure == .25])
  coh50_control <- median(control0$Basal_OCR.Av[control0$Dex_Exposure == .50])
  coh1_control <- median(control0$Basal_OCR.Av[control0$Dex_Exposure == 1])
  coh2_control <- median(control0$Basal_OCR.Av[control0$Dex_Exposure == 2])
  coh4_control <- median(control0$Basal_OCR.Av[control0$Dex_Exposure == 4])
  coh8_control <- median(control0$Basal_OCR.Av[control0$Dex_Exposure == 8])
  
  
  # Calculate median for each exposure within 10nM condition
  coh25_conc10um <- median(concentration10um$Basal_OCR.Av[concentration10um$Dex_Exposure == .25])
  coh50_conc10um <- median(concentration10um$Basal_OCR.Av[concentration10um$Dex_Exposure == .50])
  coh1_conc10um <- median(concentration10um$Basal_OCR.Av[concentration10um$Dex_Exposure == 1])
  coh2_conc10um <- median(concentration10um$Basal_OCR.Av[concentration10um$Dex_Exposure == 2])
  coh4_conc10um <- median(concentration10um$Basal_OCR.Av[concentration10um$Dex_Exposure == 4])
  coh8_conc10um <- median(concentration10um$Basal_OCR.Av[concentration10um$Dex_Exposure == 8])
  
  line0points <- c(coh25_control, coh50_control, coh1_control, coh2_control,  coh4_control,  coh8_control)
  line1points10um <- c(coh25_conc10um, coh50_conc10um, coh1_conc10um, coh2_conc10um,  coh4_conc10um,  coh8_conc10um)
  
  line0 <- approxfun(as.numeric(x_values), line0points)
  line1.10um <- approxfun(as.numeric(x_values), line1points10um)
  
  difference10um <- function(x) line0(x) - line1.10um(x)
  
  areas10um <- integrate(difference10um, lower = min(x_values), upper = max(x_values))$value
  
  tstats10um[i]<- (areas10um)
  
  
}
TestStat10um = -45.44549
pval10um <- (sum(tstats10um >= (TestStat10um)) / length(tstats10um))

