setwd("C:/Users/iberd/OneDrive/Desktop")
CompleteData = read.csv("FullData.csv")
View(CompleteData)
library(tidyverse)
library(tidyr)
library(dplyr)

#BASAL OCR, 10nM Condition at each exposure

#data wrangling to create the test stat table
ControlValues <-subset(CompleteData, Dex_Concentration == "Control")
View(ControlValues)
Concentration10data <- subset(CompleteData, CompleteData$Dex_Concentration == "10nM")
ConcentrationDataAll <- subset(CompleteData, Dex_Concentration != "Control")
BinaryData.2 <- CompleteData
BinaryData.2$Concentration <- BinaryData.2$Dex_Concentration 
BinaryData.2$Concentration<-ifelse(BinaryData.2$Concentration=="Control",0,1) 

mean(ConcentrationDataAll$Basal_OCR.Av[ConcentrationDataAll$Dex_Exposure =="8h"|
                                    ConcentrationDataAll$Dex_Concentration == "100nM"])
mean(ControlValues$Basal_OCR.Av[ControlValues$Dex_Exposure == "8h"])

mean(ConcentrationDataAll$Basal_OCR.Av[ConcentrationDataAll$Dex_Exposure =="8h"|
                                         ConcentrationDataAll$Dex_Concentration == "100nM"]) - mean(ControlValues$Basal_OCR.Av[ControlValues$Dex_Exposure == "8h"])

#Step 1: Compute Difference in Means

#mean for control at each exposure
means0_combos <- ControlValues %>%
  group_by(Dex_Exposure) %>%
  summarize(mean_OCR = mean(Basal_OCR.Av))

#mean for each concentration at each exposure
means_combos <- ConcentrationDataAll %>%
  group_by(Dex_Exposure, Dex_Concentration) %>%
  summarize(mean_OCR = mean(Basal_OCR.Av))

View(means0_combos)

#VISUALIZATION
ggplot(means_combos, aes(x = Dex_Exposure, y = Dex_Concentration, color = mean_OCR)) +
  geom_point() +
  geom_text(aes(label = sprintf("%.3f", mean_OCR)), hjust = 0, vjust = 0) + 
  labs(x = "Dex Exposure", y = "Dex Concentration", color = "Value Column") +
  ggtitle("Mean OCR for each Exp. + Conc. Combo")
ggplot(means0_combos, aes(x = Dex_Exposure, y = mean_OCR)) +
  geom_point() +  # Adding points
  geom_text(aes(label = sprintf("%.3f", mean_OCR)), hjust = 0, vjust = 0) +  # Adding text annotations
  labs(x = "Dex Exposure", y = "Control Value") +
  ggtitle("Mean OCR for each Exp. + Control Combo")

exposure_levels <- c("10nM", "25nM", "50nM", "100nM", "500nM", "1uM", "10uM")  # List of concentration levels

# Loop through each concentration level and calculate the differences

#15 minutes
DiffM.values.15min <- list()  # Create an empty list to store the differences
for (level in exposure_levels) {
  diff_name <- paste("DiffM.15min.", level, sep = "")
  
  # Calculate the difference
  DiffM.values.15min[[diff_name]] <- means_combos$mean_OCR[means_combos$Dex_Exposure == "15min" &
                                                                 means_combos$Dex_Concentration == level]-means0_combos$mean_OCR[means0_combos$Dex_Exposure == "15min"]
}

#30 minutes
DiffM.values.30min <- list()
for (level in exposure_levels) {
  diff_name <- paste("DiffM.30min.", level, sep = "")
  DiffM.values.30min[[diff_name]] <- means_combos$mean_OCR[means_combos$Dex_Exposure == "30min" &
                                                             means_combos$Dex_Concentration == level]-means0_combos$mean_OCR[means0_combos$Dex_Exposure == "30min"]
  
}

#1 hour
DiffM.values.1hr <- list()
for (level in exposure_levels) {
  diff_name <- paste("DiffM.1hr.", level, sep = "")
  DiffM.values.1hr[[diff_name]] <-     means_combos$mean_OCR[means_combos$Dex_Exposure == "1h" &
                                                               means_combos$Dex_Concentration == level]-means0_combos$mean_OCR[means0_combos$Dex_Exposure == "1h"] 

}

#2 hours 
DiffM.values.2hr <- list()
for (level in exposure_levels) {
  diff_name <- paste("DiffM.2hr.", level, sep = "")
  DiffM.values.2hr[[diff_name]] <-     means_combos$mean_OCR[means_combos$Dex_Exposure == "2h" &
                                                               means_combos$Dex_Concentration == level] - 
    means0_combos$mean_OCR[means0_combos$Dex_Exposure == "2h"] 

}

#4 hours
DiffM.values.4hr <- list()
for (level in exposure_levels) {
  diff_name <- paste("DiffM.4hr.", level, sep = "")
  DiffM.values.4hr[[diff_name]] <-     means_combos$mean_OCR[means_combos$Dex_Exposure == "4h" &
                                                               means_combos$Dex_Concentration == level] - means0_combos$mean_OCR[means0_combos$Dex_Exposure == "4h"] 

}

#8 hours
DiffM.values.8hr <- list()
for (level in exposure_levels) {
  diff_name <- paste("DiffM.8hr.", level, sep = "")
  DiffM.values.8hr[[diff_name]] <-     means_combos$mean_OCR[means_combos$Dex_Exposure == "8h" &
                                                               means_combos$Dex_Concentration == level]-means0_combos$mean_OCR[means0_combos$Dex_Exposure == "8h"] 

}

#_______________________________________________________________________________________

#creat matrix for each combination of concentration + exposure

#Creating Matrix 15min
Matriz10nm.15min{
MatrixData15min <- subset(BinaryData.2,  Dex_Exposure == "15min")
MatrixData15min.10nm <- subset(MatrixData15min, Dex_Concentration == "10nM" |
                               Dex_Concentration == "Control")
nrow(MatrixData15min.10nm)

num15.10 <- 520

num_W_rep15min.10nm <- choose(num15.10, num_treatments)
matrix_15min.10nm <- matrix(nrow = num15.10, ncol = num_loop)

for (k in 1:num_loop) {
  matrix_15min.10nm[, k] <- sample(MatrixData15min.10nm$Concentration)
}

num_unique_allocations10.15 <- ncol(unique(matrix_15min.10nm, MARGIN = 2))


T.rep10.15 <- rep(NA,length = num_unique_allocations10.15)
for(k in 1:num_unique_allocations10.15){
  
  T.rep10.15[k] <- mean(MatrixData15min.10nm$Basal_OCR.Av[matrix_15min.10nm[, k]==1]) - 
    mean(MatrixData15min.10nm$Basal_OCR.Av[matrix_15min.10nm[, k]==0])
}

k=1



#-----------------------------------------------------------------

#find pval
p.val10.15 <- sum(T.rep10.15 <= 0.06820438)/length(T.rep10.15)

}
Matriz25nm.15min{
  MatrixData15min <- subset(BinaryData.2,  Dex_Exposure == "15min")
  MatrixData15min.25nm <- subset(MatrixData15min, Dex_Concentration == "25nM" |
                                   Dex_Concentration == "Control")
  nrow(MatrixData15min.25nm)
  
  num15.25 <- 520
  
  num_W_rep15min.25nm <- choose(num15.25, num_treatments)
  matrix_15min.25nm <- matrix(nrow = num15.25, ncol = num_loop)
  
  for (k in 1:num_loop) {
    matrix_15min.25nm[, k] <- sample(MatrixData15min.25nm$Concentration)
  }
  
  num_unique_allocations25.15 <- ncol(unique(matrix_15min.25nm, MARGIN = 2))
  
  
  T.rep25.15 <- rep(NA,length = num_unique_allocations25.15)
  for(k in 1:num_unique_allocations25.15){
    
    T.rep25.15[k] <- mean(MatrixData15min.25nm$Basal_OCR.Av[matrix_15min.25nm[, k]==1]) - 
      mean(MatrixData15min.25nm$Basal_OCR.Av[matrix_15min.25nm[, k]==0])
  }

  #find pval
  p.val25.15 <- sum(T.rep25.15 <= -1.825281)/length(T.rep25.15)
 
}
Matriz50nm.15min{
  MatrixData15min <- subset(BinaryData.2,  Dex_Exposure == "15min")
  MatrixData15min.50nm <- subset(MatrixData15min, Dex_Concentration == "50nM" |
                                   Dex_Concentration == "Control")
  nrow(MatrixData15min.50nm)
  num15.50 <- 520
  
  num_W_rep15min.50nm <- choose(num15.50, num_treatments)
  matrix_15min.50nm <- matrix(nrow = num15.50, ncol = num_loop)
  
  for (k in 1:num_loop) {
    matrix_15min.50nm[, k] <- sample(MatrixData15min.50nm$Concentration)
  }
  
  num_unique_allocations50.15 <- ncol(unique(matrix_15min.50nm, MARGIN = 2))
  
  
  T.rep50.15 <- rep(NA,length = num_unique_allocations50.15)
  for(k in 1:num_unique_allocations50.15){
    
    T.rep50.15[k] <- mean(MatrixData15min.50nm$Basal_OCR.Av[matrix_15min.50nm[, k]==1]) - 
      mean(MatrixData15min.50nm$Basal_OCR.Av[matrix_15min.50nm[, k]==0])
  }
  
  #find pval
  p.val50.15 <- sum(T.rep50.15 <= -1.425102)/length(T.rep50.15)
  
}
Matriznm100.15min{
    MatrixData15min <- subset(BinaryData.2,  Dex_Exposure == "15min")
    MatrixData15min.100nm <- subset(MatrixData15min, Dex_Concentration == "100nM" |
                                     Dex_Concentration == "Control")
    nrow(MatrixData15min.100nm)
    
    num15.100 <- 520
    
    num_W_rep15min.100nm <- choose(num15.100, num_treatments)
    matrix_15min.100nm <- matrix(nrow = num15.100, ncol = num_loop)
    
    for (k in 1:num_loop) {
      matrix_15min.100nm[, k] <- sample(MatrixData15min.100nm$Concentration)
    }
    
    num_unique_allocations100.15 <- ncol(unique(matrix_15min.100nm, MARGIN = 2))
    
    
    T.rep100.15 <- rep(NA,length = num_unique_allocations100.15)
    for(k in 1:num_unique_allocations100.15){
      
      T.rep100.15[k] <- mean(MatrixData15min.100nm$Basal_OCR.Av[matrix_15min.100nm[, k]==1]) - 
        mean(MatrixData15min.100nm$Basal_OCR.Av[matrix_15min.100nm[, k]==0])
    }
    
    
    
    #-----------------------------------------------------------------
    
    #find pval
    p.val100.15 <- sum(T.rep100.15 <= -2.311091)/length(T.rep100.15)
 

}
Matrix500nm.15min{
  MatrixData15min <- subset(BinaryData.2,  Dex_Exposure == "15min")
  MatrixData15min.500nm <- subset(MatrixData15min, Dex_Concentration == "500nM" |
                                    Dex_Concentration == "Control")
  nrow(MatrixData15min.500nm)
  
  num15.500 <- 520
  
  num_W_rep15min.500nm <- choose(num15.500, num_treatments)
  matrix_15min.500nm <- matrix(nrow = num15.500, ncol = num_loop)
  
  for (k in 1:num_loop) {
    matrix_15min.500nm[, k] <- sample(MatrixData15min.500nm$Concentration)
  }
  
  num_unique_allocations500.15 <- ncol(unique(matrix_15min.500nm, MARGIN = 2))
  
  
  T.rep500.15 <- rep(NA,length = num_unique_allocations500.15)
  for(k in 1:num_unique_allocations500.15){
    
    T.rep500.15[k] <- mean(MatrixData15min.500nm$Basal_OCR.Av[matrix_15min.500nm[, k]==1]) - 
      mean(MatrixData15min.500nm$Basal_OCR.Av[matrix_15min.500nm[, k]==0])
  }
  

  p.val500.15 <- sum(T.rep500.15 <=  -0.2562645)/length(T.rep500.15)
 
  
}
Matriz1um.15min{
  MatrixData15min <- subset(BinaryData.2,  Dex_Exposure == "15min")
  MatrixData15min.1um <- subset(MatrixData15min, Dex_Concentration == "1uM" |
                                   Dex_Concentration == "Control")
  nrow(MatrixData15min.1um)
  
  num15.1um <- 520
  
  num_W_rep15min.1um <- choose(num15.1um, num_treatments)
  matrix_15min.1um <- matrix(nrow = num15.1um, ncol = num_loop)
  
  for (k in 1:num_loop) {
    matrix_15min.1um[, k] <- sample(MatrixData15min.1um$Concentration)
  }
  
  num_unique_allocations1um.15 <- ncol(unique(matrix_15min.1um, MARGIN = 2))
  
  
  T.rep1um.15 <- rep(NA,length = num_unique_allocations1um.15)
  for(k in 1:num_unique_allocations1um.15){
    
    T.rep1um.15[k] <- mean(MatrixData15min.1um$Basal_OCR.Av[matrix_15min.1um[, k]==1]) - 
      mean(MatrixData15min.1um$Basal_OCR.Av[matrix_15min.1um[, k]==0])
  }
  

  #find pval
  p.val1um.15 <- sum(T.rep1um.15 <= 1.032506)/length(T.rep1um.15)

}
Matriz10um.15min{
  MatrixData15min <- subset(BinaryData.2,  Dex_Exposure == "15min")
  MatrixData15min.10um <- subset(MatrixData15min, Dex_Concentration == "10uM" |
                                   Dex_Concentration == "Control")
  nrow(MatrixData15min.10um)
  
  num15.10um <- 520
  
  num_W_rep15min.10um <- choose(num15.10um, num_treatments)
  matrix_15min.10um <- matrix(nrow = num15.10um, ncol = num_loop)
  
  for (k in 1:num_loop) {
    matrix_15min.10um[, k] <- sample(MatrixData15min.10um$Concentration)
  }
  
  num_unique_allocations10um.15 <- ncol(unique(matrix_15min.10um, MARGIN = 2))
  
  
  T.rep10um.15 <- rep(NA,length = num_unique_allocations10um.15)
  for(k in 1:num_unique_allocations10um.15){
    
    T.rep10um.15[k] <- mean(MatrixData15min.10um$Basal_OCR.Av[matrix_15min.10um[, k]==1]) - 
      mean(MatrixData15min.10um$Basal_OCR.Av[matrix_15min.10um[, k]==0])
  }

  #find pval
  p.val10um.15 <- sum(T.rep10um.15 <= 0.08295954)/length(T.rep10um.15)

}

#Creating Matrix 30min
Matriz10nm.30min{
  MatrixData30min <- subset(BinaryData.2,  Dex_Exposure == "30min")
  MatrixData30min.10nm <- subset(MatrixData15min, Dex_Concentration == "10nM" |
                                   Dex_Concentration == "Control")
  nrow(MatrixData30min.10nm)
  
  num30.10 <- 520
  
  num_W_rep30min.10nm <- choose(num30.10, num_treatments)
  matrix_30min.10nm <- matrix(nrow = num30.10, ncol = num_loop)
  
  for (k in 1:num_loop) {
    matrix_30min.10nm[, k] <- sample(MatrixData15min.10nm$Concentration)
  }
  
  num_unique_allocations30.15 <- ncol(unique(matrix_30min.10nm, MARGIN = 2))
  
  
  T.rep10.30 <- rep(NA,length = num_unique_allocations30.15)
  for(k in 1:num_unique_allocations30.15){
    
    T.rep10.30[k] <- mean(MatrixData30min.10nm$Basal_OCR.Av[matrix_30min.10nm[, k]==1]) - 
      mean(MatrixData15min.10nm$Basal_OCR.Av[matrix_30min.10nm[, k]==0])
  }
 
  
  #find pval
  p.val10.30 <- sum(T.rep10.30 <=  2.139632)/length(T.rep10.30)

}
Matriz25nm.30min{
  MatrixData30min <- subset(BinaryData.2,  Dex_Exposure == "30min")
  MatrixData30min.25nm <- subset(MatrixData30min, Dex_Concentration == "25nM" |
                                   Dex_Concentration == "Control")
  nrow(MatrixData30min.25nm)
  
  num30.25 <- 528
  
  num_W_rep30min.25nm <- choose(num30.25, num_treatments)
  matrix_30min.25nm <- matrix(nrow = num30.25, ncol = num_loop)
  
  for (k in 1:num_loop) {
    matrix_30min.25nm[, k] <- sample(MatrixData30min.25nm$Concentration)
  }
  
  num_unique_allocations25.30 <- ncol(unique(matrix_30min.25nm, MARGIN = 2))
  
  
  T.rep25.30 <- rep(NA,length = num_unique_allocations25.30)
  for(k in 1:num_unique_allocations25.30){
    
    T.rep25.30[k] <- mean(MatrixData30min.25nm$Basal_OCR.Av[matrix_30min.25nm[, k]==1]) - 
      mean(MatrixData30min.25nm$Basal_OCR.Av[matrix_30min.25nm[, k]==0])
  }
  
  #find pval
  p.val25.30 <- sum(T.rep25.30 <=  1.432302)/length(T.rep25.30)
  
}
Matriz50nm.30min{
  MatrixData30min <- subset(BinaryData.2,  Dex_Exposure == "30min")
  MatrixData30min.50nm <- subset(MatrixData30min, Dex_Concentration == "50nM" |
                                   Dex_Concentration == "Control")
  nrow(MatrixData30min.50nm)
  num30.50 <- 528
  
  num_W_rep30min.50nm <- choose(num30.50, num_treatments)
  matrix_30min.50nm <- matrix(nrow = num30.50, ncol = num_loop)
  
  for (k in 1:num_loop) {
    matrix_30min.50nm[, k] <- sample(MatrixData30min.50nm$Concentration)
  }
  
  num_unique_allocations50.30 <- ncol(unique(matrix_30min.50nm, MARGIN = 2))
  
  
  T.rep50.30 <- rep(NA,length = num_unique_allocations50.30)
  for(k in 1:num_unique_allocations50.30){
    
    T.rep50.30[k] <- mean(MatrixData30min.50nm$Basal_OCR.Av[matrix_30min.50nm[, k]==1]) - 
      mean(MatrixData30min.50nm$Basal_OCR.Av[matrix_30min.50nm[, k]==0])
  }
  
  #find pval
  p.val50.30 <- sum(T.rep50.30 <=  0.171706)/length(T.rep50.30)
  
}
Matriznm100.30min{
  MatrixData30min <- subset(BinaryData.2,  Dex_Exposure == "30min")
  MatrixData30min.100nm <- subset(MatrixData30min, Dex_Concentration == "100nM" |
                                    Dex_Concentration == "Control")
  nrow(MatrixData30min.100nm)
  
  num30.100 <- 528
  
  num_W_rep30min.100nm <- choose(num15.100, num_treatments)
  matrix_30min.100nm <- matrix(nrow = num30.100, ncol = num_loop)
  
  for (k in 1:num_loop) {
    matrix_30min.100nm[, k] <- sample(MatrixData30min.100nm$Concentration)
  }
  
  num_unique_allocations100.30 <- ncol(unique(matrix_30min.100nm, MARGIN = 2))
  
  
  T.rep100.30 <- rep(NA,length = num_unique_allocations100.30)
  for(k in 1:num_unique_allocations100.30){
    
    T.rep100.30[k] <- mean(MatrixData30min.100nm$Basal_OCR.Av[matrix_30min.100nm[, k]==1]) - 
      mean(MatrixData30min.100nm$Basal_OCR.Av[matrix_30min.100nm[, k]==0])
  }

  
  #find pval
  p.val100.30 <- sum(T.rep100.30 <= -0.03059165)/length(T.rep100.30)
  
  
}
Matrix500nm.30min{
  MatrixData30min <- subset(BinaryData.2,  Dex_Exposure == "30min")
  MatrixData30min.500nm <- subset(MatrixData30min, Dex_Concentration == "500nM" |
                                    Dex_Concentration == "Control")
  nrow(MatrixData30min.500nm)
  
  num30.500 <- 528
  
  num_W_rep30min.500nm <- choose(num30.500, num_treatments)
  matrix_30min.500nm <- matrix(nrow = num30.500, ncol = num_loop)
  
  for (k in 1:num_loop) {
    matrix_30min.500nm[, k] <- sample(MatrixData30min.500nm$Concentration)
  }
  
  num_unique_allocations500.30 <- ncol(unique(matrix_30min.500nm, MARGIN = 2))
  
  
  T.rep500.30 <- rep(NA,length = num_unique_allocations500.30)
  for(k in 1:num_unique_allocations500.30){
    
    T.rep500.30[k] <- mean(MatrixData30min.500nm$Basal_OCR.Av[matrix_30min.500nm[, k]==1]) - 
      mean(MatrixData30min.500nm$Basal_OCR.Av[matrix_30min.500nm[, k]==0])
  }
  
  
  p.val500.30<- sum(T.rep500.30 <=  0.7246484)/length(T.rep500.30)
  
  
}
Matriz1um.30min{
  MatrixData30min <- subset(BinaryData.2,  Dex_Exposure == "30min")
  MatrixData30min.1um <- subset(MatrixData30min, Dex_Concentration == "1uM" |
                                  Dex_Concentration == "Control")
  nrow(MatrixData30min.1um)
  
  num30.1um <- 528
  
  num_W_rep30min.1um <- choose(num30.1um, num_treatments)
  matrix_30min.1um <- matrix(nrow = num30.1um, ncol = num_loop)
  
  for (k in 1:num_loop) {
    matrix_30min.1um[, k] <- sample(MatrixData30min.1um$Concentration)
  }
  
  num_unique_allocations1um.30 <- ncol(unique(matrix_30min.1um, MARGIN = 2))
  
  
  T.rep1um.30 <- rep(NA,length = num_unique_allocations1um.30)
  for(k in 1:num_unique_allocations1um.30){
    
    T.rep1um.30[k] <- mean(MatrixData30min.1um$Basal_OCR.Av[matrix_30min.1um[, k]==1]) - 
      mean(MatrixData30min.1um$Basal_OCR.Av[matrix_30min.1um[, k]==0])
  }
  
  
  #find pval
  p.val1um.30 <- sum(T.rep1um.30<= 1.639357)/length(T.rep1um.30)
  
}
Matriz10um.30min{
  MatrixData30min <- subset(BinaryData.2,  Dex_Exposure == "30min")
  MatrixData30min.10um <- subset(MatrixData30min, Dex_Concentration == "10uM" |
                                   Dex_Concentration == "Control")
  nrow(MatrixData30min.10um)
  
  num30min.10um <- 528
  
  num_W_rep30min.10um <- choose(num30min.10um, num_treatments)
  matrix_30min.10um <- matrix(nrow = num30min.10um, ncol = num_loop)
  
  for (k in 1:num_loop) {
    matrix_30min.10um[, k] <- sample(MatrixData30min.10um$Concentration)
  }
  
  num_unique_allocations10um.30min<- ncol(unique(matrix_30min.10um, MARGIN = 2))
  
  
  T.rep10um.30min <- rep(NA,length = num_unique_allocations10um.30min)
  for(k in 1:num_unique_allocations10um.30min){
    
    T.rep10um.30min[k] <- mean(MatrixData30min.10um$Basal_OCR.Av[matrix_30min.10um[, k]==1]) - 
      mean(MatrixData30min.10um$Basal_OCR.Av[matrix_30min.10um[, k]==0])
  }
  
  #find pval
  p.val10um.30min <- sum(T.rep10um.30min<= 3.664758)/length(T.rep10um.30min)
  
}

#Creating Matrix 1h
Matriz10nm.1h{
  MatrixData1h <- subset(BinaryData.2,  Dex_Exposure == "1h")
  MatrixData1h.10nm <- subset(MatrixData15min, Dex_Concentration == "10nM" |
                                   Dex_Concentration == "Control")
  nrow(MatrixData1h.10nm)
  
  num1h.10 <- 520
  
  num_W_rep1h.10nm <- choose(num1h.10, num_treatments)
  matrix_1h.10nm <- matrix(nrow = num1h.10, ncol = num_loop)
  
  for (k in 1:num_loop) {
    matrix_1h.10nm[, k] <- sample(MatrixData1h.10nm$Concentration)
  }
  
  num_unique_allocations1h.10 <- ncol(unique(matrix_1h.10nm, MARGIN = 2))
  
  
  T.rep10.1h <- rep(NA,length = num_unique_allocations1h.10)
  for(k in 1:num_unique_allocations1h.10){
    
    T.rep10.1h[k] <- mean(MatrixData1h.10nm$Basal_OCR.Av[matrix_1h.10nm[, k]==1]) - 
      mean(MatrixData1h.10nm$Basal_OCR.Av[matrix_1h.10nm[, k]==0])
  }
  
  
  #find pval
  p.val10.1h <- sum(T.rep10.1h <= -1.195562)/length(T.rep10.1h)
  
}
Matriz25nm.1h{
  MatrixData1h<- subset(BinaryData.2,  Dex_Exposure == "1h")
  MatrixData1h.25nm <- subset(MatrixData30min, Dex_Concentration == "25nM" |
                                   Dex_Concentration == "Control")
  nrow(MatrixData1h.25nm)
  
  num1h.25 <- 528
  
  num_W_rep1h.25nm <- choose(num1h.25, num_treatments)
  matrix_1h.25nm <- matrix(nrow = num1h.25, ncol = num_loop)
  
  for (k in 1:num_loop) {
    matrix_1h.25nm[, k] <- sample(MatrixData1h.25nm$Concentration)
  }
  
  num_unique_allocations25.1h <- ncol(unique(matrix_1h.25nm, MARGIN = 2))
  
  
  T.rep25.1h<- rep(NA,length = num_unique_allocations25.1h)
  for(k in 1:num_unique_allocations25.1h){
    
    T.rep25.1h[k] <- mean(MatrixData1h.25nm$Basal_OCR.Av[matrix_1h.25nm[, k]==1]) - 
      mean(MatrixData1h.25nm$Basal_OCR.Av[matrix_1h.25nm[, k]==0])
  }
  
  #find pval
  p.val25.1h <- sum(T.rep25.1h <=  -1.699022)/length(T.rep25.1h)
  
}
Matriz50nm.1h{
  MatrixData1h <- subset(BinaryData.2,  Dex_Exposure == "1h")
  MatrixData1h.50nm <- subset(MatrixData1h, Dex_Concentration == "50nM" |
                                   Dex_Concentration == "Control")
  nrow(MatrixData1h.50nm)
  num1h.50 <- 548
  
  num_W_rep1h.50nm <- choose(num30.50, num_treatments)
  matrix_1h.50nm <- matrix(nrow = num1h.50, ncol = num_loop)
  
  for (k in 1:num_loop) {
    matrix_1h.50nm[, k] <- sample(MatrixData1h.50nm$Concentration)
  }
  
  num_unique_allocations50.1h <- ncol(unique(matrix_1h.50nm, MARGIN = 2))
  
  
  T.rep50.1h<- rep(NA,length = num_unique_allocations50.1h)
  for(k in 1:num_unique_allocations50.1h){
    
    T.rep50.1h[k] <- mean(MatrixData1h.50nm$Basal_OCR.Av[matrix_1h.50nm[, k]==1]) - 
      mean(MatrixData1h.50nm$Basal_OCR.Av[matrix_1h.50nm[, k]==0])
  }
  
  #find pval
  p.val50.1h<- sum(T.rep50.1h <=-1.966834)/length(T.rep50.1h)
  
}
Matriznm100.1h{
  MatrixData1h <- subset(BinaryData.2,  Dex_Exposure == "1h")
  MatrixData1h.100nm <- subset(MatrixData1h, Dex_Concentration == "100nM" |
                                    Dex_Concentration == "Control")
  nrow(MatrixData1h.100nm)
  
  num1h.100 <- 548
  
  num_W_rep1h.100nm <- choose(num1h.100, num_treatments)
  matrix_1h.100nm <- matrix(nrow = num1h.100, ncol = num_loop)
  
  for (k in 1:num_loop) {
    matrix_1h.100nm[, k] <- sample(MatrixData1h.100nm$Concentration)
  }
  
  num_unique_allocations100.1h <- ncol(unique(matrix_1h.100nm, MARGIN = 2))
  
  
  T.rep100.1h <- rep(NA,length = num_unique_allocations100.1h)
  for(k in 1:num_unique_allocations100.1h){
    
    T.rep100.1h[k] <- mean(MatrixData1h.100nm$Basal_OCR.Av[matrix_1h.100nm[, k]==1]) - 
      mean(MatrixData1h.100nm$Basal_OCR.Av[matrix_1h.100nm[, k]==0])
  }
  
  
  #find pval
  p.val100.1h <- sum(T.rep100.1h <= -1.1492)/length(T.rep100.1h)
  
  
}
Matrix500nm.1h{
  MatrixData1h <- subset(BinaryData.2,  Dex_Exposure == "1h")
  MatrixData1h.500nm <- subset(MatrixData1h, Dex_Concentration == "500nM" |
                                    Dex_Concentration == "Control")
  nrow(MatrixData1h.500nm)
  
  num1h.500 <- 548
  
  num_W_rep1h.500nm <- choose(num1h.500, num_treatments)
  matrix_1h.500nm <- matrix(nrow = num1h.500, ncol = num_loop)
  
  for (k in 1:num_loop) {
    matrix_1h.500nm[, k] <- sample(MatrixData1h.500nm$Concentration)
  }
  
  num_unique_allocations500.1h <- ncol(unique(matrix_1h.500nm, MARGIN = 2))
  
  
  T.rep500.1h <- rep(NA,length = num_unique_allocations500.1h)
  for(k in 1:num_unique_allocations500.1h){
    
    T.rep500.1h[k] <- mean(MatrixData1h.500nm$Basal_OCR.Av[matrix_1h.500nm[, k]==1]) - 
      mean(MatrixData1h.500nm$Basal_OCR.Av[matrix_1h.500nm[, k]==0])
  }
  
  
  p.val500.1h<- sum(T.rep500.1h <=  -0.7352085)/length(T.rep500.1h)
  
  
}
Matriz1um.1h{
  MatrixData1h <- subset(BinaryData.2,  Dex_Exposure == "1h")
  MatrixData1h.1um <- subset(MatrixData1h, Dex_Concentration == "1uM" |
                                  Dex_Concentration == "Control")
  nrow(MatrixData1h.1um)
  
  num1h.1um <- 548
  
  num_W_rep1h.1um <- choose(num1h.1um, num_treatments)
  matrix_1h.1um <- matrix(nrow = num1h.1um, ncol = num_loop)
  
  for (k in 1:num_loop) {
    matrix_1h.1um[, k] <- sample(MatrixData1h.1um$Concentration)
  }
  
  num_unique_allocations1um.1h <- ncol(unique(matrix_1h.1um, MARGIN = 2))
  
  
  T.rep1um.1h <- rep(NA,length = num_unique_allocations1um.1h)
  for(k in 1:num_unique_allocations1um.1h){
    
    T.rep1um.1h[k] <- mean(MatrixData1h.1um$Basal_OCR.Av[matrix_1h.1um[, k]==1]) - 
      mean(MatrixData1h.1um$Basal_OCR.Av[matrix_1h.1um[, k]==0])
  }
  
  
  #find pval
  p.val1um.1h <- sum(T.rep1um.1h<=  -0.106165)/length(T.rep1um.1h)
  
}
Matriz10um.1h{
  MatrixData1h <- subset(BinaryData.2,  Dex_Exposure == "1h")
  MatrixData1h.10um <- subset(MatrixData15min, Dex_Concentration == "10uM" |
                                   Dex_Concentration == "Control")
  nrow(MatrixData1h.10um)
  
  num1h.10um <- 520
  
  num_W_rep1h.10um <- choose(num1h.10um, num_treatments)
  matrix_1h.10um <- matrix(nrow = num1h.10um, ncol = num_loop)
  
  for (k in 1:num_loop) {
    matrix_1h.10um[, k] <- sample(MatrixData1h.10um$Concentration)
  }
  
  num_unique_allocations10um.1h<- ncol(unique(matrix_1h.10um, MARGIN = 2))
  
  
  T.rep10um.1h<- rep(NA,length = num_unique_allocations10um.1h)
  for(k in 1:num_unique_allocations10um.1h){
    
    T.rep10um.1h[k] <- mean(MatrixData1h.10um$Basal_OCR.Av[matrix_1h.10um[, k]==1]) - 
      mean(MatrixData1h.10um$Basal_OCR.Av[matrix_1h.10um[, k]==0])
  }
  
  #find pval
  p.val10um.1h <- sum(T.rep10um.1h<=  1.670323)/length(T.rep10um.1h)
  
}

#Creating Matrix 2h
Matriz10nm.2h{
  MatrixData2h <- subset(BinaryData.2,  Dex_Exposure == "2h")
  MatrixData2h.10nm <- subset(MatrixData2h, Dex_Concentration == "10nM" |
                                Dex_Concentration == "Control")
  nrow(MatrixData2h.10nm)
  
  num2h.10 <- 588
  
  num_W_rep2h.10nm <- choose(num2h.10, num_treatments)
  matrix_2h.10nm <- matrix(nrow = num2h.10, ncol = num_loop)
  
  for (k in 1:num_loop) {
    matrix_2h.10nm[, k] <- sample(MatrixData2h.10nm$Concentration)
  }
  
  num_unique_allocations2h.10 <- ncol(unique(matrix_2h.10nm, MARGIN = 2))
  
  
  T.rep10.2h <- rep(NA,length = num_unique_allocations2h.10)
  for(k in 1:num_unique_allocations2h.10){
    
    T.rep10.2h[k] <- mean(MatrixData2h.10nm$Basal_OCR.Av[matrix_2h.10nm[, k]==1]) - 
      mean(MatrixData2h.10nm$Basal_OCR.Av[matrix_2h.10nm[, k]==0])
  }
  
  
  #find pval
  p.val10.2h <- sum(T.rep10.2h <= 1.629559)/length(T.rep10.2h)
  
}
Matriz25nm.2h{
  MatrixData2h<- subset(BinaryData.2,  Dex_Exposure == "2h")
  MatrixData2h.25nm <- subset(MatrixData2h, Dex_Concentration == "25nM" |
                                Dex_Concentration == "Control")
  nrow(MatrixData2h.25nm)
  
  num2h.25 <- 588
  
  num_W_rep2h.25nm <- choose(num2h.25, num_treatments)
  matrix_2h.25nm <- matrix(nrow = num2h.25, ncol = num_loop)
  
  for (k in 1:num_loop) {
    matrix_2h.25nm[, k] <- sample(MatrixData2h.25nm$Concentration)
  }
  
  num_unique_allocations25.2h <- ncol(unique(matrix_2h.25nm, MARGIN = 2))
  
  
  T.rep25.2h<- rep(NA,length = num_unique_allocations25.2h)
  for(k in 1:num_unique_allocations25.2h){
    
    T.rep25.2h[k] <- mean(MatrixData2h.25nm$Basal_OCR.Av[matrix_2h.25nm[, k]==1]) - 
      mean(MatrixData2h.25nm$Basal_OCR.Av[matrix_2h.25nm[, k]==0])
  }
  
  #find pval
  p.val25.2h <- sum(T.rep25.2h <=  1.118905)/length(T.rep25.2h)
  
}
Matriz50nm.2h{
  MatrixData2h <- subset(BinaryData.2,  Dex_Exposure == "2h")
  MatrixData2h.50nm <- subset(MatrixData2h, Dex_Concentration == "50nM" |
                                Dex_Concentration == "Control")
  nrow(MatrixData2h.50nm)
  num2h.50 <- 588
  
  num_W_rep1h.50nm <- choose(num2h.50, num_treatments)
  matrix_2h.50nm <- matrix(nrow = num2h.50, ncol = num_loop)
  
  for (k in 1:num_loop) {
    matrix_2h.50nm[, k] <- sample(MatrixData2h.50nm$Concentration)
  }
  
  num_unique_allocations50.2h <- ncol(unique(matrix_2h.50nm, MARGIN = 2))
  
  
  T.rep50.2h<- rep(NA,length = num_unique_allocations50.2h)
  for(k in 1:num_unique_allocations50.2h){
    
    T.rep50.2h[k] <- mean(MatrixData2h.50nm$Basal_OCR.Av[matrix_2h.50nm[, k]==1]) - 
      mean(MatrixData2h.50nm$Basal_OCR.Av[matrix_2h.50nm[, k]==0])
  }
  
  #find pval
  p.val50.2h<- sum(T.rep50.2h <=1.605882)/length(T.rep50.2h)
  
}
Matriznm100.2h{
  MatrixData2h <- subset(BinaryData.2,  Dex_Exposure == "2h")
  MatrixData2h.100nm <- subset(MatrixData2h, Dex_Concentration == "100nM" |
                                 Dex_Concentration == "Control")
  nrow(MatrixData2h.100nm)
  
  num2h.100 <- 588
  
  num_W_rep2h.100nm <- choose(num2h.100, num_treatments)
  matrix_2h.100nm <- matrix(nrow = num2h.100, ncol = num_loop)
  
  for (k in 1:num_loop) {
    matrix_2h.100nm[, k] <- sample(MatrixData2h.100nm$Concentration)
  }
  
  num_unique_allocations100.2h <- ncol(unique(matrix_2h.100nm, MARGIN = 2))
  
  
  T.rep100.2h <- rep(NA,length = num_unique_allocations100.2h)
  for(k in 1:num_unique_allocations100.2h){
    
    T.rep100.2h[k] <- mean(MatrixData2h.100nm$Basal_OCR.Av[matrix_2h.100nm[, k]==1]) - 
      mean(MatrixData2h.100nm$Basal_OCR.Av[matrix_2h.100nm[, k]==0])
  }
  
  
  #find pval
  p.val100.2h <- sum(T.rep100.2h <= 1.019617)/length(T.rep100.2h)
  
  
}
Matrix500nm.2h{
  MatrixData2h <- subset(BinaryData.2,  Dex_Exposure == "2h")
  MatrixData2h.500nm <- subset(MatrixData2h, Dex_Concentration == "500nM" |
                                 Dex_Concentration == "Control")
  nrow(MatrixData2h.500nm)
  
  num2h.500 <- 588
  
  num_W_rep2h.500nm <- choose(num2h.500, num_treatments)
  matrix_2h.500nm <- matrix(nrow = num2h.500, ncol = num_loop)
  
  for (k in 1:num_loop) {
    matrix_2h.500nm[, k] <- sample(MatrixData2h.500nm$Concentration)
  }
  
  num_unique_allocations500.2h <- ncol(unique(matrix_2h.500nm, MARGIN = 2))
  
  
  T.rep500.2h <- rep(NA,length = num_unique_allocations500.2h)
  for(k in 1:num_unique_allocations500.2h){
    
    T.rep500.2h[k] <- mean(MatrixData2h.500nm$Basal_OCR.Av[matrix_2h.500nm[, k]==1]) - 
      mean(MatrixData2h.500nm$Basal_OCR.Av[matrix_2h.500nm[, k]==0])
  }
  
  
  p.val500.2h<- sum(T.rep500.2h <=  2.664589)/length(T.rep500.2h)
  
  
}
Matriz1um.2h{
  MatrixData2h <- subset(BinaryData.2,  Dex_Exposure == "2h")
  MatrixData2h.1um <- subset(MatrixData1h, Dex_Concentration == "1uM" |
                               Dex_Concentration == "Control")
  nrow(MatrixData2h.1um)
  
  num2h.1um <- 548
  
  num_W_rep2h.1um <- choose(num2h.1um, num_treatments)
  matrix_2h.1um <- matrix(nrow = num2h.1um, ncol = num_loop)
  
  for (k in 1:num_loop) {
    matrix_2h.1um[, k] <- sample(MatrixData2h.1um$Concentration)
  }
  
  num_unique_allocations1um.2h <- ncol(unique(matrix_2h.1um, MARGIN = 2))
  
  
  T.rep1um.2h <- rep(NA,length = num_unique_allocations1um.2h)
  for(k in 1:num_unique_allocations1um.2h){
    
    T.rep1um.2h[k] <- mean(MatrixData2h.1um$Basal_OCR.Av[matrix_2h.1um[, k]==1]) - 
      mean(MatrixData2h.1um$Basal_OCR.Av[matrix_2h.1um[, k]==0])
  }
  
  
  #find pval
  p.val1um.2h <- sum(T.rep1um.2h<=  3.358541)/length(T.rep1um.2h)
  
}
Matriz10um.2h{
  MatrixData2h <- subset(BinaryData.2,  Dex_Exposure == "2h")
  MatrixData2h.10um <- subset(MatrixData15min, Dex_Concentration == "10uM" |
                                Dex_Concentration == "Control")
  nrow(MatrixData2h.10um)
  
  num2h.10um <- 520
  
  num_W_rep2h.10um <- choose(num2h.10um, num_treatments)
  matrix_2h.10um <- matrix(nrow = num2h.10um, ncol = num_loop)
  
  for (k in 1:num_loop) {
    matrix_2h.10um[, k] <- sample(MatrixData2h.10um$Concentration)
  }
  
  num_unique_allocations10um.2h<- ncol(unique(matrix_2h.10um, MARGIN = 2))
  
  
  T.rep10um.2h<- rep(NA,length = num_unique_allocations10um.2h)
  for(k in 1:num_unique_allocations10um.2h){
    
    T.rep10um.2h[k] <- mean(MatrixData2h.10um$Basal_OCR.Av[matrix_2h.10um[, k]==1]) - 
      mean(MatrixData2h.10um$Basal_OCR.Av[matrix_2h.10um[, k]==0])
  }
  
  #find pval
  p.val10um.2h <- sum(T.rep10um.2h<=  2.425812)/length(T.rep10um.2h)
  
}

#Creating Matrix 4h
Matriz10nm.4h{
  MatrixData4h <- subset(BinaryData.2,  Dex_Exposure == "4h")
  MatrixData4h.10nm <- subset(MatrixData4h, Dex_Concentration == "10nM" |
                                Dex_Concentration == "Control")
  nrow(MatrixData4h.10nm)
  
  num4h.10 <- 668
  
  num_W_rep4h.10nm <- choose(num4h.10, num_treatments)
  matrix_4h.10nm <- matrix(nrow = num4h.10, ncol = num_loop)
  
  for (k in 1:num_loop) {
    matrix_4h.10nm[, k] <- sample(MatrixData4h.10nm$Concentration)
  }
  
  num_unique_allocations4h.10 <- ncol(unique(matrix_4h.10nm, MARGIN = 2))
  
  
  T.rep10.4h <- rep(NA,length = num_unique_allocations4h.10)
  for(k in 1:num_unique_allocations4h.10){
    
    T.rep10.4h[k] <- mean(MatrixData4h.10nm$Basal_OCR.Av[matrix_4h.10nm[, k]==1]) - 
      mean(MatrixData4h.10nm$Basal_OCR.Av[matrix_4h.10nm[, k]==0])
  }
  
  
  #find pval
  p.val10.4h <- sum(T.rep10.4h <= 1.821774)/length(T.rep10.4h)
  
}
Matriz25nm.4h{
  MatrixData4h<- subset(BinaryData.2,  Dex_Exposure == "4h")
  MatrixData4h.25nm <- subset(MatrixData4h, Dex_Concentration == "25nM" |
                                Dex_Concentration == "Control")
  nrow(MatrixData4h.25nm)
  
  num4h.25 <- 668
  
  num_W_rep24h.25nm <- choose(num4h.25, num_treatments)
  matrix_4h.25nm <- matrix(nrow = num4h.25, ncol = num_loop)
  
  for (k in 1:num_loop) {
    matrix_4h.25nm[, k] <- sample(MatrixData4h.25nm$Concentration)
  }
  
  num_unique_allocations25.4h <- ncol(unique(matrix_4h.25nm, MARGIN = 2))
  
  
  T.rep25.4h<- rep(NA,length = num_unique_allocations25.4h)
  for(k in 1:num_unique_allocations25.4h){
    
    T.rep25.4h[k] <- mean(MatrixData4h.25nm$Basal_OCR.Av[matrix_4h.25nm[, k]==1]) - 
      mean(MatrixData4h.25nm$Basal_OCR.Av[matrix_4h.25nm[, k]==0])
  }
  
  #find pval
  p.val25.4h <- sum(T.rep25.4h <=  2.129578)/length(T.rep25.4h)
  
}
Matriz50nm.4h{
  MatrixData4h <- subset(BinaryData.2,  Dex_Exposure == "4h")
  MatrixData4h.50nm <- subset(MatrixData4h, Dex_Concentration == "50nM" |
                                Dex_Concentration == "Control")
  nrow(MatrixData4h.50nm)
  num4h.50 <- 668
  
  num_W_rep4h.50nm <- choose(num4h.50, num_treatments)
  matrix_4h.50nm <- matrix(nrow = num4h.50, ncol = num_loop)
  
  for (k in 1:num_loop) {
    matrix_4h.50nm[, k] <- sample(MatrixData4h.50nm$Concentration)
  }
  
  num_unique_allocations50.4h <- ncol(unique(matrix_4h.50nm, MARGIN = 2))
  
  
  T.rep50.4h<- rep(NA,length = num_unique_allocations50.4h)
  for(k in 1:num_unique_allocations50.4h){
    
    T.rep50.4h[k] <- mean(MatrixData4h.50nm$Basal_OCR.Av[matrix_4h.50nm[, k]==1]) - 
      mean(MatrixData4h.50nm$Basal_OCR.Av[matrix_4h.50nm[, k]==0])
  }
  
  #find pval
  p.val50.4h<- sum(T.rep50.4h <=2.993151)/length(T.rep50.4h)
  
}
Matriznm100.4h{
  MatrixData4h <- subset(BinaryData.2,  Dex_Exposure == "4h")
  MatrixData4h.100nm <- subset(MatrixData4h, Dex_Concentration == "100nM" |
                                 Dex_Concentration == "Control")
  nrow(MatrixData4h.100nm)
  
  num4h.100 <- 668
  
  num_W_rep4h.100nm <- choose(num4h.100, num_treatments)
  matrix_4h.100nm <- matrix(nrow = num4h.100, ncol = num_loop)
  
  for (k in 1:num_loop) {
    matrix_4h.100nm[, k] <- sample(MatrixData4h.100nm$Concentration)
  }
  
  num_unique_allocations100.4h <- ncol(unique(matrix_4h.100nm, MARGIN = 2))
  
  
  T.rep100.4h <- rep(NA,length = num_unique_allocations100.4h)
  for(k in 1:num_unique_allocations100.4h){
    
    T.rep100.4h[k] <- mean(MatrixData4h.100nm$Basal_OCR.Av[matrix_4h.100nm[, k]==1]) - 
      mean(MatrixData4h.100nm$Basal_OCR.Av[matrix_4h.100nm[, k]==0])
  }
  
  
  #find pval
  p.val100.4h <- sum(T.rep100.4h <= 3.331625)/length(T.rep100.4h)
  
  
}
Matrix500nm.4h{
  MatrixData4h <- subset(BinaryData.2,  Dex_Exposure == "4h")
  MatrixData4h.500nm <- subset(MatrixData4h, Dex_Concentration == "500nM" |
                                 Dex_Concentration == "Control")
  nrow(MatrixData4h.500nm)
  
  num4h.500 <- 668
  
  num_W_rep4h.500nm <- choose(num4h.500, num_treatments)
  matrix_4h.500nm <- matrix(nrow = num4h.500, ncol = num_loop)
  
  for (k in 1:num_loop) {
    matrix_4h.500nm[, k] <- sample(MatrixData4h.500nm$Concentration)
  }
  
  num_unique_allocations500.4h <- ncol(unique(matrix_4h.500nm, MARGIN = 2))
  
  
  T.rep500.4h <- rep(NA,length = num_unique_allocations500.4h)
  for(k in 1:num_unique_allocations500.4h){
    
    T.rep500.4h[k] <- mean(MatrixData4h.500nm$Basal_OCR.Av[matrix_4h.500nm[, k]==1]) - 
      mean(MatrixData4h.500nm$Basal_OCR.Av[matrix_4h.500nm[, k]==0])
  }
  
  
  p.val500.4h<- sum(T.rep500.4h <=  3.969358)/length(T.rep500.4h)
  
  
}
Matriz1um.4h{
  MatrixData4h <- subset(BinaryData.2,  Dex_Exposure == "4h")
  MatrixData4h.1um <- subset(MatrixData4h, Dex_Concentration == "1uM" |
                               Dex_Concentration == "Control")
  nrow(MatrixData4h.1um)
  
  num4h.1um <- 668
  
  num_W_rep4h.1um <- choose(num4h.1um, num_treatments)
  matrix_4h.1um <- matrix(nrow = num4h.1um, ncol = num_loop)
  
  for (k in 1:num_loop) {
    matrix_4h.1um[, k] <- sample(MatrixData4h.1um$Concentration)
  }
  
  num_unique_allocations1um.4h <- ncol(unique(matrix_4h.1um, MARGIN = 2))
  
  
  T.rep1um.4h <- rep(NA,length = num_unique_allocations1um.4h)
  for(k in 1:num_unique_allocations1um.4h){
    
    T.rep1um.4h[k] <- mean(MatrixData4h.1um$Basal_OCR.Av[matrix_4h.1um[, k]==1]) - 
      mean(MatrixData4h.1um$Basal_OCR.Av[matrix_4h.1um[, k]==0])
  }
  
  
  #find pval
  p.val1um.4h <- sum(T.rep1um.4h<=  4.954242)/length(T.rep1um.4h)
  
}
Matriz10um.4h{
  MatrixData4h <- subset(BinaryData.2,  Dex_Exposure == "4h")
  MatrixData4h.10um <- subset(MatrixData4h, Dex_Concentration == "10uM" |
                                Dex_Concentration == "Control")
  nrow(MatrixData4h.10um)
  
  num4h.10um <- 668
  
  num_W_rep4h.10um <- choose(num4h.10um, num_treatments)
  matrix_4h.10um <- matrix(nrow = num4h.10um, ncol = num_loop)
  
  for (k in 1:num_loop) {
    matrix_4h.10um[, k] <- sample(MatrixData4h.10um$Concentration)
  }
  
  num_unique_allocations10um.4h<- ncol(unique(matrix_4h.10um, MARGIN = 2))
  
  
  T.rep10um.4h<- rep(NA,length = num_unique_allocations10um.4h)
  for(k in 1:num_unique_allocations10um.4h){
    
    T.rep10um.4h[k] <- mean(MatrixData4h.10um$Basal_OCR.Av[matrix_4h.10um[, k]==1]) - 
      mean(MatrixData4h.10um$Basal_OCR.Av[matrix_4h.10um[, k]==0])
  }
  
  #find pval
  p.val10um.4h <- sum(T.rep10um.4h<=  3.288206)/length(T.rep10um.4h)
  
}

#Creating Matrix 8h
Matriz10nm.8h{
  MatrixData8h <- subset(BinaryData.2,  Dex_Exposure == "8h")
  MatrixData8h.10nm <- subset(MatrixData8h, Dex_Concentration == "10nM" |
                                Dex_Concentration == "Control")
  nrow(MatrixData8h.10nm)
  
  num8h.10 <- 828
  
  num_W_rep8h.10nm <- choose(num8h.10, num_treatments)
  matrix_8h.10nm <- matrix(nrow = num8h.10, ncol = num_loop)
  
  for (k in 1:num_loop) {
    matrix_8h.10nm[, k] <- sample(MatrixData8h.10nm$Concentration)
  }
  
  num_unique_allocations8h.10 <- ncol(unique(matrix_8h.10nm, MARGIN = 2))
  
  
  T.rep10.8h <- rep(NA,length = num_unique_allocations8h.10)
  for(k in 1:num_unique_allocations8h.10){
    
    T.rep10.8h[k] <- mean(MatrixData8h.10nm$Basal_OCR.Av[matrix_8h.10nm[, k]==1]) - 
      mean(MatrixData8h.10nm$Basal_OCR.Av[matrix_8h.10nm[, k]==0])
  }
  
  
  #find pval
  p.val10.8h <- sum(T.rep10.8h <= 1.818691)/length(T.rep10.8h)
  
}
Matriz25nm.8h{
  MatrixData8h<- subset(BinaryData.2,  Dex_Exposure == "8h")
  MatrixData8h.25nm <- subset(MatrixData8h, Dex_Concentration == "25nM" |
                                Dex_Concentration == "Control")
  nrow(MatrixData8h.25nm)
  
  num8h.25 <- 828
  
  num_W_rep8h.25nm <- choose(num8h.25, num_treatments)
  matrix_8h.25nm <- matrix(nrow = num8h.25, ncol = num_loop)
  
  for (k in 1:num_loop) {
    matrix_8h.25nm[, k] <- sample(MatrixData8h.25nm$Concentration)
  }
  
  num_unique_allocations25.8h <- ncol(unique(matrix_8h.25nm, MARGIN = 2))
  
  
  T.rep25.8h<- rep(NA,length = num_unique_allocations25.8h)
  for(k in 1:num_unique_allocations25.8h){
    
    T.rep25.8h[k] <- mean(MatrixData8h.25nm$Basal_OCR.Av[matrix_8h.25nm[, k]==1]) - 
      mean(MatrixData8h.25nm$Basal_OCR.Av[matrix_8h.25nm[, k]==0])
  }
  
  #find pval
  p.val25.8h <- sum(T.rep25.8h <=  1.875025)/length(T.rep25.8h)
  
}
Matriz50nm.8h{
  MatrixData8h <- subset(BinaryData.2,  Dex_Exposure == "8h")
  MatrixData8h.50nm <- subset(MatrixData8h, Dex_Concentration == "50nM" |
                                Dex_Concentration == "Control")
  nrow(MatrixData8h.50nm)
  num8h.50 <- 828
  
  num_W_rep8h.50nm <- choose(num4h.50, num_treatments)
  matrix_8h.50nm <- matrix(nrow = num8h.50, ncol = num_loop)
  
  for (k in 1:num_loop) {
    matrix_8h.50nm[, k] <- sample(MatrixData8h.50nm$Concentration)
  }
  
  num_unique_allocations50.8h <- ncol(unique(matrix_8h.50nm, MARGIN = 2))
  
  
  T.rep50.8h<- rep(NA,length = num_unique_allocations50.8h)
  for(k in 1:num_unique_allocations50.8h){
    
    T.rep50.8h[k] <- mean(MatrixData8h.50nm$Basal_OCR.Av[matrix_8h.50nm[, k]==1]) - 
      mean(MatrixData8h.50nm$Basal_OCR.Av[matrix_8h.50nm[, k]==0])
  }
  
  #find pval
  p.val50.8h<- sum(T.rep50.8h <=2.273082)/length(T.rep50.8h)
  
}
Matriznm100.8h{
  MatrixData8h <- subset(BinaryData.2,  Dex_Exposure == "8h")
  MatrixData8h.100nm <- subset(MatrixData8h, Dex_Concentration == "100nM" |
                                 Dex_Concentration == "Control")
  nrow(MatrixData8h.100nm)

  num8h.100 <- 828
  
  num_W_rep8h.100nm <- choose(num8h.100, num_treatments)
 matrix_8h.100nm <- matrix(nrow = num8h.100, ncol = num_loop)
  
  for (k in 1:num_loop) {
    matrix_8h.100nm[, k] <- sample(MatrixData8h.100nm$Concentration)
  }
  
  num_unique_allocations100.8h <- ncol(unique(matrix_8h.100nm, MARGIN = 2))
  
  T.rep100.8h <- rep(NA,length = num_unique_allocations100.8h)
  for(k in 1:num_unique_allocations100.8h){
    
    T.rep100.8h[k] <- (mean(MatrixData8h.100nm$Basal_OCR.Av[matrix_8h.100nm[, k]==1]) - 
      mean(MatrixData8h.100nm$Basal_OCR.Av[matrix_8h.100nm[, k]==0]))
  }
  
  
  #find pval
  p.val100.8h <- sum(T.rep100.8h <= 3.066767)/length(T.rep100.8h)
  
  
#Visualization
  ggplot(data.frame(T_rep = T.rep100.8h), aes(x = T_rep)) +
    geom_histogram(binwidth = 0.1, fill = "skyblue", color = "black", alpha = 0.7) +
    geom_vline(xintercept = 3.066767, linetype = "dashed", color = "red", size = 1) +
    geom_text(aes(x = 2, y = 5000, label = paste("p-value: <1/1000")), color = "black", size = 4) +
    geom_text(aes(x = 2, y = 2500, label = paste("Observed Test Statistic")), color = "red", size = 4) +
    geom_text(aes(x = 2, y = 2000, label = paste("= 3.066767")), color = "black", size = 4) +
    labs(x = "Test Statistic (T)", y = "Frequency", title = "Distribution of Test Statistic") +
    theme_minimal()
}
Matrix500nm.8h{
  MatrixData8h <- subset(BinaryData.2,  Dex_Exposure == "8h")
  MatrixData8h.500nm <- subset(MatrixData8h, Dex_Concentration == "500nM" |
                                 Dex_Concentration == "Control")
  nrow(MatrixData8h.500nm)
  
  num8h.500 <- 828
  
  num_W_rep8h.500nm <- choose(num8h.500, num_treatments)
  matrix_8h.500nm <- matrix(nrow = num8h.500, ncol = num_loop)
  
  for (k in 1:num_loop) {
    matrix_8h.500nm[, k] <- sample(MatrixData8h.500nm$Concentration)
  }
  
  num_unique_allocations500.8h <- ncol(unique(matrix_8h.500nm, MARGIN = 2))
  
  
  T.rep500.8h <- rep(NA,length = num_unique_allocations500.8h)
  for(k in 1:num_unique_allocations500.8h){
    
    T.rep500.8h[k] <- mean(MatrixData8h.500nm$Basal_OCR.Av[matrix_8h.500nm[, k]==1]) - 
      mean(MatrixData8h.500nm$Basal_OCR.Av[matrix_8h.500nm[, k]==0])
  }
  
  
  p.val500.8h<- sum(T.rep500.8h <=  2.129999)/length(T.rep500.8h)
  
  
}
Matriz1um.8h{
  MatrixData8h <- subset(BinaryData.2,  Dex_Exposure == "8h")
  MatrixData8h.1um <- subset(MatrixData8h, Dex_Concentration == "1uM" |
                               Dex_Concentration == "Control")
  nrow(MatrixData8h.1um)
  
  num8h.1um <- 828
  
  num_W_rep8h.1um <- choose(num8h.1um, num_treatments)
  matrix_8h.1um <- matrix(nrow = num8h.1um, ncol = num_loop)
  
  for (k in 1:num_loop) {
    matrix_8h.1um[, k] <- sample(MatrixData8h.1um$Concentration)
  }
  
  num_unique_allocations1um.8h <- ncol(unique(matrix_8h.1um, MARGIN = 2))
  
  
  T.rep1um.8h <- rep(NA,length = num_unique_allocations1um.8h)
  for(k in 1:num_unique_allocations1um.8h){
    
    T.rep1um.8h[k] <- mean(MatrixData8h.1um$Basal_OCR.Av[matrix_8h.1um[, k]==1]) - 
      mean(MatrixData8h.1um$Basal_OCR.Av[matrix_8h.1um[, k]==0])
  }
  
  
  #find pval
  p.val1um.8h <- sum(T.rep1um.8h <=  2.620658)/length(T.rep1um.8h)
  
}
Matriz10um.8h{
  MatrixData8h <- subset(BinaryData.2,  Dex_Exposure == "8h")
  MatrixData8h.10um <- subset(MatrixData8h, Dex_Concentration == "10uM" |
                                Dex_Concentration == "Control")
  nrow(MatrixData8h.10um)
  
  num8h.10um <- 828
  
  num_W_rep8h.10um <- choose(num8h.10um, num_treatments)
  matrix_8h.10um <- matrix(nrow = num8h.10um, ncol = num_loop)
  
  for (k in 1:num_loop) {
    matrix_8h.10um[, k] <- sample(MatrixData8h.10um$Concentration)
  }
  
  num_unique_allocations10um.8h<- ncol(unique(matrix_8h.10um, MARGIN = 2))
  
  
  T.rep10um.8h<- rep(NA,length = num_unique_allocations10um.8h)
  for(k in 1:num_unique_allocations10um.8h){
  
    T.rep10um.8h[k] <- mean(MatrixData8h.10um$Basal_OCR.Av[matrix_8h.10um[, k]==1]) - 
      mean(MatrixData8h.10um$Basal_OCR.Av[matrix_8h.10um[, k]==0])
  }
  
  #find pval
  p.val10um.8h <- sum(T.rep10um.8h<=  2.096873)/length(T.rep10um.8h)
  
}


