#Using data from Mean Difference Code - 8h & 100nM

#Spread of Means Differences across combinations of exposure and concentration

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


#Visualization of Frequency & Distribution of Test Statistics

ggplot(data.frame(T_rep = T.rep100.8h), aes(x = T_rep)) +
  geom_histogram(binwidth = 0.1, fill = "skyblue", color = "black", alpha = 0.7) +
  geom_vline(xintercept = 3.066767, linetype = "dashed", color = "red", size = 1) +
  geom_text(aes(x = 2, y = 5000, label = paste("p-value: <1/1000")), color = "black", size = 4) +
  geom_text(aes(x = 2, y = 2500, label = paste("Observed Test Statistic")), color = "red", size = 4) +
  geom_text(aes(x = 2, y = 2000, label = paste("= 3.066767")), color = "black", size = 4) +
  labs(x = "Test Statistic (T)", y = "Frequency", title = "Distribution of Test Statistic") +
  theme_minimal()


#AUC visualization using data from FisherExact AUC code

#Control
plot(exposures, controlOCRvals, type = "o", ylim = range(controlOCRvals),
     xlim = range(exposures), main = "Control treatment plot",
     xlab = "Exposure", ylab = "OCR values")
#100nM
plot(exposures, conc100nm, type = "o", ylim = range(conc100nm),
     xlim = range(exposures), main = "100nM treatment plot",
     xlab = "Exposure", ylab = "OCR values")

# Plot100nM treatment on the same graph and fills in AUC
lines(exposures, conc100nm, type = "o", col = "red")

polygon(c(exposures, rev(exposures)), c(controlOCRvals, rev(conc100nm)), col = "gray", border = NA)
