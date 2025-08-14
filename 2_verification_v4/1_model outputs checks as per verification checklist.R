# script for model verification after calibration to ensure everything add in 
# the plots are for one average hospital
# we translated the outputs for the scale of the entire England when implementing the intervention, sensitivity analysis, and scenario analysis

#==================================================#
# run the model with calibrated inputs
#==================================================#

# clear the environment
rm(list = ls())

#===============================================================================
# 1. package installations # only install the ones that are not installed yet

install_if_not_installed <- function(pkg) {
  install_packages <- pkg[!pkg %in% installed.packages()[, "Package"]]
  if (length(install_packages)) install.packages(install_packages)
}

# install packages if not installed
install_if_not_installed(c("deSolve", "ggplot2", "dplyr", "reshape2", "tidyverse", 
                           "readr", "readxl", "lubridate", "purrr", "data.table", 
                           "writexl", "git2r"))

# load multiple libraries
lapply(c("deSolve", "ggplot2", "dplyr", "reshape2", "tidyverse", "readr", "readxl", 
         "lubridate", "purrr", "data.table", "writexl", "git2r", "scales"), library, character.only = TRUE)

# load model chunk
source("2_verification_v4/verification_inputs/model1_v4_for_verification.R")

#==================================================#
# plot the basic model outputs 
#==================================================#

# initiate a blank list to save the plots
plots <- list()

# check the total population size

# total population size (community + hospital)
total_population <- rowSums(out1[, c("S_c1", "C_c1", "I_c1", "S_h1", "C_h1", "I_h1")])

# set margins for the plots
par(mar = c(5, 7, 4, 2))  

plot(out1[, "time"], total_population, 
     type = "l", 
     lwd = 2, 
     col = "#000000", 
     xlab = "Number of days", 
     ylab = "", 
     main = "Total population (community and hospital)",
     ylim = c(75000, 81000),
     yaxt = "n") # the model starts settling in after time == 3650

axis(side = 2,
     at = seq(75000, 81000, by = 1000),
     labels = format(seq(75000, 81000, by = 1000), big.mark = ","),
     las = 1)

mtext("Number of people", side = 2, line = 5, cex = 1.1)

# save as R object
plots[["total_population"]] <- recordPlot()


# total population size (community)
total_pop_comm <- rowSums(out1[, c("S_c1", "C_c1", "I_c1")])

# set margins for the plots
par(mar = c(5, 7, 4, 2))  

plot(out1[, "time"], total_pop_comm, 
     type = "l", 
     lwd = 2, 
     col = "#4D4D4D", 
     xlab = "Number of days", 
     ylab = "", 
     main = "Total population (community)",
     ylim = c(75000, 81000),
     yaxt = "n") # the model starts settling in after time == 3650

axis(side = 2,
     at = seq(75000, 81000, by = 1000),
     labels = format(seq(75000, 81000, by = 1000), big.mark = ","),
     las = 1)

mtext("Number of people", side = 2, line = 5, cex = 1.1)

# save as R object
plots[["total_population_community"]] <- recordPlot()


# total population size (hospital)
total_pop_hosp <- rowSums(out1[, c("S_h1", "C_h1", "I_h1")])

# set margins for the plots
par(mar = c(5, 7, 4, 2))  

plot(out1[, "time"], total_pop_hosp, 
     type = "l", 
     lwd = 2, 
     col = "#009E73", 
     xlab = "Number of days", 
     ylab = "", 
     main = "Total population (single hospital)",
     ylim = c(0, 450),
     yaxt = "n") # the model starts settling in after time == 3650

axis(side = 2,
     at = seq(0, 450, by = 150),
     labels = format(seq(0, 450, by = 150), big.mark = ","),
     las = 1)

mtext("Number of patients", side = 2, line = 5, cex = 1.1)

# save as R object
plots[["total_population_hospital"]] <- recordPlot()

#-----------------------------------#
# community compartments
#-----------------------------------#
# susceptible - community
par(mar = c(5, 7, 4, 2))  

plot(out1[, "time"], out1[, "S_c1"], 
     type = "l", 
     lwd = 2, 
     col = "#0072B2", 
     xlab = "Number of days", 
     ylab = "", 
     main = "Number of susceptibles in the community",
     ylim = c(78340, 78360),
     yaxt = "n") # the model starts settling in after time == 3650

axis(side = 2,
     at = seq(78340, 78360, by = 5),
     labels = format(seq(78340, 78360, by = 5), big.mark = ","),
     las = 1)

mtext("Number of people", side = 2, line = 5, cex = 1.1)

# save as R object
plots[["susceptible_community"]] <- recordPlot()



# colonised - community
par(mar = c(5, 7, 4, 2))  

plot(out1[, "time"], out1[, "C_c1"], 
     type = "l", 
     lwd = 2, 
     col = "#E69F00", 
     xlab = "Number of days", 
     ylab = "", 
     main = "Number of colonised in the community",
     ylim = c((min(out1[, "C_c1"])), 1690),
     yaxt = "n") # the model starts settling in after time == 3650

axis(side = 2,
     at = seq((min(out1[, "C_c1"])), 1690, by = 2),
     labels = format(seq(round(min(out1[, "C_c1"])), 1690, by = 2), big.mark = ","),
     las = 1)

mtext("Number of people", side = 2, line = 5, cex = 1.1)

# save as R object
plots[["colonised_community"]] <- recordPlot()


# BSI - community
par(mar = c(5, 7, 4, 2))  

plot(out1[, "time"], out1[, "I_c1"], 
     type = "l", 
     lwd = 2, 
     col = "#CC79A7", 
     xlab = "Number of days", 
     ylab = "", 
     main = "Number of BSI in the community",
     yaxt = "n") # the model starts settling in after time == 3650

axis(side = 2,
     at = seq(min(out1[, "I_c1"]), max(out1[, "I_c1"]), length.out = 5),
     labels = formatC(seq(min(out1[, "I_c1"]), max(out1[, "I_c1"]), length.out = 5),
                      format = "f", digits = 5),
     las = 1)

mtext("Number of people", side = 2, line = 5, cex = 1.1)

# save as R object
plots[["BSI_community"]] <- recordPlot()



#-----------------------------------#
# hospital compartments
#-----------------------------------#
# susceptible - hospital
par(mar = c(5, 7, 4, 2)) 

plot(out1[, "time"], out1[, "S_h1"], 
     type = "l", 
     lwd = 2, 
     col = "#009E73", 
     xlab = "Number of days", 
     ylab = "", 
     main = "Number of susceptibles in the hospital",
     ylim = c(min(out1[, "S_h1"]), max(out1[, "S_h1"])),
     yaxt = "n") # the model starts settling in after time == 3650

axis(side = 2,
     at = seq(min(out1[, "S_h1"]), max(out1[, "S_h1"]), length.out = 5),
     labels = formatC(seq(min(out1[, "S_h1"]), max(out1[, "S_h1"]), length.out = 5),
                      format = "f", digits = 2),
     las = 1)

mtext("Number of patients", side = 2, line = 5, cex = 1.1)

# save as R object
plots[["susceptible_hospital"]] <- recordPlot()


# colonised - hospital
par(mar = c(5, 7, 4, 2)) 

plot(out1[, "time"], out1[, "C_h1"], 
     type = "l", 
     lwd = 2, 
     col = "#D55E00", 
     xlab = "Number of days", 
     ylab = "", 
     main = "Number of colonised in the hospital",
     ylim = c(min(out1[, "C_h1"]), max(out1[, "C_h1"])),
     yaxt = "n") # the model starts settling in after time == 3650

axis(side = 2,
     at = seq(min(out1[, "C_h1"]), max(out1[, "C_h1"]), length.out = 5),
     labels = formatC(seq(min(out1[, "C_h1"]), max(out1[, "C_h1"]), length.out = 5),
                      format = "f", digits = 2),
     las = 1)

mtext("Number of patients", side = 2, line = 5, cex = 1.1)

# save as R object
plots[["colonised_hospital"]] <- recordPlot()


# BSI - hospital
par(mar = c(5, 7, 4, 2)) 

plot(out1[, "time"], out1[, "I_h1"], 
     type = "l", 
     lwd = 2, 
     col = "#BA55D3", 
     xlab = "Number of days", 
     ylab = "", 
     main = "Number of BSI in the hospital",
     ylim = c(min(out1[, "I_h1"]), max(out1[, "I_h1"])),
     yaxt = "n") # the model starts settling in after time == 3650


axis(side = 2,
     at = seq(min(out1[, "I_h1"]), max(out1[, "I_h1"]), length.out = 5),
     labels = formatC(seq(min(out1[, "I_h1"]), max(out1[, "I_h1"]), length.out = 5),
                      format = "f", digits = 4),
     las = 1)

mtext("Number of patients", side = 2, line = 5, cex = 1.1)

# save as R object
plots[["BSI_hospital"]] <- recordPlot()



#==================================================#
# save all plots to a file
#==================================================#
# create a directory to save the plots
if (!dir.exists("2_verification_v4/verification_outputs")) {
  dir.create("2_verification_v4/verification_outputs", recursive = TRUE)
}

# save all plots to a file
saveRDS(plots, file = "2_verification_v4/verification_outputs/basic_model_checks_plots.rds")

#==================================================#
# save plots as PNG files
#==================================================#

plots <- readRDS("2_verification_v4/verification_outputs/basic_model_checks_plots.rds")

for (name in names(plots)){
  png(file.path("2_verification_v4/verification_outputs", paste0(name, ".png")), width = 1200, height = 1500, res = 200)
  replayPlot(plots[[name]])
  dev.off()
}

