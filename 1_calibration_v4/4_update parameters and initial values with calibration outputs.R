# clear the environment
#rm(list = ls())
#===========================================================================================================
# save the rounded values of final parameter combinations with the lowest error score for verification step
#===========================================================================================================
# load multiple libraries
lapply(c("deSolve", "ggplot2", "dplyr", "reshape2", "tidyverse", "readr", "readxl", 
         "lubridate", "data.table", "writexl", "git2r", "scales", "lhs", "sensitivity", "scatterplot3d", "plotly"), 
       library, character.only = TRUE)

# no scientific notation
options(scipen = 999)

# round the values to 3 decimal places for beta_h, 5 decimal places for kappa_h and 7 decimal places for kappa_1
best_comb <- data.frame(beta_h = 0.1177, kappa_h = 0.00124, kappa_1 = 0.0000455)

# load the existing parameter and initial values table
initial_values <- read.csv("1_calibration_v4/calibration_inputs/inputs - initial_values - for calibration - model1 - v5.csv", header = TRUE, stringsAsFactors = FALSE) 
parameters <- read.csv("1_calibration_v4/calibration_inputs/inputs - parameters_v4 - for calibration - model1 - v5.csv", header = TRUE, stringsAsFactors = FALSE) 

# convert data frames to data tables
initial_values <- as.data.table(initial_values)
parameters <- as.data.table(parameters)

# ensure everything is in numeric form 
initial_values <- initial_values[, value := as.numeric(gsub("[^0-9.-]","",value))]

# update the parameters table with the calibrated values
parameters[epi_symbol == "beta_h", baseline_value := best_comb$beta_h] # use the calibrated beta_h value
parameters[epi_symbol == "kappa_h", baseline_value := best_comb$kappa_h]  # use the calibrated kappa_h value
parameters[epi_symbol == "kappa_1", baseline_value := best_comb$kappa_1]  # use the calibrated kappa_1 value

# check the updated parameters table
parameters

# save the updated parameters table to csv
write.csv(parameters, "2_verification_v4/verification_inputs/inputs - parameters_v5 - for verification - model1.csv", row.names = FALSE)
write.csv(parameters, "3_intervention_v4/intervention_inputs/inputs - parameters_v5 - for intervention - model1.csv", row.names = FALSE)

# update the initial value of each compartment

# load the model and input script
source("1_calibration_v4/calibration_inputs/model1_v4_for_calibration.R")
source("1_calibration_v4/calibration_inputs/inputs_v4_calibration_model1.R")$value

#==================================================================================================
#### 2. Setup initial values #### #ensure the order is the same as the order of the equations above
#==================================================================================================
istate1 <- c(S_c1 = S_c1, C_c1 = C_c1, I_c1 = I_c1,
S_h1 = S_h1, C_h1 = C_h1, I_h1 = I_h1, 
D_h1 = D_h1, CInc_1a = CInc_1a, CInc_1b = CInc_1b, CInc_1c = CInc_1c, CInc_1d = CInc_1d) 

#===============================================================
#### 3. Setup parameters ####
#===============================================================
parameters1 <- c(mu_1 = mu_1, #recruitment rate and natural mortality rate
sigma_1a = sigma_1a, #hospital admissions rate for the susceptible and colonised population in the community
sigma_1b = sigma_1b, #hospital admissions rate for persons with BSI in community
epsilon_1 = epsilon_1, #natural clearance rate from colonisation (only in the community)
gamma_1 = gamma_1, #hospital discharge rate 
kappa_1 =  kappa_1, #kappa_1, #the daily rate of progression from colonised state to BSI in the community
kappa_h =  kappa_h,#kappa_h,  #the daily rate of progression from colonisation to BSI in the hospital
tau_h = tau_h, #the daily rate of recovery from BSI back to the colonisation state (in the hospital)
omega_h = omega_h, #the daily mortality rate from BSI (in the hospital) #also used as a recruitment rate to maintain constant total population
beta_h =  beta_h) #transmission coefficient in the hospital

# input the final parameter combination to the model
parameters1["beta_h"] <- best_comb$beta_h
parameters1["kappa_h"] <- best_comb$kappa_h
parameters1["kappa_1"] <- best_comb$kappa_1

#===============================================================
#### 4. Setup timesteps ####
#===============================================================
start1 <- 0 #start time for the simulation
end1 <- 365*200 #end point of the timestep #for calibration, the timestep is set up to be arbitrarily long to ensure the system reaches equilibrium before we move on to the next step
tps1 <- seq(start1, end1, by = 1) #time steps for the simulation

# run the model 
out1 <- ode(y = istate1, times = tps1, func = klebsiella_1, parms = parameters1, method = "lsoda") 

# check at which time step the model starts to reach equilibrium 

# number of susceptible - community
plot(out1[, "time"], out1[, "S_c1"], 
     type = "l", 
     lwd = 2, 
     col = "#0072B2", 
     xlab = "Number of days", 
     ylab = "Number of people", 
     main = "Number of susceptibles in the community",
     ylim = c(min(out1[, "S_c1"]), max(out1[, "S_c1"]))) # the model starts settling in after time == 365*60

# number of colonised - community
plot(out1[, "time"], out1[, "C_c1"], 
     type = "l", 
     lwd = 2, 
     col = "#E69F00", 
     xlab = "Number of days", 
     ylab = "Number of people", 
     main = "Number of colonised in the community",
     ylim = c(min(out1[, "C_c1"]), max(out1[, "C_c1"]))) # the model starts settling in after time == 365*60

# number of BSI - hospital
plot(out1[, "time"], out1[, "I_c1"], 
     type = "l", 
     lwd = 2, 
     col = "#CC79A7", 
     xlab = "Number of days", 
     ylab = "Number of people", 
     main = "Number of BSI in the community",
     ylim = c(min(out1[, "I_c1"]), max(out1[, "I_c1"]))) # the model starts settling in after time == 365*60


# create the dataframe from the output
df.after_calibration <- as.data.frame(out1) 

# use the values at equilibrium as initial states to shorten the running time when applying intervention
# discard run-in period for intervention
df.after_calibration_eq <- df.after_calibration %>% filter(time == max(time)) 

# update the initial values table with the equilibrium values we get from this calibration phase
initial_values[epi_symbol == "Sc1", value := df.after_calibration_eq$S_c1]
initial_values[epi_symbol == "Cc1", value := df.after_calibration_eq$C_c1]
initial_values[epi_symbol == "Ic1", value := df.after_calibration_eq$I_c1]
initial_values[epi_symbol == "Sh1", value := df.after_calibration_eq$S_h1]
initial_values[epi_symbol == "Ch1", value := df.after_calibration_eq$C_h1]
initial_values[epi_symbol == "Ih1", value := df.after_calibration_eq$I_h1]
initial_values[epi_symbol == "Dh1", value := df.after_calibration_eq$D_h1]
initial_values[epi_symbol == "CInc_1a", value := df.after_calibration_eq$CInc_1a]
initial_values[epi_symbol == "CInc_1b", value := df.after_calibration_eq$CInc_1b]
initial_values[epi_symbol == "CInc_1c", value := df.after_calibration_eq$CInc_1c]
initial_values[epi_symbol == "CInc_1d", value := df.after_calibration_eq$CInc_1d]

# check the updated initial values table
initial_values

# save the updated initial values table to csv
# we are going to use initial values from our calibration for verifications so the initial values at equilibrium after calibration will just be saved for intervention
write.csv(initial_values, "3_intervention_v4/intervention_inputs/inputs - initial_values - for intervention - model1.csv", row.names = FALSE)


# use the same initial values for calibration in verification to show the change overtime
# but we can use shorter run-in period for verification as we know the model reaches equilibrium before time == 365*200
df.after_calibration_eq <- df.after_calibration %>% filter(time == 0) 

# update the initial values table with the equilibrium values we get from this calibration phase
initial_values[epi_symbol == "Sc1", value := df.after_calibration_eq$S_c1]
initial_values[epi_symbol == "Cc1", value := df.after_calibration_eq$C_c1]
initial_values[epi_symbol == "Ic1", value := df.after_calibration_eq$I_c1]
initial_values[epi_symbol == "Sh1", value := df.after_calibration_eq$S_h1]
initial_values[epi_symbol == "Ch1", value := df.after_calibration_eq$C_h1]
initial_values[epi_symbol == "Ih1", value := df.after_calibration_eq$I_h1]
initial_values[epi_symbol == "Dh1", value := df.after_calibration_eq$D_h1]
initial_values[epi_symbol == "CInc_1a", value := df.after_calibration_eq$CInc_1a]
initial_values[epi_symbol == "CInc_1b", value := df.after_calibration_eq$CInc_1b]
initial_values[epi_symbol == "CInc_1c", value := df.after_calibration_eq$CInc_1c]
initial_values[epi_symbol == "CInc_1d", value := df.after_calibration_eq$CInc_1d]

# check the updated initial values table
initial_values

write.csv(initial_values, "2_verification_v4/verification_inputs/inputs - initial_values - for verification - model1.csv", row.names = FALSE)







