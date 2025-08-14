#==========================================================================================================
# This script generates Latin hypercube samples for the sensitivity analysis
#==========================================================================================================
# prepare the environment

# clear the environment
rm(list = ls())

# load required packages
if(!requireNamespace("lhs", quietly = TRUE)) install.packages("lhs") #change all install.package codes to install only when needed!
if(!requireNamespace("sensitivity", quietly = TRUE)) install.packages("sensitivity") #change all install.package codes to install only when needed!

# load the lhs package
require(lhs)
set.seed(123)

# load required sources
source("4_sensitivity_analysis_v4/sensitivity_inputs/model1_v4_for_sensitivity analysis.R")

#==========================================================================================================
# define sample size and parameter counts
#==========================================================================================================
# define the number of samples per parameter
sample_size <- 2000

# define the number of parameters
parameter_count <- 6 # number of parameters to vary

# generate Latin hypercube samples
# lhs_samples <- randomLHS(sample_size, parameter_count)

#=================================================================
# load parameter tables
# Load initial values and parameters from CSV files
#=================================================================

# load initial values and parameters from CSV files
initial_values <- read.csv("4_sensitivity_analysis_v4/sensitivity_inputs/inputs - initial_values - for sensitivity analysis - model1.csv", header = TRUE, stringsAsFactors = FALSE) 
parameters <- read.csv("4_sensitivity_analysis_v4/sensitivity_inputs/inputs - parameters_v5 - for sensitivity analysis - model1.csv", header = TRUE, stringsAsFactors = FALSE) 

# convert data frames to data tables
initial_values <- as.data.table(initial_values)
parameters <- as.data.table(parameters)

# ensure everything is in numeric form 
initial_values <- initial_values[, value := as.numeric(gsub("[^0-9.-]","",value))]

#=============================================================================================
# get min and max of each parameter

# define parameter range
# min
min_beta_h <- min_beta_h
min_gamma_1 <- min_gamma_1
min_kappa_1 <- min_kappa_1
min_epsilon_1 <- min_epsilon_1
min_kappa_h <- min_kappa_h
min_tau_h <- min_tau_h

# max
max_beta_h <- max_beta_h
max_gamma_1 <- max_gamma_1
max_kappa_1 <- max_kappa_1
max_epsilon_1 <- max_epsilon_1
max_kappa_h <- max_kappa_h
max_tau_h <- max_tau_h
#=============================================================================================



#=============================================================================================
# sample from a predefined range of each parameter based on literature values
# define the range of each parameter
# is qunif the right function to use here?
# partially adapted from https://github.com/moyinNUHS/abxduration_abm/blob/main/run_models/prcc_simple3state.R

parameter_ranges_klebsiella1 <- list(
  beta_h = c("qunif", list(min = min_beta_h, max = max_beta_h)), #recruitment rate and natural mortality rate
  epsilon_1 = c("qunif", list(min = min_epsilon_1, max = max_epsilon_1)), #natural clearance rate from colonisation (only in the community)
  gamma_1 = c("qunif", list(min = min_gamma_1, max = max_gamma_1)), #hospital discharge rate 
  kappa_1 = c("qunif", list(min = min_kappa_1, max = max_kappa_1)), #the daily rate of progression from colonised state to BSI in the community
  kappa_h = c("qunif", list(min = min_kappa_h, max = max_kappa_h)), #the daily rate of progression from colonisation to BSI in the hospital
  tau_h = c("qunif", list(min = min_tau_h, max = max_tau_h)) #the daily rate of recovery from BSI back to the colonisation state (in the hospital)
)
#==============================================================================================


#==============================================================================================
# create latin hypercube samples from the parameter_ranges_klebsiella1
#==============================================================================================
# develop a matrix of Latin hypercube samples
lhs_matrix_klebsiella1 <- randomLHS(sample_size, length(parameter_ranges_klebsiella1))

# print the first few samples
#print(head(scaled_samples))

# scale the Latin hypercube samples to the parameter ranges
lhs_samples_klebsiella1 <- matrix(0, nrow = sample_size, ncol = length(parameter_ranges_klebsiella1))

# assign the names of the parameters to the columns of the scaled samples matrix
colnames(lhs_samples_klebsiella1) <- names(parameter_ranges_klebsiella1)

# iterate over each parameter to scale the samples
for (i in seq_along(parameter_ranges_klebsiella1)) {
  # get the parameter name
  param_name <- names(parameter_ranges_klebsiella1)[i]
  
  # get the distribution and min and max values for the parameter
  param_dist <- parameter_ranges_klebsiella1[[param_name]][[1]]
  param_min <- parameter_ranges_klebsiella1[[param_name]][[2]]
  param_max <- parameter_ranges_klebsiella1[[param_name]][[3]]
  
  # scale the samples to the real parameter range
  lhs_samples_klebsiella1[, i] <- lhs_matrix_klebsiella1[, i] * (param_max - param_min) + param_min
}

# check the first few scaled samples
head(lhs_samples_klebsiella1)

# transform to data frame for easier handling
lhs_samples_klebsiella1_df <- as.data.frame(lhs_samples_klebsiella1)

View(lhs_samples_klebsiella1)

# check once again
head(lhs_samples_klebsiella1_df) # ok, looks good -> the ranges assigned are obtained as expected

# save the scaled samples to a CSV file
write.csv(lhs_samples_klebsiella1_df, "4_sensitivity_analysis_v4/latin_hypercube_samples_klebsiella1.csv", row.names = FALSE)


#==============================================================================================
# run the model with the generated samples
# adapted from https://bogaotory.github.io/training-r/intro-sa/sa.nb.html
# define a function to run the model with the given parameters

# load the model
source("4_sensitivity_analysis_v4/sensitivity_inputs/model1_v4_for_sensitivity analysis.R")

model_run <- function(parameters_row, state, times, model, method) {
    # convert parameter combinations row to a list
    parameters <- as.list(parameters_row)
  
    # run the model with the given parameters and initial values
    # no event list for now, as we are not running the intervention and to avoid ugly decimals
    out <- ode(y = state, times = times, func = model, parms = parameters, method = "lsoda") 
    
    # convert output to data frame for easier handling
    out_df <- as.data.frame(out)
    
    # calculate the incidence for CInc_1c (the number of new cases in the community)
    out_df$inc_bsi_hosp <- c(0, diff(out_df$CInc_1c)) # calculate the incidence of BSI hospitalisation
    
    # calculate the incidence for CInc_1d (the number of new cases in the hospital)
    out_df$inc_hosp_onset <- c(0, diff(out_df$CInc_1d)) # calculate the incidence of hospital-onset-cases
    
    # calculate total number of cases
    out_df$total_cases <- out_df$inc_bsi_hosp  +  out_df$inc_hosp_onset
    
    return(out_df)
}

fixed_params <- c(
  mu_1 = mu_1,
  sigma_1a = sigma_1a,
  sigma_1b = sigma_1b,
  omega_h = omega_h,
  vacc_coverage = 0,
  rho_1 = 0
)

output_test <- model_run(
  parameters_row = c(fixed_params, lhs_samples_klebsiella1_df[3,]),
  state = istate1,
  times = tps1,
  model = klebsiella_1,
  method =  "lsoda"
)

# save the outputs for checking 
write.csv(output_test, "4_sensitivity_analysis_v4/sensitivity_outputs/latin_hypercube_samples_klebsiella1_propagated_check.csv", row.names = FALSE)

#==============================================================================================

function_name <- "run.varied.parameters"
results.dataframe <- NULL # initialize an empty list to store results

# Loop through each row of the scaled samples data frame
if(!file.exists(function_name)){
  start_time_stamp <- Sys.time() # start the timer
  
  for (rr in 1:nrow(lhs_samples_klebsiella1_df)) {
    
    # Extract the parameters for the current sample
    p1 <- c(lhs_samples_klebsiella1_df[rr, ]) #,vacc.on = TRUE 
    
    # Define fixed parameters
    fixed_params <- c(
      mu_1 = mu_1,
      sigma_1a = sigma_1a,
      sigma_1b = sigma_1b,
      omega_h = omega_h,
      vacc_coverage = 0,
      rho_1 = 0
    )
    
    # Run the model with the current parameters
    run <- model_run(parameters = c(fixed_params,p1),
                     state = istate1,
                     times = tps1,
                     model = klebsiella_1,
                     method =  "lsoda")
    
    # combine results
    combine.df <- c(p1, run) # combine the varied parameters and results into a single list
    
    # bind the results to the parameters
    results.dataframe <- bind_rows(results.dataframe, combine.df) # bind the parameters to the results data frame
  }
  end_time_stamp <- Sys.time() # end the timer
  end_time_stamp - start_time_stamp # print the time taken to run the model
  
  # Save the results 
  saveRDS(results.dataframe, file = paste0("4_sensitivity_analysis_v4/", function_name, ".rds")) # save the results to a file
}
#=============================================================================================
# Load the results from the saved file
results.dataframe.file <- readRDS(paste0("4_sensitivity_analysis_v4/run.varied.parameters.rds"))

options(scipen = 999)

tail(results.dataframe.file,5) # check the first few rows of the results data frame

View(results.dataframe.file %>% filter(time >= 1824))
#=============================================================================================
# Save the results to a CSV file
#write.csv(results.dataframe.file, "4_sensitivity_analysis_v4/klebsiella_1_lhs_results_propagated.csv", row.names = FALSE)
#=============================================================================================





