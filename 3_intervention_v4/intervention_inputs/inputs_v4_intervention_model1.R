# sub-script to store parameters values and initial values for the model

# load necessary library
library(data.table)

# load initial values and parameters from CSV files
initial_values <- read.csv("3_intervention_v4/intervention_inputs/inputs - initial_values - for intervention - model1.csv", header = TRUE, stringsAsFactors = FALSE) 
parameters <- read.csv("3_intervention_v4/intervention_inputs/inputs - parameters_v5 - for intervention - model1.csv", header = TRUE, stringsAsFactors = FALSE) 


# convert data frames to data tables
initial_values <- as.data.table(initial_values)
parameters <- as.data.table(parameters)

# ensure everything is in numeric form 
initial_values <- initial_values[, value := as.numeric(gsub("[^0-9.-]","",value))]

# define the initial values
# initial values for klebsiella_1
S_c1 <- initial_values[initial_values$epi_symbol == "Sc1", value] #susceptible in the community
C_c1  <- initial_values[initial_values$epi_symbol == "Cc1", value] #colonised in the community 
I_c1 <-  initial_values[initial_values$epi_symbol == "Ic1", value] #bsi in the community
VS_c1 <-  initial_values[initial_values$epi_symbol == "VS_c1", value] #vaccinated susceptible in the community
VC_c1 <-  initial_values[initial_values$epi_symbol == "VC_c1", value] #vaccinated colonised in the community
S_h1 <- initial_values[initial_values$epi_symbol == "Sh1", value] #susceptible in the hospital
C_h1 <- initial_values[initial_values$epi_symbol == "Ch1", value] #colonised in the hospital
I_h1 <- initial_values[initial_values$epi_symbol == "Ih1", value] #bsi in the hospital
VS_h1 <-  initial_values[initial_values$epi_symbol == "VS_h1", value] #vaccinated susceptible in the hospital
VC_h1 <-  initial_values[initial_values$epi_symbol == "VC_h1", value] #vaccinated colonised in the hospital
D_h1 <- initial_values[initial_values$epi_symbol == "Dh1", value] #cumulative number of deaths associated with BSI
CInc_1a <- initial_values[initial_values$epi_symbol == "CInc_1a", value] #cumulative incidence of colonisation in the hospital 
CInc_1b <- initial_values[initial_values$epi_symbol == "CInc_1b", value] #cumulative incidence of community-onset cases
CInc_1c <- initial_values[initial_values$epi_symbol == "CInc_1c", value] #cumulative incidence of BSI hospitalisations
CInc_1d <- initial_values[initial_values$epi_symbol == "CInc_1d", value] #cumulative incidence of hospital-onset cases


# define the parameters
# parameters for klebsiella_1
# make sure annual rates are converted to daily rates!!!
mu_1 <- parameters[parameters$epi_symbol == "mu_1", baseline_value] #recruitment rate and natural mortality rate
beta_h <- parameters[parameters$epi_symbol == "beta_h", baseline_value] #transmission coefficient in the hospital
sigma_1a <- parameters[parameters$epi_symbol == "sigma_1a", baseline_value] #hospital admissions rate for the susceptible and colonised population in the community
sigma_1b <- parameters[parameters$epi_symbol == "sigma_1b", baseline_value] #hospital admissions rate for persons with BSI in community
gamma_1 <- parameters[parameters$epi_symbol == "gamma_1", baseline_value] #hospital discharge rate 
kappa_1 <- parameters[parameters$epi_symbol == "kappa_1", baseline_value] #the daily rate of progression from colonised state to BSI in the community
epsilon_1 <- parameters[parameters$epi_symbol == "epsilon_1", baseline_value] #natural clearance rate from colonisation (only in the community)
rho_1 <- parameters[parameters$epi_symbol == "rho_1", baseline_value]	#vaccine efficacy level
kappa_h <- parameters[parameters$epi_symbol == "kappa_h", baseline_value] #the daily rate of progression from colonisation to BSI in the hospital
tau_h <- parameters[parameters$epi_symbol == "tau_h", baseline_value] #the daily rate of recovery from BSI back to the colonisation state (in the hospital)
omega_h <- parameters[parameters$epi_symbol == "omega_h", baseline_value] #the daily mortality rate from BSI (in the hospital) #also used as a recruitment rate to maintain constant total population
vacc_coverage <- parameters[parameters$epi_symbol == "v_1", baseline_value] #the vaccination coverage in the community