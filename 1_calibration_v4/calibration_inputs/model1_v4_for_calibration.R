# this script is the model chunk for calibration
# define ordinary differential equations here
# treated as input for "calibration function chunk" script

#===============================================================================
#### 1. Setup differential equations ####
#### all 65+ year-olds #### 
#===============================================================================
klebsiella_1 <- function(time, state, parameters){
  with(as.list(c(state, parameters)), {
    
    #-----------------------#
    #calculating populations#
    #-----------------------#
    
    # total population in the hospital
    P_hosp_c1 <- S_h1 + C_h1 + I_h1  #define total 65+ year-olds population 
    
    # total population outside the hospital (within the community)
    P_community <- S_c1 + C_c1 + I_c1  #define total 65+ year-olds population in the community
    
    # total older adults population
    P1 <- P_hosp_c1 + P_community #define total 65+ year-olds population
    
    # infectious populations     
    I_hosp_c1 <- C_h1 + I_h1 #define size of infectious population for the entire older adults population in the hospital
    
    #------------------------------#
    #calculating force of infection#
    #------------------------------#
    
    # force of infection in the hospital
    lambda_h <- beta_h*I_hosp_c1/P_hosp_c1 #define lambda for the entire older adults population within the hospital
    
    #--------------------------------------------#
    #maintaining constant hospital population#
    #--------------------------------------------#
    
    # this logic will also allow the system to admit patients according to our original admissions rate but also self-adjust the final admissions against the available beds
    # thus, if the available beds turn out to be less than the amount required to absorb admissions from colonised and susceptibles, the number of final admissions will be reduced, and vice-versa
    # the implication is also that only the composition of patients will change (the proportion of susceptibles vs. colonised vs. BSI patients) but the number of total patients in the hospital will roughly be constant
    # the main purpose of this logic is also to see more clearly the impact of the vaccination we would implement later
    
    #------------------------------------------------------------#
    #hospital admissions will reduce the number of available beds#
    #------------------------------------------------------------#
    
    # prioritised entry for the BSI patients as they are severely ill
    entry_bsi <- sigma_1b*I_c1
    
    # entry of colonised population that are as "random" (follow regular admissions rate) as the susceptible population as all of them do not have symptoms
    entry_colonised <- sigma_1a*C_c1
    
    # entry of colonised population that are "random" (follow regular admissions rate) as all of them do not have symptoms
    entry_susceptibles <- sigma_1a*S_c1
    
    #-------------------------------------------------------------------------#
    #hospital discharges and deaths will increase the number of available beds#
    #-------------------------------------------------------------------------#
    
    # discharges of colonised patients from the hospital follows the average length of stay
    discharge_colonised <- gamma_1*C_h1
    
    # discharges of susceptible patients from the hospital follows the average length of stay
    discharge_susceptibles <- gamma_1*S_h1 
    
    # deaths within the hospital that will free up beds
    deaths_hosp <- mu_1*P_hosp_c1 + omega_h*I_h1 # exits due to deaths will free up beds
    
    
    # total remaining "available beds" we must adjust for is the number of freed up beds capacity from discharges and deaths minus the number of regular admissions 
    avail_beds <- discharge_colonised + discharge_susceptibles + deaths_hosp - entry_bsi - entry_colonised - entry_susceptibles  #total net hospital exit #minus the difference between entry and exit of colonised and susceptibles so that when the entry is higher than exit, the available beds are reduced, and vice-versa
    
    # the remaining "available beds" will only be distributed among non-BSI patients as the BSI patients are ensured to be admitted as they are severely ill
    non_bsi_entry_pool <- S_c1 + C_c1#non-bsi entry pool in the community
    
    
    #----------------------------------------------------#
    #topup admissions to maintain fixed hospital capacity#
    #----------------------------------------------------#
    
    topup_Sc1 <- S_c1/non_bsi_entry_pool*avail_beds #proportion of susceptible in the community
    
    topup_Cc1 <- C_c1/non_bsi_entry_pool*avail_beds #proportion of colonised in the community
    

    #-----------------------------------------------------------------------------------------------------#
    #all 65+ year-olds population differential equations # make sure the inflows and outflows are balanced#
    #-----------------------------------------------------------------------------------------------------#
    
    #------------------#
    #community dynamics
    #------------------#
    
    #community unvaccinated
    dS_c1 <- omega_h*I_h1 + mu_1*P1 + gamma_1*S_h1 + epsilon_1*C_c1 - sigma_1a*S_c1 - topup_Sc1 - mu_1*S_c1 
    dC_c1 <- gamma_1*C_h1 - kappa_1*C_c1 - epsilon_1*C_c1 - sigma_1a*C_c1 - topup_Cc1 - mu_1*C_c1 
    dI_c1 <- kappa_1*C_c1 - sigma_1b*I_c1 + - mu_1*I_c1 
    
    #------------------#
    #hospital dynamics
    #------------------#
    
    #hospital unvaccinated
    dS_h1 <- sigma_1a*S_c1 + topup_Sc1 - gamma_1*S_h1 - lambda_h*S_h1 - mu_1*S_h1 
    dC_h1 <- sigma_1a*C_c1 + topup_Cc1 + lambda_h*S_h1 + tau_h*I_h1 - gamma_1*C_h1 - kappa_h*C_h1 - mu_1*C_h1 
    dI_h1 <- kappa_h*C_h1 + sigma_1b*I_c1 - tau_h*I_h1 - omega_h*I_h1 - mu_1*I_h1 
    
    #----------------------------#
    #deaths and incidence counter#
    #----------------------------#
    
    # deaths counter instead of including a death as a main compartment to simplify the equations
    # the death counter is only for the hospital population as we assume that all deaths happen in the hospital since the BSI patients will be immediately hospitalised
    # the deaths rates for for calculating BSI-associated deaths include both death rate attributed to the BSI and natural death rate for the BSI compartment
    dD_h1 <- (omega_h + mu_1)*I_h1  
    
    # incidence counter for community and hospital
    dCInc_1a <- lambda_h*S_h1 # incidence of colonisation in the hospital
    dCInc_1b <- kappa_1*C_c1  # progression from colonised state to BSI state in the community 
    dCInc_1c <- sigma_1b*I_c1 # BSI hospitalisations from the community 
    dCInc_1d <- kappa_h*C_h1 # hospital-onset cases 
    
    
    # list down the ode #MAKE SURE the order is the same as the order of the equations above!
    list(c(S_c1 = dS_c1, C_c1 = dC_c1, I_c1 = dI_c1, 
           S_h1 = dS_h1, C_h1 = dC_h1, I_h1 = dI_h1,
           D_h1 = dD_h1, CInc_1a = dCInc_1a, CInc_1b = dCInc_1b, CInc_1c = dCInc_1c, CInc_1d = dCInc_1d)) 
  })
}


#==================================================================================================
#load initial values and parameters values sub script
source("1_calibration_v4/calibration_inputs/inputs_v4_calibration_model1.R")$value
#==================================================================================================

#==================================================================================================
#### 2. Setup initial values #### #ensure the order is the same as the order of the equations above
#==================================================================================================
#istate1 <- c(S_c1 = S_c1, C_c1 = C_c1, I_c1 = I_c1,
             #S_h1 = S_h1, C_h1 = C_h1, I_h1 = I_h1, 
             #D_h1 = D_h1, CInc_1a = CInc_1a, CInc_1b = CInc_1b, CInc_1c = CInc_1c, CInc_1d = CInc_1d) 

#===============================================================
#### 3. Setup parameters ####
#===============================================================
#parameters1 <- c(mu_1 = mu_1, #recruitment rate and natural mortality rate
                 #sigma_1a = sigma_1a, #hospital admissions rate for the susceptible and colonised population in the community
                 #sigma_1b = sigma_1b, #hospital admissions rate for persons with BSI in community
                 #epsilon_1 = epsilon_1, #natural clearance rate from colonisation (only in the community)
                 #gamma_1 = gamma_1, #hospital discharge rate 
                 #kappa_1 =  kappa_1, #kappa_1, #the daily rate of progression from colonised state to BSI in the community
                 #kappa_h =  kappa_h,#kappa_h,  #the daily rate of progression from colonisation to BSI in the hospital
                 #tau_h = tau_h, #the daily rate of recovery from BSI back to the colonisation state (in the hospital)
                 #omega_h = omega_h, #the daily mortality rate from BSI (in the hospital) #also used as a recruitment rate to maintain constant total population
                 #beta_h =  beta_h) #transmission coefficient in the hospital

#===============================================================
#### 4. Setup timesteps ####
#===============================================================
#start1 <- 0 #start time for the simulation
#end1 <- 365*200 #end point of the timestep #for calibration, the timestep is set up to be arbitrarily long to ensure the system reaches equilibrium before we move on to the next step
#tps1 <- seq(start1, end1, by = 1) #time steps for the simulation

#===============================================================
#### 5. Run the model and check the outputs ####
#===============================================================
#run the model 
#we are using the euler method such that the timesteps are read as integers
#out1 <- ode(y = istate1, times = tps1, func = klebsiella_1, parms = parameters1, method = "lsoda") 

# check outputs

#total_cases_eq <- (tail(out1[, "CInc_1c"], 1) - tail(out1[, "CInc_1c"], 2)[1] + tail(out1[, "CInc_1d"], 1) - tail(out1[, "CInc_1d"], 2)[1])* 365 * 134

#hosp_cases_pred <- (tail(out1[, "CInc_1d"], 1) - tail(out1[, "CInc_1d"], 2)[1]) * 365 * 134 #total hospital cases per year

#bsi_hospitalisation <- (tail(out1[, "CInc_1c"], 1) - tail(out1[, "CInc_1c"], 2)[1]) * 365 * 134 #total BSI hospitalisation per year

#col_comm_pred <- round(tail(out1[, "C_c1"], 1) / (tail(out1[, "S_c1"], 1) + tail(out1[, "C_c1"], 1) + tail(out1[, "I_c1"], 1)),4)*100

# check if the outputs are at equilibrium
#if (tail(out1[, "S_c1"], 1) - tail(out1[, "S_c1"], 2)[1] == 0.0 & 
    #tail(out1[, "C_c1"], 1) - tail(out1[, "C_c1"], 2)[1] == 0.0 & 
    #tail(out1[, "I_c1"], 1) - tail(out1[, "I_c1"], 2)[1] == 0.0) {
  #print("The model has reached equilibrium.")
#} else {
  #print("The model has not reached equilibrium.")
#}

#===============================================================
#### 6.Save the pre-calibration model outputs ####
#===============================================================
# save the output to a csv file # make sure the output is in the same order as the equations above
# make sure the outputs without vaccination is the same as the output of the calibration at equilibrium
#write.csv(out1, file = "1_calibration_v4/calibration_outputs/csv/klebsiella_1_.csv") #save the output to a csv file