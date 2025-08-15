# This script is for generating outcomes with baseline vaccination scenario (70% coverage and 70% efficacy)
#===================================================================================================
# Research Q:  Potential impacts of the hypothetical vaccine on the burden of K. pneumoniae BSI
# Modelling Q: How many cases and deaths would occur without and with intervention?
#===================================================================================================

# clear the environment
rm(list = ls())


#-------------------------------------------------------------------------------------------------#
# Generate outcomes for vaccination scenario #### 
#-------------------------------------------------------------------------------------------------#
# initial checking before applying the real scenario ####

# check what happens when the vaccine has full efficacy before proceeding to the baseline efficacy of 70%
# load the model code 
source("3_intervention_v4/intervention_inputs/model1_v4_for_intervention.R")

# produce output with vaccination
parameters1["vacc_coverage"] <- 0.7 # set vaccination rate to 0.7 (baseline coverage)
parameters1["rho_1"] <- 0.7 # set vaccination efficacy to 100% (full efficacy)

run <- ode(y = istate1, times = tps1, func = klebsiella_1, parms = parameters1, method = "lsoda", events = event_list)

# Load the model output
# Add columns containing the desired values
df.after.intervention2 <- as_tibble(as.data.frame(run)) %>%  
  mutate(colonised_h1 = C_h1 + VC_h1, #calculate the total infected individuals from the hospital
         incidence_col_hosp = c(0, diff(CInc_1a)),
         incidence_comm_cases = c(0, diff(CInc_1b)),
         bsi_hospitalisation = c(0, diff(CInc_1c)),
         hosp_onset_cases = c(0, diff(CInc_1d)), 
         total_incidence = (bsi_hospitalisation + hosp_onset_cases),
         bsi_associated_deaths = c(0, diff(D_h1))) %>% pivot_longer(names_to = "variable", cols = !1)

View(df.after.intervention2)


#---------------------------------------------------------#
# total number of cases and deaths for the whole England
#---------------------------------------------------------#
# Total NHS Trust
nhs_trust <- 134

# define time bound as we only investigate the impact within 5 years and the vaccination impact will happen at 365*80
start_time <- 365*0
end_time <- 365*5



#-------------------------#
# total number of cases
#-------------------------#

# Calculate the total number of bsi cases in the baseline 
total_bsi_cases_with_intervention2_hosp <- df.after.intervention2 %>% 
  filter(variable == "total_incidence") %>% #Filter solutions 
  filter(time>start_time, time<=end_time) %>% 
  mutate(year = floor((as.integer(time) - 1) / 365) + 1)

total_bsi_cases_with_intervention2_grouped <- total_bsi_cases_with_intervention2_hosp %>% 
  filter(year > 0) %>% # filter out the first year
  group_by(year) %>% 
  summarise(total_cases = sum(value)*nhs_trust) %>% 
  mutate(total_cases = round(total_cases)) #perform rounding to get integer

total_bsi_cases_with_intervention2_grouped

# total number of cases over 5 years without intervention at the baseline 2023-2024 figure would be 6,495
total_number_of_cases_5_years <- sum(total_bsi_cases_with_intervention2_grouped$total_cases)

total_number_of_cases_5_years

# save to csv for use later
write.csv(total_bsi_cases_with_intervention2_grouped, "3_intervention_v4/intervention_outputs/csv/cases_per_year_with_intervention2.csv")

write.csv(total_number_of_cases_5_years, "3_intervention_v4/intervention_outputs/csv/cases_5_years_with_intervention2.csv")

#-------------------------#
# total number of deaths
#-------------------------#

# Calculate the total number of deaths in the baseline 
total_deaths_with_intervention2_hosp <- df.after.intervention2 %>% 
  filter(variable == "bsi_associated_deaths") %>% #Filter solutions 
  filter(time>start_time, time<=end_time) %>% 
  mutate(year = floor((as.integer(time) - 1) / 365) + 1)

total_deaths_with_intervention2_grouped <- total_deaths_with_intervention2_hosp %>% 
  filter(year > 0) %>% # filter out the first year
  group_by(year) %>% 
  summarise(total_deaths = sum(value)*nhs_trust) %>% 
  mutate(total_deaths = round(total_deaths)) #perform rounding to get integer

total_deaths_with_intervention2_grouped

# total number of cases over 5 years without intervention at the baseline 2023-2024 figure would be 6,495
total_number_of_deaths_5_years <- sum(total_deaths_with_intervention2_grouped$total_deaths) 

total_number_of_deaths_5_years


# save to csv for use later
write.csv(total_deaths_with_intervention2_grouped, "3_intervention_v4/intervention_outputs/csv/deaths_per_year_with_intervention2.csv")

write.csv(total_number_of_deaths_5_years, "3_intervention_v4/intervention_outputs/csv/deaths_5_years_with_intervention2.csv")


#-----------------------------------------#
# bsi hospitalisations
#----------------------------------------#

# Calculate the total number of bsi cases in the baseline 
bsi_hosp_with_intervention2 <- df.after.intervention2 %>% 
  filter(variable == "bsi_hospitalisation") %>% #Filter solutions 
  filter(time>start_time, time<=end_time) %>% 
  mutate(year = floor((as.integer(time) - 1) / 365) + 1)

bsi_hosp_with_intervention2_grouped <- bsi_hosp_with_intervention2 %>% 
  filter(year > 0) %>% # filter out the first year
  group_by(year) %>% 
  summarise(total_bsi_hosp = sum(value)*nhs_trust) %>% 
  mutate(total_bsi_hosp = round(total_bsi_hosp)) #perform rounding to get integer

bsi_hosp_with_intervention2_grouped

# breakdown of the number of cases over 5 years without intervention at the baseline 2023-2024
bsi_hosp_5_years <- sum(bsi_hosp_with_intervention2_grouped$total_bsi_hosp)

bsi_hosp_5_years

# save to csv for use later
write.csv(bsi_hosp_with_intervention2_grouped, "3_intervention_v4/intervention_outputs/csv/bsi_hosp_per_year_with_intervention2.csv")

write.csv(bsi_hosp_5_years, "3_intervention_v4/intervention_outputs/csv/bsi_hosp_5_years_with_intervention2.csv")


#-----------------------------------------#
# hospital-onset cases
#----------------------------------------#

# Calculate the total number of bsi cases in the baseline 
hosp_onset_with_intervention2 <- df.after.intervention2 %>% 
  filter(variable == "hosp_onset_cases") %>% #Filter solutions 
  filter(time>start_time, time<=end_time) %>% 
  mutate(year = floor((as.integer(time) - 1) / 365) + 1)

hosp_onset_with_intervention2_grouped <- hosp_onset_with_intervention2 %>% 
  filter(year > 0) %>% # filter out the first year
  group_by(year) %>% 
  summarise(total_hosp_onset = sum(value)*nhs_trust) %>% 
  mutate(total_hosp_onset = round(total_hosp_onset)) #perform rounding to get integer

hosp_onset_with_intervention2_grouped

# breakdown of the number of cases over 5 years without intervention at the baseline 2023-2024
hosp_onset_5_years <- sum(hosp_onset_with_intervention2_grouped$total_hosp_onset)

hosp_onset_5_years

# save to csv for use later
write.csv(hosp_onset_with_intervention2_grouped, "3_intervention_v4/intervention_outputs/csv/hosp_onset_cases_per_year_with_intervention2.csv")

write.csv(hosp_onset_5_years, "3_intervention_v4/intervention_outputs/csv/hosp_onset_cases_5_years_with_intervention2.csv")



#-------------------------#
# cumulative outcomes
#-------------------------#

# year by year

cum_cases_per_year <- read.csv("3_intervention_v4/intervention_outputs/csv/cases_per_year_with_intervention2.csv")

cum_deaths_per_year <- read.csv("3_intervention_v4/intervention_outputs/csv/deaths_per_year_with_intervention2.csv")

cum_bsi_hosp_per_year <- read.csv("3_intervention_v4/intervention_outputs/csv/bsi_hosp_per_year_with_intervention2.csv")

cum_hosp_onset_per_year <- read.csv("3_intervention_v4/intervention_outputs/csv/hosp_onset_cases_per_year_with_intervention2.csv")

# cumulative
cum_cases_per_year$cum_cases <- cumsum(cum_cases_per_year$total_cases)

cum_deaths_per_year$cum_deaths <- cumsum(cum_deaths_per_year$total_deaths)

cum_bsi_hosp_per_year$cum_bsi_hosp <- cumsum(cum_bsi_hosp_per_year$total_bsi_hosp)

cum_hosp_onset_per_year$cum_hosp_onset <- cumsum(cum_hosp_onset_per_year$total_hosp_onset)


# save to csv for use later
write.csv(cum_cases_per_year, "3_intervention_v4/intervention_outputs/csv/cases_per_year_with_intervention2.csv")

write.csv(cum_deaths_per_year, "3_intervention_v4/intervention_outputs/csv/deaths_per_year_with_intervention2.csv")

write.csv(cum_bsi_hosp_per_year, "3_intervention_v4/intervention_outputs/csv/bsi_hosp_per_year_with_intervention2.csv")

write.csv(cum_hosp_onset_per_year, "3_intervention_v4/intervention_outputs/csv/hosp_onset_cases_per_year_with_intervention2.csv")



# all years
cum_cases_5_years <- read.csv("3_intervention_v4/intervention_outputs/csv/cases_5_years_with_intervention2.csv")

cum_deaths_5_years <- read.csv("3_intervention_v4/intervention_outputs/csv/deaths_5_years_with_intervention2.csv")

cum_bsi_hosp_5_years <- read.csv("3_intervention_v4/intervention_outputs/csv/bsi_hosp_5_years_with_intervention2.csv")

cum_hosp_onset_5_years <- read.csv("3_intervention_v4/intervention_outputs/csv/hosp_onset_cases_5_years_with_intervention2.csv")


names(cum_cases_5_years)[1] <- "number_of_cumulative" 
names(cum_cases_5_years)[2] <- "total_cases_5_years" 

names(cum_deaths_5_years)[1] <- "number_of_cumulative" 
names(cum_deaths_5_years)[2] <- "total_deaths_5_years" 


names(cum_bsi_hosp_5_years)[1] <- "number_of_cumulative" 
names(cum_bsi_hosp_5_years)[2] <- "total_bsi_hosp_5_years" 

names(cum_hosp_onset_5_years)[1] <- "number_of_cumulative" 
names(cum_hosp_onset_5_years)[2] <- "total_hosp_onset_5_years" 


# save to csv for use later
write.csv(cum_cases_5_years, "3_intervention_v4/intervention_outputs/csv/cases_5_years_with_intervention2.csv")

write.csv(cum_deaths_5_years, "3_intervention_v4/intervention_outputs/csv/deaths_5_years_with_intervention2.csv")


write.csv(cum_bsi_hosp_5_years, "3_intervention_v4/intervention_outputs/csv/bsi_hosp_5_years_with_intervention2.csv")

write.csv(cum_hosp_onset_5_years, "3_intervention_v4/intervention_outputs/csv/hosp_onset_cases_5_years_with_intervention2.csv")
