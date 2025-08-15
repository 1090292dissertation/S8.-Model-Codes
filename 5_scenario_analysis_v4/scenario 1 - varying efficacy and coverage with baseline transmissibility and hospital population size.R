# Visualise dynamics overtime without and with vaccination

# clear the environment
rm(list = ls())

# community colonisation, hospital colonisation, and hospital BSI without vaccination
# multiply up by 134 NHS trust to get the total number for the whole England
#-------------------------------------------------------------------------------------------------#
# Generate outcomes for low prevalence of colonisation in the community scenario #### 
#-------------------------------------------------------------------------------------------------#
# load the model code 
source("5_scenario_analysis_v4/scenario_inputs/model1_v4_for_scenario analysis.R")

#proportion to colonised is also lower 
#parameters1["p"] <- 0.0013
#parameters1["beta_h"] <- 0.10991111

# define time bound as we only investigate the impact within 5 years and the vaccination impact will happen at 365*80
#===============================================================
#### 4. Setup the event list ####
#===============================================================

event_list <- list(func = vaccination_event,
                   time = c(0)) #time at which the vaccination event occurs -> the solver will implement it at t+1 

#===============================================================

#===============================================================
#### 4. Setup timesteps ####
#===============================================================
start1 <- 0 #start time for the simulation
end1 <- 365*85 #end point of the timestep #ensure the system reach equilibrium before the intervention is run
tps1 <- seq(start1, end1, by = 1) #time steps for the simulation

#===============================================================
#### 5. Run the model and check the outputs ####
#===============================================================
#run the model 
#we are using the euler method such that the timesteps are read as integers
out1 <- ode(y = istate1, times = tps1, func = klebsiella_1, parms = parameters1, method = "lsoda", events = event_list)

# extract equilibrium values
eq_values <- tail(out1, 1)

# save as csv
write.csv(eq_values, file = "5_scenario_analysis_v4/scenario_outputs/csv/eq_values_low_transmissibility.csv")

#====================================================================================
# discard non-equilibrium values to implement vaccination with shorter running time
#====================================================================================
# initial values at equilibrium
eq_values <- read.csv("5_scenario_analysis_v4/scenario_outputs/csv/eq_values_low_transmissibility.csv")

# use equilibrium values to shorten running time for intervention
istate1["S_h1"] <- eq_values$S_h1
istate1["C_h1"] <- eq_values$C_h1 
istate1["I_h1"] <- eq_values$I_h1

istate1["S_c1"] <- eq_values$S_c1
istate1["C_c1"] <- eq_values$C_c1
istate1["I_c1"] <- eq_values$I_c1

#===============================================================
#### 4. Setup the event list ####
#===============================================================

event_list <- list(func = vaccination_event,
                   time = c(0)) #time at which the vaccination event occurs -> the solver will implement it at t+1 

#===============================================================

#===============================================================
#### 4. Setup timesteps ####
#===============================================================
start1 <- 0 #start time for the simulation
end1 <- 365*5 #end point of the timestep #ensure the system reach equilibrium before the intervention is run
tps1 <- seq(start1, end1, by = 1) #time steps for the simulation

# check run outputs
out1 <- ode(y = istate1, times = tps1, func = klebsiella_1, parms = parameters1, method = "lsoda", events = event_list)


#===============================================================
#### 5. Run the simulation ####
#===============================================================
# varying vaccination coverage and efficacy
sim1 <- seq(0, 1, by = 0.1) # define the vaccination coverage 
sim2 <- seq(0, 1, by = 0.1) # define the vaccination efficacy

# create grid to check different combinations
grid <- expand.grid(vacc_coverage = sim1, rho_1 = sim2)

# Load the model output
# Add columns containing the desired values
df.scen1<-NULL # warning! the data frame stores the results so the next outputs will be amended based on the last value stored
for(i in 1:nrow(grid)){
  parameters1["vacc_coverage"] <- grid$vacc_coverage[i];
  parameters1["rho_1"] <- grid$rho_1[i]; # set the vaccination efficacy to the current value in the loop)
  #parameters1["beta_h"] <- 0.10991111
  run <- ode(y = istate1, times = tps1, func = klebsiella_1, parms = parameters1, method = "lsoda", events = event_list);
  tmp <- as_tibble(as.data.frame(run)) %>%  
    mutate(colonised_h1 = C_h1 + VC_h1, #calculate the total infected individuals from the hospital
           incidence_col_hosp = c(0, diff(CInc_1a)),
           incidence_comm_cases = c(0, diff(CInc_1b)),
           bsi_hospitalisation = c(0, diff(CInc_1c)),
           incidence_hosp_cases = c(0, diff(CInc_1d)), 
           incidence_bsi_hosp = (bsi_hospitalisation + incidence_hosp_cases),
           bsi_associated_deaths = c(0, diff(D_h1)),
           vacc_coverage =  grid$vacc_coverage[i],
           vacc_efficacy = grid$rho_1[i],
           #community_prevalence_aux =  grid$p[i],
           community_prevalence = (C_c1 + VC_c1)/(S_c1 + C_c1 + VS_c1 + VC_c1 + I_c1),
           vacc_cov_percent = paste0(grid$vacc_coverage[i]*100, "%"),
           vacc_efficacy_percent = paste0(grid$rho_1[i]*100, "%")) %>% 
    select(everything())
  df.scen1<-bind_rows(df.scen1, tmp) 
}
View(df.scen1)

View(df.scen1 %>% filter(time > 1823))

# add the number of vaccine column 
# get the number of vaccines
num_vaccines <- df.scen1 %>% filter(time == 1)

View(num_vaccines)

df.scen1_long <- df.scen1 %>% 
  pivot_longer(cols = -c(time, vacc_coverage, vacc_efficacy, vacc_cov_percent, vacc_efficacy_percent, community_prevalence),
               names_to = "variable", 
               values_to = "value")
# save as RDS
#saveRDS(df.scen1_join, file = "3_intervention_v4/intervention_outputs/df.scen1.rds")

#saveRDS(df.scen1_long, file = "3_intervention_v4/intervention_outputs/df.scen1_long.rds")

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
total_bsi_cases_with_scen1_hosp <- df.scen1_long %>% 
  filter(variable == "incidence_bsi_hosp") %>% #Filter solutions 
  filter(time>start_time, time<=end_time) %>% 
  mutate(year = floor((as.integer(time) - 1) / 365) + 1)

total_bsi_cases_with_scen1_grouped <- total_bsi_cases_with_scen1_hosp %>% 
  filter(year > 0) %>% # filter out the first year
  group_by(year, vacc_coverage, vacc_efficacy) %>% 
  summarise(total_cases = sum(value)*nhs_trust) %>% 
  mutate(total_cases = round(total_cases)) #perform rounding to get integer

total_bsi_cases_with_scen1_grouped

View(total_bsi_cases_with_scen1_grouped)

# total number of cases over 5 years 
total_number_of_cases_5_years <- total_bsi_cases_with_scen1_grouped %>% 
  group_by(vacc_coverage, vacc_efficacy) %>% 
  summarise(total_cases = sum(total_cases))

View(total_number_of_cases_5_years)

# save to csv for use later
write.csv(total_bsi_cases_with_scen1_grouped, "5_scenario_analysis_v4/scenario_outputs/csv/cases_per_year_with_scen1.csv")

write.csv(total_number_of_cases_5_years, "5_scenario_analysis_v4/scenario_outputs/csv/cases_5_years_with_scen1.csv")

#-------------------------#
# cumulative outcomes
#-------------------------#

# all years
cum_cases_5_years <- read.csv("5_scenario_analysis_v4/scenario_outputs/csv/cases_per_year_with_scen1.csv")

View(cum_cases_5_years)



# remove redundant combinations
cum_cases_5_years_no_int_clean <- cum_cases_5_years  %>% filter(vacc_coverage == 0, vacc_efficacy == 0)



# remove redundant combinations
cum_cases_5_years_with_int_clean <- cum_cases_5_years  %>% filter(vacc_coverage != 0, vacc_efficacy != 0)

View(cum_cases_5_years_with_int_clean)

# baseline cases per year without intervention
baseline_cases <- cum_cases_5_years_no_int_clean$total_cases[1]

# insert the number of cases without intervention to each vaccine coverage and efficacy scenario
cases_with_intervention_varying_coverage_efficacy <- cum_cases_5_years_with_int_clean %>% 
  mutate(no_int_cases = baseline_cases)


# calculate the outcome
# calculate baseline cumulative deaths for each coverage and efficacy combination
cum_cases_per_year <- cases_with_intervention_varying_coverage_efficacy %>% 
  arrange(vacc_coverage, vacc_efficacy, year) %>% 
  group_by(vacc_coverage, vacc_efficacy) %>% 
  mutate(no_int_cases_cum = cumsum(no_int_cases),
         cum_cases = cumsum(total_cases)) %>% 
  ungroup()



cum_cases_averted <- data.frame(year = cum_cases_per_year$year, #so the year starts at 1
                                vacc_coverage = cum_cases_per_year$vacc_coverage, 
                                vacc_efficacy = cum_cases_per_year$vacc_efficacy, 
                                no_int_cases_cum = cum_cases_per_year$no_int_cases_cum,
                                with_int_baseline_cases_cum = cum_cases_per_year$cum_cases)


View(cum_cases_averted)

# save to csv for use later
write.csv(cum_cases_averted, "5_scenario_analysis_v4/scenario_outputs/csv/cases_5_years_with_scen1.csv")


## calculate the total number of cases averted
cases_averted_scenario1 <- cum_cases_averted %>% 
  mutate(cases_averted_cum = abs(with_int_baseline_cases_cum - no_int_cases_cum))

View(cases_averted_scenario1)

# select the total number of deaths averted over 5 years
total_cases_averted_scenario1 <- cases_averted_scenario1 %>% filter(year == 5, vacc_coverage != 0, vacc_efficacy !=0)

total_cases_averted_scenario1$percent_reduction_cases <- 
  round((total_cases_averted_scenario1$cases_averted_cum/total_cases_averted_scenario1$no_int_cases_cum),3)*100

View(total_cases_averted_scenario1)

total_cases_averted_scenario1 <- total_cases_averted_scenario1 %>% mutate(vacc_cov_percent = vacc_coverage*100,
                                                                          vacc_eff_percent = vacc_efficacy*100)

View(total_cases_averted_scenario1)


# Calculate the total number of vaccines
num_vaccines$total_num_vaccination <- round((num_vaccines$VS_c1 + num_vaccines$VC_c1 + num_vaccines$VS_h1 + num_vaccines$VC_h1)*nhs_trust,0)

View(num_vaccines)

# create a lookup table to assign the number of vaccination per coverage and efficacy
lookup_num_vaccines <- num_vaccines %>% select(c(vacc_coverage, total_num_vaccination)) 

# filter out zero coverage
lookup_num_vaccines  <- lookup_num_vaccines %>% filter(vacc_coverage != 0) %>% distinct() 

lookup_num_vaccines$vacc_cov_percent<- round(lookup_num_vaccines$vacc_coverage*100,0)

View(lookup_num_vaccines)

# Join the table 

total_cases_averted_scenario1
View(total_cases_averted_scenario1)

total_cases_averted_scenario1_join <- total_cases_averted_scenario1 %>% 
  left_join(lookup_num_vaccines, by = "vacc_cov_percent")

View(total_cases_averted_scenario1_join)

total_cases_averted_scenario1_join$cases_averted_per_vaccine <- 
  round((total_cases_averted_scenario1_join$cases_averted_cum/total_cases_averted_scenario1_join$total_num_vaccination),4)

total_cases_averted_scenario1_join$cases_averted_per_1000_vaccine <- 
  total_cases_averted_scenario1_join$cases_averted_per_vaccine*1000

total_cases_averted_scenario1_join$num_vacc_per_one_case_averted <- 
  round((total_cases_averted_scenario1_join$total_num_vaccination/total_cases_averted_scenario1_join$cases_averted_cum),0)

View(total_cases_averted_scenario1_join)

total_cases_averted_scenario1_join2 <- total_cases_averted_scenario1_join #%>% filter(vacc_cov_percent %in% c(10, 40, 70, 100), vacc_eff_percent %in% c(10, 40, 70, 100))

View(total_cases_averted_scenario1_join2)

percent_reduction_quantiles <- quantile(total_cases_averted_scenario1_join2$percent_reduction_cases, probs = c(0, 0.25, 0.75, 1))


# plot the outcomes
plot_cases_scenario1_absolute_reduction <- ggplot(total_cases_averted_scenario1_join,
                                                 aes(x = vacc_cov_percent, y = vacc_eff_percent, fill = cases_averted_cum)) + 
  geom_tile(color = "white", linewidth = 0.5) + 
  #geom_text(aes(label = cases_averted_cum), color = "black", size = 4) +
  geom_text(data = filter(total_cases_averted_scenario1_join, 
                          vacc_cov_percent %in% c(10, 40, 70, 100),
                          vacc_eff_percent %in% c(10, 40, 70, 100)), 
            aes(label = round(cases_averted_cum, 1)), size = 6) +
  scale_fill_gradientn(colors = c(
    "#FDE7BD",  # soft yellow
    "#F8B3D0",  # pastel pink
    "#C19BD3",  # muted lavender
    "#8676B7"),   # soft purple (low saturation)
    limits = c(0, 30000)
  )+
  scale_x_continuous(breaks = sort(unique(total_cases_averted_scenario1_join$vacc_cov_percent)),
                     labels = paste0(sort(unique(total_cases_averted_scenario1_join$vacc_cov_percent)))) +
  
  scale_y_continuous(breaks = sort(unique(total_cases_averted_scenario1_join$vacc_eff_percent)),
                     labels = paste0(sort(unique(total_cases_averted_scenario1_join$vacc_eff_percent)))) +
  labs(
    title = "Scenario 1 (varying coverage and efficacy only)\nNumber of cases averted\n(over 5 years)",
    x = "Vaccine coverage (%)",
    y = "Vaccine efficacy (%)",
    fill = NULL
  ) + 
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18),
    legend.title = element_blank(),
    legend.text = element_text(size = 24),
    legend.key.height = unit(2.7, "cm"),  # increase legend height
    legend.key.width  = unit(0.4, "cm"), # decrease legend width
    legend.margin = margin(0, 0, 0, 0),   # shrink inside legend padding
    legend.box.margin = margin(0, 0, 0, -5) # pull legend closer to plot
  )

ggsave(plot_cases_scenario1_absolute_reduction, file = "5_scenario_analysis_v4/scenario_outputs/plots/scen1/plot_cases_scenario1c_absolute_reduction.png", width = 7, height = 7)


plot_cases_scenario1_percent_reduction <- ggplot(total_cases_averted_scenario1_join,
                                                 aes(x = vacc_cov_percent, y = vacc_eff_percent, fill = percent_reduction_cases)) + 
  geom_tile(color = "white", linewidth = 0.5) + 
  #geom_text(aes(label = percent_reduction_cases), color = "black", size = 4) +
  geom_text(data = filter(total_cases_averted_scenario1_join, 
                          vacc_cov_percent %in% c(10, 40, 70, 100),
                          vacc_eff_percent %in% c(10, 40, 70, 100)), 
            aes(label = round(percent_reduction_cases, 1)), size = 6) +
  scale_fill_gradientn(colors = c(
    "#FDE7BD",  # soft yellow
    "#F8B3D0",  # pastel pink
    "#C19BD3",  # muted lavender
    "#8676B7"),   # soft purple (low saturation)
    limits = c(0, 100)
  )+
  scale_x_continuous(breaks = sort(unique(total_cases_averted_scenario1_join$vacc_cov_percent)),
                     labels = paste0(sort(unique(total_cases_averted_scenario1_join$vacc_cov_percent)))) +
  
  scale_y_continuous(breaks = sort(unique(total_cases_averted_scenario1_join$vacc_eff_percent)),
                     labels = paste0(sort(unique(total_cases_averted_scenario1_join$vacc_eff_percent)))) +
  labs(
    title = "Scenario 1 (varying coverage and efficacy only)\nPercentage of cases averted\n(over 5 years)",
    x = "Vaccine coverage (%)",
    y = "Vaccine efficacy (%)",
    fill = NULL
  ) + 
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18),
    legend.title = element_blank(),
    legend.text = element_text(size = 24),
    legend.key.height = unit(2.7, "cm"),  # increase legend height
    legend.key.width  = unit(0.4, "cm"), # decrease legend width
    legend.margin = margin(0, 0, 0, 0),   # shrink inside legend padding
    legend.box.margin = margin(0, 0, 0, -5) # pull legend closer to plot
  )

ggsave(plot_cases_scenario1_percent_reduction, file = "5_scenario_analysis_v4/scenario_outputs/plots/scen1/plot_cases_scenario1c_percent_reduction.png", width = 7, height = 7)



plot_cases_scenario1_efficiency1 <- ggplot(total_cases_averted_scenario1_join,
                                           aes(x = vacc_cov_percent, y = vacc_eff_percent, fill = cases_averted_per_1000_vaccine)) + 
  geom_tile(color = "white", linewidth = 0.5) + 
  #geom_text(aes(label = cases_averted_per_1000_vaccine), color = "black", size = 4) +
  geom_text(data = filter(total_cases_averted_scenario1_join, 
                          vacc_cov_percent %in% c(10, 40, 70, 100),
                          vacc_eff_percent %in% c(10, 40, 70, 100)), 
            aes(label = round(cases_averted_per_1000_vaccine, 2)), size = 6) +
  scale_fill_gradientn(colors = c(
    "#FDE7BD",  # soft yellow
    "#F8B3D0",  # pastel pink
    "#C19BD3",  # muted lavender
    "#8676B7"),   # soft purple (low saturation)
    limits = c(0, 40)
  )+
  scale_x_continuous(breaks = sort(unique(total_cases_averted_scenario1_join$vacc_cov_percent)),
                     labels = paste0(sort(unique(total_cases_averted_scenario1_join$vacc_cov_percent)))) +
  
  scale_y_continuous(breaks = sort(unique(total_cases_averted_scenario1_join$vacc_eff_percent)),
                     labels = paste0(sort(unique(total_cases_averted_scenario1_join$vacc_eff_percent)))) +
  labs(
    title = "Scenario 1 (varying coverage and efficacy only)\nCases averted per 1000 vaccinated persons\n(over 5 years)",
    x = "Vaccine coverage (%)",
    y = "Vaccine efficacy (%)",
    fill = NULL
  ) + 
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18),
    legend.title = element_blank(),
    legend.text = element_text(size = 24),
    legend.key.height = unit(2.7, "cm"),  # increase legend height
    legend.key.width  = unit(0.4, "cm"), # decrease legend width
    legend.margin = margin(0, 0, 0, 0),   # shrink inside legend padding
    legend.box.margin = margin(0, 0, 0, -5) # pull legend closer to plot
  )

ggsave(plot_cases_scenario1_efficiency1, file = "5_scenario_analysis_v4/scenario_outputs/plots/scen1/plot_cases_scenario1c_efficiency1.png", width = 7, height = 7)


plot_cases_scenario1_efficiency2 <- ggplot(total_cases_averted_scenario1_join,
                                           aes(x = vacc_cov_percent, y = vacc_eff_percent, fill = num_vacc_per_one_case_averted)) + 
  geom_tile(color = "white", linewidth = 0.5) + 
  #geom_text(aes(label = num_vacc_per_one_case_averted), color = "black", size = 4) +
  geom_text(data = filter(total_cases_averted_scenario1_join, 
                          vacc_cov_percent %in% c(10, 40, 70, 100),
                          vacc_eff_percent %in% c(10, 40, 70, 100)), 
            aes(label = round(num_vacc_per_one_case_averted, 2)), size = 6) +
  scale_fill_gradientn(colors = c(
    "#FDE7BD",  # soft yellow
    "#F8B3D0",  # pastel pink
    "#C19BD3",  # muted lavender
    "#8676B7"),   # soft purple (low saturation)
    limits = c(0, 8000)
  )+
  scale_x_continuous(breaks = sort(unique(total_cases_averted_scenario1_join$vacc_cov_percent)),
                     labels = paste0(sort(unique(total_cases_averted_scenario1_join$vacc_cov_percent)))) +
  
  scale_y_continuous(breaks = sort(unique(total_cases_averted_scenario1_join$vacc_eff_percent)),
                     labels = paste0(sort(unique(total_cases_averted_scenario1_join$vacc_eff_percent)))) +
  labs(
    title = "Scenario 1 (varying coverage and efficacy only)\nNumber of vaccinated persons per one case averted\n(over 5 years)",
    x = "Vaccine coverage (%)",
    y = "Vaccine efficacy (%)",
    fill = NULL
  ) + 
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18),
    legend.title = element_blank(),
    legend.text = element_text(size = 24),
    legend.key.height = unit(2.7, "cm"),  # increase legend height
    legend.key.width  = unit(0.4, "cm"), # decrease legend width
    legend.margin = margin(0, 0, 0, 0),   # shrink inside legend padding
    legend.box.margin = margin(0, 0, 0, -5) # pull legend closer to plot
  )

ggsave(plot_cases_scenario1_efficiency2, file = "5_scenario_analysis_v4/scenario_outputs/plots/scen1/plot_cases_scenario1c_efficiency2.png", width = 7, height = 7)

