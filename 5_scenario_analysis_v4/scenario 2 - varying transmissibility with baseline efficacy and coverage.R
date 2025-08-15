# Visualise dynamics under varying transmissibility

# clear the environment
rm(list = ls())

# community colonisation, hospital colonisation, and hospital BSI without vaccination
# multiply up by 134 NHS trust to get the total number for the whole England
#-------------------------------------------------------------------------------------------------#
# Generate outcomes for low prevalence of colonisation in the community scenario #### 
#-------------------------------------------------------------------------------------------------#
# load the model code 
# use the same model code as the one we used for sensitivity analysis 
source("5_scenario_analysis_v4/scenario_inputs/model1_v4_for_scenario analysis.R")
# use initial values that are the same as the initial values during calibration so the curve displayed includes settling in period

#===============================================================
#### 4. Setup the event list ####
#===============================================================

event_list <- list(func = vaccination_event,
                   time = c(365*80-1)) #time at which the vaccination event occurs -> the solver will implement it at t+1 

#===============================================================

#===============================================================
#### 4. Setup timesteps for scenario analysis ####
#===============================================================
start1 <- 0 #start time for the simulation
end1 <- 365*85 #end point of the timestep #ensure the system reach equilibrium before the intervention is run
tps1 <- seq(start1, end1, by = 1) #time steps for the simulation

# check run outputs
#out1 <- ode(y = istate1, times = tps1, func = klebsiella_1, parms = parameters1, method = "lsoda", events = event_list)


#-------------------------------------------------------------------------------------------------#
# Generate outcomes for vaccination scenario #### 
#-------------------------------------------------------------------------------------------------#
# load the model code 
source("5_scenario_analysis_v4/scenario_inputs/model1_v4_for_scenario analysis.R")

sim1 <- seq(0, 0.7, length.out = 2) # define the vaccination coverage 
sim2 <- seq(0, 0.7, length.out = 2) # define the vaccination efficacy
sim3 <- seq(min_beta_h, max_beta_h, length.out = 10) # define the transmissibility range

# create grid to check different combinations
grid <- expand.grid(vacc_coverage = sim1, rho_1 = sim2, beta_h = sim3)

# Load the model output
# Add columns containing the desired values
df.scen2<-NULL # warning! the data frame stores the results so the next outputs will be amended based on the last value stored
for(i in 1:nrow(grid)){
  parameters1["vacc_coverage"] <- grid$vacc_coverage[i];
  parameters1["rho_1"] <- grid$rho_1[i]; # set the vaccination efficacy to the current value in the loop)
  parameters1["beta_h"] <- grid$beta_h[i];
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
           transmission_coefficient =  grid$beta_h[i],
           vacc_cov_percent = paste0(grid$vacc_coverage[i]*100, "%"),
           vacc_efficacy_percent = paste0(grid$rho_1[i]*100, "%")) %>% 
    select(everything())
  df.scen2<-bind_rows(df.scen2, tmp) 
}
View(df.scen2)

#check and make sure every output has reached equilibrium before vaccination is implemented
options(scipen = 999)

check_equilibrium <- df.scen2 %>% filter(time>=365*80-2, time<=365*80-1, vacc_coverage == 0, vacc_efficacy == 0) %>% 
  mutate(across(where(is.numeric), ~round(.x, 8))) %>% select(c(time, transmission_coefficient, S_c1, C_c1, I_c1, S_h1, C_h1, I_h1))

View(check_equilibrium) # ok #last row == the previous row

# filter equilibrium values right at a time step before vaccination
check_equilibrium2 <- df.scen2 %>% 
  mutate(difference = S_c1 - lag(S_c1)) %>%  filter(time >= 365*80-1) %>% 
  mutate(difference = S_c1 - lag(S_c1, 1)) %>%  filter(time >= 365*80-1)

# view the check equilibrium output
View(check_equilibrium2 %>% filter(vacc_coverage == 0, difference > abs(1e6))) #should be none #ok

df.scen2_long <- df.scen2 %>% 
  pivot_longer(cols = -c(time, vacc_coverage, vacc_efficacy, vacc_cov_percent, vacc_efficacy_percent, transmission_coefficient),
               names_to = "variable", 
               values_to = "value")
# save as RDS
#saveRDS(df.scen2, file = "5_scenario_analysis_v4/scenario_outputs/df.scen2.rds")

#saveRDS(df.scen2_long, file = "5_scenario_analysis_v4/scenario_outputs/df.scen2_long.rds")

#---------------------------------------------------------#
# total number of cases and deaths for the whole England
#---------------------------------------------------------#
# Total NHS Trust
nhs_trust <- 134

# define time bound as we only investigate the impact within 5 years and the vaccination impact will happen at 365*80
start_time <- 365*80
end_time <- 365*85


#-------------------------#
# total number of cases
#-------------------------#

# Calculate the total number of bsi cases in the baseline 
total_bsi_cases_with_scen2_hosp <- df.scen2_long %>% 
  filter(variable == "incidence_bsi_hosp") %>% #Filter solutions 
  filter(time>start_time, time<=end_time) %>% 
  mutate(year = floor((as.integer(time) - 1) / 365) + 1)

View(total_bsi_cases_with_scen2_hosp %>% filter(transmission_coefficient == 0.194))

total_bsi_cases_with_scen2_grouped <- total_bsi_cases_with_scen2_hosp %>% 
  filter(year > 0) %>% # filter out the first year
  group_by(year, vacc_coverage, vacc_efficacy, transmission_coefficient) %>% 
  summarise(total_cases = sum(value)*nhs_trust) %>% 
  mutate(total_cases = round(total_cases)) #perform rounding to get integer

total_bsi_cases_with_scen2_grouped

View(total_bsi_cases_with_scen2_grouped)

# total number of cases over 5 years 
total_number_of_cases_5_years <- total_bsi_cases_with_scen2_grouped %>% 
  group_by(vacc_coverage, vacc_efficacy, transmission_coefficient) %>% 
  summarise(total_cases = sum(total_cases))

View(total_number_of_cases_5_years)

# save to csv for use later
write.csv(total_bsi_cases_with_scen2_grouped, "5_scenario_analysis_v4/scenario_outputs/csv/cases_per_year_with_scen2.csv")

write.csv(total_number_of_cases_5_years, "5_scenario_analysis_v4/scenario_outputs/csv/cases_5_years_with_scen2.csv")


#-------------------------#
# cumulative outcomes
#-------------------------#

# all years
cum_cases_5_years <- read.csv("5_scenario_analysis_v4/scenario_outputs/csv/cases_5_years_with_scen2.csv")

View(cum_cases_5_years)


# remove redundant combinations
cum_cases_5_years_no_int_clean <- cum_cases_5_years  %>% filter(vacc_coverage == 0, vacc_efficacy == 0)



# remove redundant combinations
cum_cases_5_years_with_int_clean <- cum_cases_5_years  %>% filter(vacc_coverage == 0.7, vacc_efficacy == 0.7)


# calculate the outcome
outcome_scen2_cases <- data.frame(transmission_coefficient = cum_cases_5_years_no_int_clean$transmission_coefficient,
                            no_int = cum_cases_5_years_no_int_clean$total_cases,
                            with_int = cum_cases_5_years_with_int_clean$total_cases,
                            vacc_coverage = cum_cases_5_years_with_int_clean$vacc_coverage,
                            vacc_efficacy = cum_cases_5_years_with_int_clean$vacc_efficacy)

outcome_scen2_cases <- outcome_scen2_cases %>% mutate(cases_averted = abs(with_int - no_int),
                                                      percent_reduction = cases_averted/no_int)


# save to csv for use later
write.csv(outcome_scen2_cases, "5_scenario_analysis_v4/scenario_outputs/csv/cases_5_years_with_scen2.csv")


# plot the oucome
outcome_scen2_cases <- read.csv("5_scenario_analysis_v4/scenario_outputs/csv/cases_5_years_with_scen2.csv")

outcome_scen2_cases[is.na(outcome_scen2_cases)] <- 0

outcome_scen2_cases$transmission_coefficient2 <-  outcome_scen2_cases$transmission_coefficient %>% round(3)

outcome_scen2_cases$percent_reduction2 <-  round(outcome_scen2_cases$percent_reduction,3)*100 

View(outcome_scen2_cases)

# calculate the number required to vaccinate 
num_vaccines <- df.scen2 %>% filter(time == 365*80)

num_vaccines$total_num_vaccination <- round((num_vaccines$VS_c1 + num_vaccines$VC_c1 + num_vaccines$VS_h1 + num_vaccines$VC_h1)*nhs_trust,0)

View(num_vaccines)

# create a lookup table to assign the number of vaccination per coverage and efficacy
lookup_num_vaccines <- num_vaccines %>% select(c(vacc_coverage, vacc_efficacy, transmission_coefficient, total_num_vaccination)) 

# filter out zero coverage
lookup_num_vaccines  <- lookup_num_vaccines %>% filter(vacc_coverage != 0, vacc_efficacy != 0) %>% distinct() 

lookup_num_vaccines$vacc_cov_percent<- round(lookup_num_vaccines$vacc_coverage*100,0)

lookup_num_vaccines$transmission_coefficient2 <-  lookup_num_vaccines$transmission_coefficient %>% round(3)

View(lookup_num_vaccines)

# Add total number of vaccination column # the total population in the community will be slightly different across prevalence as the proportion of population in the hospital is different
outcome_scen2_cases <- left_join(outcome_scen2_cases, lookup_num_vaccines, by = c("vacc_coverage", "vacc_efficacy", "transmission_coefficient2"))

View(outcome_scen2_cases)

# outcome per vaccination
outcome_scen2_cases$case_averted_per_1000_vacc <- round((outcome_scen2_cases$cases_averted/outcome_scen2_cases$total_num_vaccination)*1000, 1)

outcome_scen2_cases$case_averted_per_1000_vacc <- ifelse(is.infinite(outcome_scen2_cases$case_averted_per_1000_vacc), 0, outcome_scen2_cases$case_averted_per_1000_vacc)

# effort to avert one case
outcome_scen2_cases$num_vacc_per_one_case_averted <- round(outcome_scen2_cases$total_num_vaccination/outcome_scen2_cases$cases_averted, 1)

outcome_scen2_cases$num_vacc_per_one_case_averted <- ifelse(is.infinite(outcome_scen2_cases$num_vacc_per_one_case_averted), 0, outcome_scen2_cases$num_vacc_per_one_case_averted)

# filter the evenly spaced out items so the plot is not too cluttered
outcome_scen2_cases_clean <- outcome_scen2_cases %>% filter(X %in% c(2, 4, 6, 8, 10))

View(outcome_scen2_cases_clean)


# plot absolute number of cases in stacked bar chart
outcome_scen2_cases_clean$cases_averted2 <- ifelse(outcome_scen2_cases_clean$cases_averted == 0, NA, outcome_scen2_cases_clean$cases_averted)

df.stacked <- pivot_longer(outcome_scen2_cases_clean, 
                           cols=c("with_int", "cases_averted"),
                           names_to = "quantity",
                           values_to = "value")
View(df.stacked)

total_cases_stacked_bars <- df.stacked %>%
  group_by(transmission_coefficient2) %>%
  summarise(total_cases = sum(value), .groups = "drop")

View(total_cases_stacked_bars)

ggplot(df.stacked, aes(x = factor(transmission_coefficient2), y = value, fill = quantity)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("with_int" = "#F9D9D3", "cases_averted" = "lightblue"),
                    labels = c("with_int" = "Number of cases with vaccination",   # rename legend label
                               "cases_averted" = "Number of cases averted")) +
  geom_text(aes(label = ifelse(value == 0, "R0<1", value)),
            position = position_stack(vjust = 0.5), # places text in middle of each bar segment
            color = "black", size = 6) +
  
  ylim(0, 165000) +
  labs(
    title = "Scenario 2\n(varying transmissibility in the hospital)",
    x = "Transmission coefficient (βh)",
    y = "Total number of cases without vaccination \n(over 5 years)",
    fill = NULL
  ) +
  theme_minimal() +
  theme (panel.grid = element_blank(),
         title = element_text(size = 18),
         axis.title = element_text(size = 18, hjust = 0.5),
         axis.title.x = element_text(size = 18),
         axis.text.y = element_text(size = 18),
         axis.text.x = element_text(size = 18),
         legend.position = "inside",
         legend.position.inside = c(0.35, 0.75),
         legend.text = element_text(size = 16)) 

ggsave("5_scenario_analysis_v4/scenario_outputs/plots/scen2/stacked_chart.png", width = 8, height = 6, units = "in", dpi = 300)


# plot absolute and percentage reduction in the number of cases with bubble plot
ggplot(outcome_scen2_cases_clean, aes(x = factor(transmission_coefficient2), y = percent_reduction2)) +
  geom_point(aes(size = cases_averted), 
             shape = 21, fill = "lightblue", color = "black") +
  geom_text(aes(label = ifelse(cases_averted ==0, "R0<1", cases_averted)), size = 6, color = "black", vjust = 0.5) +
  geom_text(data = subset(outcome_scen2_cases_clean, cases_averted == 0),
            aes(x = factor(transmission_coefficient2), label = "R0<1"),
            y = 40 + 0.2, vjust = 0, size = 6) +
  coord_cartesian(ylim = c(40, 50)) +
  scale_size_area(
    name = "Number of cases averted",
    max_size = 36,                # biggest bubble in mm
    breaks = c(7000, 45000, 65000),
    labels = scales::comma
  ) +
  annotate("text", x = 0.9, y = 41.5, label = "Bubble size: Number of cases averted (over 5 years)",
           hjust = 0, vjust = 0, size = 6, fontface = "italic") + 
  labs(
    x = "Transmission coefficient (βh)",
    y = "Percentage of cases averted\n(over 5 years)",
    title = "Scenario 2\n(varying transmissibility in the hospital)"
  ) +
  theme_minimal() +
  theme (panel.grid = element_blank(),
         title = element_text(size = 18),
         axis.title = element_text(size = 18, hjust = 0.5),
         axis.title.x = element_text(size = 18),
         axis.text.y = element_text(size = 18),
         axis.text.x = element_text(size = 18),
         legend.position = "none") 

ggsave("5_scenario_analysis_v4/scenario_outputs/plots/scen2/bubble_plot.png", width = 8, height = 6, units = "in", dpi = 300)


# plot number of cases averted per 1000 vaccinated persons
ggplot(outcome_scen2_cases_clean, aes(x = factor(transmission_coefficient2), y = case_averted_per_1000_vacc)) +
  geom_col(fill = "#FCECC2") +
  geom_text(aes(label = ifelse(case_averted_per_1000_vacc == 0, "R0<1", case_averted_per_1000_vacc)), vjust = -0.5, size = 6) +
  ylim(0, 10) +
  labs(
    x = "Transmission coefficient (βh)",
    y = "Number of cases averted\nper 1000 vaccinated persons \n(over 5 years)",
    title = "Scenario 2\n(varying transmissibility in the hospital)"
  ) +
  theme_minimal() +
  theme (panel.grid = element_blank(),
         title = element_text(size = 18),
         axis.title = element_text(size = 18, hjust = 0.5),
         axis.title.x = element_text(size = 18),
         axis.text.y = element_blank(),
         axis.text.x = element_text(size = 18),
         axis.title.y = element_text(size =18))  

ggsave("5_scenario_analysis_v4/scenario_outputs/plots/scen2/case_per_1000_vacc_varying_transmissibility.png", width = 7.5, height = 6, units = "in", dpi = 300)


# plot number of vaccinated persons per one case averted
ggplot(outcome_scen2_cases_clean, aes(x = factor(transmission_coefficient2), y = num_vacc_per_one_case_averted)) +
  geom_col(fill = "#F5E4EC") +
  geom_text(aes(label = ifelse(num_vacc_per_one_case_averted == 0, "R0<1", num_vacc_per_one_case_averted)), vjust = -0.5, size = 6) +
  ylim(0, 1200) +
  labs(
    x = "Transmission coefficient (βh)",
    y = "Number of vaccinated persons\nper one case averted \n(over 5 years)",
    title = "Scenario 2\n(varying transmissibility in the hospital)"
  ) +
  theme_minimal() +
  theme (panel.grid = element_blank(),
         title = element_text(size = 18),
         axis.title = element_text(size = 18, hjust = 0.5),
         axis.title.x = element_text(size = 18),
         axis.text.y = element_blank(),
         axis.text.x = element_text(size = 18),
         axis.title.y = element_text(size =18)) 

ggsave("5_scenario_analysis_v4/scenario_outputs/plots/scen2/num_vacc_per_case_varying_transmissibility.png", width = 7.5, height = 6, units = "in", dpi = 300)


# for checking

# plot absolute reduction in the number of cases
ggplot(outcome_scen2_cases_clean, aes(x = factor(transmission_coefficient2), y = cases_averted)) +
  geom_col(fill = "darkmagenta") +
  geom_text(aes(label = ifelse(cases_averted == 0, "no \ncases", cases_averted)), vjust = -0.5, size = 6) +
  ylim(0, 165000) +
  labs(
    x = "Transmission coefficient (βh)",
    y = "",
    title = "Scenario 2 (varying transmissibility in the hospital)\n\nAbsolute reduction in number of cases \n(over 5 years)"
  ) +
  theme_minimal() +
  theme (panel.grid = element_blank(),
         title = element_text(size = 18),
         axis.title = element_text(size = 18, hjust = 0.5),
         axis.title.x = element_text(size = 18),
         axis.text.y = element_blank(),
         axis.text.x = element_text(size = 18)) 

ggsave("5_scenario_analysis_v4/scenario_outputs/plots/scen2/absolute_reduction_varying_transmissibility.png", width = 6, height = 6, units = "in", dpi = 300)

# plot percentage reduction in the number of cases
ggplot(outcome_scen2_cases_clean, aes(x = factor(transmission_coefficient2), y = percent_reduction2)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = ifelse(percent_reduction2 == 0, "R0<1", percent_reduction2)), vjust = -0.5, size = 6) +
  ylim(0, 100) +
  labs(
    x = "Transmission coefficient (βh)",
    y = "",
    title = "Scenario 2 (varying transmissibility in the hospital)\n\nPercentage reduction in number of cases \n(over 5 years)"
  ) +
  theme_minimal() +
  theme (panel.grid = element_blank(),
         title = element_text(size = 18),
         axis.title = element_text(size = 18, hjust = 0.5),
         axis.title.x = element_text(size = 18),
         axis.text.y = element_blank(),
         axis.text.x = element_text(size = 18)) 

ggsave("5_scenario_analysis_v4/scenario_outputs/plots/scen2/percentage_reduction_varying_transmissibility.png", width = 6, height = 6, units = "in", dpi = 300)

