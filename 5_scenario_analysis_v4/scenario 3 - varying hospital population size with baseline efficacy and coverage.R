# varying community prevalence
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
# use the same model code as the one we used for sensitivity analysis 
source("5_scenario_analysis_v4/scenario_inputs/model1_v4_for_scenario analysis.R")

#==================================================================================================
#### 2. Setup initial values #### #ensure the order is the same as the order of the equations above
#==================================================================================================
istate1 <- c(S_c1 = S_c1,
             C_c1 = C_c1, # if the scenario analysis is turned off, use baseline value, otherwise follow the desired prevalence of colonised in the community
             I_c1 = I_c1,
             VS_c1 = VS_c1, 
             VC_c1 = VC_c1, 
             S_h1 = S_h1, #use baseline or modify proportion of patients in the hospital for scenario analysis 
             C_h1 = C_h1, #use baseline or modify proportion of patients in the hospital for scenario analysis 
             I_h1 = I_h1,
             VS_h1 = VS_h1, 
             VC_h1 = VC_h1,
             D_h1 = D_h1, CInc_1a = CInc_1a, CInc_1b = CInc_1b, CInc_1c = CInc_1c, CInc_1d = CInc_1d,
             col_adm = 0, col_discharges = 0,
             admissions = 0, discharges = 0) 

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
out1 <- ode(y = istate1, times = tps1, func = klebsiella_1, parms = parameters1, method = "lsoda", events = event_list)

#-------------------------------------------------------------------------------------------------#
# Generate outcomes for vaccination scenario #### 
#-------------------------------------------------------------------------------------------------#
# load the model code 
turn_on_scenario <- 1 #0 for turn off, 1 for turn on

sim1 <- seq(0, 0.7, length.out = 2) # define the vaccination coverage 
sim2 <- seq(0, 0.7, length.out = 2) # define the vaccination efficacy
sim3 <- seq(0.004, 0.02, length.out = 5) # define the proportion of population in the hospital

# create grid to check different combinations
grid <- expand.grid(vacc_coverage = sim1, rho_1 = sim2, ph = sim3)

# Load the model output
# Add columns containing the desired values
df.scen3<-NULL # warning! the data frame stores the results so the next outputs will be amended based on the last value stored
for(i in 1:nrow(grid)){
  parameters1["vacc_coverage"] <- grid$vacc_coverage[i];
  parameters1["rho_1"] <- grid$rho_1[i]; # set the vaccination efficacy to the current value in the loop)
  istate1["S_c1"] <-  80470 - S_c1 * grid$ph[i] - C_c1 * grid$ph[i]; #keep the total population the same while modifying proportion of patients in the hospital for scenario analysis 
  istate1["C_c1"] <-  80470 - C_c1 * grid$ph[i]; #keep the total population the same while modifying modify proportion of patients in the hospital for scenario analysis 
  istate1["C_c1"] <-  0; #keep the total population the same while modifying modify proportion of patients in the hospital for scenario analysis 
  istate1["S_h1"] <-  S_c1 * grid$ph[i]; #modify proportion of patients in the hospital for scenario analysis 
  istate1["C_h1"] <-  C_c1 * grid$ph[i]; #modify proportion of patients in the hospital for scenario analysis 
  istate1["I_h1"] <-  0; #modify proportion of patients in the hospital for scenario analysis 
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
           percent_hosp =  grid$ph[i],
           vacc_cov_percent = paste0(grid$vacc_coverage[i]*100, "%"),
           vacc_efficacy_percent = paste0(grid$rho_1[i]*100, "%")) %>% 
    select(everything())
  df.scen3<-bind_rows(df.scen3, tmp) 
}
View(df.scen3)

#check and make sure every output has reached equilibrium before vaccination is implemented
options(scipen = 999)

check_equilibrium <- df.scen3 %>% filter(time>=365*80-2, time<=365*80-1, vacc_coverage == 0, vacc_efficacy == 0) %>% 
  mutate(across(where(is.numeric), ~round(.x, 8))) %>% select(c(time, percent_hosp, S_c1, C_c1, I_c1, S_h1, C_h1, I_h1))

View(check_equilibrium) # ok #last row == the previous row

# filter equilibrium values right at a time step before vaccination
check_equilibrium3 <- df.scen3 %>% 
  mutate(difference = S_c1 - lag(S_c1)) %>%  filter(time >= 365*80-1) %>% 
  mutate(difference = S_c1 - lag(S_c1, 1)) %>%  filter(time >= 365*80-1)

# view the check equilibrium output
View(check_equilibrium3 %>% filter(vacc_coverage == 0, difference > abs(1e6))) #should be none #ok

# As we want the data to be grouped according to the prevalence before intervention, we need to assign label for each group
# create a lookup table based on the unique combination of auxiliary parameter and the prevalence before intervention
lookup_table <- df.scen3 %>% filter (time == 29200) %>%  # pick a random time point where the system has reached steady state
  distinct(percent_hosp)

# assign original prevalence before intervention to intervention column
# check if a row has intervention applied 
df.scen3$intervention <- ifelse(df.scen3$VS_c1 > 0, TRUE, FALSE)

# check if the label works
View(head(df.scen3))
View(tail(df.scen3))

# when a row is an intervention row, assign original prevalence of colonisation in the community
# join the dataframe 
df.scen3_join <- df.scen3 %>% 
  left_join(lookup_table, by = "percent_hosp")

# check if the join table works
View(df.scen3_join %>% filter(time >29200))
View(df.scen3_join %>% filter(time >29200))

df.scen3_long <- df.scen3_join %>% 
  pivot_longer(cols = -c(time, vacc_coverage, vacc_efficacy, vacc_cov_percent, vacc_efficacy_percent, percent_hosp),
               names_to = "variable", 
               values_to = "value")
# save as RDS
#saveRDS(df.scen3_join, file = "5_scenario_analysis_v4/scenario_outputs/df.scen3.rds")

#saveRDS(df.scen3_long, file = "5_scenario_analysis_v4/scenario_outputs/df.scen3_long.rds")

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
total_bsi_cases_with_scen3_hosp <- df.scen3_long %>% 
  filter(variable == "incidence_bsi_hosp") %>% #Filter solutions 
  filter(time>start_time, time<=end_time) %>% 
  mutate(year = floor((as.integer(time) - 1) / 365) + 1)

total_bsi_cases_with_scen3_grouped <- total_bsi_cases_with_scen3_hosp %>% 
  filter(year > 0) %>% # filter out the first year
  group_by(year, vacc_coverage, vacc_efficacy, percent_hosp) %>% 
  summarise(total_cases = sum(value)*nhs_trust) %>% 
  mutate(total_cases = round(total_cases)) #perform rounding to get integer

total_bsi_cases_with_scen3_grouped

View(total_bsi_cases_with_scen3_grouped)

# total number of cases over 5 years 
total_number_of_cases_5_years <- total_bsi_cases_with_scen3_grouped %>% 
  group_by(vacc_coverage, vacc_efficacy, percent_hosp) %>% 
  summarise(total_cases = sum(total_cases))

View(total_number_of_cases_5_years)

# save to csv for use later
write.csv(total_bsi_cases_with_scen3_grouped, "5_scenario_analysis_v4/scenario_outputs/csv/cases_per_year_with_scen3.csv")

write.csv(total_number_of_cases_5_years, "5_scenario_analysis_v4/scenario_outputs/csv/cases_5_years_with_scen3.csv")

#-------------------------#
# cumulative outcomes
#-------------------------#

# all years
cum_cases_5_years <- read.csv("5_scenario_analysis_v4/scenario_outputs/csv/cases_5_years_with_scen3.csv")

View(cum_cases_5_years)



# remove redundant combinations
cum_cases_5_years_no_int_clean <- cum_cases_5_years  %>% filter(vacc_coverage == 0, vacc_efficacy == 0)



# remove redundant combinations
cum_cases_5_years_with_int_clean <- cum_cases_5_years  %>% filter(vacc_coverage == 0.7, vacc_efficacy == 0.7)



# calculate the outcome
outcome_scen3_cases <- data.frame(no_int = cum_cases_5_years_no_int_clean$total_cases,
                                  with_int = cum_cases_5_years_with_int_clean$total_cases,
                                  percent_hosp = round(cum_cases_5_years_with_int_clean$percent_hosp*100,3),
                                  vacc_coverage = cum_cases_5_years_with_int_clean$vacc_coverage,
                                  vacc_efficacy = cum_cases_5_years_with_int_clean$vacc_efficacy)

View(outcome_scen3_cases)

outcome_scen3_cases <- outcome_scen3_cases %>% mutate(cases_averted = abs(with_int - no_int),
                                                      percent_reduction = cases_averted/no_int)


# save to csv for use later
write.csv(outcome_scen3_cases, "5_scenario_analysis_v4/scenario_outputs/csv/cases_5_years_with_scen3.csv")

View(outcome_scen3_cases)


# plot the outcome
outcome_scen3_cases <- read.csv("5_scenario_analysis_v4/scenario_outputs/csv/cases_5_years_with_scen3.csv")

outcome_scen3_cases[is.na(outcome_scen3_cases)] <-  0

outcome_scen3_cases$percent_reduction2 <-  round(outcome_scen3_cases$percent_reduction,3)*100 

View(outcome_scen3_cases)

# calculate the number required to vaccinate 
num_vaccines <- df.scen3 %>% filter(time == 365*80)

num_vaccines$total_num_vaccination <- round((num_vaccines$VS_c1 + num_vaccines$VC_c1 + num_vaccines$VS_h1 + num_vaccines$VC_h1)*nhs_trust,0)

View(num_vaccines)

# create a lookup table to assign the number of vaccination per coverage and efficacy
lookup_num_vaccines <- num_vaccines %>% select(c(vacc_coverage, vacc_efficacy, percent_hosp, total_num_vaccination)) 

# filter out zero coverage
lookup_num_vaccines  <- lookup_num_vaccines %>% filter(vacc_coverage == 0.7, vacc_efficacy ==0.7) %>% distinct() 

lookup_num_vaccines$vacc_cov_percent<- round(lookup_num_vaccines$vacc_coverage*100,0)

lookup_num_vaccines$percent_hosp<- round(lookup_num_vaccines$percent_hosp*100,3)

View(lookup_num_vaccines)

# Add total number of vaccination column # the total population in the community will be slightly different across prevalence as the proportion of population in the hospital is different
outcome_scen3_cases <- left_join(outcome_scen3_cases, lookup_num_vaccines, by = c("percent_hosp","vacc_coverage","vacc_efficacy"))

View(outcome_scen3_cases)

# outcome per vaccination
outcome_scen3_cases$case_averted_per_1000_vacc <- round((outcome_scen3_cases$cases_averted/outcome_scen3_cases$total_num_vaccination)*1000, 1)

outcome_scen3_cases$case_averted_per_1000_vacc <- ifelse(is.infinite(outcome_scen3_cases$case_averted_per_1000_vacc), 0, outcome_scen3_cases$case_averted_per_1000_vacc)


# effort to avert one case
outcome_scen3_cases$num_vacc_per_one_case_averted <- round(outcome_scen3_cases$total_num_vaccination/outcome_scen3_cases$cases_averted, 1)

outcome_scen3_cases$num_vacc_per_one_case_averted <- ifelse(is.infinite(outcome_scen3_cases$num_vacc_per_one_case_averted), 0, outcome_scen3_cases$num_vacc_per_one_case_averted)

# filter the evenly spaced out items so the plot is not too cluttered
outcome_scen3_cases_clean <- outcome_scen3_cases #%>% filter(X %in% c(4, 8, 12, 16, 20))

View(outcome_scen3_cases_clean)

# plot absolute number of cases in stacked bar chart
outcome_scen3_cases_clean$cases_averted2 <- ifelse(outcome_scen3_cases_clean$cases_averted == 0, NA, outcome_scen3_cases_clean$cases_averted)

df.stacked <- pivot_longer(outcome_scen3_cases_clean, 
                           cols=c("with_int", "cases_averted"),
                           names_to = "quantity",
                           values_to = "value")
View(df.stacked)

total_cases_stacked_bars <- df.stacked %>%
  group_by(percent_hosp) %>%
  summarise(total_cases = sum(value), .groups = "drop")

View(total_cases_stacked_bars)

ggplot(df.stacked, aes(x = factor(percent_hosp), y = value, fill = quantity)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("with_int" = "#F9D9D3", "cases_averted" = "lightblue"),
                    labels = c("with_int" = "Number of cases with vaccination",   # rename legend label
                               "cases_averted" = "Number of cases averted")) +
  geom_text(aes(label = ifelse(value == 0, "R0<1", value)),
            position = position_stack(vjust = 0.5), # places text in middle of each bar segment
            color = "black", size = 8) +
  
  ylim(0, 400000) +
  labs(
    title = "Scenario 3\n(varying hospital population size)",
    x = "Percentage of hospital population out of total population",
    y = "Total number of cases without vaccination \n(over 5 years)",
    fill = NULL
  ) +
  theme_minimal() +
  theme (panel.grid = element_blank(),
         title = element_text(size = 24),
         axis.title = element_text(size = 24, hjust = 0.5),
         axis.title.x = element_text(size = 24),
         axis.text.y = element_text(size = 24),
         axis.text.x = element_text(size = 24),
         legend.position = "inside",
         legend.position.inside = c(0.35, 0.75),
         legend.text = element_text(size = 20)) 

ggsave("5_scenario_analysis_v4/scenario_outputs/plots/scen3/stacked_chart_withlabel.png", width = 10.5, height = 8.5, units = "in", dpi = 300)


ggplot(df.stacked, aes(x = factor(percent_hosp), y = value, fill = quantity)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("with_int" = "#F9D9D3", "cases_averted" = "lightblue"),
                    labels = c("with_int" = "Number of cases with vaccination",   # rename legend label
                               "cases_averted" = "Number of cases averted")) +
  geom_text(data = subset(df.stacked, factor(percent_hosp, levels = unique(percent_hosp)) != unique(percent_hosp)[1]),
            aes(label = value),
            position = position_stack(vjust = 0.5), # places text in middle of each bar segment
            color = "black", size = 8) +
  
  ylim(0, 400000) +
  labs(
    title = "Scenario 3\n(varying hospital population size)",
    x = "Percentage of hospital population out of total population",
    y = "Total number of cases without vaccination \n(over 5 years)",
    fill = NULL
  ) +
  theme_minimal() +
  theme (panel.grid = element_blank(),
         title = element_text(size = 24),
         axis.title = element_text(size = 24, hjust = 0.5),
         axis.title.x = element_text(size = 24),
         axis.text.y = element_text(size = 24),
         axis.text.x = element_text(size = 24),
         legend.position = "inside",
         legend.position.inside = c(0.35, 0.75),
         legend.text = element_text(size = 20)) 

ggsave("5_scenario_analysis_v4/scenario_outputs/plots/scen3/stacked_chart_nofirstbarlabel.png", width = 10.5, height = 8.5, units = "in", dpi = 300)


# plot absolute and percentage reduction in the number of cases with bubble plot
ggplot(outcome_scen3_cases_clean, aes(x = factor(percent_hosp), y = percent_reduction2), size = cases_averted2) +
  geom_point(aes(size = cases_averted), 
             shape = 21, fill = "lightblue", color = "black") +
  geom_text(aes(label = ifelse(cases_averted ==0, "R0<1", cases_averted)), size = 8, color = "black", vjust = 0.5) +
  geom_text(data = subset(outcome_scen3_cases_clean, cases_averted == 0),
            aes(x = factor(percent_hosp), label = "R0<1"),
            y = 40 + 0.2, vjust = 0, size = 8) +
  coord_cartesian(ylim = c(40, 52.5)) +
  scale_size_area(
    name = "Number of cases averted",
    max_size = 36,                # biggest bubble in mm
    breaks = c(7000, 45000, 65000),
    labels = scales::comma
  ) +
  annotate("text", x = 0.9, y = 41, label = "Bubble size: Number of cases averted (over 5 years)",
           hjust = 0, vjust = 0, size = 8, fontface = "italic") + 
  labs(
    x = "Percentage of hospital population out of total population",
    y = "Percentage of cases averted\n(over 5 years)",
    title = "Scenario 3\n(varying hospital population size)"
  ) +
  theme_minimal() +
  theme (panel.grid = element_blank(),
         title = element_text(size = 24),
         axis.title = element_text(size = 24, hjust = 0.5),
         axis.title.x = element_text(size = 24),
         axis.text.y = element_text(size = 24),
         axis.text.x = element_text(size = 24),
         legend.position = "none") 

ggsave("5_scenario_analysis_v4/scenario_outputs/plots/scen3/bubble_plot.png", width = 10.5, height = 8.5, units = "in", dpi = 300)


# plot cases averted per 1000 vaccinated persons
ggplot(outcome_scen3_cases_clean, aes(x = factor(percent_hosp), y = case_averted_per_1000_vacc)) +
  geom_col(fill = "#FCECC2") +
  geom_text(aes(label = ifelse(case_averted_per_1000_vacc == 0, "no \ncases", case_averted_per_1000_vacc)), vjust = -0.5, size = 8) +
  ylim(0, 25) +
  labs(
    x = "Percentage of hospital population out of total population",
    y = "Number of cases averted\nper 1000 vaccinated persons \n(over 5 years)",
    title = "Scenario 3\n(varying hospital population size)"
  ) +
  theme_minimal() +
  theme (panel.grid = element_blank(),
         title = element_text(size = 20),
         axis.title = element_text(size = 18, hjust = 0.5),
         axis.title.x = element_text(size = 20),
         axis.text.x = element_text(size = 20),
         axis.title.y = element_text(size = 20)) 

ggsave("5_scenario_analysis_v4/scenario_outputs/plots/scen3/case_per_1000_vacc_varying_percentage_in_hospital.png", width = 8, height = 7, units = "in", dpi = 300)


# plot number of vaccinated persons per one case averted
ggplot(outcome_scen3_cases_clean, aes(x = factor(percent_hosp), y = num_vacc_per_one_case_averted)) +
  geom_col(fill = "#F5E4EC") +
  geom_text(aes(label = ifelse(num_vacc_per_one_case_averted == 0, "no \ncases", num_vacc_per_one_case_averted)), vjust = -0.5, size = 8) +
  ylim(0, 1200) +
  labs(
    x = "Percentage of hospital population out of total population",
    y = "Number of vaccinated persons\nper one case averted \n(over 5 years)",
    title = "Scenario 3\n(varying hospital population size)"
  ) +
  theme_minimal() +
  theme (panel.grid = element_blank(),
         title = element_text(size = 20),
         axis.title = element_text(size = 18, hjust = 0.5),
         axis.title.x = element_text(size = 20),
         axis.text.y = element_blank(),
         axis.text.x = element_text(size = 20),
         axis.title.y = element_text(size = 20)) 

ggsave("5_scenario_analysis_v4/scenario_outputs/plots/scen3/num_vacc_per_case_varying_percentage_in_hospital.png", width = 8, height = 7, units = "in", dpi = 300)



# for checking
# plot absolute reduction in the number of cases
ggplot(outcome_scen3_cases_clean, aes(x = factor(percent_hosp), y = cases_averted)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = ifelse(cases_averted == 0, "no \ncases", cases_averted)), vjust = -0.5, size = 6) +
  ylim(0, 165000) +
  labs(
    x = "Percentage of hospital population out of total population",
    y = "",
    title = "Scenario 3\n(varying hospital population size)\n\nAbsolute reduction in number of cases \n(over 5 years)"
  ) +
  theme_minimal() +
  theme (panel.grid = element_blank(),
         title = element_text(size = 20),
         axis.title = element_text(size = 18, hjust = 0.5),
         axis.title.x = element_text(size = 20),
         axis.text.y = element_blank(),
         axis.text.x = element_text(size = 20)) 

ggsave("5_scenario_analysis_v4/scenario_outputs/plots/scen3/absolute_reduction_varying_percentage_in_hospital.png", width = 8, height = 7, units = "in", dpi = 300)

# plot percentage reduction in the number of cases
ggplot(outcome_scen3_cases_clean, aes(x = factor(percent_hosp), y = percent_reduction2)) +
  geom_col(fill = "darkmagenta") +
  geom_text(aes(label = ifelse(percent_reduction2 == 0, "no \ncases", percent_reduction2)), vjust = -0.5, size = 6) +
  ylim(0, 100) +
  labs(
    x = "Percentage of hospital population out of total population",
    y = "",
    title = "Scenario 3\n(varying hospital population size)\n\nPercentage reduction in number of cases \n(over 5 years)"
  ) +
  theme_minimal() +
  theme (panel.grid = element_blank(),
         title = element_text(size = 20),
         axis.title = element_text(size = 18, hjust = 0.5),
         axis.title.x = element_text(size = 20),
         axis.text.y = element_blank(),
         axis.text.x = element_text(size = 20)) 

ggsave("5_scenario_analysis_v4/scenario_outputs/plots/scen3/percentage_reduction_varying_percentage_in_hospital.png", width = 8, height = 7, units = "in", dpi = 300)
