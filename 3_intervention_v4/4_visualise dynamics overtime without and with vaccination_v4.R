# Visualise dynamics overtime without and with vaccination

# clear the environment
rm(list = ls())


# community colonisation, hospital colonisation, and hospital BSI without vaccination
# multiply up by 134 NHS trust to get the total number for the whole England

#-------------------------------------------------------------------------------------------------#
# Generate outcomes for no vaccination scenario #### 
#-------------------------------------------------------------------------------------------------#
# load the model code 
source("3_intervention_v4/intervention_inputs/model1_v4_for_intervention.R")

# define time bound as we only investigate the impact within 5 years and the vaccination impact will happen at 365*80
#===============================================================
#### 4. Setup timesteps ####
#===============================================================
start1 <- 0 #start time for the simulation
end1 <- 365*5 #end point of the timestep #we used equilibrium values after calibration 
tps1 <- seq(start1, end1, by = 1) #time steps for the simulation

#===============================================================
#### 5. Run the model and check the outputs ####
#===============================================================
#run the model 
#we are using the euler method such that the timesteps are read as integers
out1 <- ode(y = istate1, times = tps1, func = klebsiella_1, parms = parameters1, method = "lsoda", events = event_list)


#-------------------------------------------------------------------------------------------------#
# Generate outcomes for with vaccination scenario #### 
#-------------------------------------------------------------------------------------------------#
# load the model code 
source("3_intervention_v4/intervention_inputs/model1_v4_for_intervention.R")

# set parameters to turn on vaccination
# produce output with vaccination
parameters1["vacc_coverage"] <- 0.7 # set vaccination rate to 0.7 (baseline coverage)
parameters1["rho_1"] <- 0.7 # set vaccination efficacy to 100% (full efficacy)

# define time bound as we only investigate the impact within 5 years and the vaccination impact will happen at 365*80
#===============================================================
#### 4. Setup timesteps ####
#===============================================================
start1 <- 0 #start time for the simulation
end1 <- 365*5 #end point of the timestep #ensure the system reach equilibrium before the intervention is run
tps1 <- seq(start1, end1, by = 1) #time steps for the simulation

#===============================================================
#### 4. Setup the event list ####
#===============================================================

event_list <- list(func = vaccination_event,
                   time = c(0)) #time at which the vaccination event occurs -> the solver will implement it at t+1 or time 29200 (365*80)

#===============================================================
#### 5. Run the model and check the outputs ####
#===============================================================
#run the model 
#we are using the euler method such that the timesteps are read as integers
out2 <- ode(y = istate1, times = tps1, func = klebsiella_1, parms = parameters1, method = "lsoda", events = event_list)

df.before.intervention <- as.data.frame(out1) %>% 
  mutate(num_col_comm = C_c1 + VC_c1,
         num_col_hosp = C_h1 + VC_h1,
         num_bsi_hosp = I_h1,
         inc_col_hosp = c(0, diff(CInc_1a)),
         bsi_hosp = c(0, diff(CInc_1c)),
         hosp_onset = c(0, diff(CInc_1d)),
         prev_bsi_hosp = num_bsi_hosp / (S_h1 + VS_h1 + C_h1 + VC_h1 + I_h1),
         prev_col_hosp = num_col_hosp / (S_h1 + VS_h1 + C_h1 + VC_h1 + I_h1),
         prev_col_comm = num_col_comm / (S_c1 + VS_c1 + C_c1 + VC_c1 + I_c1),
         col_admissions = c(0, diff(col_adm)),
         col_discharges = c(0, diff(col_discharges)),
         total_admissions = c(0, diff(admissions)),
         total_discharges = c(0, diff(discharges))) %>% 
  pivot_longer(cols = -time,
               names_to = "variable",
               values_to = "value") %>% mutate(scenario = "no vaccination")

df.after.intervention <-  as.data.frame(out2) %>% 
  mutate(num_col_comm = C_c1 + VC_c1,
         num_col_hosp = C_h1 + VC_h1,
         num_bsi_hosp = I_h1,
         inc_col_hosp = c(0, diff(CInc_1a)),
         bsi_hosp = c(0, diff(CInc_1c)),
         hosp_onset = c(0, diff(CInc_1d)),
         prev_bsi_hosp = num_bsi_hosp / (S_h1 + VS_h1 + C_h1 + VC_h1 + I_h1),
         prev_col_hosp = num_col_hosp / (S_h1 + VS_h1 + C_h1 + VC_h1 + I_h1),
         prev_col_comm = num_col_comm / (S_c1 + VS_c1 + C_c1 + VC_c1 + I_c1),
         col_admissions = c(0, diff(col_adm)),
         col_discharges = c(0, diff(col_discharges)),
         total_admissions = c(0, diff(admissions)),
         total_discharges = c(0, diff(discharges))) %>% 
  pivot_longer(cols = -time,
               names_to = "variable",
               values_to = "value") %>% mutate(scenario = "with vaccination")

df.combined <- bind_rows(df.before.intervention, df.after.intervention)

df.diff <- df.combined %>%
  pivot_wider(names_from = scenario, values_from = value) %>%
  rename(ymin = `with vaccination`, ymax = `no vaccination`)

df.diff$diff <- df.diff$ymax - df.diff$ymin

nhs_trust <- 134


# plot prevalence of bsi in the hospital - entire England
ggplot() +
  geom_line(data = df.combined %>% filter(variable == "prev_bsi_hosp"),
            aes(x = time, y = value*100, linetype =  scenario), color = "#CC79A7",
            linewidth = 1.2) +
  scale_linetype_manual(
    values = c("no vaccination" = "dotted", "with vaccination" = "solid")
  ) +
  #ylim(0, 100) +
  labs(title = NULL, 
       x = "Number of days", 
       y = "Hospital prevalence of BSI \n(% of hospital population)",
       fill = NULL,
       linetype = NULL) + 
  guides(fill = "none") +  #removes fill legend
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 16),
        legend.position = "top",
        legend.text  = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))

ggsave("3_intervention_v4/intervention_outputs/plots/with and without intervention plot overtime/prevalence of bsi patients in hospital.png",
       width = 5,
       height = 5,
       units = "in",
       dpi = 300)



# plot prevalence of colonisation in the hospital - entire England
ggplot() +
  geom_line(data = df.combined %>% filter(variable == "prev_col_hosp"),
            aes(x = time, y = value*100, linetype =  scenario), color = "#D55E00",
            linewidth = 1.2) +
  scale_linetype_manual(
    values = c("no vaccination" = "dotted", "with vaccination" = "solid")
  ) +
  #ylim(0, 100) +
  labs(title = NULL, 
       x = "Number of days", 
       y = "Hospital prevalence of colonisation \n(% of hospital population)",
       fill = NULL,
       linetype = NULL) + 
  guides(fill = "none") +  #removes fill legend
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 16),
        legend.position = "top",
        legend.text  = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))

ggsave("3_intervention_v4/intervention_outputs/plots/with and without intervention plot overtime/prevalence of colonised in hospital.png",
       width = 5,
       height = 5,
       units = "in",
       dpi = 300)



# plot prevalence of colonisation in the community - entire England
ggplot() +
  geom_line(data = df.combined %>% filter(variable == "prev_col_comm"),
            aes(x = time, y = value*100, linetype =  scenario), color = "#F0E442",
            linewidth = 1.2) +
  scale_linetype_manual(
    values = c("no vaccination" = "dotted", "with vaccination" = "solid")
  ) +
  #ylim(0, 100) +
  labs(title = NULL, 
       x = "Number of days", 
       y = "Community prevalence of colonisation \n(% of community population)",
       fill = NULL,
       linetype = NULL) + 
  guides(fill = "none") +  #removes fill legend
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 16),
        legend.position = "top",
        legend.text  = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))

ggsave("3_intervention_v4/intervention_outputs/plots/with and without intervention plot overtime/prevalence of colonised in community.png",
       width = 5,
       height = 5,
       units = "in",
       dpi = 300)


# plot number of of BSI patients in the hospital - entire England
ggplot() +
  geom_line(data = df.diff %>% filter(variable == "num_bsi_hosp"),
            aes(x = time, y = diff*nhs_trust), color = "#CC79A7",
            linewidth = 1.2) +
  labs(title = NULL, 
       x = "Number of days", 
       y = "Absolute reduction\nin the number of BSI patients\nin the hospital",
       fill = NULL,
       linetype = NULL) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 16),
        legend.position = "top",
        legend.text  = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))

ggsave("3_intervention_v4/intervention_outputs/plots/with and without intervention plot overtime/number of BSI in hospital.png",
       width = 5,
       height = 5,
       units = "in",
       dpi = 300)



# plot number of of colonisation in the hospital - entire England
ggplot() +
  geom_line(data = df.diff %>% filter(variable == "num_col_hosp"),
            aes(x = time, y = diff*nhs_trust), color = "#D55E00",
            linewidth = 1.2) +
  labs(title = NULL, 
       x = "Number of days", 
       y = "Absolute reduction\nin the number of colonised patients\nin the hospital",
       fill = NULL,
       linetype = NULL) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 16),
        legend.position = "top",
        legend.text  = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))

ggsave("3_intervention_v4/intervention_outputs/plots/with and without intervention plot overtime/number of colonised in hospital.png",
       width = 5,
       height = 5,
       units = "in",
       dpi = 300)


# plot number of of colonisation in the hospital - entire England
ggplot() +
  geom_line(data = df.diff %>% filter(variable == "num_col_comm"),
            aes(x = time, y = diff*nhs_trust), color = "#F0E442",
            linewidth = 1.2) +
  labs(title = NULL, 
       x = "Number of days", 
       y = "Absolute reduction\nin the number of colonised persons\nin the community",
       fill = NULL,
       linetype = NULL) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 16),
        legend.position = "top",
        legend.text  = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))

ggsave("3_intervention_v4/intervention_outputs/plots/with and without intervention plot overtime/number of colonised in community.png",
       width = 5,
       height = 5,
       units = "in",
       dpi = 300)


# additional table for proportion of BSI and colonised upon admissions and proportion of colonised upon discharges
bsi_adm <- df.combined %>% filter(variable == "bsi_hosp", time >0) %>% 
  group_by(scenario) %>%
  summarise(total_bsi_adm = sum(value)*nhs_trust) %>% 
  mutate(total_bsi_adm = round(total_bsi_adm))


col_adm <- df.combined %>% filter(variable == "col_admissions", time >0) %>% 
  group_by(scenario) %>%
  summarise(total_col_adm = sum(value)*nhs_trust) %>% 
  mutate(total_col_adm = round(total_col_adm))


col_discharges <- df.combined %>% filter(variable == "col_discharges", time >0) %>% 
  group_by(scenario) %>%
  summarise(total_col_discharges = sum(value)*nhs_trust) %>% 
  mutate(total_col_discharges = round(total_col_discharges))


adm <- df.combined %>% filter(variable == "total_admissions", time >0) %>% 
  group_by(scenario) %>%
  summarise(total_adm = sum(value)*nhs_trust) %>% 
  mutate(total_adm = round(total_adm))


disc <- df.combined %>% filter(variable == "total_discharges", time >0) %>% 
  group_by(scenario) %>%
  summarise(total_discharges = sum(value)*nhs_trust) %>% 
  mutate(total_discharges = round(total_discharges))

# proportion of BSI and colonised out of total admissions and proportion of colonised out of total discharges
df.adm.disc <- cbind(adm, disc, bsi_adm, col_adm, col_discharges)

df.adm.disc.clean <- df.adm.disc[, !duplicated(names(df.adm.disc))]

extra_summary_table <- df.adm.disc.clean %>% mutate(percent_bsi_adm = round(total_bsi_adm/total_adm*100,2),
                             percent_col_adm = round(total_col_adm/total_adm*100,2),
                             percent_col_disc = round(total_col_discharges/total_discharges*100,2))


extra_summary_table_transposed <- t(extra_summary_table)