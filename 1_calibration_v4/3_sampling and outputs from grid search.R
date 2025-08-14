# GRID SEARCH #
# clear the environment
rm(list = ls())


#---------------------------------------------------#
# output from the selected samples ####
#---------------------------------------------------#
# load multiple libraries
lapply(c("deSolve", "ggplot2", "dplyr", "reshape2", "tidyverse", "readr", "readxl", 
         "lubridate", "data.table", "writexl", "git2r", "scales", "lhs", "sensitivity", "scatterplot3d", "plotly", 
         "grid"), 
       library, character.only = TRUE)

# load the model and input sources
source("1_calibration_v4/calibration_inputs/model1_v4_for_calibration.R")
source("1_calibration_v4/calibration_inputs/inputs_v4_calibration_model1.R")$value

istate1 <- c(S_c1 = S_c1, C_c1 = C_c1, I_c1 = I_c1,
             S_h1 = S_h1, C_h1 = C_h1, I_h1 = I_h1, 
             D_h1 = D_h1, CInc_1a = CInc_1a, CInc_1b = CInc_1b, CInc_1c = CInc_1c, CInc_1d = CInc_1d) 

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

tps1 <- seq(0, 365*200, by = 1) 

# parameter combinations from the first optimisation

# combining the "best combinations from set 1 and set 2 -> narrow the range further
# beta_h : 0.1142 to 0.1204
# kappa_h: 0.00104 to 0.00148
# kappa_1: 0.0000451 to 0.0000456

# pick a few values from here so the samples are not too many
beta_h <- seq(0.1142, 0.1204, by  = 0.0007) 
kappa_h <- seq(0.00104, 0.00148, by  = 0.0002)
kappa_1 <- seq(0.0000451,  0.0000456, by  = 0.0000001)

# final samples for grid search
beta_h <- c(0.1142, 0.1156, 0.1177, 0.1204)
kappa_h <- seq(0.00104, 0.00148, by  = 0.0002)
kappa_1 <- c(0.0000451, 0.0000452, 0.0000455, 0.0000456)

grid <- expand.grid(beta_h = beta_h, kappa_h = kappa_h, kappa_1 = kappa_1)

df.calibration_3 <-NULL
for(i in 1:nrow(grid)){
  parameters1["kappa_1"]<- grid$kappa_1[i];
  parameters1["beta_h"] <- grid$beta_h[i]; # use the calibrated beta_h value
  parameters1["kappa_h"] <- grid$kappa_h[i]; # use the calibrated kappa_h value
  run<-ode(times=tps1, y=istate1, func=klebsiella_1, parms=parameters1, method = "euler"); 
  tmp<-as_tibble(as.data.frame(run)) %>% 
    mutate(incidence_col_hosp = c(0, diff(CInc_1a)),
           incidence_comm_cases = c(0, diff(CInc_1b)),
           bsi_hospitalisation = c(0, diff(CInc_1c)),
           incidence_hosp_cases = c(0, diff(CInc_1d)), 
           incidence_bsi_hosp = bsi_hospitalisation + incidence_hosp_cases) %>%
    mutate(number_col_comm = C_c1, 
           number_col_hosp = C_h1,
           number_bsi_comm = I_c1,
           number_bsi_hosp = I_h1) %>%
    mutate(population_comm = S_c1 + C_c1 + I_c1,
           population_hosp = S_h1 + C_h1 + I_h1,
           prevalence_col_comm = number_col_comm/population_comm,
           kappa_1 = grid$kappa_1[i],
           kappa_h = grid$kappa_h[i],
           beta_h = grid$beta_h[i])%>%
    select(time, kappa_1, kappa_h, beta_h, incidence_col_hosp, incidence_comm_cases, bsi_hospitalisation, incidence_hosp_cases, incidence_bsi_hosp,  
           population_comm, population_hosp, prevalence_col_comm)
  df.calibration_3<-bind_rows(df.calibration_3, tmp)
}

View(df.calibration_3)

# save to csv
##write.csv(df.calibration_3, "1_calibration_v4/calibration_outputs/csv/klebsiella_1_calibration_3.csv", row.names = FALSE)

# check output # plot the outputs of hospital population at equilibrium
cases_calibration_3 <- df.calibration_3 %>% 
  mutate(year = floor((as.integer(time) - 1) / 365) + 1)

View(cases_calibration_3)

nhs_trust <- 134

cases_grouped_by_year_calibration_3 <- cases_calibration_3 %>%
  group_by(beta_h, kappa_h, kappa_1, year) %>% 
  summarise(total_hosp_incidence = sum(incidence_hosp_cases),
            total_bsi_hospitalisation = sum(bsi_hospitalisation),
            prevalence_col_comm = mean(prevalence_col_comm)) %>% 
  mutate(hosp_incidence_all_trust = round(total_hosp_incidence*nhs_trust),
         bsi_hospitalised_all_trust  = round(total_bsi_hospitalisation*nhs_trust),
         total_incidence_all_trust = hosp_incidence_all_trust + bsi_hospitalised_all_trust,
         prop_bsi_hospitalised_all_trust = round((bsi_hospitalised_all_trust / total_incidence_all_trust)* 100,2),
         colonisation_prevalence_comm_percent_all_trust  = round(prevalence_col_comm*100,2)) #perform rounding to get integer

View(cases_grouped_by_year_calibration_3)

# save to csv
#write.csv(cases_grouped_by_year_calibration_3, "1_calibration_v4/calibration_outputs/csv/klebsiella_1_calibration_3_bsi_incidence_summarised_per_year.csv", row.names = FALSE)

# check if the outputs has reached equilibrium
# check that the difference between the maximum time point and the one time point before is 0
check_equilibrium_calibration_3 <- cases_grouped_by_year_calibration_3 %>% 
  mutate(difference = hosp_incidence_all_trust - lag(hosp_incidence_all_trust)) %>%  filter(year >= max(year)-1) %>% 
  mutate(difference = total_bsi_hospitalisation - lag(total_bsi_hospitalisation, 1)) %>%  filter(year >= max(year)-1)

# view the check equilibrium output
View(check_equilibrium_calibration_3) #ok


# save to csv
write.csv(check_equilibrium_calibration_3, "1_calibration_v4/calibration_outputs/csv/klebsiella_1_calibration_3_eq_values.csv", row.names = FALSE)

# turn off scientific notation
options(scipen = 999)

# load the output
calibration_output_1 <- read.csv("1_calibration_v4/calibration_outputs/csv/klebsiella_1_calibration_3_eq_values.csv")

# create faceted horizontal point plot to visualise how the 3 parameter values combinations produce different outputs
# combine parameters into a single label
df <- calibration_output_1 %>% filter(year == max(year)) %>% 
  mutate(param_comb = paste0("βh = ", round(beta_h,4), 
                             ", κh = ", round(kappa_h,5), 
                             ", κ1 = ", round(kappa_1,7)))

View(df)

# transform the data to long format for ggplot
df.long <- df %>% filter(total_incidence_all_trust > 6000, 
                             prop_bsi_hospitalised_all_trust > 56,
                             colonisation_prevalence_comm_percent_all_trust >1.5, colonisation_prevalence_comm_percent_all_trust <3)

View(df.long)

df.long <- df.long %>%
  pivot_longer(cols = c(hosp_incidence_all_trust, bsi_hospitalised_all_trust, colonisation_prevalence_comm_percent_all_trust),
               names_to = "outputs", values_to = "value")

View(df.long)

# rename output
df.long <- df.long %>%
  mutate(outputs = recode(outputs,
                          colonisation_prevalence_comm_percent_all_trust = "colonisation prevalence - community (%)",
                          hosp_incidence_all_trust = "total hospital-onset cases",
                          bsi_hospitalised_all_trust = "total number of BSI hospitalisations"
  ))

View(df.long)

# save as csv
#write.csv(df.long, "1_calibration_v4/calibration_outputs/csv/filtered_combinations1.csv")

# tidy up the order of parameter values combinations so the visual is neat
df.long <- df.long %>%
  mutate(param_comb = factor(param_comb, levels = unique(param_comb)))

View(df.long)

# observed values
obs_hosp_cases <- 2742
obs_bsi_hosp_prop <- 3754
obs_col_comm_prev <- 2


obs_values <- data.frame(
  outputs = c("total hospital-onset cases", "total number of BSI hospitalisations", "colonisation prevalence - community (%)"),
  value = c(obs_hosp_cases, obs_bsi_hosp_prop, obs_col_comm_prev))


# plot
ggplot(df.long, aes(x = value, y = param_comb))+
  geom_point(color = "steelblue", size =2) +
  geom_vline(data = obs_values, aes(xintercept = value), linetype = "dashed", color = "#D55E00") +
  facet_wrap(~ fct_rev(factor(outputs)), scales = "free_x") + 
  theme_minimal() + 
  labs(title = NULL,
       x = "Output value",
       y = "Parameter value combination"
  ) +
  theme(
    axis.text.y = element_text(hjust = 0, size = 16),        # left-align y-axis labels
    axis.text.x = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    strip.text = element_text(face = "bold", size =16),
    
    plot.margin = margin(t = 10, r = 30, b = 10, l = 10), 
    panel.spacing = unit(10, "pt"))

# final selected combinations
#beta_h = 0.1177
#kappa_h = 0.00124
#kappa_1 = 0.0000455

# save plot
ggsave("1_calibration_v4/calibration_outputs/plots/klebsiella_1_calibration_3_outputs_plot.png", width = 18.5, height = 10, dpi = 300)


