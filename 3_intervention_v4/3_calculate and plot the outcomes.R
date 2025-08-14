# This script is for summarising the vaccination outcomes
#===================================================================================================
# Research Q:  Potential impacts of the hypothetical vaccine on the burden of K. pneumoniae BSI
# Modelling Q: How many cases and deaths would occur without and with intervention?
#===================================================================================================

# clear the environment
rm(list = ls())

#===============================================================================
# 1. package installations # only install the ones that are not installed yet

install_if_not_installed <- function(pkg) {
  install_packages <- pkg[!pkg %in% installed.packages()[, "Package"]]
  if (length(install_packages)) install.packages(install_packages)
}

# install packages if not installed
install_if_not_installed(c("deSolve", "ggplot2", "dplyr", "reshape2", "tidyverse", 
                           "readr", "readxl", "lubridate", "purrr", "data.table", 
                           "writexl", "git2r"))

# load multiple libraries
lapply(c("deSolve", "ggplot2", "dplyr", "reshape2", "tidyverse", "readr", "readxl", 
         "lubridate", "purrr", "data.table", "writexl", "git2r", "scales"), library, character.only = TRUE)

#===============================================================================

#------------------------------------------------------------------------------#
# baseline scenario
# compare no intervention vs. intervention (with 70% coverage - 70% efficacy) 
#------------------------------------------------------------------------------#

# calculate the number of cases averted - cumulative per year
# number of cases without intervention
cases_no_intervention <- read.csv("3_intervention_v4/intervention_outputs/csv/cases_per_year_no_intervention.csv", row.names = NULL)

# number of cases with intervention
cases_with_intervention_baseline <- read.csv("3_intervention_v4/intervention_outputs/csv/cases_per_year_with_intervention2.csv", row.names = NULL)

# create a dataframe to calculate the number of cases averted
merge_cases <- data.frame(year = cases_no_intervention$year, #so the year starts at 1
           no_int_cases = cases_no_intervention$cum_cases,
           with_int_baseline_cases = cases_with_intervention_baseline$cum_cases)

cum_cases_averted <- merge_cases %>% 
  mutate(averted_cases_cum = abs(with_int_baseline_cases - no_int_cases))


# calculate the number of deaths averted - cumulative per year
# number of deaths without intervention
deaths_no_intervention <- read.csv("3_intervention_v4/intervention_outputs/csv/deaths_per_year_no_intervention.csv", row.names = NULL)


# number of deaths with intervention
deaths_with_intervention_baseline <- read.csv("3_intervention_v4/intervention_outputs/csv/deaths_per_year_with_intervention2.csv", row.names = NULL)

# create a dataframe to calculate the number of deaths averted
merge_deaths <- data.frame(year = deaths_no_intervention$year, #so the year starts at 1
                          no_int_deaths = deaths_no_intervention$cum_deaths,
                          with_int_baseline_deaths = deaths_with_intervention_baseline$cum_deaths)

cum_deaths_averted <- merge_deaths %>% 
  mutate(averted_deaths_cum = abs(with_int_baseline_deaths - no_int_deaths))

#-----------------------------------------------------------------------------------#
# calculate the number of cases averted - total over 5 years and average per year
#-----------------------------------------------------------------------------------#
# number of cases averted over 5 years
total_cases_averted <- cum_cases_averted$averted_cases_cum[5]

# percentage of reduction over 5 years
cases_percent_reduction <- cum_cases_averted$averted_cases_cum[5]/cum_cases_averted$no_int_cases[5]

# number of cases averted per year
cases_averted_annual_average <- round(total_cases_averted/5,0)


# calculate the number of deaths averted - total over 5 years and average per year
# number of deaths averted over 5 years
total_deaths_averted <- cum_deaths_averted$averted_deaths_cum[5]

# number of deaths averted per year
deaths_percent_reduction <- cum_deaths_averted$averted_deaths_cum[5]/cum_deaths_averted$no_int_deaths[5]

# number of deaths averted per year
deaths_averted_annual_average <- total_deaths_averted/5

#------------------------------------------------------#
# create a summary table for total cases and deaths 
#------------------------------------------------------#
summary_table <- data.frame(
  item = c("Total number of cases", "Total number of cases", "Total number of deaths", "Total number of deaths"),
  timeframe = c("Over 5 years", "Per year", "Over 5 years", "Per year"),
  without_int = c(cum_cases_averted$no_int_cases[5], cum_cases_averted$no_int_cases[5]/5, 
                  cum_deaths_averted$no_int_deaths[5], cum_deaths_averted$no_int_deaths[5]/5),
  with_int = c(cum_cases_averted$with_int_baseline_cases[5], cum_cases_averted$with_int_baseline_cases[5]/5, 
               cum_deaths_averted$with_int_baseline_deaths[5],  cum_deaths_averted$with_int_baseline_deaths[5]/5),
  abs_reduction = c(total_cases_averted, total_cases_averted/5,
                total_deaths_averted, total_deaths_averted/5)
) %>%  mutate(across(where(is.numeric), round))

summary_table$percentage_reduction <- round((summary_table$abs_reduction/summary_table$without_int)*100,1)


#---------------------------------------------------#
# bsi hospitalisations and hosp-onset cases
#---------------------------------------------------#

#------------------------------------------------------------------------------#
# baseline scenario
# compare no intervention vs. intervention (with 70% coverage - 70% efficacy) 
#------------------------------------------------------------------------------#

# calculate the number of cases averted - cumulative per year
# number of cases without intervention
bsi_hosp_no_intervention <- read.csv("3_intervention_v4/intervention_outputs/csv/bsi_hosp_per_year_no_intervention.csv", row.names = NULL)

# number of cases with intervention
bsi_hosp_with_intervention_baseline <- read.csv("3_intervention_v4/intervention_outputs/csv/bsi_hosp_per_year_with_intervention2.csv", row.names = NULL)

# create a dataframe to calculate the number of cases averted
merge_bsi_hosp <- data.frame(year = bsi_hosp_no_intervention$year, #so the year starts at 1
                          no_int_bsi_hosp = bsi_hosp_no_intervention$cum_bsi_hosp,
                          with_int_baseline_bsi_hosp = bsi_hosp_with_intervention_baseline$cum_bsi_hosp)

cum_bsi_hosp_averted <- merge_bsi_hosp %>% 
  mutate(averted_bsi_hosp_cum = abs(with_int_baseline_bsi_hosp- no_int_bsi_hosp))


# calculate the number of hospital onset cases averted - cumulative per year
# number of hospital onset cases without intervention
hosp_onset_no_intervention <- read.csv("3_intervention_v4/intervention_outputs/csv/hosp_onset_cases_per_year_no_intervention.csv", row.names = NULL)


# number of hospital onset cases with intervention
hosp_onset_with_intervention_baseline <- read.csv("3_intervention_v4/intervention_outputs/csv/hosp_onset_cases_per_year_with_intervention2.csv", row.names = NULL)

# create a dataframe to calculate the number of hospital onset cases averted
merge_hospital_onset_cases <- data.frame(year = hosp_onset_no_intervention$year, #so the year starts at 1
                           no_int_hospital_onset_cases = hosp_onset_no_intervention$cum_hosp_onset,
                           with_int_baseline_hospital_onset_cases = hosp_onset_with_intervention_baseline$cum_hosp_onset)

cum_hosp_onset_averted <- merge_hospital_onset_cases %>% 
  mutate(averted_hosp_onset_cum = abs(with_int_baseline_hospital_onset_cases - no_int_hospital_onset_cases))


#-----------------------------------------------------------------------------------#
# calculate the number of cases averted - total over 5 years and average per year
#-----------------------------------------------------------------------------------#
# number of cases averted over 5 years
total_bsi_hosp_averted <- cum_bsi_hosp_averted$averted_bsi_hosp_cum[5]

# percentage of reduction over 5 years
bsi_hosp_percent_reduction <- cum_bsi_hosp_averted$averted_bsi_hosp_cum[5]/cum_bsi_hosp_averted$no_int_bsi_hosp[5]

# number of cases averted per year
bsi_hosp_averted_annual_average <- round(total_bsi_hosp_averted/5,0)


# calculate the number of deaths averted - total over 5 years and average per year
# number of deaths averted over 5 years
total_hosp_onset_averted <- cum_hosp_onset_averted$averted_hosp_onset_cum[5]

# number of deaths averted per year
hosp_onset_percent_reduction <- cum_hosp_onset_averted$averted_hosp_onset_cum[5]/cum_hosp_onset_averted$no_int_hospital_onset_cases[5]

# number of deaths averted per year
hosp_onset_averted_annual_average <- round(total_hosp_onset_averted/5,0)


#-----------------------------------------------------------------#
# create a summary table for total bsi_hosp and hosp_onset_cases 
#-----------------------------------------------------------------#
summary_table2 <- data.frame(
  item = c("BSI hospitalisations", "BSI hospitalisations", "Hospital-onset bsi_hosp", "Hospital-onset bsi_hosp"),
  timeframe = c("Over 5 years", "Per year", "Over 5 years", "Per year"),
  without_int = c(cum_bsi_hosp_averted$no_int_bsi_hosp[5], cum_bsi_hosp_averted$no_int_bsi_hosp[5]/5, 
                  cum_hosp_onset_averted$no_int_hospital_onset_cases[5], cum_hosp_onset_averted$no_int_hospital_onset_cases[5]/5),
  with_int = c(cum_bsi_hosp_averted$with_int_baseline_bsi_hosp[5], cum_bsi_hosp_averted$with_int_baseline_bsi_hosp[5]/5, 
               cum_hosp_onset_averted$with_int_baseline_hospital_onset_cases[5],  cum_hosp_onset_averted$with_int_baseline_hospital_onset_cases[5]/5),
  abs_reduction = c(total_bsi_hosp_averted, total_bsi_hosp_averted/5,
                    total_hosp_onset_averted, total_hosp_onset_averted/5)
) %>%  mutate(across(where(is.numeric), round))

summary_table2$percentage_reduction <- round((summary_table2$abs_reduction/summary_table2$without_int)*100,1)

