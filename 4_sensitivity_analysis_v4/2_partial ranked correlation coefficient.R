#==========================================================================================================
# This script generates Latin hypercube samples for the Klebsiella care homes model parameters and runs the model with these samples.
#==========================================================================================================
# prepare the environment

# clear the environment
rm(list = ls())

# load required packages
if(!requireNamespace("lhs", quietly = TRUE)) install.packages("lhs") #change all install.package codes to install only when needed!
if(!requireNamespace("sensitivity", quietly = TRUE)) install.packages("sensitivity") #change all install.package codes to install only when needed!

# load the lhs package
require(lhs)
require(sensitivity)
#require(spartan) #the package no longer exists
set.seed(123)

# load required sources
source("4_sensitivity_analysis_v4/sensitivity_inputs/model1_v4_for_sensitivity analysis.R")


#===========================================================================================================
# Feeding in LHS samples into the PRCC function
#===========================================================================================================

#===========================================================================================================
# define inputs and outputs for the PRCC function

# load saved LHS samples from result data frame
results.dataframe.file <- readRDS("4_sensitivity_analysis_v4/run.varied.parameters.rds")
data_for_prcc <- results.dataframe.file 

View(data_for_prcc)

# filter only combinations with equilibrium 
data_for_prcc_eq <- data_for_prcc %>% 
  mutate(difference = total_cases - lag(total_cases)) %>%  filter(time >= max(time)-1) %>% 
  mutate(difference = total_cases - lag(total_cases, 1)) %>%  filter(time >= max(time)-1)

data_for_prcc_eq <- data_for_prcc_eq %>% filter(time >= max(time))

# view the check equilibrium output
View(data_for_prcc_eq) #ok

# filter equilibrium values only 
# tolerance is 4 decimal places so we do not exclude combinations that have actually reached steady state 
# 4 decimal places is considered negligible difference between the final time step and the previous time step 
# and thus if the absolute difference between the final time step and the previous time step is below 5 decimal places we can consider them as having reached equilibrium values
data_for_prcc_eq_clean <- data_for_prcc_eq %>% filter(abs(difference) < 1e-4)

# view the check equilibrium output
View(data_for_prcc_eq_clean) #ok

# get rid of invalid combinations
# increase tolerance to 6 decimal places because sometimes the solver throws a really small "negatives" that are not true negatives but actually 0
data_for_prcc_eq_clean2 <- data_for_prcc_eq_clean %>%
  filter(if_all(8:29, ~abs(.x) > -1e-6)) 

View(data_for_prcc_eq_clean2)

# inputs for PRCC
X <- data_for_prcc_eq_clean[,1:6]
View(X)


# outputs
# include outcome of interest (total number of cases)
# define the outputs for PRCC
Y <- data_for_prcc_eq_clean %>% select(total_cases)


# check that the outputs are in the correct format and order
View(Y)

# save the inputs and outputs for checking later
saveRDS(X, "4_sensitivity_analysis_v4/X_for_prcc.RData")
saveRDS(Y, "4_sensitivity_analysis_v4/Y_for_prcc.RData")

#===========================================================================================================

# run the PRCC function # a two-side PRCC is the default option
prcc_results_no_boot <- pcc(X = X, y = Y, rank = TRUE) # without bootstrapping
prcc_results_with_boot <- pcc(X = X, y = Y, rank = TRUE, nboot = 1000) # with bootstrapping

#===========================================================================================================

# save as data frame 
prcc_df_no_boot <- as.data.frame(prcc_results_no_boot$PRCC)
prcc_df_with_boot <- as.data.frame(prcc_results_with_boot$PRCC)
#===========================================================================================================

# insert column name for the parameters and rename "original" to "PRCC"
prcc_df_no_boot$Variable <- rownames(prcc_df_no_boot) # the variable name can be renamed later to description for a more straightforward visualisation
prcc_df_no_boot$PRCC <- prcc_df_no_boot$original 

prcc_df_with_boot$Variable <- rownames(prcc_df_with_boot) # the variable name can be renamed later to description for a more straightforward visualisation
prcc_df_with_boot$PRCC <- prcc_df_with_boot$original

#===========================================================================================================
# save the results # without confidence intervals #should I include the confidence intervals?
saveRDS(prcc_df_no_boot, file = paste0("4_sensitivity_analysis_v4/prcc_df_no_boot.RData")) # save the results to a file
saveRDS(prcc_df_with_boot, file = paste0("4_sensitivity_analysis_v4/prcc_df_with_boot.RData")) # save the results to a file
#===========================================================================================================

#===========================================================================================================
# Visualise the correlations coefficients in the form of tornado plots
#===========================================================================================================
library(ggplot2)

# assign name to each parameter so the reader can easily interpret the plot
prcc_df_no_boot_reorder <- prcc_df_no_boot %>%
  mutate(variable_name = recode(Variable, 
                                beta_h = "Transmission coefficient (hospital)",
                                gamma_1 = "Hospital discharge rate",
                                kappa_1 = "Progression from colonised state to BSI (community)",
                                epsilon_1 = "Natural clearance rate (community)",
                                kappa_h = "Progression from colonised state to BSI (hospital)",
                                tau_h = "Treatment rate of BSI patients")) %>% 
  mutate(var_name = as.character(variable_name))

# reorder based on absolute values
prcc_df_no_boot_reorder <- prcc_df_no_boot_reorder %>%
  mutate(parameter = forcats::fct_reorder(var_name, abs(PRCC))) #, .desc = TRUE

prcc_df_with_boot_reorder <- prcc_df_with_boot %>%
  mutate(variable_name = recode(Variable, 
                                beta_h = "Transmission coefficient (hospital)",
                                gamma_1 = "Hospital discharge rate",
                                kappa_1 = "Progression from colonised state to BSI (community)",
                                epsilon_1 = "Natural clearance rate (community)",
                                kappa_h = "Progression from colonised state to BSI (hospital)",
                                tau_h = "Treatment rate of BSI patients")) %>% 
  mutate(var_name = as.character(variable_name))

# reorder based on absolute values
prcc_df_with_boot_reorder <- prcc_df_with_boot_reorder %>%
  mutate(parameter = forcats::fct_reorder(var_name, abs(PRCC))) #, .desc = TRUE

# Function to create a tornado plot
tornado_plot <- function(data, title) {
  ggplot(data, aes(x = parameter, y = PRCC)) +
    geom_bar(stat = "identity", fill = "lightgrey") +
    coord_flip() +
    scale_y_continuous(limits = c(-1, 1)) +  
    labs(title = title, x = "", y = "Partial Rank Correlation Coefficient (PRCC)") +
    theme_minimal() + 
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_text(hjust = 1, size = 16),
      axis.text.x = element_text(size =16),
      axis.title = element_text(size =16)
    )
}

# return tornado plots for both no bootstrapping and bootstrapping results
# remember to check whether the relationships are as expected
tornado_plot_no_boot <- tornado_plot(prcc_df_no_boot_reorder, "")
tornado_plot_with_boot <- tornado_plot(prcc_df_with_boot_reorder, "PRCC Results with Bootstrapping")

#===========================================================================================================
# save the plots
# 4/6/2025: the results are not as expected, so I will need to check the inputs and outputs and their combinations again
ggsave("4_sensitivity_analysis_v4/sensitivity_outputs/klebsiella1_prcc_tornado_no_boot.png", tornado_plot_no_boot, width = 10, height = 6) 
ggsave("4_sensitivity_analysis_v4/sensitivity_outputs/klebsiella1_prcc_tornado_with_boot.png", tornado_plot_with_boot, width = 10, height = 6)
#===========================================================================================================

cor(data_for_prcc_eq_clean$gamma_1, data_for_prcc_eq_clean$total_cases, method = "spearman")
