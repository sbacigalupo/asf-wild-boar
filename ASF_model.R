library(tidyverse)

setwd("/Users/sonnybacigalupo/Documents/GitHub/asf-wild-boar")

# User define wd
if(Sys.info()["user"]=="adamkuchars" | Sys.info()["user"]=="akucharski" | Sys.info()["user"]=="adamkucharski") {
  setwd("~/Documents/GitHub/asf-wild-boar/")
}

# Set local path ----------------------------------------------
wdir <- getwd()
   
# Run simulation model ----------------------------------------------
# For each function, outputs are saved in 'dir_pick' directory.

n_run_pick <- 100 # model iterations


# Load model functions
source("SEID_model_functions.R")

# Create list to store model iterations

List = list()

# Baseline scenario

base_scenario <- c(beta= 0.75/1500,
                   beta_wbp= 1/30/1500,
                   zeta= 1/7,
                   infection_death= 1/10)


# Run models

for(ii in 1:n_run_pick){
  models <- outbreak(max_time = 365,
                     beta = 0.75/1500,
                     beta_wbp = 1/30/1500,
                     zeta = 1/7,
                     infection_death = 1/10)
  List[[length(List)+1]] = models
}


# Sub-list

# populations
list_pop <- lapply(List, function(x) x$populations) 

# outputs
list_output<- lapply(List, function(x) x$outputs)   



# Plot single outbreak

plot_outbreak(
  scenario = base_scenario,
  n_run = 24
)



 
 # Probability of infection in pigs over outbreak
 
 infect_pig(n_runs = length(list_output))
         
 # Size of epidemic over outbreak (Dead wild boar)

 size(n_runs = length(list_output))
 
 
 

# Plot median outbreak

plot_median(
  
)
