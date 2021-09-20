library(tidyverse)
#rray(NA,dim=c(n_runs,length(time.vals),length(scenario_names)),dimnames=list(NULL,NULL,scenario_names))

setwd("/Users/sonnybacigalupo/Documents/GitHub/asf-wild-boar")

# User define wd
if(Sys.info()["user"]=="adamkuchars" | Sys.info()["user"]=="akucharski" | Sys.info()["user"]=="adamkucharski") {
  setwd("~/Documents/GitHub/asf-wild-boar/")
}

# Set local path ----------------------------------------------
wdir <- getwd()
   
# Run simulation model ----------------------------------------------
# For each function, outputs are saved in 'dir_pick' directory.

n_run_pick <- 1000 # model iterations

# Create array to store model runs

patchNames=c("forest","border","outside")
stateNames=c("S","E","I","D","prob_wb_wb","prob_wb_p")
n.patch=length(patchNames)
n.state=length(stateNames)
time.val=length(1:365)

store <- array(dim=c(length(1:n_run_pick),length(1:365),n.patch,n.state),dimnames=list(NULL,NULL,patchNames,stateNames))

# Load model functions
source("SEID_model_functions.R")


# Run model and store values

for(ii in 1:n_run_pick){
  model = outbreak(max_time =365,
                   beta = 0.75/1500,
                   beta_wbp = 1/30/1500,
                   zeta = 1/7,
                   infection_death = 1/10)
  store[ii,,,] = model
}

# Exploring data

store[2,,,]


store[150,1:65,,"E"] 



# Attempt multiple scenarios
#beta_range  <- c(0.7,2.2)
#
#  for(mm in beta_range){
#    for(ii in 1:n_run_pick){
#    f <-outbreak(max_time = 365,
#                    beta = mm/1500,
#                     beta_wbp = 1/30/1500,
#                    zeta = 1/7,
#                     infection_death = 1/10)
#    List[[length(List)+1]] = f
#}
#}


# Plot wild boar numbers and infection probability in single outbreak

plot_outbreak(
  n_run = 150
)


 
 # Probability of infection in pigs over outbreak
 
 infect_pig(n_run = n_run_pick)
         
 # Size of epidemic over outbreak (Dead wild boar)

 size(n_run = n_run_pick)
 
 
 

# Plot median outbreak

plot_median(
  
)
