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

n_run_pick <- 10 # model iterations
n_scenario <- 6


# Create array to store model runs

patchNames=c("forest","border","outside")
stateNames=c("S","E","I","D","prob_wb_wb","prob_wb_p")
scenarioNames=c("beta 0.7","beta 1","beta 1.3","beta 1.6","beta 1.9","beta 2.2")
n.patch=length(patchNames)
n.state=length(stateNames)
n.scenario=length(scenarioNames)
time.val=length(1:365)

store <- array(dim=c(n.scenario,length(1:n_run_pick),length(1:365),n.patch,n.state),dimnames=list(scenarioNames,NULL, NULL,patchNames,stateNames))


# Load model functions
source("SEID_model_functions.R")

for(ii in 1:n_run_pick){
  model = outbreak(
    beta = 0.75/1500
    )
  store[ii,,,] = model
}


beta_range <- c(0.7,2) #seq(0.7,2.2,0.3)
store_beta <- array(dim=c(length(beta_range),length(1:n_run_pick),length(1:365),n.patch,n.state),dimnames=list(NULL,NULL, NULL,patchNames,stateNames))

# Loop over beta
for(mm in 1:length(beta_range)){
  # Loop over runs
  for(ii in 1:n_run_pick){
    model_ii <- outbreak(
      beta = beta_range[mm]/1500
    )
  store[mm,ii,,,] <- model_ii
  } # End loop runs
} # End loop beta


# Check debugs
model_ii[,"forest","S"]

store[1,2,1:10,,"S"]

# Plot wild boar numbers and infection probability in single outbreak


plot_outbreak(
  n_run = 6)



 
 # Probability of infection in pigs over outbreak
 
 infect_pig(n_run = n_run_pick)
         
 # Size of epidemic over outbreak (Dead wild boar)

 size(n_run = n_run_pick)
 
 
 

# Plot median outbreak

plot_median(
  
)
