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

n_run_pick <- 250 # model iterations


# Load model functions
source("SEID_model_functions.R")




#################################################################
##
## Run model looping over n_runs only
##
#################################################################


# Create array to store model runs

patchNames=c("forest","border","outside")
stateNames=c("S","E","I","D","prob_wb_wb","prob_wb_p")
n.patch=length(patchNames)
n.state=length(stateNames)
time.val=length(1:365)

store <- array(dim=c(length(1:n_run_pick),time.val,n.patch,n.state),dimnames=list(NULL, NULL,patchNames,stateNames))


# Loop over runs
for(ii in 1:n_run_pick){
  model = outbreak(beta = 0.75/1000,
                   beta_wbp = 1/30/1000,
                   base_cull_effort = 1,
                   outbreak_cull_effort = 1)
  store[ii,,,] = model
} # End loop runs




#################################################################
##
## Run model looping over scenarios and n_runs 
##
#################################################################


# Analyse and create array for different scenario values -------------------------------------------

beta_names <- c("0.7","1.6","2.2")
beta_interspecies_names <- c("1/30","7/30","14/30") # median visitation rate, maximum visitation rate, upper 95% CI visitation rate
base_cull_names <- c("none","cull_rate","2*cull_rate")
outbreak_cull_names <- c("none","cull_rate","2*cull_rate","5*cull_rate","10*cull_rate","20*cull_rate")

beta_range <- c(0.7,1.6,2.2) # range of beta 0.7-2.2
beta_interspecies_range <- c(1/30,7/30,14/30) # median, maximum and upper95%CI of visits
base_cull_range <- c(0,1,2) # cull rate (rate to maintain population at 1500) increase by factor of base_cull_range[x]
outbreak_cull_range <- c(0,1,2,5,10,20)

#store_beta <- array(dim=c(length(beta_range),length(1:n_run_pick),time.val,n.patch,n.state),dimnames=list(NULL,NULL,NULL,patchNames,stateNames))
store_scenarios <- array(dim=c(length(base_cull_range),length(outbreak_cull_range),length(beta_interspecies_range), length(beta_range),length(1:n_run_pick),time.val,n.patch,n.state),dimnames=list(base_cull_names,outbreak_cull_names, beta_interspecies_names,beta_names,NULL,NULL,patchNames,stateNames))


# Loop over scenarios
for(pp in 1:length(base_cull_range)){
  for(qq in 1:length(outbreak_cull_range)){
for(nn in 1:length(beta_interspecies_range)){
  for(mm in 1:length(beta_range)){
  # Loop over runs
    for(ii in 1:n_run_pick){
    model_ii <- outbreak(
      base_cull_effort = base_cull_range[pp],
      outbreak_cull_effort = outbreak_cull_range[qq],
      beta = beta_range[mm]/1000,
      beta_wbp = beta_interspecies_range[nn]/1000
    )
    store_scenarios[pp,qq,nn,mm,ii,,,] <- model_ii
    } # End loop runs
  } # End loop beta
} # End loop beta_interspecies
  }# End loop outbreak cull
} # End loop cull
  
# Check debugs
model_ii[,"forest","S"]

store_scenarios[2,1,2,3,1,1:6,,"S"]
store_scenarios["cull_rate","20*cull_rate","1/30","2.2",1,365,"forest","S"]



#######################################
##
## Return outputs for scenarios
##
#######################################

# Probably better to be in a function once working

# store_beta[mm, ii, time.val, patch, state]   Format of store_beta array

# Set stores for outputs

store_outputs <- NULL # Reset each time for loop modified



# array to store calculations. Add to 'dim', add loop and add extra store_beta dimensions as scenarios added

store_outbreak_size <- array(dim=c(length(base_cull_range),length(outbreak_cull_range),length(beta_interspecies_range),length(beta_range),length(1:n_run_pick)),dimnames=list(NULL,NULL,NULL,NULL,NULL))
store_spillover     <- array(dim=c(length(base_cull_range),length(outbreak_cull_range),length(beta_interspecies_range),length(beta_range),length(1:n_run_pick),n.patch),dimnames=list(NULL,NULL,NULL,NULL,NULL,patchNames))
store_survivors     <- array(dim=c(length(base_cull_range),length(outbreak_cull_range),length(beta_interspecies_range),length(beta_range),length(1:n_run_pick)),dimnames=list(NULL,NULL,NULL,NULL,NULL))
store_peak          <- array(dim=c(length(base_cull_range),length(outbreak_cull_range),length(beta_interspecies_range),length(beta_range),length(1:n_run_pick)),dimnames=list(NULL,NULL,NULL,NULL,NULL))
store_peak_time     <- array(dim=c(length(base_cull_range),length(outbreak_cull_range),length(beta_interspecies_range),length(beta_range),length(1:n_run_pick)),dimnames=list(NULL,NULL,NULL,NULL,NULL))
store_death10_time  <- array(dim=c(length(base_cull_range),length(outbreak_cull_range),length(beta_interspecies_range),length(beta_range),length(1:n_run_pick)),dimnames=list(NULL,NULL,NULL,NULL,NULL))
store_death20_time  <- array(dim=c(length(base_cull_range),length(outbreak_cull_range),length(beta_interspecies_range),length(beta_range),length(1:n_run_pick)),dimnames=list(NULL,NULL,NULL,NULL,NULL))
store_spillover_10  <- array(dim=c(length(base_cull_range),length(outbreak_cull_range),length(beta_interspecies_range),length(beta_range),length(1:n_run_pick),n.patch),dimnames=list(NULL,NULL,NULL,NULL,NULL,patchNames))
store_spillover_20  <- array(dim=c(length(base_cull_range),length(outbreak_cull_range),length(beta_interspecies_range),length(beta_range),length(1:n_run_pick),n.patch),dimnames=list(NULL,NULL,NULL,NULL,NULL,patchNames))

  
# Loop over scenario and runs to extract output data

for (ff in 1:length(base_cull_range)){
  for(ee in 1:length(outbreak_cull_range)){
for(gg in 1:length(beta_interspecies_range)){
for (hh in 1:length(beta_range)){
  for (jj in 1:n_run_pick){

    # Calculate and store outbreak sizes
    
    outbreak_size <- sum(max(store_scenarios[ff,ee,gg,hh, jj, ,"forest","D"]),
                        max(store_scenarios[ff,ee,gg,hh, jj, ,"border","D"]),
                        max(store_scenarios[ff,ee,gg,hh, jj, ,"outside","D"])
                        )
    
    store_outbreak_size[ff,ee,gg,hh,jj] <- outbreak_size
    
    # Calculate and store spillover probabilities
    
    spillover_forest <-  sum(store_scenarios[ff,ee,gg,hh, jj, ,"forest","prob_wb_p"])
    spillover_border <- sum(store_scenarios[ff,ee,gg,hh, jj, ,"border","prob_wb_p"])
    spillover_outside <- sum(store_scenarios[ff,ee,gg,hh, jj, ,"outside","prob_wb_p"])
    
    store_spillover[ff,ee,gg,hh,jj,"forest"] <- spillover_forest
    store_spillover[ff,ee,gg,hh,jj,"border"] <- spillover_border
    store_spillover[ff,ee,gg,hh,jj,"outside"] <- spillover_outside
    
    # Calculate and store survivors of outbreak (number of susceptible individuals at 365th day) 
    
    survivors <- sum(store_scenarios[ff,ee,gg,hh,jj, 365,"forest","S"],
                         store_scenarios[ff,ee,gg,hh, jj, 365,"border","S"],
                         store_scenarios[ff,ee,gg,hh, jj, 365,"outside","S"]
                         )
    
    
    store_survivors[ff,ee,gg,hh,jj] <- survivors
    
    
    
    # If outbreak takes off (sum(max(Dead)))>100 store values as above
    
    if (sum(max(store_scenarios[ff,ee,gg,hh,jj, ,"forest","D"]),
            max(store_scenarios[ff,ee,gg,hh,jj, ,"border","D"]),
            max(store_scenarios[ff,ee,gg,hh,jj, ,"outside","D"])
    )>100){
    
      # Calculate and store time at which epidemic peak is reached (time at maximum number of infected individuals))
      # (Exclude no_outbreak)
      
    max_I <- max(store_scenarios[ff,ee,gg,hh, jj, ,"forest","I"]) # maximum daily number of infected individuals
    
    t_max_I <- match(max_I,store_scenarios[ff,ee,gg,hh,jj,,"forest","I"]) # Time when max_I reached
    
    store_peak[ff,ee,gg,hh,jj] <- max_I
    store_peak_time[ff,ee,gg,hh,jj] <- t_max_I
    
    
    
    # Calculate and store time at which number of deaths in forest is 10% and 20% of starting population (time of likely detection in wild boar))
    # (Exclude no_outbreak)
    
    time_death10 <- min(which(store_scenarios[ff,ee,gg,hh,jj,,"forest","D"]>100))
    time_death20 <- min(which(store_scenarios[ff,ee,gg,hh,jj,,"forest","D"]>200))
    
    time_death10[!is.finite(time_death10)] <- 1e-20 # change inf -> 1e-20
    time_death20[!is.finite(time_death20)] <- 1e-20 # change inf -> 1e-20
    
    store_death10_time[ff,ee,gg,hh,jj] <- time_death10
    store_death20_time[ff,ee,gg,hh,jj] <- time_death20
    
    
    # Calculate and store spillover probabilities when deaths reach 10% of population
    # (Exclude no_outbreak)
    
    spillover_forest_10 <- sum(store_scenarios[ff,ee,gg,hh, jj, 1:time_death10,"forest","prob_wb_p"])
    spillover_border_10 <- sum(store_scenarios[ff,ee,gg,hh, jj, 1:time_death10,"border","prob_wb_p"])
    spillover_outside_10 <- sum(store_scenarios[ff,ee,gg,hh, jj, 1:time_death10,"outside","prob_wb_p"])
    
    store_spillover_10[ff,ee,gg,hh,jj,"forest"] <- spillover_forest_10
    store_spillover_10[ff,ee,gg,hh,jj,"border"] <- spillover_border_10
    store_spillover_10[ff,ee,gg,hh,jj,"outside"] <- spillover_outside_10
    
    # Calculate and store spillover probabilities when deaths reach 20% of population
    # (Exclude no_outbreak)
    
    spillover_forest_20 <- sum(store_scenarios[ff,ee,gg,hh, jj, 1:time_death20,"forest","prob_wb_p"])
    spillover_border_20 <- sum(store_scenarios[ff,ee,gg,hh, jj, 1:time_death20,"border","prob_wb_p"])
    spillover_outside_20 <- sum(store_scenarios[ff,ee,gg,hh, jj, 1:time_death20,"outside","prob_wb_p"])
    
    store_spillover_20[ff,ee,gg,hh,jj,"forest"] <- spillover_forest_20
    store_spillover_20[ff,ee,gg,hh,jj,"border"] <- spillover_border_20
    store_spillover_20[ff,ee,gg,hh,jj,"outside"] <- spillover_outside_20
    
    } else { # If outbreak does not take off, assign improbably small value to exclude from analysis
      store_peak[ff,ee,gg,hh,jj] <- 1e-20
      store_peak_time[ff,ee,gg,hh,jj] <- 1e-20
      store_death10_time[ff,ee,gg,hh,jj] <- 1e-20
      store_death20_time[ff,ee,gg,hh,jj] <- 1e-20
      store_spillover_10[ff,ee,gg,hh,jj,"forest"] <- 1e-20
      store_spillover_10[ff,ee,gg,hh,jj,"border"] <- 1e-20
      store_spillover_10[ff,ee,gg,hh,jj,"outside"] <- 1e-20
      store_spillover_20[ff,ee,gg,hh,jj,"forest"] <- 1e-20
      store_spillover_20[ff,ee,gg,hh,jj,"border"] <- 1e-20
      store_spillover_20[ff,ee,gg,hh,jj,"outside"] <- 1e-20
    }
    
  } # end loop run
  
  
  
  # Calculate number of non-outbreaks
  
  no_outbreak <- sum(store_outbreak_size[ff,ee,gg,hh,] < 100)
  
  
  
  # Calculate number of extinctions
  
  extinctions <- sum(store_survivors[ff,ee,gg,hh,] == 0)
  
  
  
  # Index scenarios for table
  
  beta = beta_range[hh] 
  beta_interspecies = beta_interspecies_range[gg]
  base_cull_effort = base_cull_range[ff]
  outbreak_cull_effort = outbreak_cull_range[ee]
  
  
  
#--------------------------------------------------------------------------------------
  # Store outputs
  
  # (Remove non-outbreaks from store_peak, store_peak_time, store_death_time and spillover_10/20 by insisting value must be >1e-20)

  store_outputs <- rbind(store_outputs,c(base_cull_effort,outbreak_cull_effort,beta_interspecies,beta,# scenarios
                                         median(store_outbreak_size[ff,ee,gg,hh,]), # size of outbreak
                                         quantile(store_outbreak_size[ff,ee,gg,hh,],c(0.25,0.75)),
                                         min(store_outbreak_size[ff,ee,gg,hh,]),
                                         max(store_outbreak_size[ff,ee,gg,hh,]),
                                         quantile(store_outbreak_size[ff,ee,gg,hh,],c(0.025,0.975)),
                                         median(store_spillover[ff,ee,gg,hh,,"forest"]), # probability boar-pig transmission
                                         quantile(store_spillover[ff,ee,gg,hh,,"forest"],c(0.25,0.75)),
                                         min(store_spillover[ff,ee,gg,hh,,"forest"]),
                                         max(store_spillover[ff,ee,gg,hh,,"forest"]),
                                         quantile(store_spillover[ff,ee,gg,hh,,"forest"],c(0.025,0.975)),
                                         median(store_spillover[ff,ee,gg,hh,,"border"]),
                                         quantile(store_spillover[ff,ee,gg,hh,,"border"],c(0.25,0.75)),
                                         min(store_spillover[ff,ee,gg,hh,,"border"]),
                                         max(store_spillover[ff,ee,gg,hh,,"border"]),
                                         quantile(store_spillover[ff,ee,gg,hh,,"border"],c(0.025,0.975)),
                                         median(store_spillover[ff,ee,gg,hh,,"outside"]),
                                         quantile(store_spillover[ff,ee,gg,hh,,"outside"],c(0.25,0.75)),
                                         min(store_spillover[ff,ee,gg,hh,,"outside"]),
                                         max(store_spillover[ff,ee,gg,hh,,"outside"]),
                                         quantile(store_spillover[ff,ee,gg,hh,,"outside"],c(0.025,0.975)),
                                         no_outbreak, # non-outbreaks
                                         (no_outbreak/(max(n_run_pick)))*100,
                                         extinctions, # extinctions
                                         (extinctions/(max(n_run_pick)))*100,
                                         median(store_peak[ff,ee,gg,hh,][store_peak[ff,ee,gg,hh,]>1e-20]), # peak daily infections removing non-outbreaks
                                         quantile(store_peak[ff,ee,gg,hh,][store_peak[ff,ee,gg,hh,]>1e-20],c(0.25,0.75)),
                                         min(store_peak[ff,ee,gg,hh,][store_peak[ff,ee,gg,hh,]>1e-20]),
                                         max(store_peak[ff,ee,gg,hh,][store_peak[ff,ee,gg,hh,]>1e-20]),
                                         quantile(store_peak[ff,ee,gg,hh,][store_peak[ff,ee,gg,hh,]>1e-20],c(0.025,0.975)),
                                         n_run_pick-length(store_peak[ff,ee,gg,hh,][store_peak[ff,ee,gg,hh,]>1e-20]),
                                         median(store_peak_time[ff,ee,gg,hh,][store_peak_time[ff,ee,gg,hh,]>1e-20]), # peak time of daily infections removing non-outbreaks
                                         quantile(store_peak_time[ff,ee,gg,hh,][store_peak_time[ff,ee,gg,hh,]>1e-20],c(0.25,0.75)),
                                         min(store_peak_time[ff,ee,gg,hh,][store_peak_time[ff,ee,gg,hh,]>1e-20]),
                                         max(store_peak_time[ff,ee,gg,hh,][store_peak_time[ff,ee,gg,hh,]>1e-20]),
                                         quantile(store_peak_time[ff,ee,gg,hh,][store_peak_time[ff,ee,gg,hh,]>1e-20],c(0.025,0.975)),
                                         n_run_pick-length(store_peak_time[ff,ee,gg,hh,][store_peak_time[ff,ee,gg,hh,]>1e-20]),
                                         median(store_death10_time[ff,ee,gg,hh,][store_death10_time[ff,ee,gg,hh,]>1e-20]), # day when deaths>10% removing non-outbreaks
                                         quantile(store_death10_time[ff,ee,gg,hh,][store_death10_time[ff,ee,gg,hh,]>1e-20],c(0.25,0.75)),
                                         min(store_death10_time[ff,ee,gg,hh,][store_death10_time[ff,ee,gg,hh,]>1e-20]),
                                         max(store_death10_time[ff,ee,gg,hh,][store_death10_time[ff,ee,gg,hh,]>1e-20]),
                                         quantile(store_death10_time[ff,ee,gg,hh,][store_death10_time[ff,ee,gg,hh,]>1e-20],c(0.025,0.975)),
                                         n_run_pick-length(store_death20_time[ff,ee,gg,hh,][store_death20_time[ff,ee,gg,hh,]>1e-20]),
                                         median(store_death20_time[ff,ee,gg,hh,][store_death20_time[ff,ee,gg,hh,]>1e-20]), # day when deaths>20% removing non-outbreaks
                                         quantile(store_death20_time[ff,ee,gg,hh,][store_death20_time[ff,ee,gg,hh,]>1e-20],c(0.25,0.75)),
                                         min(store_death20_time[ff,ee,gg,hh,][store_death20_time[ff,ee,gg,hh,]>1e-20]),
                                         max(store_death20_time[ff,ee,gg,hh,][store_death20_time[ff,ee,gg,hh,]>1e-20]),
                                         quantile(store_death20_time[ff,ee,gg,hh,][store_death20_time[ff,ee,gg,hh,]>1e-20],c(0.025,0.975)),
                                         n_run_pick-length(store_death20_time[ff,ee,gg,hh,][store_death20_time[ff,ee,gg,hh,]>1e-20]),
                                         median(store_spillover_10[ff,ee,gg,hh,,"forest"][store_spillover_10[ff,ee,gg,hh,,"forest"]>1e-20]), # probability boar-pig transmission before 10% deaths removing non-outbreaks
                                         quantile(store_spillover_10[ff,ee,gg,hh,,"forest"][store_spillover_10[ff,ee,gg,hh,,"forest"]>1e-20],c(0.25,0.75)),
                                         min(store_spillover_10[ff,ee,gg,hh,,"forest"][store_spillover_10[ff,ee,gg,hh,,"forest"]>1e-20]),
                                         max(store_spillover_10[ff,ee,gg,hh,,"forest"][store_spillover_10[ff,ee,gg,hh,,"forest"]>1e-20]),
                                         quantile(store_spillover_10[ff,ee,gg,hh,,"forest"][store_spillover_10[ff,ee,gg,hh,,"forest"]>1e-20],c(0.025,0.975)),
                                         n_run_pick-length(store_spillover_10[ff,ee,gg,hh,,"forest"][store_spillover_10[ff,ee,gg,hh,,"forest"]>1e-20]),
                                         median(store_spillover_10[ff,ee,gg,hh,,"border"][store_spillover_10[ff,ee,gg,hh,,"border"]>1e-20]),
                                         quantile(store_spillover_10[ff,ee,gg,hh,,"border"][store_spillover_10[ff,ee,gg,hh,,"border"]>1e-20],c(0.25,0.75)),
                                         min(store_spillover_10[ff,ee,gg,hh,,"border"][store_spillover_10[ff,ee,gg,hh,,"border"]>1e-20]),
                                         max(store_spillover_10[ff,ee,gg,hh,,"border"][store_spillover_10[ff,ee,gg,hh,,"border"]>1e-20]),
                                         quantile(store_spillover_10[ff,ee,gg,hh,,"border"][store_spillover_10[ff,ee,gg,hh,,"border"]>1e-20],c(0.025,0.975)),
                                         n_run_pick-length(store_spillover_10[ff,ee,gg,hh,,"border"][store_spillover_10[ff,ee,gg,hh,,"border"]>1e-20]),
                                         median(store_spillover_10[ff,ee,gg,hh,,"outside"][store_spillover_10[ff,ee,gg,hh,,"outside"]>1e-20]),
                                         quantile(store_spillover_10[ff,ee,gg,hh,,"outside"][store_spillover_10[ff,ee,gg,hh,,"outside"]>1e-20],c(0.25,0.75)),
                                         min(store_spillover_10[ff,ee,gg,hh,,"outside"][store_spillover_10[ff,ee,gg,hh,,"outside"]>1e-20]),
                                         max(store_spillover_10[ff,ee,gg,hh,,"outside"][store_spillover_10[ff,ee,gg,hh,,"outside"]>1e-20]),
                                         quantile(store_spillover_10[ff,ee,gg,hh,,"outside"][store_spillover_10[ff,ee,gg,hh,,"outside"]>1e-20],c(0.025,0.975)),
                                         n_run_pick-length(store_spillover_10[ff,ee,gg,hh,,"outside"][store_spillover_10[ff,ee,gg,hh,,"outside"]>1e-20]),
                                         median(store_spillover_20[ff,ee,gg,hh,,"forest"][store_spillover_20[ff,ee,gg,hh,,"forest"]>1e-20]), # probability boar-pig transmission before 20% deaths removing non-outbreaks
                                         quantile(store_spillover_20[ff,ee,gg,hh,,"forest"][store_spillover_20[ff,ee,gg,hh,,"forest"]>1e-20],c(0.25,0.75)),
                                         min(store_spillover_20[ff,ee,gg,hh,,"forest"][store_spillover_20[ff,ee,gg,hh,,"forest"]>1e-20]),
                                         max(store_spillover_20[ff,ee,gg,hh,,"forest"][store_spillover_20[ff,ee,gg,hh,,"forest"]>1e-20]),
                                         quantile(store_spillover_20[ff,ee,gg,hh,,"forest"][store_spillover_20[ff,ee,gg,hh,,"forest"]>1e-20],c(0.025,0.975)),
                                         n_run_pick-length(store_spillover_20[ff,ee,gg,hh,,"forest"][store_spillover_20[ff,ee,gg,hh,,"forest"]>1e-20]),
                                         median(store_spillover_20[ff,ee,gg,hh,,"border"][store_spillover_20[ff,ee,gg,hh,,"border"]>1e-20]),
                                         quantile(store_spillover_20[ff,ee,gg,hh,,"border"][store_spillover_20[ff,ee,gg,hh,,"border"]>1e-20],c(0.25,0.75)),
                                         min(store_spillover_20[ff,ee,gg,hh,,"border"][store_spillover_20[ff,ee,gg,hh,,"border"]>1e-20]),
                                         max(store_spillover_20[ff,ee,gg,hh,,"border"][store_spillover_20[ff,ee,gg,hh,,"border"]>1e-20]),
                                         quantile(store_spillover_20[ff,ee,gg,hh,,"border"][store_spillover_20[ff,ee,gg,hh,,"border"]>1e-20],c(0.025,0.975)),
                                         n_run_pick-length(store_spillover_20[ff,ee,gg,hh,,"border"][store_spillover_20[ff,ee,gg,hh,,"border"]>1e-20]),
                                         median(store_spillover_20[ff,ee,gg,hh,,"outside"][store_spillover_20[ff,ee,gg,hh,,"outside"]>1e-20]),
                                         quantile(store_spillover_20[ff,ee,gg,hh,,"outside"][store_spillover_20[ff,ee,gg,hh,,"outside"]>1e-20],c(0.25,0.75)),
                                         min(store_spillover_20[ff,ee,gg,hh,,"outside"][store_spillover_20[ff,ee,gg,hh,,"outside"]>1e-20]),
                                         max(store_spillover_20[ff,ee,gg,hh,,"outside"][store_spillover_20[ff,ee,gg,hh,,"outside"]>1e-20]),
                                         quantile(store_spillover_20[ff,ee,gg,hh,,"outside"][store_spillover_20[ff,ee,gg,hh,,"outside"]>1e-20],c(0.025,0.975)),
                                         n_run_pick-length(store_spillover_20[ff,ee,gg,hh,,"outside"][store_spillover_20[ff,ee,gg,hh,,"outside"]>1e-20]))
                                         
  )
                         
  
  
  # Convert
  
  store_outputs <- as_tibble(store_outputs)
  names(store_outputs) <- c("base_cull_effort",
                            "outbreak_cull_effort",
                            "beta_interspeciesx1000",
                            "betax1000",
                            
                            "median_outbreak_size",
                            "quantile0.25_outbreak_size",
                            "quantile0.75_outbreak_size",
                            "min_outbreak_size", 
                            "max_outbreak_size",
                            "quantile0.025_outbreak_size",
                            "quantile0.975_outbreak_size",
                            
                            "median_spillover_forest",
                            "quantile0.25_spillover_forest",
                            "quantile0.75_spillover_forest",
                            "min_spillover_forest", 
                            "max_spillover_forest",
                            "quantile0.025_spillover_forest",
                            "quantile0.975_spillover_forest",
                            
                            "median_spillover_border",
                            "quantile0.25_spillover_border",
                            "quantile0.75_spillover_border",
                            "min_spillover_border", 
                            "max_spillover_border",
                            "quantile0.025_spillover_border",
                            "quantile0.975_spillover_border",
                            
                            "median_spillover_outside",
                            "quantile0.25_spillover_outside",
                            "quantile0.75_spillover_outside",
                            "min_spillover_outside", 
                            "max_spillover_outside",
                            "quantile0.025_spillover_outside",
                            "quantile0.975_spillover_outside",
                            
                            "number_no_outbreak",
                            "percent_no_outbreak",
                            
                            "number_extinctions",
                            "percent_extinctions",
                            
                            "median_peak_infection",
                            "quantile0.25_peak_infection",
                            "quantile0.75_peak_infection",
                            "min_peak_infection", 
                            "max_peak_infection",
                            "quantile0.025_peak_infection",
                            "quantile0.975_peak_infection",
                            "check_peak_infection",
                            
                            "median_time_peak_infection",
                            "quantile0.25_time_peak_infection",
                            "quantile0.75_time_peak_infection",
                            "min_time_peak_infection", 
                            "max_time_peak_infection",
                            "quantile0.025_time_peak_infection",
                            "quantile0.975_time_peak_infection",
                            "check_time_peak_infection",
                            
                            "median_time_10death",
                            "quantile0.25_time_10death",
                            "quantile0.75_time_10death",
                            "min_time_10death", 
                            "max_time_10death",
                            "quantile0.025_time_10death",
                            "quantile0.975_time_10death",
                            "check_time_10death",
                            
                            "median_time_20death",
                            "quantile0.25_time_20death",
                            "quantile0.75_time_20death",
                            "min_time_20death", 
                            "max_time_20death",
                            "quantile0.025_time_20death",
                            "quantile0.975_time_20death",
                            "check_time_20death",
                            
                            "median_spillover10_forest",
                            "quantile0.25_spillover10_forest",
                            "quantile0.75_spillover10_forest",
                            "min_spillover10_forest", 
                            "max_spillover10_forest",
                            "quantile0.025_spillover10_forest",
                            "quantile0.975_spillover10_forest",
                            "check_spillover10_forest",
                            
                            "median_spillover10_border",
                            "quantile0.25_spillover10_border",
                            "quantile0.75_spillover10_border",
                            "min_spillover10_border", 
                            "max_spillover10_border",
                            "quantile0.025_spillover10_border",
                            "quantile0.975_spillover10_border",
                            "check_spillover10_border",
                            
                            "median_spillover10_outside",
                            "quantile0.25_spillover10_outside",
                            "quantile0.75_spillover10_outside",
                            "min_spillover10_outside", 
                            "max_spillover10_outside",
                            "quantile0.025_spillover10_outside",
                            "quantile0.975_spillover10_outside",
                            "check_spillover10_outside",
                            
                            "median_spillover20_forest",
                            "quantile0.25_spillover20_forest",
                            "quantile0.75_spillover20_forest",
                            "min_spillover20_forest", 
                            "max_spillover20_forest",
                            "quantile0.025_spillover20_forest",
                            "quantile0.975_spillover20_forest",
                            "check_spillover20_forest",
                            
                            "median_spillover20_border",
                            "quantile0.25_spillover20_border",
                            "quantile0.75_spillover20_border",
                            "min_spillover20_border", 
                            "max_spillover20_border",
                            "quantile0.025_spillover20_border",
                            "quantile0.975_spillover20_border",
                            "check_spillover20_border",
                            
                            "median_spillover20_outside",
                            "quantile0.25_spillover20_outside",
                            "quantile0.75_spillover20_outside",
                            "min_spillover20_outside", 
                            "max_spillover20_outside",
                            "quantile0.025_spillover20_outside",
                            "quantile0.975_spillover20_outside",
                            "check_spillover20_outside")
  
  store_outputs$base_cull_effort<-as.numeric(store_outputs$base_cull_effort)
  store_outputs$outbreak_cull_effort<-as.numeric(store_outputs$outbreak_cull_effort)
  store_outputs$beta_interspeciesx1000<-as.numeric(store_outputs$beta_interspeciesx1000)
  store_outputs$betax1000<-as.numeric(store_outputs$betax1000)
  
  store_outputs$median_outbreak_size<-as.numeric(store_outputs$median_outbreak_size)
  store_outputs$quantile0.25_outbreak_size<-as.numeric(store_outputs$quantile0.25_outbreak_size)
  store_outputs$quantile0.75_outbreak_size<-as.numeric(store_outputs$quantile0.75_outbreak_size)
  store_outputs$min_outbreak_size<-as.numeric(store_outputs$min_outbreak_size) 
  store_outputs$max_outbreak_size<-as.numeric(store_outputs$max_outbreak_size)
  store_outputs$quantile0.025_outbreak_size<-as.numeric(store_outputs$quantile0.025_outbreak_size)
  store_outputs$quantile0.975_outbreak_size<-as.numeric(store_outputs$quantile0.975_outbreak_size)
  
  store_outputs$median_spillover_forest<-as.numeric(store_outputs$median_spillover_forest)
  store_outputs$quantile0.25_spillover_forest<-as.numeric(store_outputs$quantile0.25_spillover_forest)
  store_outputs$quantile0.75_spillover_forest<-as.numeric(store_outputs$quantile0.75_spillover_forest)
  store_outputs$min_spillover_forest<-as.numeric(store_outputs$min_spillover_forest) 
  store_outputs$max_spillover_forest<-as.numeric(store_outputs$max_spillover_forest)
  store_outputs$quantile0.025_spillover_forest<-as.numeric(store_outputs$quantile0.025_spillover_forest)
  store_outputs$quantile0.975_spillover_forest<-as.numeric(store_outputs$quantile0.975_spillover_forest)
  
  store_outputs$median_spillover_border<-as.numeric(store_outputs$median_spillover_border)
  store_outputs$quantile0.25_spillover_border<-as.numeric(store_outputs$quantile0.25_spillover_border)
  store_outputs$quantile0.75_spillover_border<-as.numeric(store_outputs$quantile0.75_spillover_border)
  store_outputs$min_spillover_border<-as.numeric(store_outputs$min_spillover_border) 
  store_outputs$max_spillover_border<-as.numeric(store_outputs$max_spillover_border)
  store_outputs$quantile0.025_spillover_border<-as.numeric(store_outputs$quantile0.025_spillover_border)
  store_outputs$quantile0.975_spillover_border<-as.numeric(store_outputs$quantile0.975_spillover_border)
  
  store_outputs$median_spillover_outside<-as.numeric(store_outputs$median_spillover_outside)
  store_outputs$quantile0.25_spillover_outside<-as.numeric(store_outputs$quantile0.25_spillover_outside)
  store_outputs$quantile0.75_spillover_outside<-as.numeric(store_outputs$quantile0.75_spillover_outside)
  store_outputs$min_spillover_outside<-as.numeric(store_outputs$min_spillover_outside) 
  store_outputs$max_spillover_outside<-as.numeric(store_outputs$max_spillover_outside)
  store_outputs$quantile0.025_spillover_outside<-as.numeric(store_outputs$quantile0.025_spillover_outside)
  store_outputs$quantile0.975_spillover_outside<-as.numeric(store_outputs$quantile0.975_spillover_outside)
  
  store_outputs$number_no_outbreak<-as.numeric(store_outputs$number_no_outbreak)
  store_outputs$percent_no_outbreak<-as.numeric(store_outputs$percent_no_outbreak)
  
  store_outputs$number_extinctions<-as.numeric(store_outputs$number_extinctions)
  store_outputs$percent_extinctions<-as.numeric(store_outputs$percent_extinctions)
  
  store_outputs$median_peak_infection<-as.numeric(store_outputs$median_peak_infection)
  store_outputs$quantile0.25_peak_infection<-as.numeric(store_outputs$quantile0.25_peak_infection)
  store_outputs$quantile0.75_peak_infection<-as.numeric(store_outputs$quantile0.75_peak_infection)
  store_outputs$min_peak_infection<-as.numeric(store_outputs$min_peak_infection) 
  store_outputs$max_peak_infection<-as.numeric(store_outputs$max_peak_infection)
  store_outputs$quantile0.025_peak_infection<-as.numeric(store_outputs$quantile0.025_peak_infection)
  store_outputs$quantile0.975_peak_infection<-as.numeric(store_outputs$quantile0.975_peak_infection)
  store_outputs$check_peak_infection<-as.numeric(store_outputs$check_peak_infection)
  
  store_outputs$median_time_peak_infection<-as.numeric(store_outputs$median_time_peak_infection)
  store_outputs$quantile0.25_time_peak_infection<-as.numeric(store_outputs$quantile0.25_time_peak_infection)
  store_outputs$quantile0.75_time_peak_infection<-as.numeric(store_outputs$quantile0.75_time_peak_infection)
  store_outputs$min_time_peak_infection<-as.numeric(store_outputs$min_time_peak_infection) 
  store_outputs$max_time_peak_infection<-as.numeric(store_outputs$max_time_peak_infection)
  store_outputs$quantile0.025_time_peak_infection<-as.numeric(store_outputs$quantile0.025_time_peak_infection)
  store_outputs$quantile0.975_time_peak_infection<-as.numeric(store_outputs$quantile0.975_time_peak_infection)
  store_outputs$check_time_peak_infection<-as.numeric(store_outputs$check_time_peak_infection)
  
  store_outputs$median_time_10death<-as.numeric(store_outputs$median_time_10death)
  store_outputs$quantile0.25_time_10death<-as.numeric(store_outputs$quantile0.25_time_10death)
  store_outputs$quantile0.75_time_10death<-as.numeric(store_outputs$quantile0.75_time_10death)
  store_outputs$min_time_10death<-as.numeric(store_outputs$min_time_10death) 
  store_outputs$max_time_10death<-as.numeric(store_outputs$max_time_10death)
  store_outputs$quantile0.025_time_10death<-as.numeric(store_outputs$quantile0.025_time_10death)
  store_outputs$quantile0.975_time_10death<-as.numeric(store_outputs$quantile0.975_time_10death)
  store_outputs$check_time_10death<-as.numeric(store_outputs$check_time_10death)
  
  store_outputs$median_time_20death<-as.numeric(store_outputs$median_time_20death)
  store_outputs$quantile0.25_time_20death<-as.numeric(store_outputs$quantile0.25_time_20death)
  store_outputs$quantile0.75_time_20death<-as.numeric(store_outputs$quantile0.75_time_20death)
  store_outputs$min_time_20death<-as.numeric(store_outputs$min_time_20death) 
  store_outputs$max_time_20death<-as.numeric(store_outputs$max_time_20death)
  store_outputs$quantile0.025_time_20death<-as.numeric(store_outputs$quantile0.025_time_20death)
  store_outputs$quantile0.975_time_20death<-as.numeric(store_outputs$quantile0.975_time_20death)
  store_outputs$check_time_20death<-as.numeric(store_outputs$check_time_20death)
  
  store_outputs$median_spillover10_forest<-as.numeric(store_outputs$median_spillover10_forest)
  store_outputs$quantile0.25_spillover10_forest<-as.numeric(store_outputs$quantile0.25_spillover10_forest)
  store_outputs$quantile0.75_spillover10_forest<-as.numeric(store_outputs$quantile0.75_spillover10_forest)
  store_outputs$min_spillover10_forest<-as.numeric(store_outputs$min_spillover10_forest) 
  store_outputs$max_spillover10_forest<-as.numeric(store_outputs$max_spillover10_forest)
  store_outputs$quantile0.025_spillover10_forest<-as.numeric(store_outputs$quantile0.025_spillover10_forest)
  store_outputs$quantile0.975_spillover10_forest<-as.numeric(store_outputs$quantile0.975_spillover10_forest)
  store_outputs$check_spillover10_forest<-as.numeric(store_outputs$check_spillover10_forest)
  
  store_outputs$median_spillover10_border<-as.numeric(store_outputs$median_spillover10_border)
  store_outputs$quantile0.25_spillover10_border<-as.numeric(store_outputs$quantile0.25_spillover10_border)
  store_outputs$quantile0.75_spillover10_border<-as.numeric(store_outputs$quantile0.75_spillover10_border)
  store_outputs$min_spillover10_border<-as.numeric(store_outputs$min_spillover10_border) 
  store_outputs$max_spillover10_border<-as.numeric(store_outputs$max_spillover10_border)
  store_outputs$quantile0.025_spillover10_border<-as.numeric(store_outputs$quantile0.025_spillover10_border)
  store_outputs$quantile0.975_spillover10_border<-as.numeric(store_outputs$quantile0.975_spillover10_border)
  store_outputs$check_spillover10_border<-as.numeric(store_outputs$check_spillover10_border)
  
  store_outputs$median_spillover10_outside<-as.numeric(store_outputs$median_spillover10_outside)
  store_outputs$quantile0.25_spillover10_outside<-as.numeric(store_outputs$quantile0.25_spillover10_outside)
  store_outputs$quantile0.75_spillover10_outside<-as.numeric(store_outputs$quantile0.75_spillover10_outside)
  store_outputs$min_spillover10_outside<-as.numeric(store_outputs$min_spillover10_outside) 
  store_outputs$max_spillover10_outside<-as.numeric(store_outputs$max_spillover10_outside)
  store_outputs$quantile0.025_spillover10_outside<-as.numeric(store_outputs$quantile0.025_spillover10_outside)
  store_outputs$quantile0.975_spillover10_outside<-as.numeric(store_outputs$quantile0.975_spillover10_outside)
  store_outputs$check_spillover10_outside<-as.numeric(store_outputs$check_spillover10_outside)
  
  store_outputs$median_spillover20_forest<-as.numeric(store_outputs$median_spillover20_forest)
  store_outputs$quantile0.25_spillover20_forest<-as.numeric(store_outputs$quantile0.25_spillover20_forest)
  store_outputs$quantile0.75_spillover20_forest<-as.numeric(store_outputs$quantile0.75_spillover20_forest)
  store_outputs$min_spillover20_forest<-as.numeric(store_outputs$min_spillover20_forest) 
  store_outputs$max_spillover20_forest<-as.numeric(store_outputs$max_spillover20_forest)
  store_outputs$quantile0.025_spillover20_forest<-as.numeric(store_outputs$quantile0.025_spillover20_forest)
  store_outputs$quantile0.975_spillover20_forest<-as.numeric(store_outputs$quantile0.975_spillover20_forest)
  store_outputs$check_spillover20_forest<-as.numeric(store_outputs$check_spillover20_forest)
  
  store_outputs$median_spillover20_border<-as.numeric(store_outputs$median_spillover20_border)
  store_outputs$quantile0.25_spillover20_border<-as.numeric(store_outputs$quantile0.25_spillover20_border)
  store_outputs$quantile0.75_spillover20_border<-as.numeric(store_outputs$quantile0.75_spillover20_border)
  store_outputs$min_spillover20_border<-as.numeric(store_outputs$min_spillover20_border) 
  store_outputs$max_spillover20_border<-as.numeric(store_outputs$max_spillover20_border)
  store_outputs$quantile0.025_spillover20_border<-as.numeric(store_outputs$quantile0.025_spillover20_border)
  store_outputs$quantile0.975_spillover20_border<-as.numeric(store_outputs$quantile0.975_spillover20_border)
  store_outputs$check_spillover20_border<-as.numeric(store_outputs$check_spillover20_border)
  
  store_outputs$median_spillover20_outside<-as.numeric(store_outputs$median_spillover20_outside)
  store_outputs$quantile0.25_spillover20_outside<-as.numeric(store_outputs$quantile0.25_spillover20_outside)
  store_outputs$quantile0.75_spillover20_outside<-as.numeric(store_outputs$quantile0.75_spillover20_outside)
  store_outputs$min_spillover20_outside<-as.numeric(store_outputs$min_spillover20_outside) 
  store_outputs$max_spillover20_outside<-as.numeric(store_outputs$max_spillover20_outside)
  store_outputs$quantile0.025_spillover20_outside<-as.numeric(store_outputs$quantile0.025_spillover20_outside)
  store_outputs$quantile0.975_spillover20_outside<-as.numeric(store_outputs$quantile0.975_spillover20_outside)
  store_outputs$check_spillover20_outside<-as.numeric(store_outputs$check_spillover20_outside)
 
} # end loop beta
} # end loop beta_interspecies
  }# end loop outbreak cull range    
} # end loop base cull range


write.csv(store_outputs,paste0("store_outputs_250runs.csv"))


#################################################################
##
## Plots



# Plot single outbreak when looping over n_run only
plot_outbreak(
  n_run = 1
)



# Plot individual run for particular scenario

plot_individual_run(
  base_cull_effort = "cull_rate",
  outbreak_cull_effort = "20*cull_rate",
  beta_interspecies = "1/30",
  beta = "1.6",
  n_run = 1 # run 14 of cull rate, 1/30,0.7 is no_outbreak
)

# Plot median run for particular scenario

plot_median_outbreak(
  base_cull_effort = "cull_rate",
  outbreak_cull_effort = "20*cull_rate",
  beta_interspecies = "14/30",
  beta = "1.6"
)




#################################################################
##
## Outputs when looping over n_run only


 # Probability of infection in pigs over outbreak
 
 infect_pig(n_run = n_run_pick)
         
 # Size of epidemic over outbreak (Dead wild boar)

 outputs <- as.matrix(size(n_run = n_run_pick))
 
 

