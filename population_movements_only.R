# Version 2 - working version

library(tidyverse)

setwd("/Users/sonnybacigalupo/Documents/GitHub/asf-wild-boar")

# User define wd
if(Sys.info()["user"]=="adamkuchars" | Sys.info()["user"]=="akucharski" | Sys.info()["user"]=="adamkucharski") {
  setwd("~/Documents/GitHub/asf-wild-boar/")
}


# Set up model parameters ------------------------------------------------------------

#Set up numbers
n_population <- 3
max_time <- 1800 # in days

# Setting up an array

patchNames=c("forest","border","outside")
stateNames=c("S","E","I","D")
n.patch=length(patchNames)
n.state=length(stateNames)
c_trace_tab<- array(NA,dim=c(length(1:max_time),n.patch,n.state),dimnames=list(NULL,patchNames,stateNames))

# Initial conditions
init_pop_s <- rep(0,n.patch)
init_pop_s[1] <- 800 # Number of susceptible in forest
init_pop_s[2] <- 160 # Number of susceptible at border
init_pop_s[3] <- 40  # Number of susceptible outside forest

init_pop_e <- rep(0,n.patch)
init_pop_e[1] <- 0  # Number of pre-infectious in forest

init_pop_i <- rep(0,n.patch)
init_pop_i[1] <- 0  # Number of infectious in forest

init_pop_d <- rep(0,n.patch)
init_pop_d[1] <- 0  # Number of dead in forest

c_trace_tab[1,,"S"] <- init_pop_s # Set number of susceptible
c_trace_tab[1,,"E"] <- init_pop_e # Set number of pre-infectious 
c_trace_tab[1,,"I"] <- init_pop_i # Set number of infectious 
c_trace_tab[1,,"D"] <- init_pop_d # Set number of infectious 

c_trace_tab[1,"forest","S"] # Check numbers
c_trace_tab[1,"forest","E"] 
c_trace_tab[1,"forest","I"] 
c_trace_tab[1,"forest","D"] 

# Transmission rate  (i.e. R0*gamma)
beta <-0.00075
zeta <- 1/7 # Pre-infectious to infectious (latent period) 7 days
infection_death <- 1/5 # 5 day lifespan

# Read in movement data (per week)
move_data <- read_csv("move_matrix.csv",col_names = T) %>% as.matrix.data.frame()

# Convert data to daily rate
move_data_daily_0 <- move_data

# Make sure matrix is daily probability of movement to each other location by location
# i.e. columns sum up to 1
move_data_daily_0 <- move_data_daily_0/colSums(move_data_daily_0)

# Set up birth, death parameters per day
birth_per_capita <- 0.45*5/365 # Assumes 90% of females (50% of population) have 5 piglets per year
cull_probability <- 0.5/365
death_probability <- (0.1/365 + cull_probability) # Assumes lifespan 10 years
carrying_capacity <- c(1e4,100,1e5) # How many boar sustainable in different areas
daily_sd_movement <- 0.4 # Daily measure of variation in movement



# Run model ---------------------------------------------------------------

# Loop over time, with population movements from matrix

for(tt in 2:max_time){ # iterate over days
  
  pop_time_tt <- c_trace_tab[tt-1,,] # Population at start of this day
  new_total <- pop_time_tt # set up matrix to store new values
  
  # Loop over infection compartments
  for(kk in 1:3){
    
    # Choose compartment to move
    if(kk==1){pick_p <- "S"}
    if(kk==2){pick_p <- "E"}
    if(kk==3){pick_p <- "I"}
    
    # Add randomness to probability of movement
    move_data_daily <- matrix(rlnorm(9,mean=log(move_data_daily_0),sd=daily_sd_movement),nrow=n.patch,byrow=F) # generate random daily values for movement
    move_data_daily <- t(t(move_data_daily)/colSums(move_data_daily)) # Make sure everything sums up to 1. eg sum(move_data_daily[1])

    # Calculate movement transitions
    pop_f_to_f <- rbinom(1,  pop_time_tt[1,pick_p],move_data_daily[1,1]) # How many stay in forest
    pop_f_to_b <- rbinom(1,  pop_time_tt[1,pick_p] -pop_f_to_f,move_data_daily[2,1]/(move_data_daily[2,1]+move_data_daily[3,1])) # Of those that leave forest, how many go to border
    pop_f_to_o <- as.integer(pop_time_tt[1,pick_p] - pop_f_to_f - pop_f_to_b) # The rest must go to outside
    
    pop_b_to_f <- rbinom(1,  pop_time_tt[2,pick_p],move_data_daily[1,2]) # How many go to in forest
    pop_b_to_b <- rbinom(1,  pop_time_tt[2,pick_p]-pop_b_to_f,move_data_daily[2,2]/(move_data_daily[2,2]+move_data_daily[3,2])) # Of rest, how many stay
    pop_b_to_o <- as.integer(pop_time_tt[2,pick_p] - pop_b_to_f - pop_b_to_b) # The rest must stay outside
    
    pop_o_to_f <- rbinom(1,  pop_time_tt[3,pick_p],move_data_daily[1,3]) # How many go to forest
    pop_o_to_b <- rbinom(1,  pop_time_tt[3,pick_p]-pop_o_to_f,move_data_daily[2,3]/(move_data_daily[2,3]+move_data_daily[3,3])) # Of rest, how many go to border
    pop_o_to_o <- as.integer(pop_time_tt[3,pick_p] - pop_o_to_f - pop_o_to_b) # The rest must go to outside
  
    # Tally up new populations after movement
    new_pop_f <- pop_f_to_f+pop_b_to_f+pop_o_to_f
    new_pop_b <- pop_f_to_b+pop_b_to_b+pop_o_to_b
    new_pop_o <- pop_f_to_o+pop_b_to_o+pop_o_to_o
    
    new_total[,pick_p] <- c(new_pop_f,new_pop_b,new_pop_o)
    
  }
  
  # Add disease transitions
  # removed /rowSums(new_total) from line 120
  
  S_to_E <- rpois(3,lambda=beta*new_total[,"S"]*new_total[,"I"] ) # generate random infections 
  E_to_I <- rpois(3, lambda = zeta*new_total[,"E"])
  I_to_death <- rpois(3,lambda=infection_death*new_total[,"I"] ) # death rate 
  
  new_total[,"S"] <- pmax(new_total[,"S"] - S_to_E,0)  # (+ check S >=0)
  new_total[,"E"] <- pmax(new_total[,"E"] + S_to_E - E_to_I,0)
  new_total[,"I"] <- pmax(new_total[,"I"] + E_to_I - I_to_death,0)
  new_total[,"D"] <- new_total[,"D"] + I_to_death
  
  
  # Births and deaths into different compartments
  scaled_births_with_carrying_capacity <- birth_per_capita*pmax(0,(1-new_total[,"S"]/carrying_capacity) ) # Calculate birth rate, scaling to reduce as nearer carrying capacity
  new_births_S<- rpois(n_population,lambda = new_total[,"S"]*scaled_births_with_carrying_capacity) # Add Poisson births
  
  new_deaths_S <- rbinom(n_population,size=new_total[,"S"],prob=death_probability) # Add deaths
  new_deaths_E <- rbinom(n_population,size=new_total[,"E"],prob=death_probability) # Add deaths
  new_deaths_I <- rbinom(n_population,size=new_total[,"I"],prob=death_probability) # Add deaths

  # Store new populations
  c_trace_tab[tt,,"S"] <- new_total[,"S"] + new_births_S - new_deaths_S # Add births and subtract deaths from total
  c_trace_tab[tt,,"E"] <- new_total[,"E"] - new_deaths_E # Subtract deaths from total
  c_trace_tab[tt,,"I"] <- new_total[,"I"] - new_deaths_I # Subtract deaths from total
  c_trace_tab[tt,,"D"] <- new_total[,"D"]
  
    
  }



# Plot outplots -----------------------------------------------------------



c_trace_tab[1:5,,"S"] # Check first few iterations

c_trace_tab[1:5,,"E"]

c_trace_tab[1:5,,"I"]

c_trace_tab[1:30,,"D"]


# Plot outputs by n.patch
{
par(mfrow=c(2,2))

col_pick_states <- list("blue","orange","red","black")

{
plot(c_trace_tab[,,],col="white", 
     main = "Movements of Populations in the Forest", xlab = "Days", ylab = "Number of individuals", 
     xlim=c(0,max_time),ylim=c(1,max(c_trace_tab[,"forest",])))


for(dd in 1:n.state){
  lines(c_trace_tab[,"forest",dd],col=col_pick_states[[dd]])
}

{plot(c_trace_tab[,,],col="white", 
     main = "Movements of Populations at the Border", xlab = "Days", ylab = "Number of individuals", 
     xlim=c(0,max_time),ylim=c(1,max(c_trace_tab[,"border",])))
  
  for(dd in 1:n.state){
    lines(c_trace_tab[,"border",dd],col=col_pick_states[[dd]])
  }
}
  
  {plot(c_trace_tab[,,],col="white", 
        main = "Movements of Populations Outside", xlab = "Days", ylab = "Number of individuals", 
        xlim=c(0,max_time),ylim=c(1,max(c_trace_tab[,"outside",])))
    
    for(dd in 1:n.state){
      lines(c_trace_tab[,"outside",dd],col=col_pick_states[[dd]])
    }
    
    {plot(c_trace_tab[,,],col="white",axes = FALSE, xlab = "", ylab = "",
          main = "Legend", 
          xlim=c(0,max_time),ylim=c(1,max(c_trace_tab[,"outside",])))
      
      legend("center",legend = c("susceptible","pre-infectious","infectious","dead", beta, zeta, infection_death), 
             col = c("blue","orange","red","black","white","white","white"), lty = 1, border = "black", cex=1)
    }
  }
}
}


