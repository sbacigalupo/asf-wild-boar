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
max_time <- 365 # in days

# Setting up an array

patchNames=c("forest","border","outside")
stateNames=c("S","I")
n.patch=length(patchNames)
n.state=length(stateNames)
c_trace_tab<- array(NA,dim=c(length(1:max_time),n.patch,n.state),dimnames=list(NULL,patchNames,stateNames))

# Initial conditions
init_pop_s <- rep(0,n.patch)
init_pop_s[1] <- 1000  # Number of susceptible

init_pop_i <- rep(0,n.patch)
init_pop_i[1] <- 10  # Number of infected

c_trace_tab[1,,"S"] <- init_pop_s # Set number of susceptible in forest
c_trace_tab[1,,"I"] <- init_pop_i # Set number of infected in forest

c_trace_tab[1,"forest","S"] # Check numbers
c_trace_tab[1,"forest","I"] 

# Transmission rate  (i.e. R0*gamma)
beta <-0.22
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
  for(kk in 1:2){
    
    # Choose compartment to move
    if(kk==1){pick_p <- "S"}
    if(kk==2){pick_p <- "I"}
    
    # Add randomness to probability of movement
    move_data_daily <- matrix(rlnorm(9,mean=log(move_data_daily_0),sd=daily_sd_movement),nrow=n.patch,byrow=F) # generate random daily values for movement
    move_data_daily <- t(t(move_data_daily)/colSums(move_data_daily)) # Make sure everything sums up to 1. eg sum(move_data_daily[1])

    # Calculate susceptible transitions
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
  
  S_to_I <- rpois(3,lambda=beta*new_total[,"S"]*new_total[,"I"]/rowSums(new_total) ) # generate random infections 
  I_to_death <- rpois(3,lambda=infection_death*new_total[,"I"] ) # generate random infections 
  
  new_total[,"S"] <- pmax(new_total[,"S"] - S_to_I,0)  # (+ check S >=0)
  new_total[,"I"] <- new_total[,"I"] + S_to_I - I_to_death
  
  # Births and deaths into different compartments
  scaled_births_with_carrying_capacity <- birth_per_capita*pmax(0,(1-new_total[,"S"]/carrying_capacity) ) # Calculate birth rate, scaling to reduce as nearer carrying capacity
  new_births_S<- rpois(n_population,lambda = new_total[,"S"]*scaled_births_with_carrying_capacity) # Add Poisson births
  
  new_deaths_S <- rbinom(n_population,size=new_total[,"S"],prob=death_probability) # Add deaths
  new_deaths_I <- rbinom(n_population,size=new_total[,"I"],prob=death_probability) # Add deaths

  # Store new populations
  c_trace_tab[tt,,"S"] <- new_total[,"S"] + new_births_S - new_deaths_S # Add births and subtract deaths from total
  c_trace_tab[tt,,"I"] <- new_total[,"I"]  - new_deaths_I # Subtract deaths from total
    
  }



# Plot outplots -----------------------------------------------------------



c_trace_tab[1:5,,"S"] # Check first few iterations

c_trace_tab[1:5,,"I"]

# Plot populations

# par(mar=c(3,3,1,1))

plot(c_trace_tab[,,],col="white", 
     main = "Movements of Susceptible and Infected Populations", xlab = "Days", ylab = "Number of individuals", 
     xlim=c(0,max_time),ylim=c(1,1e2))

col_pick_s <- list("green","royalblue","red")
col_pick_i <- list("olivedrab3","skyblue","indianred1")


for(ii in 1:n.patch){
  
  lines(c_trace_tab[,ii,"S"],col=col_pick_s[[ii]])
  lines(c_trace_tab[,ii,"I"],col=col_pick_i[[ii]])
  
}

legend(1,3000,legend = c("s_forest","s_border","s_outside","i_forest","i_border","i_outside"), 
col = c("green","royalblue","red","olivedrab3","skyblue","indianred1"), lty = 1, border = "black", cex=0.5)

