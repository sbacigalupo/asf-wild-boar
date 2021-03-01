library(tidyverse)

setwd("/Users/sonnybacigalupo/Documents/RVC_PhD/R_Code/movement_model/Sonny_movement_model")


#Set up numbers
n_population <- 3
max_time <- 365 # in days
pop_store <- matrix(NA,nrow = max_time, ncol = n_population)

# Initial conditions
init_pop <- rep(0,n_population)
init_pop[1] <- 1000
pop_store[1,] <- init_pop


# Read in movement data (per week)
move_data <- read_csv("move_matrix.csv",col_names = T) %>% as.matrix.data.frame()

# Convert data to daily rate
move_data_daily_0 <- move_data

# Make sure matrix is daily probability of movement to each other location by location
# i.e. columns sum up to 1
move_data_daily_0 <- move_data_daily_0/colSums(move_data_daily_0)
# NEED TO CHECK THIS

# Remove below for now, as already checked probabilities above
# for(ii in 1:n_population){
#   outside_movements <-  sum(move_data[,ii])- move_data[ii,ii] # movements outside
#   remaining_daily <- sum(move_data[,ii])- outside_movements/7
#   
#   move_data_daily[,ii] <- move_data[,ii]/7 # convert leaving to daily
#   move_data_daily[ii,ii] <- remaining_daily # remaining daily
#   
#   move_data_daily[,ii] <- move_data_daily[,ii]/sum(move_data_daily[,ii]) # normalise to probabilty
#   
# }

# Set up birth, death parameters per day
birth_per_capita <- 0.45*5/365 # Assumes 90% of females (50% of population) have 5 piglets per year
cull_probability <- 0.5/365
death_probability <- (0.1/365 + cull_probability) # Assumes lifespan 10 years
carrying_capacity <- c(1e4,100,1e5) # How many boar sustainable in different areas
daily_sd_movement <- 0.4 # Daily measure of variation in movement

# Loop over time, with population movements from matrix

for(tt in 2:max_time){ # iterate over days
  
  pop_time_tt <- pop_store[tt-1,] # Population at start of this day
  
  #expected_on_next_day <- move_data_daily %*% pop_store[tt-1,] # Calculate expected number in each place on next day
  
  # Add randomness to probability of movement
  move_data_daily <- matrix(rlnorm(9,mean=log(move_data_daily_0),sd=daily_sd_movement),nrow=n_population,byrow=F) # generate random daily values for movement
  move_data_daily <- t(t(move_data_daily)/colSums(move_data_daily)) # Make sure everything sums up to 1

  
  # Calculate transitions
  pop_f_to_f <- rbinom(1,pop_time_tt[1],move_data_daily[1,1]) # How many stay in forest
  pop_f_to_b <- rbinom(1,pop_time_tt[1]-pop_f_to_f,move_data_daily[2,1]/(move_data_daily[2,1]+move_data_daily[3,1])) # Of those that leave forest, how many go to border
  pop_f_to_o <- pop_time_tt[1] - pop_f_to_f - pop_f_to_b # The rest must go to outside
  
  pop_b_to_f <- rbinom(1,pop_time_tt[2],move_data_daily[1,2]) # How many go to in forest
  pop_b_to_b <- rbinom(1,pop_time_tt[2]-pop_b_to_f,move_data_daily[2,2]/(move_data_daily[2,2]+move_data_daily[3,2])) # Of rest, how many stay
  pop_b_to_o <- pop_time_tt[2] - pop_b_to_f - pop_b_to_b # The rest must stay outside
  
  pop_o_to_f <- rbinom(1,pop_time_tt[3],move_data_daily[1,3]) # How many go to forest
  pop_o_to_b <- rbinom(1,pop_time_tt[3]-pop_o_to_f,move_data_daily[2,3]/(move_data_daily[2,3]+move_data_daily[3,3])) # Of rest, how many go to border
  pop_o_to_o <- pop_time_tt[3] - pop_o_to_f - pop_o_to_b # The rest must go to outside
  
  # Tally up new populations after movement
  new_pop_f <- pop_f_to_f+pop_b_to_f+pop_o_to_f
  new_pop_b <- pop_f_to_b+pop_b_to_b+pop_o_to_b
  new_pop_o <- pop_f_to_o+pop_b_to_o+pop_o_to_o
  
  new_total <- c(new_pop_f,new_pop_b,new_pop_o) # Create vector with totals
  
  # Births and deaths
  scaled_births_with_carrying_capacity <- birth_per_capita*pmax(0,(1-new_total/carrying_capacity) ) # Calculate birth rate, scaling to reduce as nearer carrying capacity
  
  new_births <- rpois(n_population,lambda = new_total*scaled_births_with_carrying_capacity) # Add Poisson births
  new_deaths <- rbinom(n_population,size=new_total,prob=death_probability) # Add deaths
  
  new_total_with_birth_and_death <- c(new_pop_f,new_pop_b,new_pop_o) + new_births - new_deaths # Add births and subtract deaths from total
  
  # Store new populations
  pop_store[tt,] <- new_total_with_birth_and_death

}



# Plot populations

par(mar=c(3,3,1,1))

plot(pop_store[,1],col="white", xlim=c(0,max_time),ylim=c(1,5e3))
col_pick <- list("red","blue","green")

for(ii in 1:n_population){
  
  lines(pop_store[,ii],col=col_pick[[ii]])
  
}

# Setting up an array

patchNames=c("forest","outside")
n.patch=length(patchNames)
c_trace_tab<- array(NA,dim=c(length(max_time),n.patch,2),dimnames=list(NULL,patchNames,c("S","I")))
c_trace_tab[1,"forest","S"]

c_trace_tab




