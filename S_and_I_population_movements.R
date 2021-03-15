library(tidyverse)

setwd("/Users/sonnybacigalupo/Documents/GitHub/asf-wild-boar")

#Set up numbers
n_population <- 6
max_time <- 365 # in days

# Setting up an array

patchNames=c("forest","border","outside")
stateNames=c("S","I")
n.patch=length(patchNames)
n.state=length(stateNames)
c_trace_tab<- array(NA,dim=c(length(1:max_time),n.patch,n.state),dimnames=list(NULL,patchNames,stateNames))


# Initial conditions
init_pop_s <- rep(0,n.patch)
init_pop_s[1] <- 950  # Number of susceptible

init_pop_i <- rep(0,n.patch)
init_pop_i[1] <- 50  # Number of infected

c_trace_tab[1,,"S"] <- init_pop_s # Set number of susceptible in forest
c_trace_tab[1,,"I"] <- init_pop_i # Set number of infected in forest


c_trace_tab[1,"forest","S"] # Check numbers
c_trace_tab[1,"forest","I"] 

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


# Loop over time, with population movements from matrix

for(tt in 2:max_time){ # iterate over days
  
  pop_time_tt <- c_trace_tab[tt-1,,] # Population at start of this day
  
  # Add randomness to probability of movement
  move_data_daily <- matrix(rlnorm(9,mean=log(move_data_daily_0),sd=daily_sd_movement),nrow=n.patch,byrow=F) # generate random daily values for movement
  move_data_daily <- t(t(move_data_daily)/colSums(move_data_daily)) # Make sure everything sums up to 1. eg sum(move_data_daily[1])
  
  # Calculate susceptible transitions
  popS_f_to_f <- rbinom(1,pop_time_tt[1],move_data_daily[1,1]) # How many stay in forest
  popS_f_to_b <- rbinom(1,pop_time_tt[1]-popS_f_to_f,move_data_daily[2,1]/(move_data_daily[2,1]+move_data_daily[3,1])) # Of those that leave forest, how many go to border
  popS_f_to_o <- as.integer(pop_time_tt[1] - popS_f_to_f - popS_f_to_b) # The rest must go to outside
  
  popS_b_to_f <- rbinom(1,pop_time_tt[2],move_data_daily[1,2]) # How many go to in forest
  popS_b_to_b <- rbinom(1,pop_time_tt[2]-popS_b_to_f,move_data_daily[2,2]/(move_data_daily[2,2]+move_data_daily[3,2])) # Of rest, how many stay
  popS_b_to_o <- as.integer(pop_time_tt[2] - popS_b_to_f - popS_b_to_b) # The rest must stay outside
  
  popS_o_to_f <- rbinom(1,pop_time_tt[3],move_data_daily[1,3]) # How many go to forest
  popS_o_to_b <- rbinom(1,pop_time_tt[3]-popS_o_to_f,move_data_daily[2,3]/(move_data_daily[2,3]+move_data_daily[3,3])) # Of rest, how many go to border
  popS_o_to_o <- as.integer(pop_time_tt[3] - popS_o_to_f - popS_o_to_b) # The rest must go to outside
  
  # Calculate infected transitions
  popI_f_to_f <- rbinom(1,pop_time_tt[4],move_data_daily[1,1]) # How many stay in forest
  popI_f_to_b <- rbinom(1,pop_time_tt[4]-popI_f_to_f,move_data_daily[2,1]/(move_data_daily[2,1]+move_data_daily[3,1])) # Of those that leave forest, how many go to border
  popI_f_to_o <- as.integer(pop_time_tt[4] - popI_f_to_f - popI_f_to_b) # The rest must go to outside
  
  popI_b_to_f <- rbinom(1,pop_time_tt[5],move_data_daily[1,2]) # How many go to in forest
  popI_b_to_b <- rbinom(1,pop_time_tt[5]-popI_b_to_f,move_data_daily[2,2]/(move_data_daily[2,2]+move_data_daily[3,2])) # Of rest, how many stay
  popI_b_to_o <- as.integer(pop_time_tt[5] - popI_b_to_f - popI_b_to_b) # The rest must stay outside
  
  popI_o_to_f <- rbinom(1,pop_time_tt[6],move_data_daily[1,3]) # How many go to forest
  popI_o_to_b <- rbinom(1,pop_time_tt[6]-popI_o_to_f,move_data_daily[2,3]/(move_data_daily[2,3]+move_data_daily[3,3])) # Of rest, how many go to border
  popI_o_to_o <- as.integer(pop_time_tt[6] - popI_o_to_f - popI_o_to_b) # The rest must go to outside
  
  # Tally up new populations after movement
  new_popS_f <- popS_f_to_f+popS_b_to_f+popS_o_to_f
  new_popS_b <- popS_f_to_b+popS_b_to_b+popS_o_to_b
  new_popS_o <- popS_f_to_o+popS_b_to_o+popS_o_to_o
  
  new_popI_f <- popI_f_to_f+popI_b_to_f+popI_o_to_f
  new_popI_b <- popI_f_to_b+popI_b_to_b+popI_o_to_b
  new_popI_o <- popI_f_to_o+popI_b_to_o+popI_o_to_o
  
  new_total <- c(new_popS_f,new_popS_b,new_popS_o,new_popI_f,new_popI_b,new_popI_o) # Create vector with totals
  
  # Births and deaths
  scaled_births_with_carrying_capacity <- birth_per_capita*pmax(0,(1-new_total/carrying_capacity) ) # Calculate birth rate, scaling to reduce as nearer carrying capacity
  
  new_births <- rpois(n_population,lambda = new_total*scaled_births_with_carrying_capacity) # Add Poisson births
  new_deaths <- rbinom(n_population,size=new_total,prob=death_probability) # Add deaths
  
  new_total_with_birth_and_death <- c(new_popS_f,new_popS_b,new_popS_o,new_popI_f,new_popI_b,new_popI_o) + new_births - new_deaths # Add births and subtract deaths from total
  
  # Store new populations

    c_trace_tab[tt,"forest","S"] <- new_total_with_birth_and_death[1]
    c_trace_tab[tt,"border","S"] <- new_total_with_birth_and_death[2]
    c_trace_tab[tt,"outside","S"] <- new_total_with_birth_and_death[3]
    c_trace_tab[tt,"forest","I"] <- new_total_with_birth_and_death[4]
    c_trace_tab[tt,"border","I"] <- new_total_with_birth_and_death[5]
    c_trace_tab[tt,"outside","I"] <- new_total_with_birth_and_death[6]
  }


c_trace_tab[1:5,,"S"] # Check first few iterations

c_trace_tab[1:5,,"I"]

# Plot populations

# par(mar=c(3,3,1,1))

plot(c_trace_tab[,,],col="white", 
     main = "Movements of Susceptible and Infected Populations", xlab = "Days", ylab = "Number of individuals", 
     xlim=c(0,max_time),ylim=c(1,3e3))

col_pick_s <- list("green","royalblue","red")
col_pick_i <- list("olivedrab3","skyblue","indianred1")


for(ii in 1:n.patch){
  
  lines(c_trace_tab[,ii,"S"],col=col_pick_s[[ii]])
  lines(c_trace_tab[,ii,"I"],col=col_pick_i[[ii]])
  
}

legend(1,3000,legend = c("s_forest","s_border","s_outside","i_forest","i_border","i_outside"), 
col = c("green","royalblue","red","olivedrab3","skyblue","indianred1"), lty = 1, border = "black", cex=0.5)

