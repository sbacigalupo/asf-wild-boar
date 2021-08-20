# Version 2 - working version
# Frequency dependent transmission

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

# Set areas

# Inside forest
area_f <- 38 # pi*3.5km(r)squared

# Border area 
area_b <- 25 # pi*4.5^2 minus the above

# Outside forest
area_o <- 190 # pi*9^2 minus the above

# Setting up an array 

patchNames=c("forest","border","outside")
stateNames=c("S","E","I","D")
n.patch=length(patchNames)
n.state=length(stateNames)
c_trace_tab<- array(NA,dim=c(length(1:max_time),n.patch,n.state),dimnames=list(NULL,patchNames,stateNames))

# Setting up an array for force of infection and other outputs

outputNames=c("i_force","boar-pig_inf")
n.output=length(outputNames)
o_trace_tab<- array(NA,dim=c(length(1:max_time),n.patch,n.output),dimnames=list(NULL,patchNames,outputNames))

# Initial conditions
init_pop_s <- rep(0,n.patch)
init_pop_s[1] <- 1200 # Number of susceptible in forest
init_pop_s[2] <- 240 # Number of susceptible at border
init_pop_s[3] <- 60  # Number of susceptible outside forest

init_pop_e <- rep(0,n.patch)
init_pop_e[1] <- 0  # Number of pre-infectious

init_pop_i <- rep(0,n.patch)
init_pop_i[1] <- 1  # Number of infectious

init_pop_d <- rep(0,n.patch)
init_pop_d[1] <- 0  # Number of dead

init_foi <- rep(0,n.patch)
o_trace_tab[1,,] <- init_foi

c_trace_tab[1,,"S"] <- init_pop_s # Set number of susceptible in forest
c_trace_tab[1,,"E"] <- init_pop_e # Set number of pre-infectious in forest
c_trace_tab[1,,"I"] <- init_pop_i # Set number of infectious in forest
c_trace_tab[1,,"D"] <- init_pop_d # Set number of infectious in forest

c_trace_tab[1,"forest","S"] # Check numbers
c_trace_tab[1,"forest","E"] 
c_trace_tab[1,"forest","I"] 
c_trace_tab[1,"forest","D"] 

# Transmission rate  (i.e. R0*gamma)

beta <-0.75/1500        # density dependent would be contact rate*transmission rate

beta_wbp <- 1/30/1500 # Contact rate between wild boar and farm (1 visit per 30 days per 1500 wild boar)
                     # assumes 100% transmission at the moment

zeta <- 1/7 # Pre-infectious to infectious (latent period) 7 days

infection_death <- 1/10 # 5 day lifespan

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
death_probability <- (1/365 + cull_probability) # Assumes lifespan 10 years
carrying_capacity <- c(1e4,300,1e5) # How many boar sustainable in different areas
daily_sd_movement <- 0.4 # Daily measure of variation in movement
a <- 0.6 # constant for density dependent growth
growth_rate = birth_per_capita - death_probability





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
  # removed /rowSums(new_total) from line 120
  
  S_to_E <- rpois(3,lambda=beta*new_total[,"S"]*new_total[,"I"] ) # generate random infections 
  E_to_I <- rpois(3, lambda = zeta*new_total[,"E"])
  I_to_death <- rpois(3,lambda=infection_death*new_total[,"I"] ) # death rate
  
  new_total[,"S"] <- pmax(new_total[,"S"] - S_to_E,0)  # (+ check S >=0)
  new_total[,"E"] <- pmax(new_total[,"E"] + S_to_E - E_to_I,0)
  new_total[,"I"] <- pmax(new_total[,"I"] + E_to_I - I_to_death,0)
  new_total[,"D"] <- new_total[,"D"] + I_to_death
  
  # Births and deaths into different compartments
#  scaled_births_with_carrying_capacity <- birth_per_capita*pmax(0,(1-new_total[,"S"]/carrying_capacity) ) # Calculate birth rate, scaling to reduce as nearer carrying capacity
#  new_births_S<- rpois(n_population,lambda = new_total[,"S"]*scaled_births_with_carrying_capacity) # Add Poisson births
  
  # Add deaths
#  new_deaths_S <- rbinom(n_population,size=new_total[,"S"],prob=death_probability) # Add deaths
#  new_deaths_E <- rbinom(n_population,size=new_total[,"E"],prob=death_probability) # Add deaths
#  new_deaths_I <- rbinom(n_population,size=new_total[,"I"],prob=death_probability) # Add deaths
  
  
  # Density dependent coefficients
  
  density_dependent_birth_coeff <- (a*growth_rate*new_total/carrying_capacity)# *new_total  # if doesn't run, remove hash here and line below. remove *total from new_

  density_dependent_death_coeff <- ((1-a)*growth_rate*new_total/carrying_capacity)#*new_total 
  
  new_births_S <- round((birth_per_capita - density_dependent_birth_coeff[,"S"])*new_total[,"S"], digits = 0)
  
  new_deaths_S <- round((rbinom(n_population,size=new_total[,"S"],prob=death_probability + density_dependent_death_coeff[,"S"])), digits=0) # Add deaths
  new_deaths_E <- round((rbinom(n_population,size=new_total[,"E"],prob=death_probability + density_dependent_death_coeff[,"E"])), digits=0) # Add deaths
  new_deaths_I <- round((rbinom(n_population,size=new_total[,"I"],prob=death_probability + density_dependent_death_coeff[,"I"])), digits=0) # Add deaths

  # might need to add new growth rate for next iteration based on actual difference in rates
  
  # Store new populations
  c_trace_tab[tt,,"S"] <- new_total[,"S"] + new_births_S - new_deaths_S # Add births and subtract deaths from total
  c_trace_tab[tt,,"E"] <- new_total[,"E"]  - new_deaths_E # Subtract deaths from total
  c_trace_tab[tt,,"I"] <- new_total[,"I"]  - new_deaths_I # Subtract deaths from total
  c_trace_tab[tt,,"D"] <- new_total[,"D"]
  
  # Calculate force of infection
 
 # for(tt in 2:max_time){ # iterate over days
    outputs_tt <- o_trace_tab[tt-1,,] # Values at the start of the day
    new_outputs <- outputs_tt         # Matrix to store new values
  
    foi <- beta*new_total[,"I"]       # Force of infection between wild boar in each area
    wb_p <- beta_wbp*new_total[,"I"]  # Force of infection from wild boar on farm
    
    new_outputs[,"i_force"] <- foi
    new_outputs[,"boar-pig_inf"] <- wb_p
    
    o_trace_tab[tt,,"i_force"] <- new_outputs[,"i_force"]
    o_trace_tab[tt,,"boar-pig_inf"] <- new_outputs[,"boar-pig_inf"]
    
  }


# Plot outplots -----------------------------------------------------------



c_trace_tab[1:5,,"S"] # Check first few iterations

c_trace_tab[1:5,,"E"]

c_trace_tab[1:5,,"I"]

c_trace_tab[1:5,,"D"]

#c_trace_tab[1:30,,"FOI"]


# Plot outputs by n.patch
{
par(mfrow=c(3,2))

col_pick_states <- list("blue","orange","red","black")
col_pick_patch <- list("darkgreen","brown","orange")

{plot(c_trace_tab[,,],col="white", 
     main = "Wild boar numbers in the forest", xlab = "Days", ylab = "Number of individuals", 
     xlim=c(0,max_time),ylim=c(1,max(c_trace_tab[,"forest",])))


for(dd in 1:n.state){
  lines(c_trace_tab[,"forest",dd],col=col_pick_states[[dd]])
  }
}

{plot(c_trace_tab[,,],col="white", 
     main = "Wild boar numbers at the border", xlab = "Days", ylab = "Number of individuals", 
     xlim=c(0,max_time),ylim=c(1,max(c_trace_tab[,"border",])))
  
  for(dd in 1:n.state){
    lines(c_trace_tab[,"border",dd],col=col_pick_states[[dd]])
  }
}
  
  {plot(c_trace_tab[,,],col="white", 
        main = "Wild boar numbers outside the forest", xlab = "Days", ylab = "Number of individuals", 
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

  {plot(o_trace_tab[,,],col="white", 
       main = "Force of infection (beta*I) over time", xlab = "Days", ylab = "Force of infection", 
       xlim=c(0,max_time),
       ylim=c(0,max(o_trace_tab[,,"i_force"])))
       legend("topright", legend = c("inside","border","outside"), 
       col = c("darkgreen","brown","orange"), lty = 1, border = "black", cex=1)
  
  
  for(dd in 1:n.patch){
    lines(o_trace_tab[,dd,"i_force"],col=col_pick_patch[[dd]])
  }
  }
{plot(o_trace_tab[,,],col="white", 
      main = "Force of infection (beta*I) on one farm over time\nassuming contact with farm leads to transmission\n(same as number of new infections)", xlab = "Days", ylab = "Force of infection", 
      xlim=c(0,max_time),
      ylim=c(0,max(o_trace_tab[,,"boar-pig_inf"])))
  legend("topright", legend = c("inside","border","outside"), 
         col = c("darkgreen","brown","orange"), lty = 1, border = "black", cex=1)
  
  
  for(dd in 1:n.patch){
    lines(o_trace_tab[,dd,"boar-pig_inf"],col=col_pick_patch[[dd]])
  }
}
}


