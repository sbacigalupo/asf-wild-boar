outbreak <- function(max_time,beta,beta_wbp,zeta,infection_death){
  # Set up model parameters ------------------------------------------------------------
  
  #Set up numbers
  n_population <- 3
  
  # Setting up an array 
  
  
  patchNames=c("forest","border","outside")
  stateNames=c("S","E","I","D")
  n.patch=length(patchNames)
  n.state=length(stateNames)
  c_trace_tab<- array(NA,dim=c(length(1:max_time),n.patch,n.state),dimnames=list(NULL,patchNames,stateNames))
  
  # Setting up an array for force of infection and other outputs
  
  outputNames=c("i_force","prob_wb_wb", "boar-pig_inf" ,"prob_wb_p")
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
  
  # Read in movement data (per week)
  
  # Suppress parse message
  col_types <- cols(
    forest = col_double(),
    border = col_double(),
    outside = col_double()
  )
  
  move_data <- read_csv("move_matrix.csv",col_types = col_types, col_names = T) %>% as.matrix.data.frame()
  #move_data <- read.csv("move_matrix.csv",header = T)
  
  # Convert data to daily rate
  move_data_daily_0 <- move_data
  
  # Make sure matrix is daily probability of movement to each other location by location
  # i.e. columns sum up to 1
  move_data_daily_0 <- move_data_daily_0/colSums(move_data_daily_0)
  
  # Set up birth, death parameters per day
  births_per_capita <- 0.45*5/365 # Assumes 90% of females (50% of population) have 5 piglets per year
  cull_probability <- 0.5/365
  deaths_per_capita <- 1/3650   # Assumes lifespan 10 years
  carrying_capacity <- c(1200,240,60) # How many boar sustainable in different areas
  daily_sd_movement <- 0.4 # Daily measure of variation in movement
  a <- 0.5                 # constant for density dependent growth (between 0 and 1)
  growth_rate <- births_per_capita - deaths_per_capita
  
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
    new_total[,"D"] <- pmax(new_total[,"D"] + I_to_death,0)
    
    ################### SG EDITED FROM HERE
    # Births and deaths into different compartments
    
    # Density dependent birth and death rates
    birth_rate <- pmax(births_per_capita - a*growth_rate*colSums(new_total[,1:3])/carrying_capacity,0) #SB amended to ensure >0
    death_rate <- deaths_per_capita + (1-a)*growth_rate*colSums(new_total[,1:3])/carrying_capacity
    
    new_births_S <- rpois(n_population, lambda = birth_rate*colSums(new_total[,1:3])) #SB amended
    
    new_deaths_S <- rpois(n_population, lambda = death_rate*new_total[,"S"]) #SB amended
    new_deaths_E <- rpois(n_population, lambda = death_rate*new_total[,"E"]) #SB amended
    new_deaths_I <- rpois(n_population, lambda = death_rate*new_total[,"I"]) #SB amended
    
    ################### SG EDITED TO HERE
    # might need to add new growth rate for next iteration based on actual difference in rates
    
    # Store new populations
    c_trace_tab[tt,,"S"] <- new_total[,"S"] + new_births_S - new_deaths_S # Add births and subtract deaths from total
    c_trace_tab[tt,,"E"] <- new_total[,"E"]  - new_deaths_E # Subtract deaths from total
    c_trace_tab[tt,,"I"] <- new_total[,"I"]  - new_deaths_I # Subtract deaths from total
    c_trace_tab[tt,,"D"] <- new_total[,"D"]
    
    # Calculate force of infection
    
    outputs_tt <- o_trace_tab[tt-1,,] # Values at the start of the day
    new_outputs <- outputs_tt         # Matrix to store new values
    
    foi <- beta*new_total[,"I"]       # Force of infection between wild boar in each area
    
    prob_wb <- 1 - (exp(-foi*1))
    
    wb_p <- beta_wbp*new_total[,"I"]  # Force of infection from wild boar on farm
    
    prob_wbp <- (wb_p*1)
    
    new_outputs[,"i_force"] <- foi
    new_outputs[,"prob_wb_wb"] <- prob_wb
    new_outputs[,"boar-pig_inf"] <- wb_p
    new_outputs[,"prob_wb_p"] <- prob_wbp
    
    o_trace_tab[tt,,"i_force"] <- new_outputs[,"i_force"]
    o_trace_tab[tt,,"boar-pig_inf"] <- new_outputs[,"boar-pig_inf"]
    o_trace_tab[tt,,"prob_wb_wb"] <- new_outputs[,"prob_wb_wb"]
    o_trace_tab[tt,,"prob_wb_p"] <- new_outputs[,"prob_wb_p"]
  }
  
  trace_tabs <- list("populations" = c_trace_tab,"outputs" = o_trace_tab)
  return(trace_tabs)
}




###########################################################
##
## Plot individual run
##
###########################################################




  # Plot outputs by n.patch
  plot_outbreak <- function(n_run,scenario){
    
    pop_trace <- list_pop[[n_run]]
    output_trace <- list_output[[n_run]]
    
    max_time <- length(pop_trace[,1,"S"])
    n.state <- length(pop_trace[1,"forest",])
    n.patch <- length(pop_trace[1,,"S"])
    
    par(mfrow=c(3,2), oma = c(0, 0, 1, 0))
    
  
    col_pick_states <- list("blue","orange","red","black")
    col_pick_patch <- list("darkgreen","brown","orange")
    
    {plot(pop_trace[,,],col="white", 
          main = "Wild boar numbers in the forest", xlab = "Days", ylab = "Number of individuals", 
          xlim=c(0,max_time),ylim=c(1,max(pop_trace[,"forest",])))
      
      
      for(dd in 1:n.state){
        lines(pop_trace[,"forest",dd],col=col_pick_states[[dd]])
      }
    }
    
    {plot(pop_trace[,,],col="white", 
          main = "Wild boar numbers at the border", xlab = "Days", ylab = "Number of individuals", 
          xlim=c(0,max_time),ylim=c(1,max(pop_trace[,"border",])))
      
      for(dd in 1:n.state){
        lines(pop_trace[,"border",dd],col=col_pick_states[[dd]])
      }
    }
    
    {plot(pop_trace[,,],col="white", 
          main = "Wild boar numbers outside the forest", xlab = "Days", ylab = "Number of individuals", 
          xlim=c(0,max_time),ylim=c(1,max(pop_trace[,"outside",])))
      
      for(dd in 1:n.state){
        lines(pop_trace[,"outside",dd],col=col_pick_states[[dd]])
      }
      
      {plot(pop_trace[,,],col="white",axes = FALSE, xlab = "", ylab = "",
            main = "Legend", 
            xlim=c(0,max_time),ylim=c(1,max(pop_trace[,"outside",])))
        
        legend("center",legend = c("susceptible","pre-infectious","infectious","dead"), 
               col = c("blue","orange","red","black","white","white","white"), lty = 1, border = "black", cex=1)
      }
    }
    
    {plot(output_trace[,,],col="white", 
          main = "Probability of infection between \nwild boar over time", xlab = "Days", ylab = "Likelihood of infection", 
          xlim=c(0,max_time),
          ylim=c(0,max(output_trace[,,"prob_wb_wb"])))
      legend("topright", legend = c("inside","border","outside"), 
             col = c("darkgreen","brown","orange"), lty = 1, border = "black", cex=1)
      
      
      for(dd in 1:n.patch){
        lines(output_trace[,dd,"prob_wb_wb"],col=col_pick_patch[[dd]])
      }
    }
    {plot(output_trace[,,],col="white", 
          main = "Probability of infection between \nwild boar and pigs over time", xlab = "Days", ylab = "Force of infection", 
          xlim=c(0,max_time),
          ylim=c(0,max(output_trace[,,"prob_wb_p"])))
      legend("topright", legend = c("inside","border","outside"), 
             col = c("darkgreen","brown","orange"), lty = 1, border = "black", cex=1)
      
      
      for(dd in 1:n.patch){
        lines(output_trace[,dd,"prob_wb_p"],col=col_pick_patch[[dd]])
      }
      mtext(paste0("Run ",n_run,"  (B_wb=",scenario[["beta"]],";  B_p=",signif(scenario[["beta_wbp"]], digits = 3),";  z(tEI)=",signif(scenario[["zeta"]], digits = 3),";  y(tID)=",signif(scenario[["infection_death"]], digits = 3),")"), line = -1, cex = 0.85, outer = TRUE)
    }
  }
  
  
  ###########################################################
  ##
  ## Return total probability of transmission to pigs "prob_wb_p"
  ##
  ###########################################################
  
infect_pig <- function(n_runs){
  store_table <- NULL
  store_r <- NULL
  for(jj in 1:n_runs){
  output_trace <- list_output[[jj]]
  infectpig_F <- sum(output_trace[,"forest","prob_wb_p"])
  infectpig_B <- sum(output_trace[,"border","prob_wb_p"])
  infectpig_O <- sum(output_trace[,"outside","prob_wb_p"])
  store_r<-rbind(store_r,c(jj,infectpig_F,infectpig_B,infectpig_O))
  }
  store_r <- as_tibble(store_r)
  names(store_r) <- c("run","infectpig_F","infectpig_B","infectpig_O")
  
#  write_csv(store_r,paste0("rstore.csv"))
  
  store_table <- rbind(store_table,c(median(store_r$infectpig_F),quantile(store_r$infectpig_F,c(0.025,0.975)),
                                     median(store_r$infectpig_B),quantile(store_r$infectpig_B,c(0.025,0.975)),
                                     median(store_r$infectpig_O),quantile(store_r$infectpig_O,c(0.025,0.975))))
  
  # Convert
  store_table <- as_tibble(store_table)
  names(store_table) <- c("infectpig_F_med","infectpig_F_95_1","infectpig_F_95_2",
                                   "infectpig_B_med","infectpig_B_95_1","infectpig_B_95_2",
                                   "infectpig_O_med","infectpig_O_95_1","infectpig_O_95_2")
  store_table$infectpig_F_med <- as.numeric(store_table$infectpig_F_med)
  store_table$infectpig_F_95_1 <- as.numeric(store_table$infectpig_F_95_1)
  store_table$infectpig_F_95_2 <- as.numeric(store_table$infectpig_F_95_2)
  store_table$infectpig_B_med <- as.numeric(store_table$infectpig_B_med)
  store_table$infectpig_B_95_1 <- as.numeric(store_table$infectpig_B_95_1)
  store_table$infectpig_B_95_2 <- as.numeric(store_table$infectpig_B_95_2)
  store_table$infectpig_O_med <- as.numeric(store_table$infectpig_O_med)
  store_table$infectpig_O_95_1 <- as.numeric(store_table$infectpig_O_95_1)
  store_table$infectpig_O_95_2 <- as.numeric(store_table$infectpig_O_95_2)
  
  # write_csv(store_table,paste0("testresults1.csv"))
}

###########################################################
##
## Return number of wild boar dead (outbreak size)
##
###########################################################

size <- function(n_runs){
  store_table <- NULL
  store_r <- NULL
  for(kk in 1:n_runs){
    pop_trace <- list_pop[[kk]]
    outbreak_size <- sum(max(pop_trace[,"forest","D"]),
                         max(pop_trace[,"border","D"]),
                         max(pop_trace[,"outside","D"]))
    
    store_r<-rbind(store_r,c(kk,outbreak_size))
  }
  store_r <- as_tibble(store_r)
  names(store_r) <- c("run","outbreak_size")
  
    write_csv(store_r,paste0("rstore.csv"))
  
  store_table <- rbind(store_table,c(median(store_r$outbreak_size),quantile(store_r$outbreak_size,c(0.025,0.975))))
  
  # Convert
  store_table <- as_tibble(store_table)
  names(store_table) <- c("outbreak_size_med","outbreak_size_95_1","outbreak_size_95_2")
  store_table$outbreak_size_med <- as.numeric(store_table$outbreak_size_med)
  store_table$outbreak_size_95_1 <- as.numeric(store_table$outbreak_size_95_1)
  store_table$outbreak_size_95_2 <- as.numeric(store_table$outbreak_size_95_2)

  
   write_csv(store_table,paste0("testresults2.csv"))
}
  
    
  