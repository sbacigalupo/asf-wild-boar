library(tidyverse)
library(kableExtra)

setwd("/Users/sonnybacigalupo/Documents/GitHub/asf-wild-boar")

# User define wd
if(Sys.info()["user"]=="adamkuchars" | Sys.info()["user"]=="akucharski" | Sys.info()["user"]=="adamkucharski") {
  setwd("~/Documents/GitHub/asf-wild-boar/")
}


# Load 'store_outputs' into Global Environment from ASF_model


#################################################################
##
## Prepare tibble for manipulation to visualise data and create tables
##
#################################################################

mod_outputs<-store_outputs # new tibble

is.num <- sapply(mod_outputs, is.numeric) # Round numeric columns to 3 decimal places
mod_outputs[is.num] <- lapply(mod_outputs[is.num], round, 3)



#################################################################
##
## Spillovers over course of outbreak
##
#################################################################



# Collate medians and ranges for spillovers

mod_outputs$spillover_forest <- paste0(mod_outputs$median_spillover_forest,
                                       "(",mod_outputs$quantile0.25_spillover_forest,",",mod_outputs$quantile0.75_spillover_forest,")",
                                       " Range: ",mod_outputs$min_spillover_forest,"-",mod_outputs$max_spillover_forest)

mod_outputs$spillover_border <- paste0(mod_outputs$median_spillover_border,
                                       "(",mod_outputs$quantile0.25_spillover_border,",",mod_outputs$quantile0.75_spillover_border,")",
                                       " Range: ",mod_outputs$min_spillover_border,"-",mod_outputs$max_spillover_border)

mod_outputs$spillover_outside <- paste0(mod_outputs$median_spillover_outside,
                                        "(",mod_outputs$quantile0.25_spillover_outside,",",mod_outputs$quantile0.75_spillover_outside,")",
                                        " Range: ",mod_outputs$min_spillover_outside,"-",mod_outputs$max_spillover_outside)

# subset by min and max beta

range_beta_outputs <- subset(mod_outputs, betax1000 == 0.7|betax1000 == 2.2)

range_beta_outputs<-arrange(range_beta_outputs, cull_effort, betax1000, beta_interspeciesx1000)

spillovers <- matrix(c(range_beta_outputs$spillover_forest,range_beta_outputs$spillover_border,range_beta_outputs$spillover_outside), ncol=6, byrow=TRUE)
colnames(spillovers) <- c("1 visit/30 days", "7 visits/30 days", "14 visits/30 days","1 visit/30 days", "7 visits/30 days", "14 visits/30 days")
rownames(spillovers) <- c("No Cull", "Current Cull", "150% Cull", "Double Cull", "No Cull", "Current Cull", "150% Cull", "Double Cull", "No Cull", "Current Cull", "150% Cull", "Double Cull")
spillovers_table <- as.table(spillovers)

kbl(spillovers_table, caption = "Expected number of spillover events during outbreak", escape = F) %>%
  kable_styling(bootstrap_options = c("striped", "hover", font_size = 1,full_width = F, position = "centre"))%>%
  pack_rows(index = c("Forest" = 4, "Border" =4, "Outside" = 4))%>%
  add_header_above(c("Beta" = 1, "0.7" = 3, "2.2" = 3))


#################################################################
##
## Spillovers before disease detected when x% of population dead (x=10 or 20)
##
#################################################################


#---------------------------------------------------------------------
# Spillovers before disease detected when 10% of population dead

# Collate medians and ranges for spillovers

mod_outputs$spillover10_forest <- paste0(mod_outputs$median_spillover10_forest,
                                       "(",mod_outputs$quantile0.25_spillover10_forest,",",mod_outputs$quantile0.75_spillover10_forest,")",
                                       " Range: ",mod_outputs$min_spillover10_forest,"-",mod_outputs$max_spillover10_forest)

mod_outputs$spillover10_border <- paste0(mod_outputs$median_spillover10_border,
                                       "(",mod_outputs$quantile0.25_spillover10_border,",",mod_outputs$quantile0.75_spillover10_border,")",
                                       " Range: ",mod_outputs$min_spillover10_border,"-",mod_outputs$max_spillover10_border)

mod_outputs$spillover10_outside <- paste0(mod_outputs$median_spillover10_outside,
                                        "(",mod_outputs$quantile0.25_spillover10_outside,",",mod_outputs$quantile0.75_spillover10_outside,")",
                                        " Range: ",mod_outputs$min_spillover10_outside,"-",mod_outputs$max_spillover10_outside)

# subset by min and max beta

range_beta_outputs <- subset(mod_outputs, betax1000 == 0.7|betax1000 == 2.2)

range_beta_outputs<-arrange(range_beta_outputs, cull_effort, betax1000, beta_interspeciesx1000)


spillovers10 <- matrix(c(range_beta_outputs$spillover10_forest,range_beta_outputs$spillover10_border,range_beta_outputs$spillover10_outside), ncol=6, byrow=TRUE)
colnames(spillovers10) <- c("1 visit/30 days", "7 visits/30 days", "14 visits/30 days","1 visit/30 days", "7 visits/30 days", "14 visits/30 days")
rownames(spillovers10) <- c("No Cull", "Current Cull", "150% Cull", "Double Cull", "No Cull", "Current Cull", "150% Cull", "Double Cull", "No Cull", "Current Cull", "150% Cull", "Double Cull")
spillovers10_table <- as.table(spillovers10)

kbl(spillovers10_table, caption = "Expected number of spillover events before disease detected in wild boar (10% population death)", escape = F) %>%
  kable_styling(bootstrap_options = c("striped", "hover", font_size = 1,full_width = F, position = "centre"))%>%
  pack_rows(index = c("Forest" = 4, "Border" =4, "Outside" = 4))%>%
  add_header_above(c("Beta" = 1, "0.7" = 3, "2.2" = 3))


#---------------------------------------------------------------------
# Spillovers before disease detected when 20% of population dead

# Collate medians and ranges for spillovers

mod_outputs$spillover20_forest <- paste0(mod_outputs$median_spillover20_forest,
                                         "(",mod_outputs$quantile0.25_spillover20_forest,",",mod_outputs$quantile0.75_spillover20_forest,")",
                                         " Range: ",mod_outputs$min_spillover20_forest,"-",mod_outputs$max_spillover20_forest)

mod_outputs$spillover20_border <- paste0(mod_outputs$median_spillover20_border,
                                         "(",mod_outputs$quantile0.25_spillover20_border,",",mod_outputs$quantile0.75_spillover20_border,")",
                                         " Range: ",mod_outputs$min_spillover20_border,"-",mod_outputs$max_spillover20_border)

mod_outputs$spillover20_outside <- paste0(mod_outputs$median_spillover20_outside,
                                          "(",mod_outputs$quantile0.25_spillover20_outside,",",mod_outputs$quantile0.75_spillover20_outside,")",
                                          " Range: ",mod_outputs$min_spillover20_outside,"-",mod_outputs$max_spillover20_outside)

# subset by min and max beta

range_beta_outputs <- subset(mod_outputs, betax1000 == 0.7|betax1000 == 2.2)

range_beta_outputs<-arrange(range_beta_outputs, cull_effort, betax1000, beta_interspeciesx1000)


spillovers20 <- matrix(c(range_beta_outputs$spillover20_forest,range_beta_outputs$spillover20_border,range_beta_outputs$spillover20_outside), ncol=6, byrow=TRUE)
colnames(spillovers20) <- c("1 visit/30 days", "7 visits/30 days", "14 visits/30 days","1 visit/30 days", "7 visits/30 days", "14 visits/30 days")
rownames(spillovers20) <- c("No Cull", "Current Cull", "150% Cull", "Double Cull", "No Cull", "Current Cull", "150% Cull", "Double Cull", "No Cull", "Current Cull", "150% Cull", "Double Cull")
spillovers20_table <- as.table(spillovers20)

kbl(spillovers20_table, caption = "Expected number of spillover events before disease detected in wild boar (20% population death)", escape = F) %>%
  kable_styling(bootstrap_options = c("striped", "hover", font_size = 1,full_width = F, position = "centre"))%>%
  pack_rows(index = c("Forest" = 4, "Border" =4, "Outside" = 4))%>%
  add_header_above(c("Beta" = 1, "0.7" = 3, "2.2" = 3))



#---------------------------------------------------------------------
# Numbers and percentages of unsuccessful outbreaks

# Collate medians and ranges for spillovers

mod_outputs$avoided_outbreaks <- paste0(mod_outputs$number_no_outbreak ,
                                         "(",mod_outputs$percent_no_outbreak,")")

# subset by one interspecies beta value

interspecies1 <- subset(mod_outputs, beta_interspeciesx1000 == 0.033)

interspecies1<-arrange(interspecies1, cull_effort, betax1000)


avoided_outbreaks <- matrix(c(interspecies1$avoided_outbreaks), ncol=6, byrow=TRUE)
colnames(avoided_outbreaks) <- c("0.7", "1.0", "1.3","1.6", "1.9", "2.2")
rownames(avoided_outbreaks) <- c("No Cull", "Current Cull", "150% Cull", "Double Cull")
avoided_outbreaks_table <- as.table(avoided_outbreaks)

kbl(avoided_outbreaks, caption = "Numbers and percentages of outbreaks avoided out of 250 simulated outbreaks", escape = F) %>%
  kable_styling(bootstrap_options = c("striped", "hover", font_size = 1,full_width = F, position = "centre"))#%>%
 # pack_rows(index = c("Forest" = 4, "Border" =4, "Outside" = 4))%>%
  #add_header_above(c("Beta" = 1, "0.7" = 3, "2.2" = 3))

#---------------------------------------------------------------------
# Numbers and percentages of unsuccessful outbreaks

# Collate medians and ranges for spillovers

mod_outputs$outbreak_size <- paste0(mod_outputs$median_outbreak_size,
                                         "(",mod_outputs$quantile0.25_outbreak_size,",",mod_outputs$quantile0.75_spillover20_forest,")",
                                         " Range: ",mod_outputs$min_outbreak_size,"-",mod_outputs$max_outbreak_size)

# subset by min and max beta

range_beta_outputs <- subset(mod_outputs, betax1000 == 0.7|betax1000 == 2.2)

range_beta_outputs<-arrange(range_beta_outputs, cull_effort, betax1000, beta_interspeciesx1000)



outbreak_size <- matrix(c(range_beta_outputs$outbreak_size), ncol=6, byrow=TRUE)
colnames(outbreak_size) <- c("0.7", "1.0", "1.3","1.6", "1.9", "2.2")
rownames(outbreak_size) <- c("No Cull", "Current Cull", "150% Cull", "Double Cull")
outbreak_size_table <- as.table(outbreak_size)

kbl(outbreak_size, caption = "Outbreak size with 25th and 75th quantiles and ranges", escape = F) %>%
  kable_styling(bootstrap_options = c("striped", "hover", font_size = 1,full_width = F, position = "centre"))#%>%
# pack_rows(index = c("Forest" = 4, "Border" =4, "Outside" = 4))%>%
#add_header_above(c("Beta" = 1, "0.7" = 3, "2.2" = 3))