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

# subset by min and max beta

range_beta_outputs <- subset(mod_outputs, betax1500 == 0.7|betax1500 == 2.2)

range_beta_outputs<-arrange(range_beta_outputs, cull_effort, betax1500, beta_interspeciesx1500)



# Collate medians and ranges for spillovers

mod_outputs$spillover_forest <- paste0(mod_outputs$median_spillover_forest,
                                       "(",mod_outputs$Lwr95CI_spillover_forest,",",mod_outputs$Upr95CI_spillover_forest,")",
                                       " Range: ",mod_outputs$min_spillover_forest,"-",mod_outputs$max_spillover_forest)

mod_outputs$spillover_border <- paste0(mod_outputs$median_spillover_border,
                                       "(",mod_outputs$Lwr95CI_spillover_border,",",mod_outputs$Upr95CI_spillover_border,")",
                                       " Range: ",mod_outputs$min_spillover_border,"-",mod_outputs$max_spillover_border)

mod_outputs$spillover_outside <- paste0(mod_outputs$median_spillover_outside,
                                        "(",mod_outputs$Lwr95CI_spillover_outside,",",mod_outputs$Upr95CI_spillover_outside,")",
                                        " Range: ",mod_outputs$min_spillover_outside,"-",mod_outputs$max_spillover_outside)


spillovers <- matrix(c(range_beta_outputs$spillover_forest,range_beta_outputs$spillover_border,range_beta_outputs$spillover_outside), ncol=6, byrow=TRUE)
colnames(spillovers) <- c("1 visit/30 days", "7 visits/30 days", "14 visits/30 days","1 visit/30 days", "7 visits/30 days", "14 visits/30 days")
rownames(spillovers) <- c("No Cull", "Current Cull", "150% Cull", "Double Cull", "No Cull", "Current Cull", "150% Cull", "Double Cull", "No Cull", "Current Cull", "150% Cull", "Double Cull")
spillovers_table <- as.table(spillovers)
spillovers_table

kbl(spillovers_table, caption = "Number of spillover events per outbreak", escape = F) %>%
  kable_styling(bootstrap_options = c("striped", "hover", font_size = 1,full_width = F, position = "centre"))%>%
  pack_rows(index = c("Forest" = 4, "Border" =4, "Outside" = 4))%>%
  add_header_above(c("Beta" = 1, "0.7" = 3, "2.2" = 3))

