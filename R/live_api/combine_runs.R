#################################
#### Combining yearwise data ####
#################################

rm(list = ls())
library(dplyr)
library(purrr)

#############
#### OPR ####
#############

# Get all .rda file paths
rda_files <- list.files("data/pridge_vs_opr/",
                        pattern = "^pct_improvement_\\d{4}\\.rda$",
                        full.names = TRUE)

# Load each file and combine the results
result <- map_dfr(rda_files, function(file) {
    load(file)  # This loads the 'result' dataframe
    result
})

save(result, file = "data/pridge_vs_opr/pct_imp_combined.rda")

#############
#### EPA ####
#############

# Get all .rda file paths
rda_files <- list.files("data/pridge_vs_epa/",
                        pattern = "^pct_improvement_\\d{4}\\.rda$",
                        full.names = TRUE)

# Load each file and combine the results
result <- map_dfr(rda_files, function(file) {
    load(file)  # This loads the 'result' dataframe
    result
})

save(result, file = "data/pridge_vs_epa/pct_imp_combined.rda")
