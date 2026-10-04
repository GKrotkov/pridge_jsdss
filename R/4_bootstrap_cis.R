# Computing bootstrap CIs for mean pct improvement on a per-year basis
rm(list = ls())
library(tidyverse)
library(gt)
library(glue)
library(parallel)

B <- 1000000
n_cores <- max(1, floor(detectCores() / 2))

######################
#### OPR Baseline ####
######################

load("data/pridge_vs_opr/pct_imp_2016_to_2026.rda")
years <- unique(result$year)

boot_means <- replicate(B, {
    mean(sample(result$pct_imp, replace = TRUE), na.rm = TRUE)
})
overall_opr_bounds <- round(quantile(boot_means, c(0.025, 0.5, 0.975)), 2)

start <- Sys.time()

cl <- makeCluster(n_cores)
clusterSetRNGStream(cl, iseed = 449)
clusterExport(cl, c("result", "B"))

opr_results <- parLapply(cl, years, function(yr) {
    year_df <- result[result$year == yr, ]
    boot_means <- replicate(B, {
        mean(sample(year_df$pct_imp, replace = TRUE), na.rm = TRUE)
    })
    quantile(boot_means, c(0.025, 0.5, 0.975))
})

stopCluster(cl)

end <- Sys.time()
execution_time_opr <- end - start

######################
#### EPA Baseline ####
######################

load("data/pridge_vs_epa/pct_imp_2016_to_2026.rda")

boot_means <- replicate(B, {
    mean(sample(result$pct_imp, replace = TRUE), na.rm = TRUE)
})
overall_epa_bounds <- round(quantile(boot_means, c(0.025, 0.5, 0.975)), 2)

years <- unique(result$year)

start <- Sys.time()

cl <- makeCluster(n_cores)
clusterSetRNGStream(cl, iseed = 449)
clusterExport(cl, c("result", "B"))

epa_results <- parLapply(cl, years, function(yr) {
    year_df <- result[result$year == yr, ]
    boot_means <- replicate(B, {
        mean(sample(year_df$pct_imp, replace = TRUE), na.rm = TRUE)
    })
    quantile(boot_means, c(0.025, 0.5, 0.975))
})

stopCluster(cl)

end <- Sys.time()
execution_time_epa <- end - start

####################
#### Tabulation ####
####################

n_events <- result |>
    group_by(year) |>
    summarize(count = n()) |>
    pull(count)

# attach the overall bounds to the front of the results
opr_lower <- c(overall_opr_bounds[1], sapply(opr_results, `[`, 1))
opr_median <- c(overall_opr_bounds[2], sapply(opr_results, `[`, 2))
opr_upper <- c(overall_opr_bounds[3], sapply(opr_results, `[`, 3))

epa_lower <- c(overall_epa_bounds[1], sapply(epa_results, `[`, 1))
epa_median <- c(overall_epa_bounds[2], sapply(epa_results, `[`, 2))
epa_upper <- c(overall_epa_bounds[3], sapply(epa_results, `[`, 3))

n_events <- c(nrow(result), n_events)

year <- c("Overall", unique(result$year)) # implicitly coerces to chr

bounds_opr <- tibble(year, n_events, lower = opr_lower,
                     median = opr_median, upper = opr_upper)

bounds_epa <- tibble(year, n_events, lower = epa_lower,
                     median = epa_median, upper = epa_upper)

bounds_tbl <- tibble(
    `Year` = year,
    `# Events` = n_events,
    `OPR Baseline` = glue("({round(opr_lower, 1)}%, {round(opr_median, 1)}%, {round(opr_upper, 1)}%)"),
    `EPA Baseline` = glue("({round(epa_lower, 1)}%, {round(epa_median, 1)}%, {round(epa_upper, 1)}%)")
)

save(bounds_opr, bounds_epa, bounds_tbl,
     file = "data/bootstrap_results.rda")
