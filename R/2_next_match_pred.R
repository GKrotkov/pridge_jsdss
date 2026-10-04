# Compare pRidge to EPA and OPR using locally cached data
# Uses the next-match-prediction framework to make a fair comparison between
# Elo ratings and regression methods like OPR/pRidge
# Uses data from data/raw/ instead of Statbotics API calls
# Last update: 5/5/2026

rm(list = ls())
library(scoutR)
library(tidyverse)
library(parallel)
library(doParallel)

# Configuration
data_dir <- "data/raw/"

# Given a vector of coefs and a response, compute the prediction error
coef_error <- function(coefs, design_row, response){
    return(response - drop(sum(design_row * coefs)))
}

# Given a vector of coefs and a blue/red design row pair, return 1 if the
# predicted winner matches the actual winner, 0 otherwise
coef_acc <- function(coefs, blue_design_row, red_design_row,
                     blue_response, red_response){
    blue_pred <- drop(sum(blue_design_row * coefs))
    red_pred  <- drop(sum(red_design_row  * coefs))
    return(as.integer((blue_pred > red_pred) == (blue_response > red_response)))
}

# Retrieve EPA progression for all teams across the event
# Now uses cached team_matches data instead of API call
get_epa_progression <- function(team_matches_data, event_key){
    result <- map(team_matches_data, ~{
        data.frame(
            team = .x$team,
            match = .x$match,
            time = .x$time,
            pre_epa = .x$epa$total_points
        )}) |>
        list_rbind() |>
        arrange(time)
    return(result)
}

# Retrieve priors from an epa progression
get_priors <- function(epa_progression, team_list) {
    priors <- epa_progression |>
        group_by(team) |>
        arrange(time) |>
        summarize(team = first(team), initial_epas = first(pre_epa))

    result <- priors$initial_epas
    names(result) <- priors$team
    return(result)
}

# Retrieve OPR coefficients given a design matrix and a response vector
get_opr_coefs <- function(design, response){
    fit <- lm(response ~ 0 + design)
    coefs <- coef(fit)
    names(coefs) <- colnames(design)
    return(coefs)
}

# using a relatively short grid for computational considerations
# default to n_cores = 1 to avoid nested parallelization
get_pridge_coefs <- function(design, response, priors, n_cores = 1){
    # avoid having exactly 0 to reduce matrix singularity
    grid <- exp(seq(log(0.01), log(20), length.out = 100))

    mses <- pridge_lambda_cv(
        design, response, priors, grid, plot_mses = FALSE, n_cores = n_cores
    )
    lambda_star <- grid[which.min(mses)]
    pridge_coefs <- scoutR:::prior_ridge(
        X = design, y = response, lambda = lambda_star, beta_0 = priors
    )
    pridge_coefs <- round(pridge_coefs, 2)
    return(list(coefs = pridge_coefs, lambda = lambda_star))
}

get_epa_coefs <- function(epa_progression, event_key, i) {
    epa_coefs <- epa_progression |>
        mutate(match = as.numeric(substr(
            match, nchar(paste0(event_key, "_qm")) + 1, nchar(match)))) |>
        filter(match <= i) |> # inclusive because we're using pre_epa
        arrange(team, match) |>
        group_by(team) |>
        summarize(coef = last(pre_epa))
    result <- epa_coefs$coef
    names(result) <- epa_coefs$team

    return(result)
}

# Helper to build a blank NA result row, now including OPR columns
blank_result <- function(){
    data.frame(
        pridge_mse = NA, epa_mse = NA, opr_mse = NA,
        pridge_epa_pct_imp = NA, pridge_opr_pct_imp = NA,
        pridge_acc = NA, epa_acc = NA, opr_acc = NA,
        lambda_mean = NA, lambda_median = NA, lambda_sd = NA,
        lambda_max = NA, lambda_min = NA
    )
}

pred_mat <- function(lo, hi){
    result <- matrix(NA, nrow = length(lo:hi), ncol = 7)
    colnames(result) <- c("match_number",
                          "pridge_error", "epa_error", "opr_error",
                          "pridge_correct", "epa_correct", "opr_correct")
    return(result)
}

summarize_results <- function(result, lambdas){
    pridge_mse <- mean(result[, "pridge_error"] ^ 2)
    epa_mse    <- mean(result[, "epa_error"]    ^ 2)
    opr_mse    <- mean(result[, "opr_error"]    ^ 2)

    # pRidge improvement as a pct of each baseline's MSE
    pridge_epa_pct_imp <- ((epa_mse - pridge_mse) / epa_mse) * 100
    pridge_opr_pct_imp <- ((opr_mse - pridge_mse) / opr_mse) * 100

    pridge_acc <- mean(result[, "pridge_correct"])
    epa_acc    <- mean(result[, "epa_correct"])
    opr_acc    <- mean(result[, "opr_correct"])

    result <- data.frame(
        pridge_mse, epa_mse, opr_mse,
        pridge_epa_pct_imp, pridge_opr_pct_imp,
        pridge_acc, epa_acc, opr_acc,
        lambda_mean   = mean(lambdas,   na.rm = TRUE),
        lambda_median = median(lambdas, na.rm = TRUE),
        lambda_sd     = sd(lambdas,     na.rm = TRUE),
        lambda_max    = max(lambdas,    na.rm = TRUE),
        lambda_min    = min(lambdas,    na.rm = TRUE)
    )
    return(scoutR:::round_numerics(result))
}

# Main comparison function - now uses pre-loaded event data
next_match_pred <- function(event_key, event_data_entry){
    # Extract required data from cached structure
    matches <- event_data_entry$matches
    team_matches <- event_data_entry$team_matches
    teams <- event_data_entry$teams

    # Check if data is valid
    if (is.null(matches) || is.null(team_matches) || is.null(teams)) {
        return(blank_result())
    }

    team_list <- scoutR:::id2int(teams)

    epa_progression <- get_epa_progression(team_matches, event_key)

    full_design <- as.matrix(lineup_design_matrix(matches))
    full_response <- c(matches$blue_score, matches$red_score)

    priors <- get_priors(epa_progression, team_list)

    # OPR is first calculable once we have one row per col (team)
    lo <- floor(length(team_list) / 2) + 1
    # start two matches later so we have enough data to fit pridge and predict
    lo <- lo + 1
    hi <- nrow(matches)

    # Check if we have enough matches
    if (hi < lo) {
        return(blank_result())
    }

    result <- pred_mat(lo, hi)
    lambdas <- numeric(length(lo:hi))

    # fit on matches up until (i - 1), predict on match (i)
    # note - will never fit OPR on the full dataset; this is WAI
    for (i in lo:hi){
        # matrix indexing trick to interleave red and blue matches
        # R reads matrices column-first, so we can rbind and then read them off
        ridx <- c(rbind(1:(i - 1), (nrow(matches) + 1):(nrow(matches) + i - 1)))
        design <- full_design[ridx, ]
        response <- full_response[ridx]

        # compute coefs for pridge, EPA, and OPR on training data
        pridge_result <- get_pridge_coefs(design, response, priors)
        pridge_coefs <- pridge_result$coefs
        lambdas[i - lo + 1] <- pridge_result$lambda
        epa_coefs <- get_epa_coefs(epa_progression, event_key, i)
        opr_coefs <- get_opr_coefs(design, response)

        # predict on match (i)
        idx <- i - lo + 1
        result[idx, ] <- cbind(
            i,
            coef_error(pridge_coefs, full_design[i, ], full_response[i]),
            coef_error(epa_coefs,   full_design[i, ], full_response[i]),
            coef_error(opr_coefs,   full_design[i, ], full_response[i]),
            coef_acc(pridge_coefs, full_design[i, ],
                     full_design[i + nrow(matches), ],
                     full_response[i], full_response[i + nrow(matches)]),
            coef_acc(epa_coefs,   full_design[i, ],
                     full_design[i + nrow(matches), ],
                     full_response[i], full_response[i + nrow(matches)]),
            coef_acc(opr_coefs,   full_design[i, ],
                     full_design[i + nrow(matches), ],
                     full_response[i], full_response[i + nrow(matches)])
        )
    }

    return(summarize_results(result, lambdas))
}

#' Load year data from file
#' @param year Integer year
#' @return year_data list or NULL if file doesn't exist
load_year_data <- function(year) {
    filename <- file.path(data_dir, paste0("year_", year, ".rda"))

    if (!file.exists(filename)) {
        stop("No data file found for year ", year,
             ". Run pull_data.R first.")
    }

    load(filename)  # Loads 'year_data' object
    return(year_data)
}

##############
#### Data ####
##############

years <- setdiff(2016:2026, 2020:2021)

cat("=============================================\n")
cat("NEXT MATCH PREDICTION COMPARISON (LOCAL DATA)\n")
cat("=============================================\n")
cat("Loading cached data from:", data_dir, "\n")
cat("Years to process:", paste(years, collapse = ", "), "\n")
cat("=============================================\n")

# Load all year data
all_year_data <- list()
for (year in years) {
    cat("Loading data for year", year, "...")
    year_data <- load_year_data(year)
    all_year_data[[as.character(year)]] <- year_data
    cat(" Done\n")
}

# Build event_keys list and event_data lookup
event_keys <- c()
event_data_lookup <- list()

for (year in years) {
    year_data <- all_year_data[[as.character(year)]]
    year_event_keys <- names(year_data$event_data)
    event_keys <- c(event_keys, year_event_keys)

    # Add to lookup
    for (key in year_event_keys) {
        event_data_lookup[[key]] <- year_data$event_data[[key]]
    }
}

# Build qualifier_events metadata for final output
qualifier_events <- lapply(years, function(year) {
    all_year_data[[as.character(year)]]$events_metadata
}) |>
    bind_rows()

cat("\nTotal events to process:", length(event_keys), "\n")
cat("Starting analysis...\n\n")

######################
#### Parallelized ####
######################

n_cores <- parallel::detectCores() %/% 2
cl <- makeCluster(n_cores)
registerDoParallel(cl)

start <- Sys.time()

results_list <- foreach(
    key = event_keys,
    .packages = c("scoutR", "dplyr", "purrr"),
    .errorhandling = "pass"
) %dopar% {
    tryCatch(
        {
            event_data_entry <- event_data_lookup[[key]]
            next_match_pred(key, event_data_entry)
        },
        error = function(e){
            blank_result()
        }
    )
}

stopCluster(cl)

finish <- Sys.time()
execution_time <- finish - start

result <- results_list |>
    bind_rows() |>
    # coerce MSEs & pct_imp columns to numeric to replace strings with NA
    mutate(
        key = event_keys,
        pridge_mse         = as.numeric(pridge_mse),
        epa_mse            = as.numeric(epa_mse),
        opr_mse            = as.numeric(opr_mse),
        pridge_epa_pct_imp = as.numeric(pridge_epa_pct_imp),
        pridge_opr_pct_imp = as.numeric(pridge_opr_pct_imp)
    ) |>
    full_join(qualifier_events, by = "key") |>
    select(key, year, week,
           pridge_epa_pct_imp, pridge_opr_pct_imp,
           pridge_mse, epa_mse, opr_mse,
           pridge_acc, epa_acc, opr_acc,
           everything())

# Save results
output_dir <- "data/next_match_pred/"
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

save(result, execution_time, n_cores,
     file = paste0(output_dir, "pct_imp_",
                   years[1], "_to_", tail(years, 1), ".rda"))

cat("\n=====================================\n")
cat("COMPLETE\n")
cat("=====================================\n")
cat("Execution time:", round(execution_time, 2),
    attr(execution_time, "units"), "\n")
cat("Events processed:", nrow(result), "\n")
cat("Events with valid results:", sum(!is.na(result$pridge_epa_pct_imp)), "\n")
cat("Results saved to:", output_dir, "\n")
cat("=====================================\n")
