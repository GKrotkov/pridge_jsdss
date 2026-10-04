# pridge_opr_comparison.R
# Computing Prior Ridge percent improvement over OPR in terms of nested CV MSE
# Uses locally cached data rather than live API calls
# Last update: 5/3/2026

rm(list = ls())
library(scoutR)
library(foreach)
library(doParallel)
library(dplyr)

# Configuration
data_dir <- "data/raw/"

# given a numeric vector of coefficients, compute the mse on the test set
coefs_mse <- function(design, response, coefs){
    preds <- design %*% coefs
    error <- drop(preds - response)
    return(mean(error ^ 2))
}

# given a numeric vector of coefficients, compute the winner prediction accuracy
# on the test set. Assumes response and design rows are ordered blue then red.
coefs_acc <- function(design, response, coefs){
    n <- nrow(design) / 2
    preds <- drop(design %*% coefs)
    correct <- (preds[1:n] > preds[(n + 1):(2 * n)]) ==
        (response[1:n] > response[(n + 1):(2 * n)])
    return(mean(correct))
}

# Extract priors from cached team_events data
get_priors <- function(team_events_data){
    priors <- sapply(team_events_data, function(te){te$epa$stats$start})
    names(priors) <- sapply(team_events_data, function(te){te$team})
    return(priors)
}

# compute CV error metrics for a given fold
cv_fold <- function(fold, fold_ids, matches, priors){
    test  <- matches[fold_ids == fold, ]
    train <- matches[fold_ids != fold, ]
    response_train <- c(train$blue_score, train$red_score)
    response_test <- c(test$blue_score, test$red_score)
    design <- as.matrix(lineup_design_matrix(matches))
    design_train <- design[fold_ids != fold, ]
    design_test <- design[fold_ids == fold, ]

    # select lambda via LOOCV on the training fold
    # use n_cores = 1 to avoid nested parallelization
    pridge_cv <- scoutR::pridge_lambda_cv(
        design_train, response_train, priors,
        grid = seq(0, 20, length.out = 100), plot_mses = FALSE, n_cores = 1
    )
    lambda_opt <- as.numeric(names(which.min(pridge_cv)))

    # Evaluate pridge on test fold
    pridge_coefs <- scoutR:::prior_ridge(
        design_train, response_train, lambda_opt, priors
    )
    pridge_mse <- coefs_mse(design_test, response_test, pridge_coefs)

    # Fit OPR on training fold and evaluate on test fold
    opr_fit <- fit_lineup_lm(train, list(red = train[, "red_score"][[1]],
                                         blue = train[, "blue_score"][[1]]))
    opr_mse <- coefs_mse(design_test, response_test, coef(opr_fit))

    # Compute winner prediction accuracy on test fold
    pridge_acc <- coefs_acc(design_test, response_test, pridge_coefs)
    opr_acc    <- coefs_acc(design_test, response_test, coef(opr_fit))

    result <- list(
        pridge_mse = pridge_mse, opr_mse = opr_mse, lambda_opt = lambda_opt,
        pridge_acc = pridge_acc, opr_acc = opr_acc
    )
    return(result)
}

# Compute Pridge improvement % over OPR in terms of MSE
# Now uses pre-loaded event data instead of API calls
pridge_opr_comparison <- function(event_key, event_data_entry){
    # Setting seed here so it takes hold inside the parallelized runs
    set.seed(449)

    # Extract required data from cached structure
    matches <- event_data_entry$matches
    team_events <- event_data_entry$team_events

    # Check if data is valid
    if (is.null(matches) || is.null(team_events)) {
        return(data.frame(pct_imp = NA, pridge_mse = NA, opr_mse = NA,
                          lambda_opt = NA, pridge_acc = NA, opr_acc = NA))
    }

    # Check if we have enough matches
    if (nrow(matches) < 4) {
        return(data.frame(pct_imp = NA, pridge_mse = NA, opr_mse = NA,
                          lambda_opt = NA, pridge_acc = NA, opr_acc = NA))
    }

    priors <- get_priors(team_events)

    # Assign matches to folds
    k <- 4
    fold_ids <- sample(rep(1:k, length.out = nrow(matches)))

    fold_results <- lapply(1:k, cv_fold, fold_ids, matches, priors)

    # Aggregate metrics (MSE + Accuracy) across folds
    pridge_mse <- mean(sapply(fold_results, function(x) x$pridge_mse))
    opr_mse    <- mean(sapply(fold_results, function(x) x$opr_mse))
    lambda_opt <- mean(sapply(fold_results, function(x) x$lambda_opt))
    pridge_acc <- mean(sapply(fold_results, function(x) x$pridge_acc))
    opr_acc    <- mean(sapply(fold_results, function(x) x$opr_acc))

    pct_imp <- ((opr_mse - pridge_mse) / opr_mse) * 100
    return(data.frame(pct_imp, pridge_mse, opr_mse,
                      lambda_opt, pridge_acc, opr_acc))
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

cat("=====================================\n")
cat("PRIDGE vs OPR COMPARISON (LOCAL DATA)\n")
cat("=====================================\n")
cat("Loading cached data from:", data_dir, "\n")
cat("Years to process:", paste(years, collapse = ", "), "\n")
cat("=====================================\n\n")

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
    dplyr::bind_rows()

cat("\nTotal events to process:", length(event_keys), "\n")
cat("Starting analysis...\n\n")

######################
#### Parallelized ####
######################

n_cores <- parallel::detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)

start <- Sys.time()

results_list <- foreach(
    key = event_keys,
    .packages = c("scoutR"),
    .errorhandling = "pass"
) %dopar% {
    tryCatch(
        {
            event_data_entry <- event_data_lookup[[key]]
            pridge_opr_comparison(key, event_data_entry)
        },
        error = function(e){
            data.frame(pct_imp = NA, pridge_mse = NA, opr_mse = NA,
                       lambda_opt = NA, pridge_acc = NA, opr_acc = NA)
        }    )
}

stopCluster(cl)

finish <- Sys.time()
execution_time <- finish - start

result <- results_list |>
    bind_rows() |>
    mutate(key = event_keys) |>
    full_join(qualifier_events, by = "key") |>
    select(key, year, week, pct_imp, pridge_mse,
           lambda_opt, opr_mse, everything())

# Save results
output_dir <- "data/pridge_vs_opr/"
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
cat("Events with valid results:", sum(!is.na(result$pct_imp)), "\n")
cat("Results saved to:", output_dir, "\n")
cat("=====================================\n")
