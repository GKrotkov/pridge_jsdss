# Computing Prior Ridge percent improvement over OPR in terms of nested CV MSE
# Compares for a range of years of interest
rm(list = ls())
library(scoutR)
library(foreach)
library(doParallel)
library(dplyr)

# given a numeric vector of coefficients, compute the mse on the test set
coefs_mse <- function(design, response, coefs){
    preds <- design %*% coefs
    error <- drop(preds - response)
    return(mean(error ^ 2))
}

get_priors <- function(event_key){
    epa <- team_events_sb(event = event_key)
    priors <- sapply(epa, function(te){te$epa$stats$start})
    names(priors) <- sapply(epa, function(te){te$team})
    return(priors)
}

# compute CV error for a given fold
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

    list(pridge_mse = pridge_mse, opr_mse = opr_mse, lambda_opt = lambda_opt)
}

# Compute Pridge improvement % over OPR in terms of MSE
# Selecting a default of k = 4 to reduce the chance of a singular matrix and
# for computational cost
pridge_opr_pct_improvement <- function(event_key, k = 4){
    # Setting seed here so it takes hold inside the parallelized runs
    set.seed(449)
    # Collect required data
    matches <- event_matches(event_key, match_type = "quals")
    epa <- team_events_sb(event = event_key)
    priors <- sapply(epa, function(te){te$epa$stats$start})
    names(priors) <- sapply(epa, function(te){te$team})

    # Assign matches to folds
    fold_ids <- sample(rep(1:k, length.out = nrow(matches)))

    fold_results <- lapply(1:k, cv_fold, fold_ids, matches, priors)

    # Aggregate MSEs across folds
    pridge_mse <- mean(sapply(fold_results, function(x) x$pridge_mse))
    opr_mse    <- mean(sapply(fold_results, function(x) x$opr_mse))
    lambda_opt <- mean(sapply(fold_results, function(x) x$lambda_opt))

    pct_imp <- ((opr_mse - pridge_mse) / opr_mse) * 100
    return(data.frame(pct_imp, pridge_mse, opr_mse, lambda_opt))
}

##############
#### Data ####
##############

years <- setdiff(2016:2025, 2020:2021)

qualifier_events <- lapply(years, events, official = TRUE) |>
    dplyr::bind_rows() |>
    dplyr::filter(event_type %in% c(0, 1))

event_keys <- qualifier_events |>
    dplyr::pull(key)

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
        {pridge_opr_pct_improvement(key)},
        error = function(e){
            data.frame(pct_imp = NA, pridge_mse = NA,
                       opr_mse = NA, lambda_opt = NA)
        }
    )
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

save(result, execution_time, n_cores,
     file = paste0("data/pridge_vs_opr/", "pct_imp_",
                   years[1], "_to_", tail(years, 1),
                   ".rda"))
