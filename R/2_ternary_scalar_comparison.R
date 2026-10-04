rm(list = ls())
library(scoutR)
library(parallel)
library(doParallel)
library(dplyr)

# ── Helper functions ───────────────────────────────────────────────────────────

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

ternary_loocv <- function(j, X, y, beta_0, lambda_vec, options) {
    mse_options <- sapply(options, function(lam) {
        lambda_trial <- lambda_vec
        lambda_trial[j] <- lam
        scoutR:::pridge_loocv(X, y, lambda_trial, beta_0, mse = TRUE)
    })
    list(best_lam = options[which.min(mse_options)],
         best_mse = min(mse_options))
}

ternary_descent_pass <- function(X, y, beta_0, lambda_vec, options,
                                 option_names, current_mse, verbose) {
    n_changed <- 0
    for (j in seq_len(ncol(X))) {
        trial <- ternary_loocv(j, X, y, beta_0, lambda_vec, options)
        if (trial$best_lam != lambda_vec[j]) {
            if (verbose) cat(sprintf(
                "  Team %s: %s -> %s (MSE: %.4f -> %.4f)\n",
                colnames(X)[j],
                option_names[which(options == lambda_vec[j])],
                option_names[which(options == trial$best_lam)],
                current_mse, trial$best_mse
            ))
            lambda_vec[j] <- trial$best_lam
            n_changed <- n_changed + 1
        }
        current_mse <- trial$best_mse
    }
    list(lambda_vec = lambda_vec, current_mse = current_mse,
         n_changed = n_changed)
}

reoptimize_lambda_star <- function(X, y, beta_0, lambda_vec,
                                   lambda_star, grid, verbose) {
    star_mask <- lambda_vec == lambda_star
    if (all(star_mask)) {
        if (verbose) cat("\nPost-hoc: all teams at lambda_star, skipping.\n")
        return(list(lambda_vec = lambda_vec, lambda_star = lambda_star))
    }
    if (!any(star_mask)) {
        if (verbose) cat("\nPost-hoc: no teams at lambda_star, skipping.\n")
        return(list(lambda_vec = lambda_vec, lambda_star = lambda_star))
    }
    if (verbose) cat("\nPost-hoc: re-optimizing lambda_star for non-pinned teams...\n")
    mses_refit <- sapply(grid, function(lam) {
        lambda_trial <- lambda_vec
        lambda_trial[star_mask] <- lam
        scoutR:::pridge_loocv(X, y, lambda_trial, beta_0, mse = TRUE)
    })
    lambda_star_refit <- grid[which.min(mses_refit)]
    lambda_vec[star_mask] <- lambda_star_refit
    if (verbose) cat(sprintf(
        "  lambda_star updated: %.4f -> %.4f\n", lambda_star, lambda_star_refit
    ))
    list(lambda_vec = lambda_vec, lambda_star = lambda_star_refit)
}

#' Ternary Vector Lambda Optimization for Prior Ridge
#'
#' Optimizes a per-team lambda vector where each team's lambda is restricted
#' to one of three values: the grid minimum, the scalar LOOCV-optimal lambda,
#' or the grid maximum. Uses greedy coordinate descent to make the search
#' tractable.
#'
#' @param X (matrix) 1-hot encoded design matrix, matches on rows, teams on cols
#' @param y (vector) response vector (typically alliance scores)
#' @param beta_0 (numeric) named prior vector aligned to colnames(X)
#' @param grid (numeric) grid of lambda values; min and max are used as the
#'   boundary options, and the LOOCV-optimal scalar is the middle option.
#' @param n_cores (int) cores for parallelization in lambda_cv step;
#'   NULL = max-1
#' @param max_iter (int) maximum number of coordinate descent passes
#' @param tol (numeric) convergence tolerance on LOOCV MSE improvement
#' @param verbose (bool) if TRUE, print progress during coordinate descent
#' @return list with:
#'   \item{coefs}{fitted prior ridge coefficients}
#'   \item{lambda_vec}{final per-team lambda vector}
#'   \item{lambda_scalar}{the scalar LOOCV-optimal lambda used as the
#'   middle option}
#'   \item{mse_history}{LOOCV MSE at the end of each coordinate descent pass}
#' @export
fit_pridge_ternary_lambda <- function(
        X, y, beta_0,
        grid = exp(seq(log(0.01), log(20), length.out = 100)),
        n_cores = 1,
        max_iter = 10,
        tol = 1e-4,
        verbose = TRUE
) {
    stopifnot("X must be a matrix" = is.matrix(X))
    stopifnot("beta_0 must match ncol(X)" = length(beta_0) == ncol(X))

    # ── Step 1: Scalar LOOCV to find lambda_star ──────────────────────────────
    if (verbose) cat("Step 1: Running scalar LOOCV over grid...\n")
    mses_scalar <- pridge_lambda_cv(
        design = X, response = y, priors = beta_0,
        grid = grid, plot_mses = FALSE, n_cores = n_cores
    )
    lambda_star <- grid[which.min(mses_scalar)]
    if (verbose) cat(sprintf("  lambda_star = %.4f\n", lambda_star))

    # ── Step 2: Initialize all teams to lambda_star ───────────────────────────
    lambda_vec <- rep(lambda_star, ncol(X))
    names(lambda_vec) <- colnames(X)
    current_mse <- scoutR:::pridge_loocv(X, y, lambda_vec, beta_0, mse = TRUE)
    if (verbose) cat(sprintf(
        "  Initial LOOCV MSE (all lambda_star): %.4f\n\n", current_mse
    ))

    # ── Steps 3 & 4: Greedy coordinate descent ────────────────────────────────
    options <- c(min(grid), lambda_star, max(grid))
    option_names <- c("lambda_min", "lambda_star", "lambda_max")
    mse_history <- numeric(max_iter)

    for (iter in seq_len(max_iter)) {
        if (verbose) cat(sprintf("Iteration %d:\n", iter))
        pass <- ternary_descent_pass(
            X, y, beta_0, lambda_vec, options, option_names, current_mse, verbose
        )
        lambda_vec  <- pass$lambda_vec
        prev_mse <- current_mse
        current_mse <- pass$current_mse
        mse_history[iter] <- current_mse
        if (verbose) cat(sprintf(
            "  End of iteration %d: LOOCV MSE = %.4f, teams changed = %d\n\n",
            iter, current_mse, pass$n_changed
        ))
        if (pass$n_changed == 0) {
            if (verbose) cat("Converged: no teams changed assignment.\n")
            break
        }
        if ((prev_mse - current_mse) < tol) {
            if (verbose) cat(sprintf(
                "Converged: MSE improvement below tolerance (%.6f).\n", tol
            ))
            break
        }
    }

    # ── Step 5: Post-hoc re-optimization of lambda_star ───────────────────────
    reopt <- reoptimize_lambda_star(
        X, y, beta_0, lambda_vec, lambda_star, grid, verbose
    )
    lambda_vec  <- reopt$lambda_vec
    lambda_star <- reopt$lambda_star

    # ── Fit final model ────────────────────────────────────────────────────────
    coefs <- prior_ridge(X, y, lambda_vec, beta_0)
    return(list(
        coefs         = round(coefs, 2),
        lambda_vec    = lambda_vec,
        lambda_scalar = lambda_star,
        mse_history   = mse_history[mse_history != 0]
    ))
}

# ── CV helpers ─────────────────────────────────────────────────────────────────

coefs_mse <- function(design_test, response_test, coefs) {
    preds <- design_test %*% coefs
    mean((drop(preds) - response_test) ^ 2)
}

get_priors <- function(team_events_data) {
    priors <- sapply(team_events_data, function(te) te$epa$stats$start)
    names(priors) <- sapply(team_events_data, function(te) te$team)
    return(priors)
}

cv_fold <- function(fold, fold_ids, design, response, priors, grid) {
    train_idx <- fold_ids != fold
    test_idx  <- fold_ids == fold

    design_train   <- design[train_idx, ]
    design_test    <- design[test_idx,  ]
    response_train <- response[train_idx]
    response_test  <- response[test_idx]

    # ── Scalar pRidge ──────────────────────────────────────────────────────────
    scalar_mses   <- pridge_lambda_cv(
        design = design_train, response = response_train, priors = priors,
        grid = grid, plot_mses = FALSE, n_cores = 1
    )
    lambda_scalar <- grid[which.min(scalar_mses)]
    scalar_coefs  <- prior_ridge(design_train, response_train,
                                 lambda_scalar, priors)
    scalar_mse    <- coefs_mse(design_test, response_test, scalar_coefs)

    # ── Ternary pRidge ─────────────────────────────────────────────────────────
    ternary_fit   <- fit_pridge_ternary_lambda(
        X = design_train, y = response_train, beta_0 = priors,
        grid = grid, n_cores = 1, verbose = FALSE
    )
    ternary_coefs <- prior_ridge(
        design_train, response_train, ternary_fit$lambda_vec, priors
    )
    ternary_mse   <- coefs_mse(design_test, response_test, ternary_coefs)

    list(
        scalar_mse    = scalar_mse,
        ternary_mse   = ternary_mse,
        lambda_scalar = lambda_scalar
    )
}

blank_result <- function() {
    data.frame(pct_imp = NA, scalar_mse = NA,
               ternary_mse = NA, lambda_opt = NA)
}

#' Ternary vs Scalar pRidge Comparison
#'
#' Runs 4-fold cross validation to compare the predictive performance of
#' ternary vector lambda pRidge against scalar lambda pRidge for a single
#' event. Mirrors the pattern of pridge_opr_comparison() from
#' 2_pridge_opr_comparison.R.
#'
#' @param event_key (chr) TBA-legal event key (e.g. "2025chcmp")
#' @param event_data_entry (list) cached event data entry with elements
#'   \code{matches} (qual match dataframe) and \code{team_events}
#'   (team_events_sb output for the event)
#' @param grid (numeric) grid of lambda values to search
#' @return single-row data.frame with columns:
#'   \item{pct_imp}{percent improvement of ternary over scalar pRidge MSE}
#'   \item{scalar_mse}{4-fold CV MSE for scalar lambda pRidge}
#'   \item{ternary_mse}{4-fold CV MSE for ternary lambda pRidge}
#'   \item{lambda_opt}{mean scalar lambda_opt across folds}
pridge_ternary_comparison <- function(
        event_key,
        event_data_entry,
        grid = exp(seq(log(0.01), log(20), length.out = 100))
) {
    matches     <- event_data_entry$matches
    team_events <- event_data_entry$team_events

    if (is.null(matches) || is.null(team_events)) return(blank_result())
    if (nrow(matches) < 4)                         return(blank_result())

    set.seed(449)

    priors <- get_priors(team_events)

    design   <- as.matrix(lineup_design_matrix(matches))
    response <- c(matches$blue_score, matches$red_score)

    # align priors to design matrix column order
    priors <- priors[match(colnames(design), scoutR:::tf(names(priors)))]

    k        <- 4
    fold_ids <- sample(rep(1:k, length.out = nrow(matches)))

    fold_results <- lapply(1:k, cv_fold,
                           fold_ids = fold_ids,
                           design   = design,
                           response = response,
                           priors   = priors,
                           grid     = grid)

    scalar_mse  <- mean(sapply(fold_results, function(x) x$scalar_mse))
    ternary_mse <- mean(sapply(fold_results, function(x) x$ternary_mse))
    lambda_opt  <- mean(sapply(fold_results, function(x) x$lambda_scalar))
    pct_imp     <- (scalar_mse - ternary_mse) / scalar_mse * 100

    return(data.frame(pct_imp, scalar_mse, ternary_mse, lambda_opt))
}

data_dir <- "data/raw/"

years <- setdiff(2016:2026, 2020:2021)

cat("================================================\n")
cat("TERNARY VS SCALAR LAMBDA COMPARISON (LOCAL DATA)\n")
cat("================================================\n")
cat("Loading cached data from:", data_dir, "\n")
cat("Years to process:", paste(years, collapse = ", "), "\n")
cat("================================================\n")

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

# Combined iterator over keys with their index for error attribution
results_list <- foreach(
    key = event_keys,
    .packages = c("scoutR"),
    .export = c(
        # user-defined functions
        "pridge_ternary_comparison",
        "cv_fold",
        "fit_pridge_ternary_lambda",
        "ternary_descent_pass",
        "ternary_loocv",
        "reoptimize_lambda_star",
        "coefs_mse",
        "get_priors",
        "blank_result"
    ),
    .errorhandling = "pass"
) %dopar% {
    tryCatch(
        {
            event_data_entry <- event_data_lookup[[key]]

            # Explicit NULL guard with an informative error rather than a
            # silent blank_result() so the logger below can distinguish
            # "data was missing" from "computation failed"
            if (is.null(event_data_entry)) {
                stop("event_data_entry is NULL: key not found in event_data_lookup")
            }

            pridge_ternary_comparison(key, event_data_entry)
        },
        error = function(e) {
            # Build a structured error log entry and attach it as an attribute
            # on blank_result() so the caller can separate logs from data
            result <- blank_result()
            attr(result, "error_log") <- list(
                key       = key,
                message   = conditionMessage(e),
                class     = class(e),
                timestamp = Sys.time()
            )
            result
        }
    )
}

stopCluster(cl)

finish <- Sys.time()
execution_time <- finish - start

# ── Separate clean results from errored ones ───────────────────────────────────

error_log <- Filter(
    Negate(is.null),
    lapply(results_list, function(x) attr(x, "error_log"))
)

n_errors <- length(error_log)
n_ok     <- length(results_list) - n_errors

if (n_errors > 0) {
    cat(sprintf(
        "\nWARNING: %d / %d cases returned errors:\n",
        n_errors, length(results_list)
    ))
    for (entry in error_log) {
        cat(sprintf(
            "  [%s] key='%s'  msg='%s'\n",
            format(entry$timestamp, "%H:%M:%S"),
            entry$key,
            entry$message
        ))
    }
} else {
    cat("\nAll cases completed without errors.\n")
}

# Save the error log alongside results so it can be inspected later
# (errors are stripped from results_list before bind_rows so the
#  data-frame structure stays clean)
for (i in seq_along(results_list)) {
    attr(results_list[[i]], "error_log") <- NULL
}

result <- results_list |>
    bind_rows() |>
    mutate(key = event_keys) |>
    full_join(qualifier_events, by = "key") |>
    select(key, year, week, pct_imp, scalar_mse,
           lambda_opt, ternary_mse, everything())

# Save results
output_dir <- "data/ternary_vs_scalar/"
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

save(result, execution_time, n_cores, error_log,
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
