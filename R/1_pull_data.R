# pull_data.R
# Fetch and store event data locally for pridge analysis pipelines
# Data is stored by year in separate .rda files

rm(list = ls())
library(scoutR)
library(dplyr)
library(purrr)

# Configuration
years <- setdiff(2016:2026, 2020:2021)  # Skip COVID years

data_dir <- "data/raw/"

if (!dir.exists(data_dir)) {
    dir.create(data_dir, recursive = TRUE)
    cat("Created directory:", data_dir, "\n")
}

#' Fetch all data for a single event
#' @param event_key String event key (e.g., "2024casj")
#' @return List containing matches, team_events, and team_matches data
fetch_event_data <- function(event_key) {
    cat("  Fetching data for event:", event_key, "\n")

    result <- tryCatch({
        # Fetch qualification matches
        matches <- event_matches(event_key, match_type = "quals")

        # Fetch team event data (for priors)
        team_events <- team_events_sb(event = event_key)

        # Fetch team match progressions (for EPA)
        team_matches <- team_matches_sb(event = event_key, elim = FALSE)

        # Get team list
        teams <- event_teams(event_key, keys = TRUE)

        list(
            matches = matches,
            team_events = team_events,
            team_matches = team_matches,
            teams = teams,
            fetch_success = TRUE,
            fetch_timestamp = Sys.time(),
            error_message = NA
        )
    }, error = function(e) {
        warning("Failed to fetch data for ", event_key, ": ", e$message)
        list(
            matches = NULL,
            team_events = NULL,
            team_matches = NULL,
            teams = NULL,
            fetch_success = FALSE,
            fetch_timestamp = Sys.time(),
            error_message = e$message
        )
    })

    return(result)
}

#' Fetch and store all data for a single year
#' @param year Integer year (e.g., 2024)
#' @return List containing all event data and metadata
fetch_year_data <- function(year) {
    cat("\n========================================\n")
    cat("Processing year:", year, "\n")
    cat("========================================\n")

    start_time <- Sys.time()

    # Get all qualifier events for this year
    qualifier_events <- events(year, official = TRUE) |>
        dplyr::filter(event_type %in% c(0, 1))  # Regional (0) and District (1)

    cat("Found", nrow(qualifier_events), "qualifier events\n\n")

    # Fetch data for each event
    event_keys <- qualifier_events$key
    event_data <- list()

    for (i in seq_along(event_keys)) {
        key <- event_keys[i]
        cat(sprintf("[%d/%d] ", i, length(event_keys)))
        event_data[[key]] <- fetch_event_data(key)

        # Be nice to the API - add a small delay
        Sys.sleep(0.1)
    }

    # Count successes and failures
    n_success <- sum(sapply(event_data, function(x) x$fetch_success))
    n_failed <- length(event_data) - n_success

    end_time <- Sys.time()
    execution_time <- end_time - start_time

    cat("\n----------------------------------------\n")
    cat("Year", year, "summary:\n")
    cat("  Successfully fetched:", n_success, "events\n")
    cat("  Failed:", n_failed, "events\n")
    cat("  Execution time:", round(execution_time, 2),
        attr(execution_time, "units"), "\n")
    cat("----------------------------------------\n")

    # Package everything together
    year_data <- list(
        year = year,
        events_metadata = qualifier_events,
        event_data = event_data,
        fetch_metadata = list(
            fetch_timestamp = end_time,
            execution_time = execution_time,
            n_events_total = nrow(qualifier_events),
            n_events_success = n_success,
            n_events_failed = n_failed,
            failed_events =
                event_keys[!sapply(event_data, function(x) x$fetch_success)]
        )
    )

    return(year_data)
}

#' Save year data to file
#' @param year_data List returned by fetch_year_data()
save_year_data <- function(year_data) {
    year <- year_data$year
    filename <- file.path(data_dir, paste0("year_", year, ".rda"))

    save(year_data, file = filename)
    cat("\nSaved data for year", year, "to:", filename, "\n")
    cat("File size:", round(file.size(filename) / 1024^2, 2), "MB\n")
}

#' Load year data from file
#' @param year Integer year
#' @return year_data list or NULL if file doesn't exist
load_year_data <- function(year) {
    filename <- file.path(data_dir, paste0("year_", year, ".rda"))

    if (!file.exists(filename)) {
        warning("No data file found for year ", year)
        return(NULL)
    }

    load(filename)  # Loads 'year_data' object
    return(year_data)
}

#' Check if data exists and is recent for a given year
#' @param year Integer year
#' @param max_age_days Maximum age in days before data is considered stale
#' @return Logical indicating if fresh data exists
check_data_freshness <- function(year, max_age_days = 30) {
    filename <- file.path(data_dir, paste0("year_", year, ".rda"))

    if (!file.exists(filename)) {
        return(FALSE)
    }

    year_data <- load_year_data(year)
    fetch_time <- year_data$fetch_metadata$fetch_timestamp
    age_days <- as.numeric(difftime(Sys.time(), fetch_time, units = "days"))

    return(age_days <= max_age_days)
}

################################################################################
# Main Execution
################################################################################

cat("=====================================\n")
cat("PULL DATA SCRIPT\n")
cat("=====================================\n")
cat("Years to process:", paste(years, collapse = ", "), "\n")
cat("Data directory:", data_dir, "\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("=====================================\n")

overall_start <- Sys.time()

# Process each year
for (year in years) {
    # Check if we already have recent data
    if (check_data_freshness(year, max_age_days = 7)) {
        cat("\nSkipping year", year, "- fresh data already exists\n")
        cat("(Delete the file or adjust max_age_days to re-fetch)\n")
        next
    }

    # Fetch and save data
    year_data <- fetch_year_data(year)
    save_year_data(year_data)
}

overall_end <- Sys.time()
overall_time <- overall_end - overall_start
cat("\n=====================================\n")
cat("COMPLETE\n")
cat("=====================================\n")
cat("Total execution time:", round(overall_time, 2),
    attr(overall_time, "units"), "\n")
cat("Data stored in:", data_dir, "\n")

# Print summary of all stored files
cat("\nStored files:\n")
files <- list.files(data_dir, pattern = "^year_[0-9]{4}\\.rda$",
                    full.names = TRUE)
for (f in files) {
    cat("  ", basename(f), "-", round(file.size(f) / 1024^2, 2), "MB\n")
}

cat("=====================================\n")
