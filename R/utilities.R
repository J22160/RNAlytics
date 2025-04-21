# TranscriptomicSuite/R/utilities.R
# Helper functions for common tasks

# Validate uploaded data
validateData <- function(df) {
  if (ncol(df) < 2 || !all(sapply(df, is.numeric))) {
    stop("File must have at least 2 numeric sample columns.")
  }
}

# Show custom notification
showNotification <- function(message, type) {
  shinyalert(title = ifelse(type == "success", "Success", "Error"), text = message, type = type)
}