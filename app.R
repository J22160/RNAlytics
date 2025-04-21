# TranscriptomicSuite/app.R
# Main entry point for the Shiny application


# Function to determine the script's directory
get_script_dir <- function() {
  # Check if running via Rscript
  args <- commandArgs(trailingOnly = FALSE)
  script_path <- sub("--file=", "", args[grep("--file=", args)])
  
  if (length(script_path) > 0) {
    # Rscript provides the full path to the script
    return(dirname(normalizePath(script_path)))
  } else {
    # Fallback for interactive sessions (e.g., RStudio or source())
    return(dirname(normalizePath(sys.frame(1)$ofile)))
  }
}

# Print working directory for debugging
print(paste("Working directory:", getwd()))

# Get the directory of this script
script_dir <- tryCatch(
  get_script_dir(),
  error = function(e) {
    # If both methods fail, assume current working directory
    message("Could not determine script directory dynamically. Using current working directory.")
    getwd()
  }
)

# Source all modular R scripts with absolute paths
source(file.path(script_dir, "R/global_config.R"))
source(file.path(script_dir, "R/ui_components.R"))
source(file.path(script_dir, "R/server_logic.R"))

# Launch the Shiny app
shinyApp(ui = appUI(), server = appServer)