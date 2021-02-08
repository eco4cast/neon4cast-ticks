# ========================================================== #
#       Workflow script for the EFI Tick Challenge           #
# ========================================================== #

# install dependencies 
devtools::install_deps()

Sys.setenv("NEONSTORE_HOME" = "/efi_neon_challenge/neonstore")
Sys.setenv("NEONSTORE_DB" = "/efi_neon_challenge/neonstore")

library(neonstore)
library(tidyverse)
library(lubridate)

run_full_workflow <- TRUE
generate_null <- TRUE



# tick data product, target sites, and dates
product <- "DP1.10093.001" 
target.sites <- c("BLAN", "ORNL", "SCBI", "SERC", "KONZ", "TALL", "UKFS")
start.date <- NA
end.date <- "2019-12-31"


message("Downloading: Tick data (DP1.10093.001)")
new_data1 <- neonstore::neon_download(product = product,
                                      site = target.sites, 
                                      type = "basic", 
                                      start_date = start.date, 
                                      end_date = end.date,
                                      .token = Sys.getenv("NEON_TOKEN"))

if(!is.null(new_data1) | run_full_workflow){
  
  message(paste0("Running Creating Tick Targets at ", Sys.time()))
  source("02_ticks_targets.R")
  message(paste0("Completed Tick Target at ", Sys.time()))
  
  if(generate_null){
    
    message(paste0("Running Tick Null Forecast at ", Sys.time()))
    source("03_nullFitAndForecast.R")
    message(paste0("Completed Tick Null Forecast at ", Sys.time()))
    
  }
}