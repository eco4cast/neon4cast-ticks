#=====================================================#
# This script creates the null model for the Tick 
# Population challenge
#
# The model is historcial means for each week at each site
#=====================================================#

library(tidyverse)
library(lubridate)
library(MMWRweek)

efi_server <- TRUE

# first load the target data set
data <- read_csv("https://data.ecoforecast.org/neon4cast-targets/ticks/ticks-targets.csv.gz", guess_max = 1e6)
sites <- data %>% 
  pull(siteID) %>% 
  unique()


# the weeks we need to forecast
sundays <- seq.Date(ymd("2021-01-03"), by = 7, length.out = 52) # all Sundays in 2021
run.date <- today()
run.month <- month(run.date)

# forecast time (the sunday date stamp that maps to mmwrWeek)
fx.time <- sundays[month(sundays) == run.month]

# forecast mmwrWeeks
forecast.weeks <- MMWRweek(fx.time)$MMWRweek

# Forecast is weekly mean and sd by site
hist_means <- function(df, target.weeks){
  
  weekly.means <- df %>% 
    group_by(siteID, mmwrWeek) %>% 
    summarise(mean = mean(amblyomma_americanum),
              sd = sd(amblyomma_americanum))
  
  
  # need to fill in NAs and need to do on all data before sub-setting to target weeks
  filler.tb <- tibble(siteID = rep(sites, length(10:44)),
                      mmwrWeek = rep(10:44, each = length(sites)))
  
  all.weeks <- left_join(filler.tb, weekly.means, by = c("siteID", "mmwrWeek")) %>% 
    arrange(siteID, mmwrWeek)
  
  # not all weeks will be accounted for - linearly interpolate if missing
  gap.filled <- all.weeks %>%
    group_by(siteID) %>%
    mutate(mean = approx(x = mmwrWeek, y = mean, xout = mmwrWeek, rule = 2)$y,
           sd = replace_na(sd, mean(sd, na.rm = TRUE)),
           sd = replace(sd, sd == 0, mean(sd))) %>% 
    filter(mmwrWeek %in% target.weeks)
  
  # now each week has a mean - sort of
  return(gap.filled)
}



create_ensembles <- function(df, nmc = 500, forecast.year = 2021) {
  
  # simulate error from the log normal (to keep the zero bound)
  # need to calculate meanLog and sdLog from the normal mean and sd
  log_norm_sim <- function(df, i){
    mu <- df$mean[[i]] + 1
    sd <- df$sd[[i]]
    meanLog <- log(mu^2 / sqrt(sd^2 + mu^2))
    sdLog <- sqrt(log(1 + (sd^2 / mu^2)))
    sim <- data.frame(
      siteID = df$siteID[[i]],
      mmwrWeek = df$mmwrWeek[[i]],
      ensemble = 1:nmc,
      density = rlnorm(nmc, meanLog, sdLog)
    )
    return(sim)
  }
  
  # make ensemble data frame, create/rename columns as needed
  ens <- 1:nrow(df) %>% 
    map_dfr(~ log_norm_sim(df, .x)) %>%
    as_tibble() %>%
    mutate(year = forecast.year,
           time = MMWRweek2Date(year, mmwrWeek)) %>%
    rename(amblyomma_americanum = density)
  
  return(ens)
  
}

# get the forecasts we want
forecast <- hist_means(data, forecast.weeks)

# generate ensembles for uncertainty and scoring
ensembles <- create_ensembles(forecast)

# finalize for EFI submission
forecast.submit <- ensembles %>% 
  select(-year, -mmwrWeek) %>% 
  mutate(forecast = 1,
         data_assimilation = 0)


# Save file as CSV in the EFI format
# [theme_name]-[time]-[team_name].csv
theme_name <- "ticks"
time <- as.character(min(forecast.submit$time))
team_name <- "EFI_avg_null"
file.name <- paste0(theme_name, "-", time, "-", team_name, ".csv.gz")

write_csv(forecast.submit, file.name)



neon4cast::submit(forecast_file = file.name, 
                  metadata = NULL, 
                  ask = FALSE)


