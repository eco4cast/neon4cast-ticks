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
  pull(site_id) %>% 
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
    group_by(site_id, mmwr_week) %>% 
    summarise(mean = mean(observed),
              sd = sd(observed))
  
  
  # need to fill in NAs and need to do on all data before sub-setting to target weeks
  filler.tb <- tibble(site_id = rep(sites, length(10:44)),
                      mmwr_week = rep(10:44, each = length(sites)))
  
  all.weeks <- left_join(filler.tb, weekly.means, by = c("site_id", "mmwr_week")) %>% 
    arrange(site_id, mmwr_week)
  
  # not all weeks will be accounted for - linearly interpolate if missing
  gap.filled <- all.weeks %>%
    group_by(site_id) %>%
    mutate(mean = approx(x = mmwr_week, y = mean, xout = mmwr_week, rule = 2)$y,
           sd = replace_na(sd, mean(sd, na.rm = TRUE)),
           sd = replace(sd, sd == 0, mean(sd))) %>% 
    filter(mmwr_week %in% target.weeks)
  
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
      site_id = df$site_id[[i]],
      mmwr_week = df$mmwr_week[[i]],
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
           time = MMWRweek2Date(year, mmwr_week)) %>%
    rename(amblyomma_americanum = density)
  
  return(ens)
  
}

# get the forecasts we want
forecast <- hist_means(data, forecast.weeks)

# generate ensembles for uncertainty and scoring
ensembles <- create_ensembles(forecast)

# finalize for EFI submission
forecast.submit <- ensembles %>% 
  select(-year, -mmwr_week) %>% 
  rename(predicted = amblyomma_americanum) |> 
  mutate(variable = "amblyomma_americanum") |> 
  mutate(start_time = lubridate::as_date(min(time)) - lubridate::weeks(1)) |> 
  select(time, start_time, site_id, variable, ensemble, predicted)

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


