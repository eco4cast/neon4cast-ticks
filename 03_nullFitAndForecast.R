#=====================================================#
# This script calibrates the null model on the tick
# training data for both species at all plots. The 
# null model is a Bayesian state-space random walk.
# The data model is Poisson as we are modeling counts
# and the process error is assumed to be Gaussian. 
#
# The model is fit using JAGS. The forecast period (the
# 34 target weeks of 2019) are added to the end of each
# data set as NAs. Therefore, JAGS is fitting the model
# AND forecasting. 
#
# After calling JAGS (fit and forecast) the script saves
# the forecast period only in the three possible formats
# recommended by the EFI forecast standards team:
# netCDF, ensemble CSV, and summary CSV. 
#=====================================================#


library(tidyverse)
library(rjags)
library(ncdf4)
library(EML)
library(uuid)
library(emld)
library(lubridate)
emld::eml_version("eml-2.2.0")

# first load the target data set
data <- read.csv("ticks-targets.csv.gz")

# for the random walk all we need are the targets and yearWeek
data <- data %>% 
  select(all_of(c("yearWeek", "targetSpecies", "targetCount", "targetPlotID", "specificTarget")))

# next, we need to extract the targets into their respective groups
# so we need data frames for each species x plot combination, retain
# the NA rows (weeks without observations), and make sure that the 
# plots that have both species present are separated into differnt
# target data sets

specific.targets <- data %>% 
  pull(specificTarget) %>% 
  unique()

data.list.master <- list()
for(i in seq_along(specific.targets)){
  data.list.master[[i]] <- data %>% 
    filter(specificTarget == specific.targets[i])
}

# names list elments
data.list.master <- set_names(data.list.master, specific.targets)

# jags model code
random_walk <- " model {
  
  # process error precision prior - uninformative
  tau.process ~ dgamma(0.001, 0.001)
  
  # first latent state prior
  x[1] ~ dpois(x.ic)
  
  # data model
  for(t in 2:n){
    y[t] ~ dpois(x[t])
  }
  
  # process model
  for(t in 2:n){
    x[t] ~ dnorm(x[t-1], tau.process) T(0,)
  }
}" 

# random walk code for forecast
  # model: the jags model to fit
  # df: data for model
  # target.weeks: sequence of epi weeks to forecast 
  # thin: the number of samples to save
run_forecast <- function(model, df, target.weeks, thin){
  
  # check we are only modelling one species at one plot
  spp <- df %>% pull(targetSpecies) %>% unique() 
  spp <- spp[!is.na(spp)]
  if(length(spp) > 1) stop("Both species in data! \n", call. = FALSE)
  
  plot <- df %>% pull(targetPlotID) %>% unique()
  if(length(plot) > 1) stop("More than one plot in data! \n", call. = FALSE)
  
  # at this point we have passed both checks
  cat("Fitting", spp, "at", plot, "\n")
  
  # n.weeks the number of weeks to forecast into the future
  # need to pad df to make the forecast in jags
  na.pad <- rep(NA, length(target.weeks))
  
  jags.data <- list() # jags takes data in a list
  jags.data$n <- nrow(df) + length(target.weeks)   # number of observations + number of forecasts
  jags.data$y <- c(pull(df, targetCount), na.pad)  # counts, NAs at end are forecast weeks
  jags.data$x.ic <- rpois(1, 5)                    # randomly set the first latent state
  
  # initialize model 
  j.model <- jags.model(
    file = textConnection(model),
    data = jags.data,
    n.chains = 3
  )
  
  # run mcmc
  loops <- 1
  jags.out <- coda.samples(
    model = j.model,
    variable.names = c("x", "tau.process"),
    n.iter = 50000
  )
  
  # check convergence
  GBR <- gelman.plot(jags.out, ask = FALSE)
  burnin <- GBR$last.iter[tail(which(apply(GBR$shrink[, , 2] > 1.1, 1, any)), 1) + 1]
  # if not burned in, keep sampling
  while(is.na(burnin)){
    loops <- loops + 1
    cat("coda.samples call:", loops, "\n")
    jags.out <- coda.samples(
      model = j.model,
      variable.names = c("x", "tau.process"),
      n.iter = 25000
    )
    GBR <- gelman.plot(jags.out, ask = FALSE)
    burnin <- GBR$last.iter[tail(which(apply(GBR$shrink[, , 2] > 1.1, 1, any)), 1) + 1]
  }
  cat("Burnin after", burnin, "iterations\n")
  out <- window(jags.out, start = burnin) # remove burnin
  
  out.mat <- as.matrix(out) # convert to matrix
  
  # we are going to save a thinned matrix (raw mcmc objects can get big)
  save.seq <- round(seq(1, nrow(out.mat), length.out = thin)) # sequence of samples to keep
  out.mat <- out.mat[save.seq,]
  
  return(out.mat)
  
}

### set forecast variables ###

# Forecasts are submitted at the end of each month March - October.
# Therefore the start week for the forecast changes depending on when
# the forecast is submitted. 
start.epi.weeks <- c(10, 14, 19, 23, 28, 32, 36, 41) # all possible start weeks
day.run <- lubridate::today() # the day the script is called

# anytime we run this script before the start of the challenge we want to forecast all 2019 target weeks
if(day.run < "2021-03-31"){ 
  start.week <- start.epi.weeks[1]
} else { # otherwise use the appropriate starting week (months are 2 ahead)
  start.week <- start.epi.weeks[month(day.run) - 2]
}

end.week <- 44 # does not change
target.weeks <- start.week:end.week
n.weeks <- length(target.weeks)       # the number of weeks forecasted, same across all targets
obs.dim <- 1                          # 1 = latent state; 2 = latent state + observation error
n.ens <- 100                          # how many samples do we want to save? 
n.targets <- length(data.list.master) # the number of unique species_plot combinations 
targetName <- names(data.list.master) # the names of specific targets
forecast <- rep(0, n.weeks)           # is this a hindcast or forecast? hindcast = 0

# forecast identifiers used in EFI standards
forecast.issue.time <- as.numeric(paste0(2019, start.week))    # start date of the forecast YYYYWW
Forecast.id <- gsub("\\D", "", day.run)    # ISO datetime should make a valid Forecast_id
ForecastProject.id <- "tickGlobalNull_RandomWalk"     # Some ID that applies to a set of forecasts

# vector mapping which week data was assimilated for the forecast.
# if a week has a data point, that index = 1, otherwise 0
data.assimilation <- rep(0, n.weeks)

# initialize array for storage
fx.array <- array(data = NA, 
                  dim = c(n.ens,   # the number of ensembles we want to save      
                          n.weeks,   # the number of weeks forecasted
                          obs.dim,   # observation dimension
                          n.targets)) # the number of targets (species_site)

# now, run fit and forecast over target data sets individually
for(i in seq_len(n.targets)){ 
  cat(i, "/", n.targets, "\n")
  
  df <- pluck(data.list.master, i) # grab data set
  
  fit <- run_forecast(
    model = random_walk,    # fit random walk
    df = df,                # subsetted data
    target.weeks = target.weeks,            # forecast 35 weeks
    thin = n.ens  
  ) 
  
  # extract forecast period and save in list
  states <- fit[,-grep("tau.process", colnames(fit))] # remove process error column
  forecast.cols <- (nrow(df)+1):ncol(states) # columns that represent the forecast
  
  fx.array[ , , obs.dim, i] <- states[,forecast.cols] # store

}

### The first EFI Standard is netCDF ###

# each forecast will get its own dir 
dir.ncfname <- file.path("Random_Walk_Fits", as.character(forecast.issue.time)) 

if(!dir.exists(dir.ncfname)) dir.create(dir.ncfname, recursive = TRUE)

source("create_netCDF.R")
create_netCDF(
  fx = fx.array,
  ncfname = file.path(dir.ncfname, "random-walk-forecast-jags.nc"),
  data.assimilation = data.assimilation, 
  obs.dim = obs.dim, 
  forecast = forecast, 
  targetName = targetName, 
  forecast.issue.time = forecast.issue.time, 
  Forecast.id = Forecast.id, 
  ForecastProject.id = ForecastProject.id
  )


### The second EFI Standard is ensemble CSV ###

# need to create a flat file

fx.df <- fx.array[,,obs.dim,1] # pull out the first forecast to set everything up
colnames(fx.df) <- paste0(2019, target.weeks)
fx.df <- fx.df %>% 
  as_tibble() %>% 
  mutate(target = targetName[1],
         ensemble = 1:nrow(.),
         data_assimilation = 0,
         forecast.issue.time_YYYYWW = as.character(forecast.issue.time),
         ForecastProject.id = ForecastProject.id,
         Forecast.id = Forecast.id) %>% 
  pivot_longer(paste0(2019, target.weeks), 
               names_to = "time",
               values_to = "individuals")

# loop through rest of targets and bind rows
for(i in 2:length(targetName)){
  fx.df.i <- fx.array[,,obs.dim,i] # pull out the first forecast to set everything up
  colnames(fx.df.i) <- paste0(2019, target.weeks)
  fx.df.i <- fx.df.i %>% 
    as_tibble() %>% 
    mutate(target = targetName[i],
           ensemble = 1:nrow(.),
           data_assimilation = 0,
           forecast.issue.time_YYYYWW = as.character(forecast.issue.time),
           ForecastProject.id = ForecastProject.id,
           Forecast.id = Forecast.id) %>% 
    pivot_longer(paste0(2019, target.weeks), 
                 names_to = "time",
                 values_to = "individuals")
  
  fx.df <- bind_rows(fx.df, fx.df.i)
}

write.csv(fx.df,
          file = file.path(dir.ncfname, "random-walk-forecast-jags.csv"))

## Publish the forecast automatically. (EFI-only)
# source("../neon4cast-shared-utilities/publish.R")
# publish(code = "03_nullFitAndForecast.R",
#         data_in = "ticks-targets.csv.gz",
#         data_out = "random-walk-forecast-jags.csv.gz",
#         # meta = "meta/eml.xml", # haven't done this yet
#         prefix = "ticks/",
#         bucket = "targets")


### The third EFI Standard is summary CSV ###
fx.df.summary <- fx.df %>% 
  group_by(target, time, data_assimilation, 
           forecast.issue.time_YYYYWW, ForecastProject.id, Forecast.id) %>% 
  summarize(mean = mean(individuals),
            Conf_interv_02.5 = quantile(individuals, 0.025),
            Conf_interv_97.5 = quantile(individuals, 0.975)) %>% 
  pivot_longer(cols = c("mean", "Conf_interv_02.5", "Conf_interv_97.5"),
               names_to = "Statistic",
               values_to = "individuals")

write.csv(fx.df,
          file = file.path(dir.ncfname, "random-walk-forecast-summary-jags.csv"))


# publish(code = "03_nullFitAndForecast.R",
#         data_in = "ticks-targets.csv.gz",
#         data_out = "random-walk-forecast-summary-jags.csv.gz",
#         # meta = "meta/eml.xml", # haven't done this yet
#         prefix = "ticks/",
#         bucket = "targets")

