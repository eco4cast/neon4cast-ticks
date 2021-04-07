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
# library(EML)
library(uuid)
# library(emld)
library(lubridate)
library(MMWRweek)
# emld::eml_version("eml-2.2.0")

efi_server <- TRUE

# first load the target data set
data <- read.csv("ticks-targets.csv.gz")

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
filter.week <- paste0("2019", start.week) %>% as.integer()

# next, we need to extract the targets into their respective groups
# so we need data frames for each species x plot combination, retain
# the NA rows (weeks without observations), and make sure that the 
# plots that have both species present are separated into differnt
# target data sets

ixodes.plots <- c(
  "BLAN_012",
  "BLAN_005",
  "SCBI_013",
  "SCBI_002",
  "SERC_001",
  "SERC_005",
  "SERC_006",
  "SERC_012",
  "ORNL_007"  
)

amblyomma.plots <- c(
  "SCBI_013",
  "SERC_001",
  "SERC_005",
  "SERC_006",
  "SERC_002",
  "SERC_012",
  "KONZ_025",
  "UKFS_001",
  "UKFS_004",
  "UKFS_003",
  "ORNL_002",
  "ORNL_040",
  "ORNL_008",
  "ORNL_007",
  "ORNL_009",
  "ORNL_003",
  "TALL_001",
  "TALL_008",
  "TALL_002"
)

# for the random walk all we need are the targets and yearWeek
data.list.ix <- list()
for(i in seq_along(ixodes.plots)){
  data.ix <- data %>% 
    select(all_of(c("yearWeek", "plotID", "ixodes_scapularis", "time"))) %>% 
    filter(plotID == ixodes.plots[i]) %>% 
    filter(yearWeek < filter.week)
  
  if(ixodes.plots[i] %in% amblyomma.plots){
    data.ix <- distinct(data.ix)
    counts <- data.ix[which(!is.na(data.ix$ixodes_scapularis)), "yearWeek"]
    data.ix <- data.ix %>% 
      filter(!(yearWeek %in% counts & is.na(ixodes_scapularis)))
  }
  
  data.list.ix[[i]] <- data.ix  
}
data.list.ix <- set_names(data.list.ix, paste0("ixodes_scapularis_", ixodes.plots))

data.list.aa <- list()
for(i in seq_along(amblyomma.plots)){
  data.aa <- data %>% 
    select(all_of(c("yearWeek", "plotID", "amblyomma_americanum", "time"))) %>% 
    filter(plotID == amblyomma.plots[i]) %>% 
    filter(yearWeek < filter.week)
  
  if(amblyomma.plots[i] %in% ixodes.plots){
    data.aa <- distinct(data.aa)
    counts <- data.aa[which(!is.na(data.aa$amblyomma_americanum)), "yearWeek"]
    data.aa <- data.aa %>% 
      filter(!(yearWeek %in% counts & is.na(amblyomma_americanum)))
  }
  
  data.list.aa[[i]] <- data.aa
}
data.list.aa <- set_names(data.list.aa, paste0("amblyomma_americanum_", amblyomma.plots))


# combine
data.list.master <- prepend(data.list.aa, data.list.ix)


  
# jags model code
random_walk <- " model {
  
  # process error precision prior - uninformative
  tau.process ~ dgamma(0.001, 0.001)
  
  # first latent state prior
  x[1] ~ dpois(x.ic)
  
  # data model
  for(t in 1:n){
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
  
  # n.weeks the number of weeks to forecast into the future
  # need to pad df to make the forecast in jags
  na.pad <- rep(NA, length(target.weeks))
  
  jags.data <- list() # jags takes data in a list
  jags.data$n <- nrow(df) + length(target.weeks)   # number of observations + number of forecasts
  jags.data$y <- c(df[,3], na.pad)  # counts, NAs at end are forecast weeks
  jags.data$x.ic <- rpois(1, 5)                    # randomly set the first latent state
  
  # initialize model 
  j.model <- jags.model(
    file = textConnection(model),
    data = jags.data,
    inits = list(x = jags.data$y),
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
end.week <- 44 # does not change
target.weeks <- start.week:end.week
n.weeks <- length(target.weeks)       # the number of weeks forecasted, same across all targets
obs.dim <- 1                          # 1 = latent state; 2 = latent state + observation error
n.ens <- 100                          # how many samples do we want to save? 
n.targets <- length(data.list.master) # the number of unique species_plot combinations 
targetName <- names(data.list.master) # the names of specific targets
forecast <- rep(0, n.weeks)           # is this a hindcast or forecast? hindcast = 0

# forecast identifiers used in EFI standards
forecast.issue.time <- paste0(2019,"-", start.week)    # start date of the forecast YYYYWW
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
  # for(i in 1:2){ 
  cat(i, "/", n.targets, "\n")
  cat("Fitting", names(data.list.master)[i], "\n")
  
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

#source("create_netCDF.R")
#create_netCDF(
#  fx = fx.array,
#  ncfname = file.path(dir.ncfname, "random-walk-forecast-jags.nc"),
#  data.assimilation = data.assimilation, 
#  obs.dim = obs.dim, 
#  forecast = forecast, 
#  targetName = targetName, 
#  forecast.issue.time = forecast.issue.time, 
#  Forecast.id = Forecast.id, 
#  ForecastProject.id = ForecastProject.id
#  )


### The second EFI Standard is ensemble CSV ###

# need to create a flat file

# get vector of species names, first extract Ixodes_scapularis
species.name <- str_extract(targetName, "ixodes_scapularis")

# NAs are Amblyomma_americanum
species.name[is.na(species.name)] <- "amblyomma_americanum"

# extract plotIDs
plot.id <- str_extract(targetName, "[[:upper:]]{4}_\\d{3}")


# convert week to date mapping to the first day of the week
date.col <- MMWRweek::MMWRweek2Date(rep(2019, length(target.weeks)), target.weeks) %>% as.character()

plot.id.unique <- plot.id %>% unique()
fx.df <- tibble()
for(p in seq_along(plot.id.unique)){
  plot.subset <- plot.id.unique[p]
  site.subset <- str_extract(plot.subset, "[[:upper:]]{4}")
  fx.index <- grep(plot.subset, targetName)
  
  # if only one species present at the plot
  if(length(fx.index) == 1){
    fx <- fx.array[,,obs.dim,fx.index]
    colnames(fx) <- date.col
    fx <- fx %>% 
      as_tibble() %>% 
      pivot_longer(all_of(date.col), 
                   names_to = "time",
                   values_to = species.name[fx.index]) %>% 
      mutate(plotID = plot.subset,
             siteID = site.subset,
             ensemble = rep(1:n.ens, each = length(target.weeks)),
             data_assimilation = 0,
             forecast = 0,
             obs_flag = obs.dim) 
    
    # if both species are present at the plot
  } else if (length(fx.index == 2)){
    
    fx.1 <- fx.array[,,obs.dim,fx.index[1]] 
    fx.2 <- fx.array[,,obs.dim,fx.index[2]] 
    colnames(fx.1) <- colnames(fx.2) <- date.col
    
    fx.1 <- fx.1 %>% 
      as_tibble() %>% 
      pivot_longer(all_of(date.col), 
                   names_to = "time",
                   values_to = species.name[fx.index[1]]) %>% 
      mutate(plotID = plot.subset,
             siteID = site.subset,
             ensemble = rep(1:n.ens, each = length(target.weeks)),
             data_assimilation = 0,
             forecast = 0,
             obs_flag = obs.dim)
    
    fx.2 <- fx.2 %>% 
      as_tibble() %>% 
      pivot_longer(all_of(date.col), 
                   names_to = "time",
                   values_to = species.name[fx.index[2]]) %>% 
      mutate(plotID = plot.subset,
             ensemble = rep(1:n.ens, each = length(target.weeks)),
             data_assimilation = 0,
             forecast = 0,
             obs_flag = obs.dim)
    
    fx <- left_join(fx.1,
                    fx.2,
                    by = c("time", "plotID", "ensemble", "data_assimilation", "forecast", "obs_flag"))
    
  }
  fx.df <- bind_rows(fx.df, fx)
}

#fx.df$siteID <- str_split_fixed(fx.df$plotID, "_", 2)[, 1]

# Save file as CSV in the
# [theme_name]-[yearWeek]-[team_name].csv
fx.file.name <- paste0("ticks-", 
                       as.character(date.col[1]), 
                       "-", 
                       ForecastProject.id, 
                       ".csv.gz")

write.csv(fx.df,
          file = file.path(dir.ncfname, fx.file.name))

#'Publish the forecast automatically.  Run only on EFI Challenge server
if(efi_server){
  source("../neon4cast-shared-utilities/publish.R")
  publish(code = "03_nullFitAndForecast.R",
          data_in = "ticks-targets.csv.gz",
          data_out = file.path(dir.ncfname, fx.file.name),
          prefix = "ticks/",
          bucket = "forecasts",
          registries = "https://hash-archive.carlboettiger.info")
}


### The third EFI Standard is summary CSV ###
#dfs <- fx.df %>% 
#  pivot_longer(cols = all_of(c("Ixodes_scapularis", "Amblyomma_americanum")), 
#               names_to = "species") %>% # make a species column
#  group_by(time, plot, obs_flag, species, data_assimilation, forecast) %>% 
#  summarize(mean = mean(value, na.rm = TRUE), # need the na.rm=TRUE for plots without both spp present
#            sd = sd(value, na.rm = TRUE),
#            Conf_interv_02.5 = quantile(value, 0.025, na.rm = TRUE),
#            Conf_interv_97.5 = quantile(value, 0.975, na.rm = TRUE)) %>%
#  pivot_longer(all_of(c("mean", "sd", "Conf_interv_02.5", "Conf_interv_97.5")),
#               names_to = "statistic") %>% # statistic column
#  pivot_wider(names_from = species,
#              values_from = value) # back to species wide

# [theme_name]-[yearWeek]-[team_name]-summary.csv
#fx.file.name <- paste0("ticks-", 
#                       as.character(date.col[1]), 
#                       "-", ForecastProject.id, 
#                       "-summary.csv.gz")

#write.csv(dfs,
#          file = file.path(dir.ncfname, fx.file.name))

#'Publish the forecast automatically.  Run only on EFI Challenge server
# if(efi_server){
#   source("../neon4cast-shared-utilities/publish.R")
#   publish(code = "03_nullFitAndForecast.R",
#           data_in = "ticks-targets.csv.gz",
#           data_out = file.path(dir.ncfname, fx.file.name),
#           prefix = "ticks/",
#           bucket = "forecasts")
# }

