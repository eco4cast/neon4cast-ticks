library(tidyverse)
library(neon4cast)
library(neonstore)
library(lubridate)
library(tools)

# functions for:
# downloading forecasts
# compile data by start date
# plotting forecasts 
# scoring
# plotting scores?



#' Download tick forecasts and scores from EFI server
#' 
#' @param dates vector of forecast start dates, default to all forecasts
#' @param dir the directory to store forecasts and scores
download_forecasts_and_scores <- function(dates = "all", dir = tempdir()){
  
  theme <- "ticks"
  
  if(dates == "all"){
    
    # need to update these as time goes on
    dates <- c(
      "2019-02-25",
      "2019-03-03",
      "2019-03-04",
      "2019-03-31",
      "2019-04-01",
      "2019-05-05",
      "2019-06-02",
      "2019-07-07",
      # "2021-03-24",
      "2021-03-31"
      # "2021-05-05"
    )
  }

  for(t in seq_along(dates)){
    message("Downloading ", dates[t])
    download_forecast(theme = theme,
                      date = dates[t],
                      dir = dir)
      
    download_scores(theme = theme,
                    date = dates[t],
                    dir = dir)
  }
}


#' Combine forecasts, scores, and observations into one data set
#' 
#' @param dest.file where the combined csv is written to
#' @param dir the directory the individual forecasts are stored
#' @param write.out do you want to write the combined data set to a csv?
#' @param get.forecasts do you want to download all forecasts
  
combine_forecasts <- function(
  dest.file = paste0("combinedForecasts_", lubridate::today()), 
  dir = tempdir(), 
  write.out = TRUE, 
  get.forecasts = FALSE
  ){
  
  base.dir <- dir
  
  if(get.forecasts){
    download_forecasts_and_scores(dir = base.dir)
  }

  fnames <- tibble(files = list.files(path = base.dir, recursive = TRUE, full.names = TRUE)) %>% 
    filter(file_ext(files) %in% c("nc","csv","gz")) %>% 
    filter(!str_detect(files, "not_in_standard")) %>% 
    mutate(basename = basename(files))
  
  d <-  unlist(str_split_fixed(tools::file_path_sans_ext(fnames$basename), pattern = "-", 5)) 
  
  data <- unlist(str_split_fixed(tools::file_path_sans_ext(fnames$basename), pattern = "-", 5)) %>% 
    as_tibble() %>% 
    unite("date", V2:V4, sep = "-") %>% 
    rename("theme" = V1,
           "team" = V5) %>% 
    bind_cols(fnames) %>% 
    select(-basename)
  
  data <- data %>% 
    filter(theme == "ticks",
           year(date) == 2019 | team == "SynergisTicks") 
  
  tick.data <- readr::read_csv("https://data.ecoforecast.org/targets/ticks/ticks-targets.csv.gz", guess_max = 1e6)
  w <- as.character(tick.data$epiWeek)
  w <- ifelse(nchar(w) == 2, w, paste0("0", w)) # add leading 0 to single digit weeks
  w <- paste0(tick.data$Year, "-W", w, "-1")
  tick.dataISO <- tick.data %>% 
    mutate(time = ISOweek::ISOweek2date(w)) 
  
  start.epi.weeks <- c(10, 14, 19, 23, 28, 32, 36, 41) # all possible start weeks
  y.vec <- rep(2019, length.out = length(start.epi.weeks))
  mmwr.start <- MMWRweek::MMWRweek2Date(y.vec, start.epi.weeks)
  iso.start <- ISOweek::ISOweek2date(paste0(as.character(y.vec), "-W", start.epi.weeks, "-1"))
  
  combined <- NULL
  for(i in 1:nrow(data)){
    message("Forecast ", i, " of ", nrow(data))
    
    d <- neon4cast:::read_forecast(data$files[i])
    d <- d %>% 
      filter(time >= "2019-03-01")
    
    if(data$team[i] == "Ticks_288.csv" & data$date[i] == "2019-02-25"){
      d <- d %>% group_by(plotID, time, statistic)
    }
    
    # which target file to use
    if(any(unique(d$time) %in% mmwr.start)){
      target.file <- tick.data
    } else if(any(unique(d$time) %in% iso.start)){
      target.file <- tick.dataISO
    }
    
    score <- neon4cast:::crps_score(d, target.file, 
                                    grouping_variables = c("plotID", "time"))
    
    for(s in 1:2){
      if(s == 1){
        species <- "ixodes_scapularis"
        target.plots <- c("BLAN_012","BLAN_005","SCBI_013","SCBI_002","SERC_001",
                          "SERC_005","SERC_006","SERC_012","ORNL_007")
      } else if (s == 2){
        species <- "amblyomma_americanum"
        target.plots <- c("SCBI_013","SERC_001","SERC_005","SERC_006","SERC_002",
                          "SERC_012","KONZ_025","UKFS_001","UKFS_004","UKFS_003",
                          "ORNL_002","ORNL_040","ORNL_008","ORNL_007","ORNL_009",
                          "ORNL_003","TALL_001","TALL_008","TALL_002")
      }
      
      score.spp <- score %>% 
        filter(plotID %in% target.plots,
               target == species)  
      
      d.spp <- d %>%
        filter(plotID %in% target.plots) %>%
        rename(fx.y = all_of(species))
      
      if("ensemble" %in% colnames(d)){
        d.spp <- d.spp %>% 
          group_by(time, plotID, forecast_start_time, horizon, team, theme) %>% 
          summarise(mean = mean(fx.y, na.rm = TRUE),
                    sd = sd(fx.y, na.rm = TRUE),
                    upper95 = quantile(fx.y, 0.975, na.rm = TRUE),
                    lower95 = quantile(fx.y, 0.025, na.rm = TRUE)) %>% 
          pivot_longer(cols = c("mean","sd", "upper95","lower95"), 
                       names_to = "statistic", 
                       values_to = "fx.y") %>% 
          select(time, plotID, forecast_start_time, horizon, team, theme, statistic, fx.y)
      }else{
        d.spp <- d.spp %>% 
          select(time, plotID, forecast_start_time, horizon, team, theme, statistic, fx.y) %>% 
          pivot_wider(names_from = statistic, 
                      values_from = fx.y, 
                      values_fn = base::mean) %>% 
          mutate(upper95 = mean + 1.96 * sd,
                 lower95 = mean - 1.96 * sd) %>% 
          pivot_longer(cols = c("mean","sd", "upper95","lower95"), 
                       names_to = "statistic", values_to = "fx.y")
      }
      
      d.spp <- d.spp %>% 
        pivot_wider(names_from = statistic, values_from = fx.y)
      
      d.score.spp <- left_join(d.spp, score.spp, 
                               by = c("team", "theme", "plotID", "time", "horizon"))
      
      target.file.spp <- target.file %>% 
        select(plotID, all_of(species), time) %>% 
        filter(time %in% unique(d.score.spp$time)) %>% 
        rename(actual = all_of(species))
      
      d.score.obs <- left_join(d.score.spp, target.file.spp,
                               by = c("plotID", "time"))
      
      d.score.obs <- d.score.obs %>% mutate(target = species)
      
      combined <- bind_rows(combined, d.score.obs)
    }
    
  }
  
  if(write.out){
    message("Writing csv...")
    write_csv(combined, 
              file = dest.file)
  }
  message("Done")
  return(combined)
}




