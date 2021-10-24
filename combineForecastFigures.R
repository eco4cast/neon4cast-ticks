library(tidyverse)
library(tools)

source("combine_forecast_functions.R")

# download_forecasts_and_scores()

combined <- combine_forecasts(get.forecasts = FALSE)
# combined <- read_csv("combinedForecasts_2021-07-13")

ticks <- read_csv("https://data.ecoforecast.org/targets/ticks/ticks-targets.csv.gz")
lat.lon <- ticks %>% 
  select(c(plotID, decimalLatitude, decimalLongitude)) %>% 
  mutate(decimalLatitude = decimalLatitude, 
         decimalLongitude = decimalLongitude) %>% 
  distinct()

nlcd <- ticks %>% 
  select(plotID, nlcdClass) %>% 
  distinct()


# combined <- left_join(combined, lat.lon, by = c("plotID"))

combined <- combined %>% 
  rename(Team = team) %>% 
  mutate(across(Team, ~ if_else(.x == "tickGlobalNull_RandomWalk", "Random Walk", .x)))
  
str(combined)

ix <- "ixodes_scapularis"
aa <- "amblyomma_americanum"

start <- c(
  "2019-02-18", # ix & aa only Ticks_288
  "2019-02-24", # ix & aa NJC
  "2019-02-25", # ix & aa 5 teams
  "2019-03-24", # ix & aa BU and null
  "2019-04-01", # ix & aa only NJC
  "2019-03-25", # ix & aa only null
  "2019-04-28", # ix & aa BU, NJC, null
  "2019-05-26", # ix & aa BU, NJC_PF, NJC, null
  "2019-06-30"  # ix & aa NJC_PF, NJC, null
)

ixodes_scapularis.plots <- c("BLAN_012","BLAN_005","SCBI_013","SCBI_002","SERC_001",
                  "SERC_005","SERC_006","SERC_012","ORNL_007")


amblyomma_americanum.plots <- c("SCBI_013","SERC_001","SERC_005","SERC_006","SERC_002",
                    "SERC_012","KONZ_025","UKFS_001","UKFS_004","UKFS_003",
                    "ORNL_002","ORNL_040","ORNL_008","ORNL_007","ORNL_009",
                    "ORNL_003","TALL_001","TALL_008","TALL_002")

both.plots <- intersect(ixodes_scapularis.plots, amblyomma_americanum.plots)

time_series <- function(combined, start.time, species2plot, team.colors, teams = "all"){
  
  df <- combined %>% 
    filter(forecast_start_time == lubridate::as_date(start.time),
           target == species2plot)
  
  if(species2plot == "amblyomma_americanum"){
    df <- df %>%
      filter(plotID %in% c("KONZ_025", "ORNL_007", "ORNL_008",
                           "SCBI_013", "SERC_005", "SERC_006",
                           "TALL_008", "UKFS_003", "UKFS_004"))
    main <- "A. americanum Nymph Forecasts"
  } else {
    main <- "I. scapularis Nymph Forecasts"
  }
  
  if(teams != "all"){
    df <- df %>% 
      filter(Team %in% teams)
  }
  
  teams <- pull(df, Team) %>% unique()
  col <- team.colors[teams]
  
  g <- df %>% 
    mutate(time = as.POSIXct(time)) %>% 
    ggplot(aes(x = time, color = Team)) +
    # geom_line(aes(y = mean), linetype = "dashed") +
    geom_ribbon(aes(x = time, ymin = lower95, ymax = upper95, fill = Team), alpha = 0.4) +
    geom_point(aes(y = actual), color = "black", shape = 20) +
    facet_wrap(vars(plotID)) +
    # facet_wrap(vars(plotID), scales = "free_y") +
    # scale_x_datetime(labels = scales::date_format("%b")) +
    theme_bw() +
    scale_fill_manual(values = col) +
    scale_color_manual(values = col) +
    labs(title = main,
         x = "2019",
         y = "Nymphs per dragged area") 
  
  return(g)
}

time_series_crps <- function(combined, start.time, species2plot, team.colors){
  df <- combined %>% 
    filter(forecast_start_time == lubridate::as_date(start.time),
           target == species2plot,
           time <= ymd("2019-08-01"))
  
  if(species2plot == "amblyomma_americanum"){
    df <- df %>%
      filter(plotID %in% c("KONZ_025", "ORNL_007", "ORNL_008",
                           "SCBI_013", "SERC_005", "SERC_006",
                           "TALL_008", "UKFS_003", "UKFS_004"))
    
    main <- "A. americanum Forecast scores"
  } else {
    main <- "I. scapularis Forecast scores"
  }
  
  teams <- pull(df, Team) %>% unique()
  col <- team.colors[teams]
  
  g <- df %>% 
    mutate(time = as.POSIXct(time)) %>% 
    filter(!is.na(score)) %>% 
    ggplot(aes(x = time)) +
    geom_jitter(aes(y = score, color = Team), alpha = 1, size = 2) +
    facet_wrap(vars(plotID)) +
    # facet_wrap(vars(plotID), scales = "free_y") +
    scale_color_manual(values = col) +
    theme_bw() +
    # scale_x_datetime(labels = scales::date_format("%b")) +
    labs(title = main,
         x = "2019",
         y = "CRPS") 
  
  return(g)
}

crps_horizon <- function(combined, species2plot, team.colors){
  df <- combined %>% 
    filter(target == species2plot,
           !is.na(score))
  
  if(species2plot == "amblyomma_americanum"){
    df <- df %>%
      filter(plotID %in% c("KONZ_025", "ORNL_007", "ORNL_008",
                           "SCBI_013", "SERC_005", "SERC_006",
                           "TALL_008", "UKFS_003", "UKFS_004"))
    
    main <- "A. americanum Forecast scores"
  } else {
    main <- "I. scapularis Forecast scores"
  }
  
  teams <- pull(df, Team) %>% unique()
  col <- team.colors[teams]
  
  g <- df  %>% 
    ggplot(aes(x = horizon)) +
    geom_point(aes(y = score, color = Team), alpha = 0.8) +
    facet_wrap(vars(plotID), scales = "free_y") +
    geom_smooth(
      data = . %>% 
        group_by(Team, plotID),
      mapping = aes(y = score, color = Team),
      method = "lm",
      se = FALSE,
      size = 0.5
    ) +
    scale_color_manual(values = col) +
    theme_bw() +
    labs(title = main,
         x = "Forecast Horizon (weeks)",
         y = "CRPS")  
  return(g)
}
# crps_horizon(combined, targets[spp], team.colors)

# theme for power point, all sizes adjust from axis size
theme_ppt <- function(axis.size){
  theme(title = element_text(size = axis.size + 10), # fig title
        axis.title = element_text(size = axis.size + 6), # main axis
        axis.text = element_text(size = axis.size), # facet axis
        legend.text = element_text(size = axis.size + 2), # legend
        strip.text = element_text(size = axis.size + 2)) # facet labels
}


team.colors <- c(
  "BU_Dem" = "#a6cee3",
  "NJC_Ticks" = "#1f78b4",
  "NJC_ETS_PF" = "#b2df8a",
  "Random Walk" = "#33a02c",
  "SynergisTicks" = "#fb9a99",
  "Ticks_288" = "#e31a1c",
  "VTicks" = "#fdbf6f"
)

axis.size <- 14
targets <- c(ix, aa)
start.dates <- pull(combined, forecast_start_time) %>% unique()

for(spp in seq_along(targets)){
  
  # scores by horizon, each species gets own plot
  g <- crps_horizon(combined, targets[spp], team.colors) +
    theme_ppt(axis.size)
  
  ggsave(filename = paste0(targets[spp], "_ScoresByHorizon"),
         plot = g, 
         device = "tiff", 
         path = "Figures")
  
  for(t in seq_along(start.dates)){
    
    axis.size <- ifelse(t >= 7, 12, 14)
    
    # forecast time series with 95% CI, each start.date_species gets own plot
    g <- time_series(combined, start.dates[t], targets[spp], team.colors, "all") +
      theme_ppt(axis.size) +
      coord_cartesian(ylim = c(0, 40))
    
    ggsave(filename = paste0(targets[spp], "_TimeSeriesFIXED_", start.dates[t]),
           plot = g, 
           device = "tiff", 
           path = "Figures")
    
    # scores through time, each start.date_species gets own plot
    # skip last start date as there are no scores
    if(t == length(start.dates)) next
    g <- time_series_crps(combined, start.dates[t], targets[spp], team.colors) +
      theme_ppt(axis.size)
    
    ggsave(filename = paste0(targets[spp], "_ScoresByTimeFIXED_", start.dates[t]),
           plot = g, 
           device = "tiff", 
           path = "Figures")
  }
}

# scores by team for both species
spp.colors <- c(
"I. scapularis" = "#fc8d59",
"A. americanum" = "#91bfdb"
)

g <- combined %>% 
  rename(Target = target) %>% 
  mutate(across(Target, ~ if_else(.x == "ixodes_scapularis", "I. scapularis", .x)),
         across(Target, ~ if_else(.x == "amblyomma_americanum", "A. americanum", .x))) %>% 
  filter(!is.na(score),
         time <= "2019-08-01") %>% 
  ggplot(aes(x = time, y = score)) +
  # geom_path(group = "Target") +
  geom_point(aes(y = score, color = Target), alpha = 0.75) +
  geom_hline(
    data = . %>% 
      group_by(Team, Target) %>% 
      summarise(mu = mean(score)),
    mapping = aes(yintercept = mu, color = Target)
  ) +
  scale_color_manual(values = spp.colors) +
  facet_wrap(vars(Team), scales = "free_y") +
  theme_bw() +
  theme_ppt(axis.size) +
  # theme(legend.position = "none") +
  labs(title = "Forecast Scores By Team",
       x = "2019",
       y = "CRPS")

g <- ggsave(filename = "ScoresByTeam_BothSpecies_withLegend",
       plot = g, 
       device = "tiff", 
       path = "Figures")


  
  
  






