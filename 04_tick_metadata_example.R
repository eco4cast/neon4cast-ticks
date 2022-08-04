#### create forecast metadata ####
# this example is for the ensemble csv forecast

library(EML)
library(tidyverse)
library(lubridate)

fx.dir <- "Random_Walk_Fits/2019-10/" # where forecast is stored
fx.file <- "ticks-2019-03-04-tickGlobalNull_RandomWalk" # file name
file.dest <- file.path(fx.dir, paste0(fx.file, ".csv.gz"))
fx <- read.csv(file.dest)

# pull out time
time <- fx %>% 
  pull(time) %>% 
  ymd() %>% 
  unique() %>% 
  sort()

# get number of ensembles
n.ens <- fx %>% 
  pull(ensemble) %>% 
  max()

# get lat.lon
ticks <- read_csv("https://data.ecoforecast.org/neon4cast-targets/ticks/ticks-targets.csv.gz")
lat.lon <- ticks %>% 
  select(c(siteID, decimalLatitude, decimalLongitude)) %>% 
  mutate(decimalLatitude = decimalLatitude, 
         decimalLongitude = decimalLongitude) %>% 
  distinct()

forecast_project_id <- "tickGlobalNull_RandomWalk" # project or team name
forecast_model_id <- "7e9c476fb483e7a77b4286db6bee34542aea6626" # last commit from model development
forecast_iteration_id <- time[1] # forecasts are submitted monthly, iteration id is the first time in each forecast
forecast_issue_time <- "2021-03-30" # the day the forecast was run

attributes <- tibble::tribble(
  ~attributeName,           ~attributeDefinition,                          ~unit,                  ~formatString, ~numberType, ~definition,
  "time",                   "[dimension]{time}",                          "year",                 "YYYY-MM-DD",  "numberType", NA,
  "siteID",                 "[dimension]{location label}",                "dimensionless",         "SSSS",       "character",  NA,
  "plotID",                 "[dimension]{location label}",                "dimensionless",         "SSSS_XXX",   "character",  NA,
  "ensemble",               "[dimension]{index of ensemble member}",      "dimensionless",         NA,           "integer",    NA,
  "obs_flag",               "[dimension]{observation error}",             "dimensionless",         NA,           "integer",    NA,
  "ixodes_scapularis",      "[variable]{Pop. abundace of I. scap.}",      "number",                NA,           "integer",    NA,
  "amblyomma_americanum",   "[variable]{Pop. abundace of A. amer.}",      "number",                NA,           "integer",    NA,
  "forecast",               "[flag]{whether time step assimilated data}", "dimensionless",         NA,           "integer",    NA,
  "data_assimilation",      "[flag]{whether time step assimilated data}", "dimensionless",         NA,           "integer",    NA
) 

## note: EML uses a different unit standard than UDUNITS. For now use EML. EFI needs to provide a custom unitList.
attributes
attrList <- set_attributes(attributes, 
                           col_classes = c("Date", "character", "character", "numeric",
                                           "numeric","numeric", "numeric","numeric", "numeric"))

## sets metadata about the file itself (name, file type, size, MD5, etc)
physical <- set_physical(file.dest,
                         recordDelimiter='\n')

## set metadata for the file as a whole
dataTable <- eml$dataTable(
  entityName = "forecast",  ## this is a standard name to allow us to distinguish this entity from 
  entityDescription = "Forecast of nymphal tick population using a state-space random walk model",
  physical = physical,
  attributeList = attrList)

me <- list(individualName = list(givenName = "John", 
                                 surName = "Foster"),
           electronicMailAddress = "fosterj@bu.edu")

taxa <- tibble::tribble(
  ~Genus,      ~Species,
  "Ixodes",    "scapularis",
  "Amblyomma", "americanum")

coverage <- 
  set_coverage(begin = first(time), 
               end = last(time),
               sci_names = taxa,
               geographicDescription = "NEON sites BLAN, ORNL, SCBI, SERC, KONZ, TALL, UKFS",
               west = min(lat.lon$decimalLongitude), 
               east = max(lat.lon$decimalLongitude), 
               north = max(lat.lon$decimalLatitude), 
               south = min(lat.lon$decimalLatitude))

keywordSet <- list(
  list(
    keywordThesaurus = "EFI controlled vocabulary",
    keyword = list("forecast",
                   "population",
                   "timeseries",
                   "ticks")
  ))

abstract <- "Random walk state-spave model. Process error is Gaussian and observation error is Poisson.
Each site by species target is fit independently."

dataset = eml$dataset(
  title = "Random Walk Null Model",
  creator = me,
  contact = list(individualName = list(givenName = "John", 
                                       surName = "Foster")),
  pubDate = forecast_issue_time,
  intellectualRights = "http://www.lternet.edu/data/netpolicy.html.",
  abstract = abstract,
  dataTable = dataTable,
  keywordSet = keywordSet,
  coverage = coverage
)

additionalMetadata <- eml$additionalMetadata(
  metadata = list(
    forecast = list(
      ## Basic elements
      timestep = "1 week", ## should be udunits parsable; already in coverage -> temporalCoverage?
      forecast_horizon = paste0(length(time), " weeks"),
      forecast_issue_time = forecast_issue_time,
      forecast_iteration_id = forecast_iteration_id,
      forecast_project_id = forecast_project_id,
      metadata_standard_version = "0.3",
      model_description = list(
        forecast_model_id = forecast_model_id,
        name = "discrete random walk",
        type = "process-based",
        repository = "https://github.com/eco4cast/neon4cast-ticks.git"
      ),
      ## MODEL STRUCTURE & UNCERTAINTY CLASSES
      initial_conditions = list(
        # Possible values: absent, present, data_driven, propagates, assimilates
        status = "propagates",
        # Number of parameters / dimensionality
        complexity = 28  ## one process error term per target (species x plot) 
      ),
      drivers = list(
        status = "absent" # no drivers in null model
      ),
      parameters = list(
        status = "absent", # no parameters in null model
      ),
      random_effects = list(
        status = "absent" # no random effects in null model
      ),
      process_error = list(
        status = "propagates",
        propagation = list(
          type = "ensemble", # ensemble vs analytic
          size = n.ens          # required if ensemble
        ),
        complexity = 1,   
        covariance = FALSE
      ),
      obs_error = list(
        status = "absent"
      )
    ) # forecast
  ) # metadata
) # eml$additionalMetadata

my_eml <- eml$eml(dataset = dataset,
                  additionalMetadata = additionalMetadata,
                  packageId = forecast_iteration_id , 
                  system = "datetime"  ## system used to generate packageId
)

eml.file <- file.path(fx.dir, paste0(fx.file, "-eml.xml"))
write_eml(my_eml, eml.file)

## check base EML
eml_validate(my_eml)

## check that the EML is also a valid EFI forecast
EFIstandards::forecast_validator(metadata)



