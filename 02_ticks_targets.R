# Download and clean tick data 


### Script to create tidy dataset of NEON tick abundances 
### To be used for the RCN Tick Forecasting Challenge 
# Tick Data
# link data from lab and field
# filter out poor quality data
# only included training data 
# through end of 2020
# A. amercanum nymphs
# standardize nymphs across forested plots
# export a .csv

### Load packages
# devtools::install_github("NEONScience/NEON-utilities/neonUtilities")
# library(neonUtilities) # for downloading data (use GitHub version)

# library(usethis)

library(tidyverse) # for data wrangling and piping (dplyr probably ok)
library(lubridate) # for finding year from dates
library(stringr) # for searching within character strings 
library(here) # for working within subdirectories
library(parallel) # for using more than one core in download
library(uuid) # for unique IDs
library(MMWRweek) # for converting from date to MMWR week

# select target species and life stage
target.species <- "Amblyomma americanum"
# target.species <- "Ixodes scapularis"
target.lifestage <- "Nymph"

# subset to the sites we know there to be "good" observations
if(target.species == "Ixodes scapularis"){
  target.sites <- c("BLAN","HARV","ORNL","SCBI","SERC","TREE")  
} else if(target.species == "Amblyomma americanum"){
  target.sites <- c("BLAN","KONZ","LENO","ORNL","OSBS","SCBI","SERC","TALL","UKFS")  
}

if(!"neonstore" %in% installed.packages()){
  library(remotes)
  remotes::install_github("cboettig/neonstore", ref = "patch/api-updates")
}

# Sys.setenv("NEONSTORE_HOME" = "/efi_neon_challenge/neonstore")
# Sys.setenv("NEONSTORE_DB" = "/efi_neon_challenge/neonstore")
library(neonstore)

efi_server <- TRUE

# get data from neon
product <- "DP1.10093.001"
sites <- c("BLAN","HARV","KONZ","LENO","ORNL","OSBS","SCBI","SERC","STEI","TALL","TREE","UKFS")
neon_download(product = product,
              site = sites)

tick.field.raw <- neon_read("tck_fielddata-basic", keep_filename = TRUE)
tick.taxon.raw <- neon_read("tck_taxonomyProcessed-basic")

# there are lots of reasons why sampling didn't occur (logistics, too wet, too cold, etc.)
# so, keep records when sampling occurred
tick.field <- tick.field.raw %>% 
  filter(totalSampledArea > 0) %>% 
  mutate(time = floor_date(collectDate, unit = "day")) %>% 
  unite(namedLocation, time, col = "occasionID", sep = "_")


tick.taxon.wide <- tick.taxon.raw %>% 
  filter(sampleCondition == "OK") %>% # remove taxonomy samples with quality issues
  mutate(sexOrAge = if_else(sexOrAge == "Female" | sexOrAge == "Male", 
                            "Adult",     # convert to Adult
                            sexOrAge),
         time = floor_date(collectDate, unit = "day")) %>% 
  unite(namedLocation, time, col = "occasionID", sep = "_") %>% 
  pivot_wider(id_cols = occasionID, # make wide by species and life stage
              names_from = c(acceptedTaxonID, sexOrAge),
              values_from = individualCount, 
              names_sep = "_",
              values_fn = {sum}, # duplicates occur because of Adults that where F/M - add them 
              values_fill = 0)

tick.joined <- left_join(tick.taxon.wide, tick.field, by = "occasionID") %>% 
  select(-NA_NA, -geodeticDatum, -samplingImpractical, -targetTaxaPresent,
         -adultCount, -nymphCount, -larvaCount, -samplingProtocolVersion, 
         -measuredBy, -sampleCode, -biophysicalCriteria, -plotType)

# all the species column names
spp.cols <- tick.joined %>% 
  select(contains("Larva"), contains("Nymph"), contains("Adult")) %>% 
  colnames()

# get matching taxon ids
taxon.ids <- tick.taxon.raw %>%
  filter(!is.na(acceptedTaxonID)) %>% 
  select(acceptedTaxonID, scientificName, taxonRank) %>% 
  distinct() 

# make longer
tick.long <- tick.joined %>% 
  pivot_longer(cols = all_of(spp.cols), 
               names_to = "taxonAge",
               values_to = "processedCount",
               values_drop_na = TRUE) %>% 
  separate(col = taxonAge, into = c("acceptedTaxonID", "lifeStage"), sep = "_")

# add taxon ids
tick.long <- left_join(tick.long, taxon.ids, by = "acceptedTaxonID") 

tick.out <- tick.long %>% 
  filter(lifeStage == target.lifestage, 
         scientificName == target.species,
         grepl("Forest", nlcdClass)) %>% 
  mutate(time = floor_date(collectDate, unit = "day"),
         time = ymd(time)) %>% 
  select(time, processedCount, totalSampledArea, scientificName, siteID, samplingMethod)

# table of sites, get lat lon
neon.sites <- neon_sites() 
neon.sites.ll <- neon.sites %>% 
  filter(siteCode %in% sites) %>% 
  select(siteCode, siteLatitude, siteLongitude) %>% 
  rename("siteID" = siteCode)

tick.standard <- tick.out %>%
  mutate(totalSampledArea = as.numeric(totalSampledArea)) %>% 
  group_by(siteID, time, scientificName) %>%
  summarise(totalCount = sum(processedCount),
            totalArea = sum(totalSampledArea),
            numPlots = n(),
            standardCount = totalCount / totalArea * 1600) %>% 
  mutate(year = year(time),
         mmwrWeek = MMWRweek(time)$MMWRweek,
         time = MMWRweek2Date(year, mmwrWeek)) %>% 
  rename("Amblyomma americanum" = standardCount) %>% 
  arrange(siteID, time) %>% 
  filter(siteID %in% target.sites) %>% 
  select(time, mmwrWeek, siteID, `Amblyomma americanum`)

# in case NEON makes a provisional data release during the challenge
# need to make sure that we filter out 2021 data that is in the "future"
run.date <- today()
year(run.date) <- year(run.date) - 1
run.mmwrWeek <- MMWRweek(run.date)$MMWRweek
challenge.time <- MMWRweek2Date(year(run.date), run.mmwrWeek)

tick.targets <- tick.standard %>% 
  filter(time < challenge.time)

# write targets to csv
write_csv(tick.targets,
          file = "ticks-targets.csv.gz")

if(efi_server){
  source("../neon4cast-shared-utilities/publish.R")
  publish(code = c("02_ticks_targets.R"),
          data_out = c("ticks-targets.csv.gz"),
          prefix = "ticks/",
          bucket = "targets",
          registries = "https://hash-archive.carlboettiger.info")
}


