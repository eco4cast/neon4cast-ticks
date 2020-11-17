#' Create a netCDF file for forecast output 
#'
#' This function takes an array of forecast output and writes a netcdf file that adheres to the EFI Forecast Standards.
#'
#' @param fx Array of forecasts. Dimensions must be number of ensembles x time x obs.dim x number of targets. (see obs.dim below)
#' @param ncfname The name of the nc file you wish to create
#' @param data.assimilation Vector that flags which time points have observation data that was assimilated. Length = n.weeks.
#' @param obs.dim Is the forecast for the latent state (obs.dim = 1) or the latent state and latent state with observation uncertainty (obs.dim = 2)
#' @param forecast Vector indicating a forecast or hindcast for each time point. Length = n.weeks
#' @param target.name Vector of target names that maps to the ensemble x time x observation of fx. Must include species name and plot name separated by "_" (i.e. "Ixodes_scapularis_BLAN_005")
#' @param forecast.issue.time start date of the forecast YYYYWW
#' @param Forecast.id Represents each forecast cycle within a ForecastProject.id. ISO datetime should make a valid Forecast.id
#' @param ForecastProject.id Some ID that applies to a set of forecasts. Represents the launch of an automated, iterative forecast. It is created each time a human modifies the forecast code
#' @param check.units Check if units are valid? Default TRUE

create_netCDF <- function(fx, ncfname, data.assimilation, 
                          obs.dim, forecast, targetName,
                          forecast.issue.time, Forecast.id, ForecastProject.id, 
                          check.units = TRUE){
  # make sequences
  ensembles <- 1:dim(fx)[1] # number of ensembles (draws from joint posterior)
  times <- seq(from = forecast.issue.time, by = 1, length = dim(fx)[2])

  ### Set dimensions ###
  # size of timestep, with units and start date
  # sequence of values that defines the dimension.
  # netCDF expresses time dimensions relative to a
  # specified start date, in this case forecast.issue.time
  timedim <- ncdim_def("time",
                       units = paste('weeks since', forecast.issue.time),
                       vals = as.numeric(times - forecast.issue.time),
                       longname = 'timestep')                                           
  ensdim <- ncdim_def("ensemble",             
                      units = "",
                      vals = ensembles,       
                      longname = 'ensemble member') 
  obsdim <- ncdim_def("obs_flag",
                      units = "",
                      vals = 1:obs.dim,
                      longname = "observation error flag. 1 = latent state, 2 = w/ obs error")
  
  if(check.units){
    library(udunits2)
    print(paste("Time units valid:", ud.is.parseable(timedim$units)))
    print(paste("Ensemble units valid:", ud.is.parseable(ensdim$units)))
    print(paste("Observation units valid:", ud.is.parseable(obsdim$units)))
  }
  
  ### Define variables ###
  fillvalue <- -1e20 
  def_list <- list()
  for(i in seq_along(targetName)){ 
    def_list[[i]] <- ncvar_def(
      name =  targetName[i],
      units = "number of individuals",
      dim = list(ensdim, timedim, obsdim),
      missval = fillvalue,
      longname = paste0('forecast for ', targetName[i]),
      prec = "single"
    )
  }
  
  def_list[[length(targetName)+1]] <- ncvar_def(
    name =  "forecast",
    units = "integer",
    dim = list(timedim),
    missval = fillvalue,
    longname = 'EFI standard forecast code. 0 = hindcast',
    prec = "single"
  )
  
  def_list[[length(targetName)+2]] <- ncvar_def(
    name =  "data_assimilation",
    units = "integer",
    dim = list(timedim),
    missval = fillvalue,
    longname = 'EFI standard data assimilation code. 0 = no data',
    prec = "single"
  )
  
  # open netCDF file
  ncout <- nc_create(ncfname, def_list, force_v4 = TRUE)

  # put output data into ncdf file
  for(i in seq_along(targetName)){ 
    for(y in 1:obs.dim){
      ncvar_put(ncout, def_list[[i]], fx[,,y,i])        ## forecasts
    }
  }
  
  ncvar_put(ncout, def_list[[length(targetName)+1]], forecast) ## forecast flag
  ncvar_put(ncout, def_list[[length(targetName)+2]], data.assimilation) ## data_assimilation flag
  
  #Global file metadata
  ncatt_put(ncout,0,"ForecastProject_id", as.character(ForecastProject.id),
            prec =  "text")
  ncatt_put(ncout,0,"Forecast_id",as.character(Forecast.id),
            prec =  "text")
  ncatt_put(ncout,0,"forecast_issue_time",as.character(forecast.issue.time),
            prec =  "text")
  
  # close (always remember to close!)
  nc_close(ncout)

}