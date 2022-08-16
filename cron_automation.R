library(cronR)

home_dir <- "/home/rstudio/"
log_dir <- "/home/rstudio/log/cron"

repo <- "neon4cast-ticks"

#Go to healthchecks.io. Create a project.  Add a check. Copy the url and add here.  
health_checks_url <- "https://hc-ping.com/9fd0a2f4-13b3-4d11-880b-1487a1d801ca"

cmd <- cronR::cron_rscript(rscript = file.path(home_dir, repo,"03_nullFitAndForecast.R"),
                           rscript_log = file.path(log_dir, "ticks-null.log"),
                           log_append = FALSE,
                           workdir = file.path(home_dir, repo),
                           trailing_arg = paste0("curl -fsS -m 10 --retry 5 -o /dev/null ", health_checks_url))
cronR::cron_add(command = cmd, frequency = "0 13 1 * *", id = 'ticks-null')