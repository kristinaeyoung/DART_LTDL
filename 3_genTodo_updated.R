###################################################################
#
# Generate a list of pixel x polygons for which to run synthetic control
#
# Input: name of response variable to use
# Output: todo.csv file
#####################################################################
rm(list = ls())

# Set your working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(dplyr)
library(tidyverse)

## Choose the var name(s) that you want to use
varnames <- "afg"
# varnames <- c("afg", "bgr", "pfg", "shr")

useCovs = FALSE # use other variables in synthetic control prediction?

library(data.table)
library(RSQLite)

source('params_updated.R')

dbFile <- dpar$dbfile

CON <- dbConnect(RSQLite::SQLite(), dbFile)

qry <- "   SELECT
               file, pixelID 
             FROM 
               meta 
           "

d <- data.table( dbGetQuery(CON, qry))

out <- list()
for( var in varnames){
  k <- d
  k$varname <- var
  k$useCovs <- useCovs
  out[[var]] <- k
}
# 
D <- rbindlist(out) %>% as_tibble()
D

write.csv(D, paste0(getwd(), '/outdir/todo.csv'), row.names=F)

cat(' there are ', nrow(D), ' sythetic controls to run \n')
cat(' assuming 10 seconds per run, total time would be ', nrow(D) *10 / (60*60), ' hours\n')
cat(' for a job lasting 1 hour, there should be ', nrow(D) / (60*60/10) , ' jobs \n')
cat(' assuming 5 seconds per run, there would be', nrow(D) / (60*60/5), ' total jobs\n')

dbDisconnect(CON)
cat('finished\n')
