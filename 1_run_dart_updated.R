#####################################
## R implementation of DART 2.0
## polygon edition
####################################

rm(list=ls())		# start with a clean workspace

# Set your working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

stime <-Sys.time()
source('dart_functions_updated.R')
source('params_updated.R')

## Load packages
required.packages <- c("raster", "rgdal", "rgeos", "sf", "stats", "gower", 'data.table', 'terra')
lapply(required.packages, require, character.only=T)
rm(required.packages)

## Import polygons
polygons <- st_read(dpar$polygons)
## Get vector of polygon IDs to loop through
pids <- as.numeric(polygons$polyID)

## Create an output folder named 'outdir' to store the DART results
dir.create(dpar$outdir, recursive = TRUE)
dir.create(file.path(dpar$outdir,'metadata'), recursive = TRUE)
file.copy('dart_functions_updated.R', file.path(dpar$outdir,'metadata')) 
save(dpar, file = file.path(dpar$outdir, 'metadata','dpar.RData'))

## Run DART using lapply
lapply(pids, function(id){try(dart_fn(id))})
