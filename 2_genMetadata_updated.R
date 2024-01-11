#
#
#  R function for populating the time-invariant data fields for each pixel
#  
#  Input: chunk #, total chunks
#  Output: writes data from chunk to SQLite database
#
###########################################################################

# Start with a clean workspace
rm(list=ls())

# Set your working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source('params_updated.R')
library(terra)
library(sf)
library(data.table)
library(RSQLite)

## raster covariates ##
## Add in other soil/climate/etc variables (raster format) that you want to use to give more context to your polygon data.

sgu <- rast('V:/PROJECTS/GAYLE_TYREE/DART/Data_GIS/Covariates/sgu_1stClass.tif')
esg <- rast('V:/PROJECTS/GAYLE_TYREE/DART/Data_GIS/Covariates/ESGs_final.tif')

arg <- commandArgs(TRUE)
chunk <- as.numeric( arg[1] )
nchunks <- as.numeric( arg[2] )
vari = arg[3] #DEPRECATED
useCovs = TRUE

ff <- list.files(dpar$outdir, full = TRUE, pattern = 'RData')
ids <- as.numeric(1:length(ff))

OUTS <- list()

for( i in ids ){
  
  f <- ff[i]
  
  try({
    
    cat(f, '\n')
    load(f)
    
    ## get metadata  ##
    meta <- padpixels
    dontKeep <- c("refrast", "othersites", "exclosures")
    meta <- meta[ , which(!names(meta) %in% dontKeep)]
    
    xy <- st_coordinates(st_centroid(padpixels))
    meta$x <- xy[,1]
    meta$y <- xy[,2]
    
    # Split out pixel-level data vs treatment level data
    meta$polyID <- padpoly$polyID
    # meta$yearStart <- padpoly$DOI    # Intervention time (date format)
    meta$yearStart <- padpoly$YOI     # Intervention time (year format)
    meta$pixelID <- 1:nrow(meta)
    meta$file <- f
    
    
    for(p in dpar$metaCols){
      try({
        meta[[p]] <- padpoly[[p]]
      })
    }
    
    ######################################################
    ######################################################  
    
    ## extract additional covariates
    
    #Soil Geomorphic units
    meta$SGU <- terra::extract(sgu, xy)[,1]
    
    #Ecological Site Groups
    meta$ESG <- terra::extract(esg, xy)[,1]
    
    ## Feel free to add in what you want!
    
    ####################################################
    ####################################################
    
    # get pre-treatment averages per pixel for all response variables
    for( variable in names(dpar$respVars)){
      
      e <- extraction$extractedTarget[[variable]]
      
      # identify years since treatment
      yy <- (2021 - ncol(e)):2020      # aka (currentYear - number of years of data in e) to previous year       
      library(lubridate)
      # y <- 1:length(yy) - grep(as.integer(year(padpoly$DOI)), yy)    # Use for date format of intervention time
      y <- 1:length(yy) - grep(as.integer(padpoly$YOI), yy)    # Use for year format of intervention time
      
      
      meta[[ paste0('preMean.',variable)]] <- rowMeans(e[, y<0], na.rm = TRUE)
      meta[[ paste0('preSD.', variable) ]] <- apply(e[, y<0], 1, sd, na.rm = TRUE)
      meta[[ paste0('preMed.', variable)]] <- apply(e[, y<0], 1, median, na.rm = TRUE)
    }
    
    
    OUTS[[as.character(i)]] <- meta
    
  })
}

OUTS
# OUTS[[1]]
M <- do.call(rbind, OUTS)
M <- M %>% st_drop_geometry()


# connect to db
dbFile <- dpar$dbfile
CON <- dbConnect(RSQLite::SQLite(), dbFile)

# write to tables
library(DBI)
cat('writing table "meta" in ', dbFile, '\n')
retry <- T
while( retry ){
  
  ## dbWriteTable wants me to drop the geometry from the sf - if it turns out that we need the geom later, then save it in a new col as a character vector first.
  res <- try({ dbWriteTable(CON, 'meta', value = M, overwrite = F, append = TRUE, row.names=F) } )
  if( class(res) != 'try-error' ) {
    retry <- FALSE
    cat('   success\n')
  } else {
    cat('    write try failed\n')
    Sys.sleep(5)
  }
}

# Doublecheck that the table wrote to DB
as.data.frame(dbListTables(CON))

dbDisconnect(CON)
cat('finished\n')