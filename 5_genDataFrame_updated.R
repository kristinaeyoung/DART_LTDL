## Script to organize the DART outcomes and Synthetic Control outcomes into two cleaner data frames.

rm(list=ls())

# Set your working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(data.table)
library(RSQLite)

## pull out data from DART analysis from sqlite db
## Input: variable to extract
## Output: RData Long data.table with 1 row per pixel-timepoint 

source('params_updated.R')

datafile =  file.path(dpar$outdir, 'metadata', paste0('all.RData'))

# extract actual values 
dbFile <- file.path(dpar$outdir,'..', paste0( basename(dpar$outdir), '.sqlite') )
dbFile <- dpar$dbfile

CON <- dbConnect(RSQLite::SQLite(), dbFile)


tables <- dbGetQuery(CON, "select name from sqlite_master WHERE type='table'")

dat <- list()

for( item in tables$name){
  
  qry <- paste0( " SELECT * FROM ", item, ";")
  dat[[item]] <- data.table( dbGetQuery( CON, qry))
  
}

#n


##   qry <- paste0(                      

##       "   SELECT                      
##              * 
##            FROM ",
##              var, " as res
##              -- npp as res
##            INNER JOIN
##              meta
##            ON res.id = meta.polyID 
##              AND
##               res.pixelID = meta.pixelID
##            WHERE 
##               res.years_after > 0
##            "
##        )

## d <- data.table( dbGetQuery(CON, qry))


save(dat, file = datafile)

as.data.frame(dbListTables(CON))

dbDisconnect(CON)

## Look at the data oucomes
getwd()
load('outdir/metadata/all.RData')

# DART pixel metadata
head(dat$meta)

# DART quantiles and Synthetic Control results
head(dat$sc)
