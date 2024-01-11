# Generate Synthetic Control timeseries for selected pixels
#
# Inputs: chunk #, N chunks
# Ouptuts: timeseries to sqlite db
#

rm(list=ls())

# Set your working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

######################################################
######################################################  
##                     Setup                        ##
######################################################
######################################################

source('dart_functions_updated.R')
source('params_updated.R')
library(data.table)
library(raster)
library(terra)
library(sf)
library(RSQLite)
library(lubridate)

args <- commandArgs(TRUE)
chunk <- as.numeric(args[1])
chunks <- as.numeric(args[2])
extract <- T


cat(args)
todo <- fread('outdir/todo.csv')

w <- as.numeric(seq(1:nrow(todo)))  

cat(' processing pixels ', min(w), ' to ', max(w), '  ', length(w), ' total\n')

tstart <- Sys.time()

######################################################
######################################################  
##            Functions                             ##
######################################################
######################################################


#CausalImpact
CI <- function(s, pre, post, covs = FALSE){
  require(CausalImpact)
  s <- data.table(s)
  
  # note CausalImpact won't work if bfast is loaded :
  #  conflict with as.zoo.data.frame :( https://github.com/joshuaulrich/quantmod/issues/168
  
  X = dcast(s, dayt ~ id, value.var = 'y')
  
  if(covs){
    cv <- setdiff(names(s), c('id', 'dayt', 'time', 'y', 'D'))
    
    for ( V in cv) {
      
      X0 <- dcast(s, dayt ~ id, value.var = V)
      X0 <- X0[,-grepl('trt', names(X0)), with = F]
      names(X0) <- paste0(V, names(X0))
      X <- cbind(X, X0)
    }
    
  }
  
  
  #rearrange
  tt = X$dayt
  ddate <- max( s[ grepl( 'trt', id ) & D == 0, dayt] )
  #pre = c(1, max(which(tt <= ddate)))
  #post = c( pre[2] + 1, length(tt) )
  
  X = X[, c('trt', grep('c', names(X), value = TRUE) ), with = F]
  
  y <- CausalImpact::CausalImpact(as.matrix(X), pre.period = pre, post.period  =post)
  
  list( effect = as.data.table(y$series[, c('point.effect', 'point.effect.upper', 'point.effect.lower',
                                            'cum.effect', 'cum.effect.lower', 'cum.effect.upper')]),  model  = y)
}

GS <- function(oi){
  require(gsynth)
  g <- gsynth(y ~ D, data = na.omit(oi), index = c("id","time"), force = "two-way", CV = TRUE,
              r = c(0, 2), se = FALSE, inference = "parametric", nboots = 200, parallel = FALSE)
  
  input <- sort(unique(oi$time))
  outpt <- rep(NA, length(input))
  ef <- g$att
  outpt[ match(names(ef), input) ] <- ef 
  
  
  cumul <- outpt
  cumul[input < 1] <- 0
  cumul <- cumsum(cumul)
  
  
  list( effect = data.table( 'point.effect' =outpt, 'point.effect.upper' = NA, point.effect.lower= NA,
                             cum.effect = cumul, cum.effect.upper = NA,
                             cum.effect.lower = NA))
}


eval_dart <- function(pix, varname, useCovs = FALSE){
  
  yr0 <- as.numeric(padpoly$yearStart)
  yr1 <- as.numeric(padpoly$yearEnd)
  
  yy <- (2021 - ncol(extraction$extractedTarget[[varname]])):2020
  pre <- c(1, grep(yr0, yy) - 1)
  post <- c(grep(yr1, yy),  length(yy))
  
  o <- longPanel(extraction, yr0, varname = varname, FALSE, FALSE)
  
  ids <- c( paste0('trt', pix), paste0('ctr', toposim$index[,pix]))
  
  oi <- o[which(o$id %in% ids),]
  oi$id <- ifelse(grepl('trt', oi$id), 'trt', oi$id)
  oi$dayt <- oi$time
  oi$time <- oi$dayt - yr0
  oi$time <- ifelse(oi$time < 1, oi$time, pmax(0, oi$dayt - yr1)) 
  
  # generate synthetic control
  ev <- try({ CI( data.table( oi ), pre, post, useCovs )$effect }, silent = TRUE)
  
  if(class(ev)[1] == 'try-error' | varname == 'desi'){     ### hard-coded to DESI - change or is it ok to keep?
    ev <- GS(oi)$effect
  }
  
  # add Dart quantile
  oi <- data.table(oi)
  
  ecSafe <- function(a, b){
    ec = try({ecdf(b)(a)}, silent = TRUE)
    if(class(ec) == 'try-error'){
      ec <- as.numeric(NA)
    }
    ec
  }
  
  
  drt = oi[, list(
    dartQuantile = ecSafe(y[id == 'trt'], y[id != 'trt']),
    dartMedianValue = median(y[id != 'trt'], na.rm = TRUE)
  ),
  by = list(dayt, time)]
  
  ev$dartQuantile <- drt[['dartQuantile']]
  ev$dartMedianValue <- drt[['dartMedianValue']]
  
  # metadata
  ev$pixelID <- pix
  ev$id <- padpoly$polyID
  meta <- oi[grep('trt', oi$id), -c(1,3)]
  ev <- cbind(ev, meta)
  
  # add climate anomaly
  yrs <- paste0('y', ev$dayt)
  
  ev <- data.table(ev)
  setnames(ev, 'time', 'years_after')
  ev <- ev[ years_after >= -1,]
  
  ev
}

######################################################
######################################################  
##        Generate Synthetic Controls               ##
######################################################
######################################################  

out <- list()
currentFile <- ''
times <- c()
for( i in w){
  
  try({
    
    st <- Sys.time()   
    cat( '    ', i, '    \n'); flush.console()
    f <- todo$file[i]
    varname <- todo$varname[i]
    useCovs <- todo$useCovs[i]
    pixel <- todo$pixelID[i]
    
    if (!currentFile == f){
      load(f)
      source('params_updated.R')
      currentFile <- f
      
      # prep climate anomaly
      xy <- st_coordinates(st_centroid(padpixels))
      
      if(extract == 'TRUE'){
        # re-extract the ts before synthetic controls
        cat('extracting ', varname, '\n')
        ras <- list()
        ras[[ varname ]] <- dpar$respVars[[varname]]
        extraction <- extract_TS(ras, padpixels, chosenCandidates, toposim) 
        
      }
    }
    
    eval_dart <- function(pix, varname, useCovs = FALSE){
      
      pix <- pixel 
      yr0 <- as.numeric(year(padpoly$EstblsD))
      yr1 <- as.numeric(year(padpoly$EstblsD))
      
      yy <- (2021 - ncol(extraction$extractedTarget[[varname]])):2020
      pre <- c(1, grep(pattern=yr0, x=yy) - 1)
      post <-c(grep(yr1, yy),  length(yy))
      
      o <- longPanel(extraction, yr0, varname = varname, FALSE, FALSE)
      
      ids <- c( paste0('trt', pix), paste0('ctr', toposim$index[,pix]))
      
      oi <- o[which(o$id %in% ids),]
      oi$id <- ifelse(grepl('trt', oi$id), 'trt', oi$id)
      oi$dayt <- oi$time
      oi$time <- oi$dayt - yr0
      oi$time <- ifelse(oi$time < 1, oi$time, pmax(0, oi$dayt - yr1)) 
      # generate synthetic control
      ev <- try({ CI( data.table( oi ), pre, post, useCovs )$effect }, silent = TRUE)
      
      if(class(ev)[1] == 'try-error' | varname == 'desi'){
        ev <- GS(oi)$effect
      }
      
      # add Dart quantile
      oi <- data.table(oi)
      
      ecSafe <- function(a, b){
        ec = try({ecdf(b)(a)}, silent = TRUE)
        if(class(ec) == 'try-error'){
          ec <- as.numeric(NA)
        }
        ec
      }
      
      
      drt = oi[, list(
        dartQuantile = ecSafe(y[id == 'trt'], y[id != 'trt']),
        dartMedianValue = median(y[id != 'trt'], na.rm = TRUE)
      ),
      by = list(dayt, time)]
      
      ev$dartQuantile <- drt[['dartQuantile']]
      ev$dartMedianValue <- drt[['dartMedianValue']]
      
      # metadata
      ev$pixelID <- pix
      ev$id <- padpoly$polyID
      meta <- oi[grep('trt', oi$id), -c(1,3)]
      ev <- cbind(ev, meta)
      
      # add climate anomaly
      yrs <- paste0('y', ev$dayt)
      
      ev <- data.table(ev)
      setnames(ev, 'time', 'years_after')
      ev <- ev[ years_after >= -1,]
      
      ev
    }
    
    ev <- eval_dart(pixel, varname, useCovs)
    ev$variable = varname
    ev$insertDate <- as.Date(Sys.time())
    setnames(ev, 'y', 'value')
    
    out[[ as.character(i) ]] <- ev
    ed <- Sys.time()
    tdiff <- ed - st
    times <- c(times, tdiff)
    cat(' that took ',tdiff , '\t\t Average: ', mean(times),' \t max: ', max(times),'\n')
    
  })
  # if running out of time, save progress
  #        if( difftime(Sys.time() ,  tstart ,units = 'secs' ) > 1.5 * 60 * 60 * .9 ){
  #          break()
  
  #     } 
}



######################################################
######################################################  
##            Add results to database               ##
######################################################
######################################################  
source('params_updated.R')
dbFile <- dpar$dbfile

response <- do.call(rbind, out)
if(nrow(response) == 0){ stop('No data')}

# Write to CSV as a backup if something goes wrong with saving to SQL
write_csv(response, paste0(getwd(), "/outdir/sc.csv"))

# connect to db
CON <- dbConnect(RSQLite::SQLite(), dbFile)

# write to tables
cat('writing table "sc" in ', dbFile, '\n')
retry <- T       # retry MUST be set to T or the file will not save!!!
# Note: 

while( retry){
  
  res <- try({ dbWriteTable(CON, "sc" , value = response, append = TRUE, overwrite = F) })
  
  if( class(res) != 'try-error'){
    retry <- FALSE
    cat('   success\n')
  }else {
    cat('   write try failed\n') 
    Sys.sleep(5)
  }
  
}

dbDisconnect(CON)
cat('\nfinished\n')


### Check the file
cat('\ncheck the saved file\n')
CON <- dbConnect(RSQLite::SQLite(), dbFile)
dbListTables(CON)
head(dbReadTable(CON, "sc"))
dbDisconnect(CON)
cat('\nfinished checking file\n')