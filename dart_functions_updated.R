## DART for large polygons - updated to use sf and terra
## Original code by Stephen Fick
## Updated by G. Tyree - gtyree@usgs.gov 
## Updated during Sept/Oct 2023

##################################################
## DART function for parallel processing         #
##################################################
#
# 1. Generate dart pixels with GOWER. 
#    Output file : x.dart  -- padpoly,bufrast,toposim
#
# 2. Extract timeseries rasters from polygons
#    Output file : x.extract
#

dart_fn <- function(id) {    
  
  cat(' ** running dart for id ', id, '\n')   
  filename <- file.path( dpar$outdir, paste0(dpar$prefix, id, '.RData'))
  if(file.exists(filename)){ return() }
  
  #internal timer
  tstart <- Sys.time()
  cat(Sys.time(), ' ', id, '\n', file = '~/snowlog.txt', append = TRUE)
  
  # get polygon
  padpoly <- get_poly(dpar, id)
  
  # get raster data
  zone <- getRasterData(padpoly, dpar, dpar$rad)
  
  # get target pixels
  pad <- get_treatment(padpoly, dpar, zone)   
  padpixels <- pad$padpixels
  
  # get table of unique combos of ec and soilps in target region
  uni <- pad$uni
  
  ## Get padstk as an object to use in the while loop below
  padstk <- pad$padstk
  
  # get candidates
  tries <- dpar$tries
  lastTry <- FALSE
  searchRadius <- dpar$rad
  
  
  # Function to find edaphic matches
  threshFunction <- function(candidates, x)  which( candidates$soilps == x[1] & 
                                                      candidates$ec > (x[2] * .95) & 
                                                      candidates$ec < (x[2] * 1.05)
  )
  
  # Placeholder to see which pixels are omitted due to lack of data
  ignoredTreatedPixels <- 0
  
  # Identify candidate pixels
  while( tries > 0){
    
    if ( tries == 1 ) lastTry = TRUE
    cat(' ** looking for candidates\n')
    cat('\t searchRadius:', searchRadius, '\n')       
    
    # get candidates
    candidates <- try({ get_candidates(dpar, searchRadius, padpoly, zone) })
    
    # check to see if any candidates were found
    if('try-error' %in% class(candidates)){
      if(grepl('ERROR1', candidates)){
        cat(' \tno non NA candidates. Doubling search radius\n');
        tries <- tries - 1
        searchRadius <- searchRadius  + dpar$rad
        zone <- getRasterData(padpoly, dpar, searchRadius)
        next()
      } else {
        stop(candidates)
      }
    } 
    
    # check to see if every edaphic class has at least 100 matching pixels
    avail  = apply(st_drop_geometry(uni[,1:2]), 1, function(x) length(threshFunction(candidates, x)) )
    cat(dim(avail),'\n')
    if( any(avail < 100 ) ) {
      
      # if not last try, fail and expand search radius
      if( !lastTry){
        cat( ' \t not enough edaphic candidates for some treated pixels \n');
        
        cat( '\t EC candidates: ' , quantile( candidates$ec, na.rm=TRUE),'\n')
        cat( '\t soilps candidates: ' , quantile( candidates$soilps, na.rm = TRUE),'\n')
        
        for(w in which(avail < 100)){
          cat('\t',unlist(uni[w,]), ': avail = ', length(avail[w]),'\n')
        }
        tries <- tries - 1
        searchRadius <- searchRadius  + dpar$rad
        zone <- getRasterData(padpoly, dpar, searchRadius)
        next()
      }
      
      # if last try, ignore pixels without enough edaphic support
      if ( lastTry) {
        
        toRemove <- uni[which(avail < 100),]$uni
        w <- which(padpixels$uni %in% toRemove)
        ignoredTreatedPixels <- w
        if(length(w) == nrow(padpixels)){
          stop('ERROR: No edaphic support for any treated pixel')
        }
        cat('\t removing ', length(w),' treated pixels due to lack of support\n')
        padpixels <- padpixels[-w,]
        uni <- uni[-which(uni$uni %in% toRemove),]
      }
    }
    cat( '\t finished search\n')      
    break() 
    
  }
  
  
  # Identify factors for classification
  candidates$LFELEMS <- as.character(candidates$LFELEMS)
  padpixels$LFELEMS <- as.character(padpixels$LFELEMS)
  
  # get similarity matrix for top 100 ref pixels for each target pixel
  cat('getting similarity index \n')
  toposim <- getGowerByGroup(uni, padpixels, candidates, threshFunction)
  chosenCandidates <- candidates[unique(c(toposim$index)),]
  
  
  ## Extract response variables
  if (length(dpar$respVars) > 0){
    
    extraction <- extract_TS(dpar$respVars, padpixels, chosenCandidates, toposim)
    
  } else {
    
    extraction <- 'NO Extraction performed because no response included in parameters object'
    
  }
  
  timeElapsed <- Sys.time() - tstart
  cat('saving data to ', filename, '\n')
  
  candidates <- candidates[,1]
  
  save( 
    padpoly, padpixels, searchRadius, candidates, 
    chosenCandidates, toposim, 
    timeElapsed, dpar, extraction, ignoredTreatedPixels, 
    file = filename) 
  
  # generate metadata and log to SQL
  cat('saving metadata\n')
  meta <- gen_metadata()
  insert_table(meta, 'meta')
  
  # save raw response timeseries
  cat('saving timeseries\n')
  ts <- gen_ts()
  ts$insertDate <- as.Date(Sys.time())
  insert_table(ts, 'ts')
  
  cat(' time elapsed: ', timeElapsed, '\n')
  cat(' number of treatment pixels: ', nrow(padpixels), '\n')
  cat(' number of pixels removed: ', length(ignoredTreatedPixels), '\n')
  cat(' final search radius: ', searchRadius, '\n')
  cat(' #COMPLETED# \n')    
  gc()
}

## Supporting functions ###
gen_ts <- function( env = parent.frame() ){
  
  # input : extraction
  # output : long table - polyID, pixelID,years_after, year, variable, value
  extraction <- get('extraction', env)
  padpoly <- get('padpoly', env)
  
  D = list()
  
  for( variable in names(extraction$extractedTarget)){
    
    e <- extraction$extractedTarget[[variable]]
    
    # identify years since treatment
    yy <- (2021 - ncol(e)):2020      # aka (currentYear - number of years of data in e) to previous year       
    library(lubridate)
    years_after <- 1:length(yy) - grep(as.integer(year(padpoly$EstblsD)), yy)
    
    
    D[[variable]] <- data.table(polyID = padpoly$polyID,
                                pixelID = rep( 1:nrow(e), ncol(e)),
                                years_after = rep(years_after, each = nrow(e)),
                                year = rep(yy, each = nrow(e)),
                                variable = variable,
                                value = c(e)) 
  }
  
  d.out <- do.call(rbind, D)
  
  d.out
  
  
}

gen_metadata <- function(env = parent.frame()){
  
  # format metadata
  
  require(data.table)
  
  ## load env variables
  dpar <- get('dpar', env)
  padpoly <- get('padpoly', env)
  extraction <- get('extraction', env)
  padpixels <- get('padpixels', env)
  filename <- get('filename', env)
  
  ## get metadata  ##
  meta <- padpixels
  
  # remove columns
  meta <- meta[ , which(!names(meta) %in% dpar$dontKeep)]
  
  # add coordinates
  xy <- st_coordinates(st_centroid(padpixels) )
  meta$x <- xy[,1]
  meta$y <- xy[,2]
  meta$pixelID <- 1:nrow(meta)
  meta$file <- filename
  
  # add metadata
  for(p in dpar$metaCols){
    try({
      meta[[p]] <- padpoly[[p]]
    })
  }
  
  #################################
  #################################
  
  ## static covariates
  cat('extracting static variables\n')
  for( v in names(dpar$staticVars)){
    cat('\t',v,'\n')
    R <- rast(dpar$staticVars[[v]])
    meta[[v]] <- terra::extract(R, padpixels)[,1]
  }
  
  
  ################################
  ################################
  
  ## get pre-treatment mean and sd for each extraction
  
  for( variable in names(extraction$extractedTarget)){
    
    variable <- "pfg"
    e <- extraction$extractedTarget[[variable]]
    
    # identify years since treatment
    yy <- (2019 - ncol(e)):2018    
    ## e could be an extraction from a multi-layer raster with nlyrs corresponding to number of years
    ## probably not supposed to have the ID col in the extraction...
    y <- 1:length(yy) - grep(padpoly$EstblsD, yy)
    # y0 <- 1:length(yy) - grep(padpoly$yearStart, yy)
    
    meta[[ paste0('preMean.',variable)]] <-rowMeans(e[,y<0], na.rm = TRUE)
    meta[[ paste0('preSD.', variable) ]] <- apply(e[,y<0],1, sd, na.rm = TRUE)
    meta[[ paste0('preMed.', variable)]] <- apply(e[,y<0],1,median, na.rm = TRUE)
  }
  
  meta$insertDate = as.Date(Sys.time())
  
  return ( meta )
}

insert_table <- function(data, tableName){
  
  # insert data into table 'tableName' for dpar$dbfile
  require(RSQLite)
  
  dbFile <- dpar$dbfile
  CON <- dbConnect(RSQLite::SQLite(), dbFile)
  
  # write to tables
  cat('writing table ', tableName, ' in ', dbFile, '\n')
  retry <- F
  while( retry ){
    
    res <- try({ dbWriteTable(CON, tableName, value = data, append = TRUE) } )
    if( class(res) != 'try-error' ) {
      retry <- FALSE
      cat('\tsuccess\n')
    } else {
      cat('\twrite try failed -- sleeping for 5 seconds\n')
      Sys.sleep(5)
    }
  }
  
  dbDisconnect(CON)
  cat('finished\n')
  
}

getGowerByGroup <- function(uni, padpixels, candidates, threshFunction){
  
  #Loop through unique combos of soilps and EC and get Gower output
  gower.out <- list( index = matrix(NA, 100, nrow(padpixels)), distance =  matrix(NA, 100, nrow(padpixels)) )
  toponames <- gsub(".tif","", names(dpar$topoVars))
  
  for(i in 1:nrow(uni)){
    
    cat('    running gower similarity for group ', i , ' of ', nrow(uni), '\n')
    
    j <- which( padpixels$uni == uni$uni[i])
    
    k <- threshFunction(candidates, unlist(uni[i,-grep('uni', names(uni))]))
    
    # select comparison vars in pad
    dat1 <- padpixels[j, c(toponames)]
    
    # select buffered candidate references
    dat2 <- candidates[k,c(toponames)]
    
    # get similarity matrix for top 100 ref pixels for each target pixel
    cat('getting similarity index \n')
    toposim <- gower_topn(x=dat1, y=dat2, n=100, nthread = 1)
    gower.out$index[,j] <- k[toposim$index]
    gower.out$distance[,j] <- toposim$distance
    
    
  }
  gower.out
}

get_poly <- function(dpar, id){
  cat(' ** loading polygon \n')
  polygons <- st_read(dpar$polygons)
  polygons$id <- polygons$polyID
  padpoly <- polygons[ which(polygons$id == id),]
  return(padpoly)
}  # Updated

# Get spatrasters of masks, topolayers, filtervars
get_fun <- function(x){
  # For all mask variables x, apply a function to their path names y:
  s <- lapply( names(x), function(y) { 
    cat('getting ', y, '\n')
    rpath <- x[[y]]
    r <- rast(rpath)
  })
  rast(s)
}    # Updated

## GLT - This function actually brings in the rasters
getRasterData <- function(padpoly, dpar, radius){ 
  
  cat(' ** loading raster data\n')   
  toponames <- gsub(".tif","", names(dpar$topoVars))
  
  masks <- get_fun(dpar$maskVars)
  names(masks) <- names(dpar$maskVars)
  topovars <- get_fun(dpar$topoVars)
  names(topovars) <- names(dpar$topoVars)
  filtervars <- get_fun(dpar$filterVars)
  names(filtervars) <- names(dpar$filterVars)
  
  rast.proj <- crs(masks[['refrast']])
  
  b <- c(masks, topovars, filtervars)
  
  # Reproject the polys if needed:
  padpoly <- st_transform(padpoly, rast.proj)
  
  ## Prepare neighborhood
  padpolybuf <- st_buffer(padpoly, dist = radius)
  padpolybuf$rastval <- 1
  
  ## crop
  cat('\tcropping predictors\n')
  zone <- crop(b, padpolybuf)
  
  return( zone )
}   # Updated

get_treatment <- function(padpoly, dpar, zone){
  
  cat(' ** getting treated pixels\n')
  
  # return raster of padpixels and unique combos of EC / soil PS
  # padpoly = pad polygon
  # dpar = list of parameters
  # zone = cropped covariate brick
  
  # apply interior buffer
  padpolyi <- st_buffer(padpoly, dist = dpar$innerRad)  
  
  if(is.null(padpolyi)){
    
    stop('ERROR not enough pixels in polygon\n')
    
  }
  
  ## Create rasters 
  padpolyi$rastval <- 1
  padrast <- rasterize(padpolyi, crop(zone[[1]], padpolyi),field=padpolyi$rastval, datatype='INT1U')
  
  # make raster stack
  padstk <- crop(zone, padpolyi)
  
  # mask out unwanted  pixels in treated area
  # Function to mask NLCD
  nlcd_fn <- function(x, y) {
    ind <- ifel(x!=21 & x!=22 & x!=24 & x!=81 & x!=82 & x!=11 & x!=12 & y==1, 1, NA)
    return(ind)
  }
  ###
  padMask <- padrast
  for( variable in dpar$interiorMaskVars){
    if(variable == 'nlcd'){
      
      padMask  <- nlcd_fn(padstk[[variable]], padMask)
    } else {
      padMask  <- padstk[[variable]] * padMask
    }  
  }
  
  padMask[padMask[] == 0] <- NA
  padstk <- padMask * padstk
  names(padstk) <-  names(zone)  
  
  # check that there are enough padpixels
  cat('\tchecking for enough non-na pixels in treatment area\n')
  if(all(is.na(padrast[]))) stop( 'ERROR2: not enough pixels in treatment polygon;')  
  
  if(length(which(!is.na(padrast[]))) == 1){
    # do something to avoid error
  }
  
  cat('\tconverting treatment raster to sf\n')
  
  padpixels <- as.polygons(x = padstk, aggregate = F, trunc=F, dissolve=F, values=T, na.rm=T, na.all=F, extent=F) %>% st_as_sf()
  padpixels$ec <- padpixels$soilec / 100 # Rescale back to actual ec units from scaled integer used to store raster
  
  # ## subSample Here
  ## GLT - Turning this off for now
  # if( dpar$subSample ){
  # 
  #   inout <- subSample(padpixels)
  #   padpixels <- padpixels[which(inout == 1),]
  # 
  # }
  
  # find unique combinations of soil ps and ec
  uni  = unique(padpixels[,c('soilps', 'ec')])      ## 
  uni$uni <- paste0(uni$soilps, '-', uni$ec)
  padpixels$uni <- paste0(padpixels$soilps, '-', padpixels$ec)
  
  return( list(padpixels = padpixels, uni = uni, padpoly = padpoly, padstk = padstk))
  
}     # Updated

get_candidates <- function(dpar, radius, padpoly, zone){
  
  #Generate raster brick of masked candidate pixels
  padpolybuf <- st_buffer(padpoly, dist = dpar$rad)
  padpolybuf$rastval <- 1
  padbufrast <- rasterize(padpolybuf, zone[[1]], field=padpolybuf$rastval, datatype='INT1U')
  
  ## Get the mask, topo, and filter variables
  # toponames <- names(dpar$topoVars)
  masks <- get_fun(dpar$maskVars)
  names(masks) <- names(dpar$maskVars)
  topovars <- get_fun(dpar$topoVars)
  names(topovars) <- names(dpar$topoVars)
  filtervars <- get_fun(dpar$filterVars)
  names(filtervars) <- names(dpar$filterVars)
  
  ## Screen out unwanted disturbances in buffer zone
  # mask out unwanted  pixels in treated area
  # Function to mask NLCD
  nlcd_fn <- function(x, y) {
    ind <- ifel(x!=21 & x!=22 & x!=24 & x!=81 & x!=82 & x!=11 & x!=12 & y==1, 1, NA)
    return(ind)
  }
  
  for( msk in names(masks)){
    if( msk == 'refrast') next()
    cat('\tmasking ', msk,'\n')
    if( msk == 'nlcd'){
      padbufrast <- nlcd_fn(zone[[msk]], padbufrast)
    } else {
      padbufrast <- padbufrast * zone[[msk]]
    }
  }
  
  # mask everything - check this...
  cat('\tmasking buffer\n')
  padbufrast[which(padbufrast[] == 0)] <- NA
  padbufstk <- padbufrast * zone
  
  # make donut
  cat('\tmaking donut\n')
  padpolyb <- st_buffer(padpoly, dist = dpar$buffer)
  padbufstk <- mask(padbufstk, padpolyb, inverse = TRUE)
  
  # propagate names
  names(padbufstk) <- names(zone)  
  
  # checks
  cat('\tchecking number of candidates in buffer\n')
  if(all(is.na(padbufrast[]))) stop( 'ERROR1: no non-na candidate pixels in buffer (every pixel masked)')    
  
  cat('\tconverting buffer to sf \n')
  padbufpixels <- as.polygons(x = padbufstk, aggregate = F, trunc=F, dissolve=F, values=T, na.rm=T, na.all=F, extent=F) %>% st_as_sf()
  
  padbufpixels$ec <- padbufpixels$soilec / 100 # Rescale back to actual ec units from scaled integer used to store raster
  
  return( padbufpixels )
}

get_shmandidates <- function(dpar, id){
  
  # get data
  cat('loading polygons\n')
  # load(dpar$polygons)
  # polygons$id <- polygons[[dpar$polygon_id_column]]
  polygons <- st_read(dpar$polygons)
  # polygons$id <- polygons$polyID
  # padpoly <- polygons[ which(polygons$id == id),]
  padpoly <- polygons[ which(polygons$polyID == id),]
  toponames <- gsub(".tif","", names(dpar$topoVars))
  
  ## GLT - update raster import method from brick() to c(rast(), ...)
  get_fun <- function(x) { 
    
    s <- lapply( names(x), function(y) { 
      cat('getting ', y, '\n')  
      r <- brick(x[[y]]) 
      # r <- rast(r)
      names(r) <- y
    })
    stack(s)
  }
  
  masks <- get_fun(dpar$maskVars)
  topovars <- get_fun(dpar$topoVars)
  filtervars <- get_fun(dpar$filterVars)
  b <- stack(masks, topovars, filtervars)
  rast.proj <- projection(masks[['refrast']])
  
  ## Prepare neighborhood
  padpolybuf <- st_buffer(padpoly, dist = dpar$rad)
  padpolybuf$rastval <- 1
  
  ## crop
  cat('cropping predictors\n')
  zone <- crop(b, padpolybuf)
  padbufrast <- rasterize(padpolybuf, zone[[1]], field=padpolybuf$rastval, datatype='INT1U')
  
  
  ## Screen out unwanted disturbances in buffer zone
  ## GLT - update masking method
  f_mask <- function(a,b) a*b
  
  nlcd_fn <- function(nlcd,padbufrast) {
    ind <- ifelse(nlcd!=21&nlcd!=22&nlcd!=24&nlcd!=81&nlcd!=82&nlcd!=11&nlcd!=12&padbufrast==1,1,NA)
    return(ind)
  }
  
  
  for( msk in names(masks)){
    if( msk == 'refrast') next()
    cat('masking ', msk,'\n')
    if( msk == 'nlcd'){
      padbufrast <- overlay( zone[[msk]], padbufrast , fun=nlcd_fn)
    } else {
      padbufrast <- overlay(padbufrast, zone[[msk]], fun = f_mask)
    }
  }
  
  # mask everything
  ## GLT - update masking method
  cat('masking buffer\n')
  padbufrast[which(padbufrast[] == 0)] <- NA
  padbufstk <- overlay(padbufrast, zone, fun = f_mask)
  
  # make donut
  cat('making donut\n')
  padpolyb <- st_buffer(padpoly, dist = dpar$buffer)
  padbufstk <- mask(padbufstk, padpolyb, inverse = TRUE)
  
  # propagate names
  names(padbufstk) <-  names(zone)  
  
  # checks
  cat('checking number of candidates in buffer\n')
  if(all(is.na(padbufrast[]))) stop( 'ERROR1: no non-na candidate pixels in buffer (every pixel masked)')    
  
  cat('converting buffer to Spatial Pixels\n')
  # padbufpixels <- as(padbufstk, "SpatialPixelsDataFrame")
  padpixels <- as.polygons(x = rast(padstk), trunc=F, dissolve=F, values=T, na.rm=T, na.all=F, extent=F) %>% st_as_sf()
  
  padbufpixels$ec <- padbufpixels$soilec / 100 # Rescale back to actual ec units from scaled integer used to store raster
  
  
  ################################################################
  
  # prepare treatment area   
  cat('preparing treatment area by inner buffer of ', dpar$innerRad, '\n')
  padpolyi <- st_buffer(padpoly, dist = dpar$innerRad)  
  if(is.null(padpolyi)){
    stop('ERROR not enough pixels in polygon\n')
    padpolyi <- polygons[ which(polygons$id == id),]
  }
  
  # for rasterization
  padpolyi$rastval <- 1
  
  ## Create rasters 
  cat('rasterizing treatment polygon\n')
  padrast <- rasterize(padpolyi, crop(zone[[1]], padpolyi),field=padpolyi$rastval, datatype='INT1U')
  
  ## if area too small, remove buffer -- should be redundant
  if(all(is.na(values(padrast)))){
    padpolyi <- polygons[ which(polygons$id == id),]
    padpolyi$rastval <- 1
    padrast <- rasterize(padpolyi, crop(zone[[1]], padpolyi),field=padpolyi$rastval, datatype='INT1U')
  }
  
  
  # make donut hole
  cat('cropping and masking treatment raster stack\n')  
  padstk <- crop(zone, padpolyi)
  
  # mask out unwanted  pixels in treated area
  padMask  <- overlay( padstk[['roadrast']], padstk[[1]], fun = f_mask)
  padMask  <- overlay( padstk[['oilgas']], padMask, fun = f_mask)
  padMask  <- overlay( padstk[['oilgas4corners']], padMask, fun = f_mask)
  padMask  <- overlay( padstk[['nlcd']], padMask, fun = nlcd_fn)
  padMask  <- overlay( padrast, padMask, fun = f_mask)
  padMask[padMask[] == 0] <- NA
  padstk <- overlay(padMask, padstk, fun = f_mask)
  names(padstk) <-  names(zone)  
  
  
  # check that there are enough padpixels
  cat('checking for enough non-na pixels in treatment area\n')
  if(all(is.na(padrast[]))) stop( 'ERROR2: not enough pixels in treatment polygon; area = ', area(padpoly)/900, ' pixels')  
  
  
  if(length(which(!is.na(padrast[]))) == 1){
    # do something to avoid error
  }
  
  
  cat('
	  Convert treatment raster to Spatial Pixels DFs 
  \n')
  
  # padpixels <- as(padstk, "SpatialPixelsDataFrame")
  padpixels <- as.polygons(x = rast(padstk), trunc=F, dissolve=F, values=T, na.rm=T, na.all=F, extent=F) %>% st_as_sf()
  
  padpixels$ec <- padpixels$soilec / 100 # Rescale back to actual ec units from scaled integer used to store raster
  
  
  ###############################################################
  
  cat('
  ## Summarizing pad pixels
  \n')
  
  padpscclasses <- as.numeric(names(summary(as.factor(padpixels$soilps))))
  getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  } # Fn to calculate mode if that is desirable
  
  padpscmode <- getmode(padpixels$psc)
  ecave <- mean(padpixels$soilec, na.rm=TRUE)
  
  # find unique combinations of soil ps and ec
  uni  = unique(padpixels[,c('soilps', 'soilec')])
  uni$uni <- paste0(uni$soilps, '-', uni$soilec)
  padpixels$uni <- match(paste0(padpixels$soilps, '-', padpixels$soilec), uni$uni)
  
  # check to see that there are at least 100 candidates per combo (using a +/- thresh value for ec)
  threshFunction <- function(x) length( which( padbufpixels$soilps == x[1] & 
                                                 padbufpixels$soilec > (x[2] * .95) & 
                                                 padbufpixels$soilec < (x[2] * 1.05)
  ))
  
  avail  = apply(uni[,1:2], 1, threshFunction)
  
  if( any(avail < 100 ) ) {
    
    # if not last try, fail and expand search radius
    cat( 'not enough pixels within threshold for some treated, \n'); 
    uni[which (avail < 100),] ; 
    stop('\nnot enough quality pixels\n')
  }
  
  
  #extend EC range by .05
  ecmin <- min(padpixels$soilec, na.rm=TRUE) - min(padpixels$soilec, na.rm=TRUE) -dpar$ecBuf
  ecmax <- max(padpixels$soilec, na.rm=TRUE) + max(padpixels$soilec, na.rm=TRUE) +dpar$ecBuf
  
  ## Select pad buffer pixels with same PSC and EC
  padbufsoilpixels <- subset(padbufpixels, padbufpixels$soilps %in% padpscclasses &    padbufpixels$soilec < ecmax & padbufpixels$soilec > ecmin)
  
  if(nrow(padbufsoilpixels) == 0){
    
    stop( 'ERROR3: not enough candidates w/correct EC/PSD in buffer')
    
  }
  
  return(list(padbufsoilpixels = padbufsoilpixels,padbufpixels=padbufpixels, padpixels = padpixels, padpoly= padpoly, padpolyi = padpolyi, padrast = padrast, padbufrast = padbufrast, padstk = padstk, padbufstk = padbufstk, ecave = ecave, ecmin = ecmin, ecmax = ecmax, zone = zone, toponames = toponames, uni = uni))
  
}
## Not including the get_shmandidates function bc it isn't actually called anywhere

# select all pixels chosen in top 100 gower pixels
get_chosen <- function(toposim, padbufsoilpixels){
  
  a <- vv <- toposim$index
  aindex <- unique(c(a))
  allpixels <- padbufsoilpixels[aindex,]
  return(allpixels)
  
}

extract_TS <- function(responseVars, padpixels, allpixels, toposim){
  
  # extract response data from dart output
  # bfile = path to brick of things to extract
  # f = filename
  # rename = replace original file with 'renamed' file
  
  ext = ext(rbind(padpixels[,1], allpixels[,1]))
  ext[1] <- ext[1] - 200
  ext[2] <- ext[2] + 200
  ext[3] <- ext[3] - 200
  ext[4] <- ext[4] + 200
  
  # box <- as(ext, 'SpatialPolygons')    ## basically a bounding box
  box <- as.polygons(ext)
  # projection(box) <- projection(padpixels)
  crs(box) <- crs(padpixels)
  
  cat('loading response variables\n'); flush.console()
  loadCrop <- function(x){
    s <- rast(x)
    b <- box
    if(crs(s) != crs(padpixels)){
      b <- project(b, crs(s))
    }
    s <- crop(s, b)
    # names(s) <- paste0(extension (basename(x), ''), 1:(dim(s)[3]))
    s
  }
  
  BB <- lapply(responseVars, loadCrop)
  extractedNames <- sapply(BB, names)
  
  
  st <- Sys.time()
  
  cat('extracting padvals\n') ; flush.console()
  # padvals <- lapply(BB, function(x) {extract(x, padpixels, progress = 'text')})
  padvals <- lapply(BB, function(x) {terra::extract(x, padpixels, ID=F)})
  cat('extracting reference vals\n') ; flush.console()
  # allvals <- lapply(BB, function(x) {extract(x, allpixels, progress = 'text')})
  allvals <- lapply(BB, function(x) {terra::extract(x, allpixels, ID=F)})
  
  get_stats <- function(x) { c( 'avg' = mean(x, na.rm= TRUE), 'med' = median(x, na.rm = TRUE), 'sd' = sd(x, na.rm = TRUE),
                                'N' = sum(!is.na(x))) }
  
  
  a <- toposim$index
  aindex <- unique(c(a)) # unique reference pixels
  
  
  timeElapsed <- Sys.time() - st
  
  return ( list( refIndex = aindex, extractedTarget =padvals, extractedReference=allvals, ExtractedTimeElapsed = timeElapsed, extractedNames = extractedNames) )
  
}


# check status of dart process
dartStatus <- function(dir){
  
  require(raster)
  ff <- list.files(dir, full = TRUE, pattern = 'RData')
  load(ff[1])
  s <- shapefile(dpar$polygons)
  pids <- s[ , dpar$polygon_id_column]
  has <- gsub('[^0-9]', '', ff)
  d <- setdiff(pids, has)
  cat('missing ', length(d), ' out of ', length(pids),' polygons\n')
  d
  
}

# split polygons by grid
split <- function(p, cellsz = 1000){
  
  x <- st_as_sf(p)
  grid <- st_make_grid(p, cellsize = cellsz)
  int <- st_intersection(x,grid)
  as(int, 'Spatial')
}




evalRun <- function(dir){
  
  setwd(dir)
  
  load('results.RData')
  source('dart_functions.R')
  
  ds <- dartStatus('..')
  k <- grep('Error', result)
  
  error <- data.frame(id = ds, error = as.character(result[k]))
  
  error$err <- gsub('area.*','', error$error)
  
  
  print(table(error$err))
  save(error, file = 'errors.RData')
}


# Sample 

subSample <- function( p ){
  
  # p = spatialPixelsDataFrame
  # create a column of 'yes' or 'no' based on raster grid w/ at least one space around every point
  
  r <- raster(p)
  nc = dim(r)[2] # rasters filled by row
  nr = dim(r)[1] # rasters filled by row
  
  r[] <- 0
  
  # make sure at least one pixel is selected
  # check if first pixel in odd column
  oddCol <- colFromX(r, coordinates(p)[1,1]) %% 2 == 1
  if(oddCol) { dotVec <-  c(1,0) } else { dotVec <- c(0,1)}
  dots <- rep(dotVec, length.out = nc)
  
  #check if first pixel in odd row
  oddRow <- rowFromY(r, coordinates(p)[1,2]) %% 2 == 1
  blank <- rep(c(0), length.out = nc)
  
  if( oddRow){
    r[] <- rep( c( dots , blank ), length.out = nc * nr)
  } else {
    r[] <- rep( c( blank , dots), length.out = nc * nr)
  }
  extract(r, p)
  
}



# Visualization

basicPlot <- function(datafile){
  
  load(datafile)
  require(raster)
  ## allpixels : SPDF of candidate pixel population
  ## artpixels : SPDF of art pixels
  ## avedist : average distance between dart pixels and target
  ## dpar : Parameter object for the extraction run
  ## padbufrast : raster of the extent of the pad buffer
  ## padpixels :  SPDF of treatment area
  ## padpoly : SPolygonsDF of treatment area
  ## padpolyi : SPolygonsDF of treatment area - buffer (the business region)
  ## padstk : raster stack of treatment area covariates
  ## timeElapsed : time for dart run
  ## toposim : distance matrix where rows = top 100 candidates, columns = target pixels
  ##     index = index in allpixels of candidate
  ##     distance = distance of candidate pixel from column id target
  ##     aindex <- unique(c(toposim$index)) # unique list of candidate reference pixels
  
  ## Basic plot of pad + buffers + candidates
  layout(matrix(1:2, 1,2))
  opar <- par()
  par(mar = c(5.1,4.1, 4.1,5.7))
  plot(padbufrast, legend = FALSE, col = 'light grey', main = 'Frequency')
  plot(padpolyi, add = TRUE)
  
  aindex <- unique(c(toposim$index))
  tb <- table(toposim$index[])
  allpixels$freq <- NA
  allpixels$freq[ match(as.numeric( names(tb) ) , aindex) ] <- tb
  plot( raster(allpixels['freq']), add = TRUE, col = rev(topo.colors(10)))
  
  ## average distance of each candidate pixel to reference
  distances = c(toposim$distance[])
  ids = c(toposim$index[])
  avedists <- aggregate( distances, list(id = ids), FUN = mean)
  allpixels$avedist <- NA
  allpixels$avedist[ match(avedists$id, aindex) ] <- avedists$x
  
  
  plot(padbufrast, legend = FALSE, col = 'light grey', main = '[Environmental] Distance')
  plot(padpolyi, add = TRUE)
  plot( raster( allpixels['avedist'] ), add = TRUE, col = rev(topo.colors(10)))
  par(opar) 
}


# plot timeseries
tsPlot <- function(datafile, xnames = 1984:2017){
  
  load(datafile)
  require(raster)
  
  target <- extraction$extractedTarget
  controls <- extraction$extractedReference
  
  yl = range(c(target[], controls[]), na.rm = TRUE)
  x <- 1:ncol(target)
  
  #png(fname, height = 7.5, width = 14, units = 'in', res = 300)
  layout(1)
  plot(x, x, ylim = yl, type = 'n', ylab = '', xaxt = 'n')
  axis(1,at = x, labels = xnames)
  
  apply(controls, 1, lines, col = rgb(0,0,0,.1))
  apply(target, 1, lines, col = rgb(1,0,0, .1))
  legend('topleft', col = c('red', 'black'), lty = 1, legend = c('treated', 'control'), bty = 'n')
  #dev.off()
  
}


# Format extracted data to long panel conventions
# 

longPanel <- function(extraction, yr, varname, covs=FALSE, testZeros=TRUE){
  
  # produce `long` data.frame with the following columns
  # id = name of pixel
  # y = response value
  # D = treatment occurred
  # time = timestep
  
  trt <- extraction$extractedTarget[[varname]]
  
  if(class(trt) == 'numeric') trt <- t(matrix(trt))
  tn <- paste0('trt', 1:nrow(trt))
  
  if(testZeros){
    if( any ( apply( trt, 1, function(x) length(which(x == 0)) > 2))) stop('too many zeros')
  }
  # to do -- add extracted year names somehow in DART process
  tm <- rev((2018:1984)[1:ncol(trt)])
  tD <- rep(0, length(tm))
  tD[grep(yr, tm):length(tm)] <- 1
  
  g1 <- data.frame(id = rep(tn, each = ncol(trt)), y = c(t(trt)), D = rep(tD, nrow(trt)), time = rep(tm, nrow(trt)), stringsAsFactors = FALSE)
  
  if(covs){
    
    toAdd <- setdiff(names(extraction$extractedTarget),varname)
    for( v in toAdd){
      
      if(v == 'npp') next()
      
      tx <- t(matrix(extraction$extractedTarget[[v]]))
      g1[[v]] <- c(t(tx))
    }
  }
  
  ctr <- extraction$extractedReference[[varname]]
  cn <- paste0('ctr', extraction$refIndex) # use spatialPixels index
  
  g2 <- data.frame(id = rep(cn, each = ncol(ctr)), y = c(t(ctr)), D = 0, time = rep(tm, nrow(ctr)), stringsAsFactors = FALSE)
  
  if(covs){
    
    toAdd <- setdiff(names(extraction$extractedReference),varname)
    for( v in toAdd){
      
      if(v == 'npp') next()
      
      tx <- t(matrix(extraction$extractedReference[[v]]))
      g2[[v]] <- c(t(tx))
    }
  }
  
  
  rbind(g1, g2)
  
  
}



get_art <- function( padbufsoilpixels, padpixels, toposim , toponames){
  
  # if there are enough pixels ...
  if(length(padbufsoilpixels) > 100){
    
    # get list of top choices per target pixel
    toposimf1 <- as.data.frame(factor(toposim$index[1,], levels = unique(toposim$index[1,])))
    colnames(toposimf1) <- "rows"
    
    # tabulate frequency
    toposimf1 <- count(toposimf1, "rows")
    
    #sort
    toposimf1$rows<-factor(toposimf1$rows, levels=toposimf1$rows[order(toposimf1$freq, decreasing=TRUE)])
    
    # create pool of candidates to draw from
    toposimrows <- as.character(toposimf1$rows)
    toposimrows <- unique(toposimrows)
    toposimrowstogo <- 100-length(toposimrows)
    
    # make list of top 25 candidates from each pixel
    toposimf25 <- as.data.frame(as.factor(toposim$index[2:25,]))
    colnames(toposimf25) <- "rows"
    toposimf25 <- count(toposimf25, "rows")
    toposimf25$rows<-factor(toposimf25$rows, levels=toposimf25$rows[order(toposimf25$freq, decreasing=TRUE)])
    
    # if still need more references
    if(toposimrowstogo > 0){
      if(length(toposimf25$rows)<toposimrowstogo){
        # add next 25 rows to selected rows
        toposimrows <- append(toposimrows, as.character(toposimf25$rows))
        toposimrows <- unique(toposimrows)
        toposimrowstogo <- 100-length(toposimrows)
      } else {
        #prepare candidates
        toposimrows <- append(toposimrows, names(summary(toposimf25$rows)[1:toposimrowstogo]))
        toposimrows <- unique(toposimrows)
        toposimrowstogo <- 100-length(toposimrows)
      }
      
      if (length(toposimrowstogo)>0){
        #look among top 100 for each pixel
        toposimf100 <- as.data.frame(as.factor(toposim$index[26:100,]))
        colnames(toposimf100) <- "rows"
        toposimf100 <- count(toposimf100, "rows")
        toposimf100$rows<-factor(toposimf100$rows, levels=toposimf100$rows[order(toposimf100$freq, decreasing=TRUE)])
        toposimrows <- append(toposimrows, names(summary(toposimf100$rows)[1:toposimrowstogo]))
        toposimrows <- unique(toposimrows)
        toposimrowstogo <- 100-length(toposimrows)
        while(toposimrowstogo > 0){
          toposimf100 <- subset(toposimf100,!(as.character(rows) %in% toposimrows))
          colnames(toposimf100) <- "rows"
          toposimf100 <- count(toposimf100, "rows")
          toposimf100$rows<-factor(toposimf100$rows, levels=toposimf100$rows[order(toposimf100$freq, decreasing=TRUE)])
          toposimrows <- append(toposimrows, names(summary(toposimf100$rows)[1:toposimrowstogo]))
          toposimrows <- unique(toposimrows)
          toposimrowstogo <- 100-length(toposimrows)
        }
      }
    } 
    artpixels <- padbufsoilpixels[as.numeric(toposimrows),][1:100,]
    rowindex <- which(toposim$index %in% toposimrows, arr.ind = T)
    distances <- toposim$distance[rowindex]
    avedist <- mean(distances)
  } else {
    artpixels <- padbufsoilpixels
    dat1 <- padpixels[,c(toponames)]
    dat2 <- artpixels[,c(toponames)]
    toposim <- gower_topn(x=dat1, y=dat2, n=length(padbufsoilpixels), nthread = 1)
    avedist <- mean(toposim$distance)
  }
  
  return(list(artpixels = artpixels, avedist = avedist))
  
}




