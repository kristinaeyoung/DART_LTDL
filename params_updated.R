# run Art

######################################################
# Parameters template

# Set your working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# creating an output directory
# dir.create(paste0(getwd(), "/outdir"))

#R list object, passes along with the code to take the parameters out as needed
# need to make sure all of these files exist ahead of running 1_run_dart_updated
dpar <- list(
  
  ## Number of cores you want to use for the parallel processing
  #ncores = 20,
  ncores = 4,
  
  ## Output directory -- now must be changed here
  # Set working directory to where this R script is saved
  
  outdir = paste0(getwd(), "/outdir"),
  
  # create a data folder to store shapefiles in (/data)
  
  polygons = list.files(path=paste0(getwd(), "/data"), pattern=".shp", full=T),
  polygon_id_column = 'polyID',
  
  ##log file
  logfile = paste0('~/dart_logfile_', format(Sys.time(), '%Y-%m-%d-%H-%M-%S'), '.txt'),
  
  ## Output Prefix
  prefix = 'ID_',
  
  ## Neighborhood radius (m)
  # rad = 3000,
  rad = 2000,
  
  ## Reduce number of treatment pixels?
  subSample = TRUE,
  
  ## Number of times to expand search radius if no candidates found
  tries = 4,
  
  ## Target variable name
  varname = 'bare',
  
  ## Inner buffer radius (shrink treatment area to avoid edge effects)
  innerRad = 30,
  buffer = 90,
  
  ## Predictions
  n_x = 20, # n columns of tiles
  n_y = 20, # n rows of tiles
  
  
  # Filter variables
  filterVars = list(
    soilec = "ec_0to60cm_100xInt_ucrb.tif",
    soilps = "UCRB_mPSC_RFE_10plus.tif"
  ),
  
  # masking variables: 1 = ok, 0 = mask
  maskVars = list(
    refrast = "refrast.tif",
    # roadrast =  "TIGER_2018_ucrb_mask.tif",
    # otherpads =  "WYCOriv_Nopads.tif",
    # oilgas =      "oil_gas_buf_ucrb_mask.tif",
    # oilgas =      "oil_gas_buf_ucrb_mask.tif",
    # oilgas4corners = "fourCorners_oilgas_mask.tif"
    # fires= "mtbs_mask.tif",
    # utblmfires= 'utfire_mask.tif',
    nlcd = "nlcd.tif"
    
  ),
  # masking variables for inside treated area
  # interiorMaskVars = c('roadrast', 'oilgas', 'oilgas4corners', 'nlcdBuf'),
  interiorMaskVars = 'nlcd',
  
  
  ## Variables for the Gower distance matrix (Topo variables)
  # Gower distance does the actual pixel matching
  topoVars = list(
    ELEVm = "ELEVm.tif",
    PCURV = "PCURV.tif",
    TCURV = "TCURV.tif",
    # RELHT1 ="RELHT1.tif",
    # RELHT32 = "RELHT32.tif",
    # RELHT128 ="RELHT128.tif",
    # RELMNHT1 ="RELMNHT1.tif",
    # RELMNHT32 ="RELMNHT32.tif",
    # RELMNHT128 ="RELMNHT128.tif",
    # MRRTF = "MRRTF.tif",
    # MRVBF = "MRVBF.tif",
    SLOPE = "SLOPE.tif",
    SOUTHNESS = "SOUTHNESS.tif",
    # EASTNESS = "EASTNESS.tif",
    TWI_TOPMODEL = "TWI_TOPMODEL.tif",
    # CAlog_10="CAlog_10.tif",
    LFELEMS = "LFELEMS.tif"
  ),
  
  ## Response variables
  ## For multi-year response variables, it's best to format them into multilayer rasters with each layer representing consecutive years of the same variable (e.g., one raster stack for annual forbs/grasses for years 2001-2020)  
  ## Provided example for RAP variables:
  respVars = list(
    afg = "",        
    bgr = "",
    pfg = "",
    shr = ""
  ),
  
  ## Metadata columns to use for downstream dataframes
  ## Artifacts from Steve's PJ paper - may be useful for other users, but they should modify as needed.
  metaCols = c('src', 'srcID', 'id', 'name', 'yearEnd', 'yearStart', 'monthEnd',
               'monthStart', 'treat', 'method', 'seeded', 'burned', 'otherids', 'target',
               'polyID', 'yearStart'),
  
  ## Metadata columns to exclude
  dontKeep =  c("refrast",  "othersites","exclosures")
  
)

## database file
dpar$dbfile = file.path( dirname( dpar$outdir ), paste0(basename(dpar$outdir), '.sqlite'))
