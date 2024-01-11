[1] "C:/Users/Kristina/OneDrive - New Mexico State University/Desktop/GIT REPOs/DART_LTDL/DART_LTDL"
> # Set your working directory
  > setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
> getwd()
[1] "C:/Users/Kristina/OneDrive - New Mexico State University/Desktop/GIT REPOs/DART_LTDL/DART_LTDL"
> # Set your working directory
  > setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
> list.files()
[1] "1_run_dart_updated.R"              "2_genMetadata_updated.R"          
[3] "3_genTodo_updated.R"               "4_genSynthCtr_updated.R"          
[5] "5_genDataFrame_updated.R"          "colorado_plateau_clip.cpg"        
[7] "colorado_plateau_clip.prj"         "colorado_plateau_clip.qmd"        
[9] "colorado_plateau_clip.shx"         "dart_functions_updated.R"         
[11] "DART_LTDL.Rproj"                   "DATA_creating_polygons_from_gdb.R"
[13] "params_updated.R"                  "README.md"                        
> dir.create(paste0(getwd(), "/outdir"))
> install.packages("terra")
WARNING: Rtools is required to build R packages but is not currently installed. Please download and install the appropriate version of Rtools before proceeding:
  
  https://cran.rstudio.com/bin/windows/Rtools/
  trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.3/terra_1.7-65.zip'
Content type 'application/zip' length 39066427 bytes (37.3 MB)
downloaded 37.3 MB

package ‘terra’ successfully unpacked and MD5 sums checked

The downloaded binary packages are in
C:\Users\Kristina\AppData\Local\Temp\RtmpY32N1x\downloaded_packages
> sessionInfo()
R version 4.3.0 (2023-04-21 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19045)

Matrix products: default


locale:
  [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8   
[3] LC_MONETARY=English_United States.utf8 LC_NUMERIC=C                          
[5] LC_TIME=English_United States.utf8    

time zone: America/Denver
tzcode source: internal

attached base packages:
  [1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
  [1] beepr_1.3       lubridate_1.9.2 forcats_1.0.0   stringr_1.5.0   dplyr_1.1.2     purrr_1.0.1    
[7] readr_2.1.4     tidyr_1.3.0     tibble_3.2.1    ggplot2_3.4.4   tidyverse_2.0.0 sf_1.0-12      

loaded via a namespace (and not attached):
  [1] gtable_0.3.3       compiler_4.3.0     tidyselect_1.2.0   Rcpp_1.0.10        scales_1.2.1      
[6] R6_2.5.1           generics_0.1.3     classInt_0.4-9     units_0.8-2        munsell_0.5.0     
[11] DBI_1.1.3          pillar_1.9.0       tzdb_0.4.0         rlang_1.1.1        utf8_1.2.3        
[16] stringi_1.7.12     audio_0.1-10       timechange_0.2.0   cli_3.6.1          withr_2.5.0       
[21] magrittr_2.0.3     class_7.3-21       grid_4.3.0         rstudioapi_0.14    hms_1.1.3         
[26] lifecycle_1.0.3    vctrs_0.6.5        KernSmooth_2.23-20 proxy_0.4-27       glue_1.6.2        
[31] fansi_1.0.4        e1071_1.7-13       colorspace_2.1-0   tools_4.3.0        pkgconfig_2.0.3   
> install.packages("R")
Warning in install.packages :
  package ‘R’ is not available for this version of R

A version of this package for your version of R might be available elsewhere,
see the ideas at
https://cran.r-project.org/doc/manuals/r-patched/R-admin.html#Installing-packages
> ?terra
No documentation for ‘terra’ in specified packages and libraries:
  you could try ‘??terra’
> library("terra")
terra 1.7.65

Attaching package: ‘terra’

The following object is masked from ‘package:tidyr’:
  
  extract

Warning message:
  package ‘terra’ was built under R version 4.3.2 
> ?terra
> download.file("http://rangeland.ntsg.umt.edu/data/rap/rap-vegetation-cover/v3/vegetation-cover-v3-1986.tif")
Error in download.file("http://rangeland.ntsg.umt.edu/data/rap/rap-vegetation-cover/v3/vegetation-cover-v3-1986.tif") : 
  argument "destfile" is missing, with no default
> download.file("http://rangeland.ntsg.umt.edu/data/rap/rap-vegetation-cover/v3/vegetation-cover-v3-1986.tif", destfile = vegetation-cover-v3-1986.tif)
Error: unexpected symbol in "download.file("http://rangeland.ntsg.umt.edu/data/rap/rap-vegetation-cover/v3/vegetation-cover-v3-1986.tif", destfile = vegetation-cover-v3-1986.tif"
> download.file("http://rangeland.ntsg.umt.edu/data/rap/rap-vegetation-cover/v3/vegetation-cover-v3-1986.tif", destfile = "vegetation-cover-v3-1986.tif")
trying URL 'http://rangeland.ntsg.umt.edu/data/rap/rap-vegetation-cover/v3/vegetation-cover-v3-1986.tif'
Content type 'image/tiff' length 41742837061 bytes (39809.1 MB)
downloaded 1259.8 MB

Error in download.file("http://rangeland.ntsg.umt.edu/data/rap/rap-vegetation-cover/v3/vegetation-cover-v3-1986.tif",  : 
                         download from 'http://rangeland.ntsg.umt.edu/data/rap/rap-vegetation-cover/v3/vegetation-cover-v3-1986.tif' failed
                       In addition: Warning messages:
                         1: In download.file("http://rangeland.ntsg.umt.edu/data/rap/rap-vegetation-cover/v3/vegetation-cover-v3-1986.tif",  :
                                               downloaded length 0 != reported length 0
                                             2: In download.file("http://rangeland.ntsg.umt.edu/data/rap/rap-vegetation-cover/v3/vegetation-cover-v3-1986.tif",  :
                                                                   URL 'http://rangeland.ntsg.umt.edu/data/rap/rap-vegetation-cover/v3/vegetation-cover-v3-1986.tif': Timeout of 60 seconds was reached
                                                                 > polygons <- vect("data/colorado_plateau.shp")
                                                                 Error: [vect] Cannot open this file as a SpatVector: C:\Users\Kristina\OneDrive - New Mexico State University\Desktop\GIT REPOs\DART_LTDL\DART_LTDL\data\colorado_plateau.shp
                                                                 In addition: Warning message:
                                                                   Unable to open C:\Users\Kristina\OneDrive - New Mexico State University\Desktop\GIT REPOs\DART_LTDL\DART_LTDL\data\colorado_plateau.shx or C:\Users\Kristina\OneDrive - New Mexico State University\Desktop\GIT REPOs\DART_LTDL\DART_LTDL\data\colorado_plateau.SHX. Set SHAPE_RESTORE_SHX config option to YES to restore or create it. (GDAL error 4) 
                                                                 > polygons <- vect("data/colorado_plateau.shp")
                                                                 Error: [vect] file does not exist: data/colorado_plateau.shp
                                                                 > polygons <- vect("data/colorado_plateau_clip.shp")
                                                                 > plot(polygons)
                                                                 > class(polygons)
                                                                 [1] "SpatVector"
                                                                 attr(,"package")
                                                                 [1] "terra"
                                                                 > head(polygons)
                                                                 data frame with 0 columns and 0 rows
                                                                 > names(polygons)
                                                                 character(0)
                                                                 > View(polygons)
                                                                 > print(polygons)
                                                                 class       : SpatVector 
                                                                 geometry    : polygons 
                                                                 dimensions  : 112, 0  (geometries, attributes)
                                                                 extent      : -112.5029, -106.599, 35.99993, 41.1252  (xmin, xmax, ymin, ymax)
                                                                 source      : colorado_plateau_clip.shp
                                                                 coord. ref. : lon/lat GCS_unknown 
                                                                 Warning message:
                                                                   In class(object) <- "environment" :
                                                                   Setting class(x) to "environment" sets attribute to NULL; result will no longer be an S4 object
                                                                 > length(polygons)
                                                                 [1] 112
                                                                 > plot(polygons, 14)
                                                                 Error: [plot] x only has 0  columns
                                                                 > plot(polygons[14])
                                                                 > plot(polygons[112])
                                                                 > plot(polgyons[c(14, 112)])
                                                                 Error in h(simpleError(msg, call)) : 
                                                                   error in evaluating the argument 'x' in selecting a method for function 'plot': object 'polgyons' not found
                                                                 > plot(polygons[c(14, 112)])
                                                                 > names(polygons)
                                                                 character(0)
                                                                 > polygons$polyID
                                                                 NULL
                                                                 > polygons$polyID <- seq(length(polygons))
                                                                 > names(polygons)
                                                                 [1] "polyID"
                                                                 Warning message:
                                                                   In class(object) <- "environment" :
                                                                   Setting class(x) to "environment" sets attribute to NULL; result will no longer be an S4 object
                                                                 > polygons$polyID
                                                                 [1]   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24
                                                                 [25]  25  26  27  28  29  30  31  32  33  34  35  36  37  38  39  40  41  42  43  44  45  46  47  48
                                                                 [49]  49  50  51  52  53  54  55  56  57  58  59  60  61  62  63  64  65  66  67  68  69  70  71  72
                                                                 [73]  73  74  75  76  77  78  79  80  81  82  83  84  85  86  87  88  89  90  91  92  93  94  95  96
                                                                 [97]  97  98  99 100 101 102 103 104 105 106 107 108 109 110 111 112
                                                                 > ?writeVector
                                                                 > 
                                                                   > writeVector(polygons, "polygons_ID.shp", filetype = "ESRI Shapefile", layer=NULL, insert=FALSE,
                                                                                 +             overwrite=FALSE, options="ENCODING=UTF-8")
                                                                 > list.files
                                                                 function (path = ".", pattern = NULL, all.files = FALSE, full.names = FALSE, 
                                                                           recursive = FALSE, ignore.case = FALSE, include.dirs = FALSE, 
                                                                           no.. = FALSE) 
                                                                   .Internal(list.files(path, pattern, all.files, full.names, recursive, 
                                                                                        ignore.case, include.dirs, no..))
                                                                 <bytecode: 0x00000236fbe0bb00>
                                                                   <environment: namespace:base>
                                                                   > list.files()
                                                                 [1] "1_run_dart_updated.R"              "2_genMetadata_updated.R"          
                                                                 [3] "3_genTodo_updated.R"               "4_genSynthCtr_updated.R"          
                                                                 [5] "5_genDataFrame_updated.R"          "Covariates"                       
                                                                 [7] "dart_functions_updated.R"          "DART_LTDL.Rproj"                  
                                                                 [9] "data"                              "DATA_creating_polygons_from_gdb.R"
                                                                 [11] "Filter_Vars"                       "Masking_Vars"                     
                                                                 [13] "outdir"                            "params_updated.R"                 
                                                                 [15] "polygons_ID.cpg"                   "polygons_ID.dbf"                  
                                                                 [17] "polygons_ID.prj"                   "polygons_ID.shp"                  
                                                                 