##spatial mapping and plotting packages
library(ggplot2)
library(rasterVis)
library(maptools)
library(maps)
library(dplyr)
library(raster)
library(rgeos)
library(rgdal)
library(sf)
library(tmap)
##spatial modeling packages
library(seegSDM)      #generation of background points biased on population density
library(dismo)
library(rJava)        #MaxENT
library(SDMtune)
library(zeallot)      # data splitting 
library(gbm)          # basic implementation
library(xgboost)      # a faster implementation of gbm
library(caret)        # an aggregator package for performing many machine learning models
library(kableExtra)   # Compile ROC reports
library(blockCV)      # spatial cross-validation
library(randomForest)
library(embarcadero)  # modeling BART
library(hSDM)
library(coda)         # summarizing and plotting mcmc outputs
library(doParallel)   # running parallel chains
library(parallel)
library(ecospat)