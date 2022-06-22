## functions file
# function file

firstDay <- function (year, month) {
  # given a year and month, return a Date object of the first day of that month
  date_string <- paste(year,
                       month,
                       '01',
                       sep = '-')
  
  date <- as.Date (date_string)
  
  return (date)
  
}

lastDay <- function (year, month) {
  # given a year and month, return a Aate object of the last day of that month
  next_month <- ifelse(month == 12,
                       1,
                       month + 1)
  
  next_year <- ifelse(month == 12,
                      year + 1,
                      year)
  
  next_date_string <- paste(next_year,
                            next_month,
                            '01',
                            sep = '-')
  next_date <- as.Date(next_date_string)
  date <- next_date - 1
  
  return (date)
}

sentenceCase <- function (text) {
  # given a vector of text strings `text`, convert to sentence case
  
  # convert all to lower case
  text <- tolower(text)
  
  # split at spaces
  text_list <- strsplit(text,
                        ' ')
  
  text_list <- lapply(text_list,
                      function(x) {
                        x[1] <- paste(toupper(substring(x[1], 1, 1)),
                                      substring(x[1], 2),
                                      sep = "")
                        x <- paste(x, collapse = ' ')
                        return(x)
                      })
  
  text_vector <- unlist(text_list)
  
  return (text_vector)
  
}

# define functions
firstTwo <- function (text) {
  # given a vector of text strings `text` subset each to only the first two
  # words (bits separated by spaces) and return this as a vector.
  text <- as.character(text)
  text_list <- strsplit(text, ' ')
  text_list <- lapply(text_list, '[', 1:2)
  text_list <- lapply(text_list, paste, collapse = ' ')
  text_vector <- unlist(text_list)
  return (text_vector)
}

rasterizeSpecies <- function(species,
                             shape,
                             raster,
                             buffer = NULL,
                             folder = 'bats/range_buff/') {
  
  # first buffer, then rasterize IUCN range map for species
  shp <- shape[shape$BINOMIAL == species, ]
  
  if (!is.null(buffer)) {
    # convert buffer from kilometres to decimal degrees, assuming at the equator
    buffer <- buffer / 111.32
    
    # buffer by this amount
    shp <- gBuffer(shp,
                   width = buffer)
  }
  
  # rasterize the shapefile 
  tmp <- rasterize(shp,
                   raster,
                   field = 1,
                   background = 0,
                   fun = 'first')
  
  writeRaster(tmp,
              filename = paste0('~/Z/zhi/ebola/',
                                folder,
                                '/',
                                gsub(' ', '_', species)),
              format = 'GTiff',
              overwrite = TRUE)
  
  rm(tmp)
  
  return (NULL)
}

tidySpecies <- function (filename, template) {
  # load a raster if it contains any of the species' range,
  # mask and resave it, else delete it
  tmp <- raster(filename)
  if (!is.na(maxValue(tmp)) && maxValue(tmp) == 1) {
    tmp <- mask(tmp,
                template)
    
    writeRaster(tmp,
                file = filename,
                overwrite = TRUE)
    
  } else {
    
    rm(tmp)
    
    file.remove(filename)
    
  }
  
  return (NULL)
  
}

subsamplePolys <- function (data, ...) {
  # given a presence-background dataset, with multiple rows for some of the
  # occurrence records, subset it so that there's only one randomly selected
  # point from each polygon and then take a bootstrap it using `subsample`.
  # Dots argument is passed to subsample.
  
  # index for background records (outbreak id = 0)
  bg_idx <- data$outbreak_id == 0
  
  # subset to get occurrence section only
  occ <- data[!bg_idx, ]
  
  # get the different outbreaks
  u <- unique(occ$outbreak_id)
  
  # loop through, picking an index for each based on the number available
  occ_idx <- sapply(u,
                    function (id, occ) {
                      idx <- which(occ$outbreak_id == id)
                      sample(idx, 1)
                    },
                    occ)
  
  # get the subsetted dataset
  dat <- rbind(occ[occ_idx, ],
               data[bg_idx, ])
  
  # randomly subsample the dataset
  ans <- subsample(dat,
                   n = nrow(dat),
                   ...)
  
  # remove the outbreak ID column
  ans <- ans[, -which(names(ans) == 'outbreak_id')]
  
  return (ans)
}

# change the polygon IDs of an SPDF so it can be rbinded to something else
# nicked from
# http://gis.stackexchange.com/questions/32732/proper-way-to-rbind-spatialpolygonsdataframes-with-identical-polygon-ids
makeUniform <- function (SPDF) {
  pref <- substitute(SPDF)  #just putting the file name in front.
  newSPDF <- spChFIDs(SPDF,
                      as.character(paste(pref,
                                         rownames(as(SPDF,
                                                     "data.frame")),
                                         sep = "_")))
  return (newSPDF)
}

summarizeStats <- function (path) {
  
  # load validation stats
  stats <- read.csv(paste0(path, 'stats.csv'),
                    row.names = 1)
  
  auc  <- c(as.character(round(mean(stats$auc,
                                    na.rm = TRUE),
                               2)),
            as.character(round(sd(stats$auc,
                                  na.rm = TRUE),
                               2)))
  
  # load relative influence stats
  relinf <- read.csv(paste0(path,
                            'relative_influence.csv'),
                     stringsAsFactors = FALSE)
  
  ans <- c(auc_mean = auc[1],
           auc_sd = auc[2],
           cov1 = relinf[1, 1],
           relinf1 = as.character(round(relinf[1, 2],
                                        1)),
           cov2 = relinf[2, 1],
           relinf2 = as.character(round(relinf[2, 2],
                                        1)),
           cov3 = relinf[3, 1],
           relinf3 = as.character(round(relinf[3, 2],
                                        1)),
           cov4 = relinf[4, 1],
           relinf4 = as.character(round(relinf[4, 2],
                                        1)),
           cov5 = relinf[5, 1],
           relinf5 = as.character(round(relinf[5, 2],
                                        1)))
  
  return (ans)
  
}

thresholdRisk <- function (risk_raster,
                           occ,
                           proportion = 1) {
  
  # given a raster `risk_raster` giving risk on the (0,1] level and a 2d
  # dataframe `occ` with columns named 'lat' and 'long' giving the latitudes
  # and longitudes of known occurrence records,
  # find the threshold value so that `proportion` fraction of the records
  # fall in areas classified as 'at risk' and return the thresholded map
  
  # extract risk values for the occurrence data
  occ_risk <- extract(risk_raster[[1]],
                      occ[, c('long', 'lat')])
  
  # remove any missing data
  occ_risk <- na.omit(occ_risk)
  
  # get the relevant quantile
  thresh <- quantile(occ_risk,
                     1 - proportion,
                     na.rm = TRUE)
  
  # classify the raster
  at_risk_raster <- risk_raster > thresh
  
  # return this
  return (at_risk_raster)
  
}

# fit an ebola niche model

# load packages
library(seegSDM)
library(snowfall)

# load functions 


# ~~~~~~~~~~
# load the data

# raster covariates
covariates5km <- load('covariates5km.RData')
e <- raster::extent(covariates5km)
plot(covariates5km)

bat<- raster("primary_mean.tif")
bat<- crop(bat,e)
bat <- resample(bat,covariates5km)
covariates5km <- addLayer(covariates5km, bat)

BA <- raster("/Volumes/Soushie/results_clim/binomial_icar_pred.tif")
BA<- crop(BA,e)
bushmeat <- resample(BA,covariates5km, method = "bilinear")
covariates5km <- addLayer(covariates5km, bushmeat)

#bushmeat_mask <- calc(bushmeat, fun=function(x){ x[x > 0] <- -999; return(x)} )


# occurrence data
occ <- read.csv("occurrence.csv")


# generate pseudo-absence data according to the population surface
# (assume perfect detection of human cases, so reported risk only
# a function of population density)
pop <- raster("pop.tif")
pop<- crop(pop, e)
set.seed(10)
bg <- bgSample(pop,
               n = 10000,
               prob = TRUE,
               replace = TRUE,
               spatial = FALSE)

colnames(bg) <- c('long', 'lat')
bg <- data.frame(bg)

# add an outbreak id to this
bg$outbreak_id <- 0

# combine the occurrence and background records
dat <- rbind(cbind(PA = rep(1, nrow(occ)),
                   occ[, c('long', 'lat', 'outbreak_id')]),
             cbind(PA = rep(0, nrow(bg)),
                   bg))

# get the covariate values
dat_covs <- raster::extract(covariates5kmBA, dat[, 2:3])

# and add them
dat_all<- cbind(dat, dat_covs)

# remove NAs
dat_all <- na.omit(dat_all)


ncpu <- 50
nboot <- ncpu * 1


# create a list with random permutations of dat_all, sampling one occurrence
# from each polygon in each iteration, the bootstrapping as usual.
# This way there's more jitter and it avoids the horrible normalization in
# gbm which otherwise gives some points a weight of ~250 (!).
data_list <- replicate(nboot,
                       subsamplePolys(dat_all,
                                      minimum = c(5, 5)),
                       simplify = FALSE)


# initialize the cluster
sfInit(parallel = TRUE, cpus = ncpu)
sfLibrary(seegSDM)

model_list <- sfLapply(data_list,
                       runBRT,
                       gbm.x = 4:ncol(data_list[[1]]),
                       gbm.y = 1,
                       pred.raster = covariates5km,
                       gbm.coords = 2:3,
                       wt = function(PA) ifelse(PA == 1, 1, sum(PA) / sum(1 - PA)))

# get cv statistics in parallel
stat_lis <- sfLapply(model_list, getStats)

# stop the cluster
sfStop()

# convert the stats list into a matrix using the do.call function
stats <- do.call("rbind", stat_lis)

# save them
write.csv(stats,"stats.csv")

# save the relative influence scores
relinf <- getRelInf(model_list)
write.csv(relinf,'relative_influence.csv')