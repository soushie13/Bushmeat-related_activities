
source("setup.R")
### open covariates data##
covars <- load("corvariates.Rdata")

### open presence data##
bm_markets <- read.csv("bushmeat_ruralxy.csv", header = TRUE)
obs.data <- bm_markets[, c("x", "y")]
plot(obs.data, pch = 20)

### define study extent ##
e <- extent(covars)
##create the background##
data(wrld_simpl)
bw <- crop(wrld_simpl, e)
plot(bw)

##===============================================
##
## 1. Random Forest
##
##===============================================
### randon forest###
library(randomForest)
###RF training
#extract data
coordinates(obs.data)<-~x+y
p <- obs.data
presence = p@coords
##create the pseudo-absences##
library(seegSDM)
pop <- raster("pop_density.tif")
pop<- crop(pop, e)
set.seed(10)
bg <- bgSample(pop,
               n = 1000,
               prob = TRUE,
               replace = TRUE,
               spatial = FALSE)

colnames(bg) <- c('x', 'y')
absence<- data.frame(bg)

### random forest  with spatial cross-validation
library(SDMtune)
covars <- stack(covars)
# Create SWD object
swd_data <- prepareSWD(species = "rural locations", p = presence, a = absence,
                   env = covars)
# Split presence locations in training (80%) and testing (20%) datasets
datasets <- trainValTest(swd_data, test = 0.2, seed = 25)
train <- datasets[[1]]
test <- datasets[[2]]

rf_model <- train(method = "RF", data = train, ntree=1000, nodesize=10, importance =T)
auc(rf_model)
tss(rf_model)
library(plotROC)
plotROC(rf_model)
map <- predict(rf_model, data = covars, type = "cloglog")
plot(map)

library(zeallot)  # For unpacking assignment
c(train, test) %<-% trainValTest(swd_data, test = 0.2, only_presence = FALSE, seed = 25)
rf_model <- train(method = "RF", data = train, ntree=1000, nodesize=10, importance =T)
cat("Training auc: ", auc(rf_model))
cat("Testing auc: ", auc(rf_model, test = test))

output <- data.frame(matrix(NA, nrow = 10, ncol = 3)) # Create an empty data.frame
colnames(output) <- c("seed", "trainAUC", "testAUC")
set.seed(25)
seeds <- sample.int(1000, 10) # Create 10 different random seeds
for (i in 1:length(seeds)) { # Loop through the seeds
  c(train, test) %<-% trainValTest(swd_data, test = 0.2, seed = seeds[i]) # Make the train/test split
  m <- train("RF", data = train) # train the model
  # Populate the output data.frame
  output[i, 1] <- seeds[i]
  output[i, 2] <- auc(m)
  output[i, 3] <- auc(m, test = test)
}
# Print the output
output
# compute the range of the testing AUC
range(output[, 3])
#Cross Validation
folds <- randomFolds(swd_data, k = 4,  seed = 25)
cv_model <- train("RF", data = swd_data, folds = folds)
cv_model
cat("Training AUC: ", auc(cv_model))
cat("Testing AUC: ", auc(cv_model, test = TRUE))
cat("Training TSS: ", tss(cv_model))
cat("Testing TSS: ", tss(cv_model, test = TRUE))
#Spatial Cross Validation
# Create spatial points data frame
library(blockCV)
sp_df <- SpatialPointsDataFrame(swd_data@coords, data = as.data.frame(swd_data@pa), proj4string = crs(covars))
e_folds <- envBlock(rasterLayer = covars, speciesData = sp_df, species = "swd_data@pa", k = 4, standardization = "standard", rasterBlock = FALSE, numLimit = 100)
scv_model <- train(method = "RF", data = swd_data, fc = "l", reg = 0.8, folds = e_folds)
cat("Training AUC: ", auc(scv_model))
cat("Testing AUC: ", auc(scv_model, test = TRUE))
cat("Training TSS: ", tss(scv_model))
cat("Testing TSS: ", tss(scv_model, test = TRUE))
plotResponse(cv_model, var = "def_1k", type = "cloglog", marginal = TRUE, fun = mean, rug = TRUE)
#tune the model hyperparameters
c(train, val, test) %<-% trainValTest(swd_data, val = 0.2, test = 0.2,  seed = 61516)
cat("# Training  : ", nrow(train@data))
cat("# Validation: ", nrow(val@data))
cat("# Testing   : ", nrow(test@data))
rf_model <- train("RF", data = train)

getTunableArgs(rf_model)
#[1] "mtry"     "ntree"    "nodesize"
# Define the values for bg
h <- list(ntree=seq(500))
# Call the gridSearch function
exp_1 <- gridSearch(rf_model, hypers = h, metric = "auc", test = val)
exp_1 <- gridSearch(rf_model, hypers = h, metric = "tss", test = val)
#Random search
h <- list(nodesize = seq(10), ntree = seq(500))
exp_6 <- randomSearch(rf_model, hypers = h, metric = "auc", test = val, pop = 10, seed = 65466)
exp_6 <- randomSearch(rf_model, hypers = h, metric = "tss", test = val, pop = 10, seed = 65466)
exp_6@results
exp_7 <- optimizeModel(rf_model, hypers = h, metric = "tss", test = val, pop = 15, gen = 2, keep_best = 0.4, keep_random = 0.2, mutation_chance = 0.4, seed = 798)

## fit final random forest model
#merge training and validation data
merged_data <- mergeSWD(train, val)
set.seed(2022)
RFfinal_model <- train("RF", data = merged_data, ntree = 500, nodesize = 5)
auc(RFfinal_model, test = test)
tss(RFfinal_model, test = test)

#predict distribution
pr <- predict(RFfinal_model, data = covars, type = "cloglog")
##plot
par(mfrow=c(1,1))
plot(pr, main='Random Forest, regression')
plot(wrld_simpl, add=TRUE, border='dark grey')
##save prediction raster
writeRaster(pr, "rf_pred.tif",format = 'GTiff', overwrite = T)
plotROC(RFfinal_model)
library(kableExtra)
modelReport(RFfinal_model,type = "cloglog", test = test,
             jk = TRUE, permut = 10, folder = "RF")

varImp(RFfinal_model)
##===============================================
##
## 2. Boosted regression trees - Gradient Boosted Machine (GBM)
##
##===============================================
library(zeallot)      # data splitting 
library(gbm)          # basic implementation
library(xgboost)      # a faster implementation of gbm
library(caret)        # an aggregator package for performing many machine learning models


# Split presence locations in training (80%) and testing (20%) datasets
datasets <- trainValTest(swd_data, test = 0.2,  seed = 345)
train <- datasets[[1]]
test <- datasets[[2]]

brt_model <- train(method = "BRT", data = train, family = "bernoulli", tree.complexity = 5,
                   learning.rate = 0.01, bag.fraction = 0.5)
auc(brt_model)
tss(brt_model)


c(train, test) %<-% trainValTest(swd_data, test = 0.2, seed = 345)
brt_model <- train(method = "BRT", data = train, family = "bernoulli", tree.complexity = 5,
                   learning.rate = 0.01, bag.fraction = 0.5)
cat("Training auc: ", auc(brt_model))
cat("Testing auc: ", auc(brt_model, test = test))

output <- data.frame(matrix(NA, nrow = 10, ncol = 3)) # Create an empty data.frame
colnames(output) <- c("seed", "trainAUC", "testAUC")
set.seed(345)
seeds <- sample.int(1000, 10) # Create 10 different random seeds
for (i in 1:length(seeds)) { # Loop through the seeds
  c(train, test) %<-% trainValTest(swd_data, test = 0.2, seed = seeds[i]) # Make the train/test split
  m <- train("BRT", data = train) # train the model
  # Populate the output data.frame
  output[i, 1] <- seeds[i]
  output[i, 2] <- auc(m)
  output[i, 3] <- auc(m, test = test)
}
# Print the output
output
# compute the range of the testing AUC
range(output[, 3])
#Cross Validation
folds <- randomFolds(swd_data, k = 4,  seed = 345)
cv_model <- train("BRT", data = swd_data, folds = folds, family = "bernoulli", tree.complexity = 5,
                  learning.rate = 0.01, bag.fraction = 0.5)
cv_model
cat("Training AUC: ", auc(cv_model))
cat("Testing AUC: ", auc(cv_model, test = TRUE))
cat("Training TSS: ", tss(cv_model))
cat("Testing TSS: ", tss(cv_model, test = TRUE))
#Spatial Cross Validation
# Create spatial points data frame
sp_df <- SpatialPointsDataFrame(swd_data@coords, data = as.data.frame(swd_data@pa), proj4string = crs(covars))
e_folds <- envBlock(rasterLayer = covars, speciesData = sp_df, species = "swd_data@pa", k = 4, standardization = "standard", rasterBlock = FALSE, numLimit = 100)
scv_model <- train(method = "BRT", data = swd_data, folds = e_folds,family = "bernoulli", tree.complexity = 5,
                   learning.rate = 0.01, bag.fraction = 0.5)
cat("Training AUC: ", auc(scv_model))
cat("Testing AUC: ", auc(scv_model, test = TRUE))
cat("Training TSS: ", tss(scv_model))
cat("Testing TSS: ", tss(scv_model, test = TRUE))
#tune the model hyperparameters
c(train, val, test) %<-% trainValTest(swd_data, val = 0.2, test = 0.2,  seed = 345)
cat("# Training  : ", nrow(train@data))
cat("# Validation: ", nrow(val@data))
cat("# Testing   : ", nrow(test@data))
brt_model <- train("BRT", data = train,family = "bernoulli", tree.complexity = 5,
                  learning.rate = 0.0025, bag.fraction = 0.5)

getTunableArgs(brt_model)
# "distribution"      "n.trees"           "interaction.depth" "shrinkage"   "bag.fraction"
# Define the values for bg
h <- list(n.trees=seq(500))
# Call the gridSearch function
exp_1 <- gridSearch(brt_model, hypers = h, metric = "auc", test = val)
exp_1 <- gridSearch(brt_model, hypers = h, metric = "tss", test = val)
exp_1@results
#Random search
h <- list(n.trees=seq(500), shrinkage=seq(0.5))
exp_6 <- randomSearch(brt_model, hypers = h, metric = "auc", test = val, pop = 10, seed = 345)
exp_6 <- randomSearch(brt_model, hypers = h, metric = "tss", test = val, pop = 10, seed = 345)
exp_6@results

## fit final BRT model
#merge training and validation data
merged_data <- mergeSWD(train, val)
set.seed(2022)
brtfinal_model <- train("BRT", data = merged_data, n.trees=500, learning.rate = 0.0025, bag.fraction = 0.5)
auc(brtfinal_model, test = test)
tss(brtfinal_model, test = test)
#predict distribution
p <- predict(brtfinal_model, data = covars, type = "cloglog")

##plot
plot(p, main='BRT prediction')
##save prediction raster
writeRaster(p, "brt_pred.tif",format = 'GTiff', overwrite = T)
library(kableExtra)
plotROC(brtfinal_model)
modelReport(brtfinal_model,folder = "BRT",type = "cloglog", test = test,
            response_curves = TRUE,jk = TRUE, permut = 10)

varImp(brtfinal_model)
##===============================================
##
## 3. Maxent
##
##===============================================
library(SDMtune)
library(rJava)
dismo::maxent()

covars <- stack(covars)
# Create SWD object
swd_data <- prepareSWD(species = "rural locations", p = presence, a = absence,
                       env = covars)

# Split presence locations in training (80%) and testing (20%) datasets
datasets <- trainValTest(swd_data, test = 0.2,  seed = 111)
train <- datasets[[1]]
test <- datasets[[2]]
maxent_model <- train(method = "Maxent", data = train, fc = "lqph", reg = 1, iter = 500)
auc(maxent_model)
tss(maxent_model)

library(zeallot)  # For unpacking assignment
c(train, test) %<-% trainValTest(swd_data, test = 0.2, seed = 111)
maxent_model <- train(method = "Maxent", data = train, fc = "lqph", reg = 1, iter = 500)
cat("Training auc: ", auc(maxent_model))
cat("Testing tss: ", auc(maxent_model, test = test))

output <- data.frame(matrix(NA, nrow = 10, ncol = 3)) # Create an empty data.frame
colnames(output) <- c("seed", "trainAUC", "testAUC")
set.seed(111)
seeds <- sample.int(1000, 10) # Create 10 different random seeds
for (i in 1:length(seeds)) { # Loop through the seeds
  c(train, test) %<-% trainValTest(swd_data, test = 0.2, seed = seeds[i]) # Make the train/test split
  m <- train("Maxent", data = train) # train the model
  # Populate the output data.frame
  output[i, 1] <- seeds[i]
  output[i, 2] <- auc(m)
  output[i, 3] <- auc(m, test = test)
}
# Print the output
output
# compute the range of the testing AUC
range(output[, 3])
#Cross Validation
folds <- randomFolds(swd_data, k = 4,  seed = 111)
cv_model <- train("Maxent", data = swd_data, folds = folds,  
                   fc = "lqph", reg = 1, iter = 500)
cv_model
cat("Training AUC: ", auc(cv_model))
cat("Testing AUC: ", auc(cv_model, test = TRUE))
cat("Training TSS: ", tss(cv_model))
cat("Testing TSS: ", tss(cv_model, test = TRUE))
#Spatial Cross Validation
# Create spatial points data frame
library(blockCV)
sp_df <- SpatialPointsDataFrame(swd_data@coords, data = as.data.frame(swd_data@pa), proj4string = crs(covars))
e_folds <- envBlock(rasterLayer = covars, speciesData = sp_df, species = "swd_data@pa", k = 4, standardization = "standard", rasterBlock = FALSE, numLimit = 100)
scv_model <- train(method = "Maxent", data = swd_data, folds = e_folds,fc = "lp", reg = 1, iter = 1000)
cat("Training AUC: ", auc(scv_model))
cat("Testing AUC: ", auc(scv_model, test = TRUE))
cat("Training TSS: ", tss(scv_model))
cat("Testing TSS: ", tss(scv_model, test = TRUE))
plotResponse(cv_model, var = "def_1k", type = "cloglog", marginal = TRUE, fun = mean, rug = TRUE)
#tune the model hyperparameters
library(zeallot)
c(train, val, test) %<-% trainValTest(swd_data, val = 0.2, test = 0.2,  seed = 111)
cat("# Training  : ", nrow(train@data))
cat("# Validation: ", nrow(val@data))
cat("# Testing   : ", nrow(test@data))
maxent_model <- train("Maxent", data = train)

getTunableArgs(maxent_model)
# "fc"   "reg"  "iter"
# Define the values for bg
h <- list(reg = seq(0.2, 1, 0.1))
# Call the gridSearch function
exp_1 <- gridSearch(maxent_model,hypers = h, metric = "auc", test = val)
exp_1 <- gridSearch(maxent_model,hypers = h, metric = "tss", test = val)
h <- list(reg = seq(0.2, 2, 0.2), fc = c("l", "lq", "lh", "lqp", "lqph", "lqpht"))
exp_5 <- gridSearch(maxent_model, hypers = h, metric = "auc", test = val)
exp_5 <- gridSearch(maxent_model, hypers = h, metric = "tss", test = val)
#Random search
h <- list(reg = seq(0.2, 5, 0.2), fc = c("l", "lq", "lh", "lp", "lqp", "lqph"))
exp_6 <- randomSearch(maxent_model, hypers = h, metric = "auc", test = val, pop = 10, seed = 111)
exp_6 <- randomSearch(maxent_model, hypers = h, metric = "tss", test = val, pop = 10, seed = 111)
exp_6@results
exp_7 <- optimizeModel(maxent_model, hypers = h, metric = "auc", test = val, pop = 15, gen = 2, keep_best = 0.4, keep_random = 0.2, mutation_chance = 0.4, seed = 798)
## fit final Maxent model
#merge training and validation data
merged_data <- mergeSWD(train, val)
set.seed(2022)
maxentfinal_model <- train("Maxent", data = merged_data, fc ="lp", reg = 4.6)
auc(maxentfinal_model, test = test)
tss(maxentfinal_model, test = test)
#predict distribution
pm  <- predict(maxentfinal_model, data = covars, type = "cloglog")

##plot
plot(pm, main='Maxent prediction')
##save prediction raster
writeRaster(pm, "maxent_pred.tif",format = 'GTiff', overwrite = T)

maxentVarImp(maxentfinal_model)
library(kableExtra)
modelReport(maxentfinal_model,type = "cloglog", folder = "MaxENT", test = test,
            response_curves = TRUE, jk = TRUE, permut = 10)
##===============================================
##
## 4. BART
##
##===============================================

library(embarcadero) 
set.seed(12345)
pb_xy <- rbind(presence, absence)
pb_pred <- raster::extract(covars, pb_xy)
pb <- c(rep(1, nrow(presence)), rep(0, nrow(absence)))
data <- data.frame( cbind(pa=pb, pb_pred) )
data <- na.omit(data)
xvars <- names(data)[!names(data)=="pa"]

BART.model <- bart.step(x.data= data[,xvars], 
                        y.data= data[,"pa"], 
                        full=TRUE,
                        quiet = TRUE)
# Saving in RData format
save(BART.model, file = "bart.RData")
summary(BART.model)


bart.map <- predict2.bart(object= BART.model,
                          x.layers = covars,
                          quantiles = c(0.025, 0.975),
                          splitby=20)
plot(bart.map[[1]])
plot(bart.map[[3]])
q <- quantile(values(bart.map[[3]]- bart.map[[2]]), 0.75, na.rm=TRUE)
uncertinity <- (bart.map[[3]]- bart.map[[2]]) > q
sp <- spartial(bart.map, preds_10, x.vars="bushmeat", equal=TRUE)
plot(uncertinity)
writeRaster(bart.map[[1]], "bart.tif" ,format = 'GTiff', overwrite = T)
writeRaster(bart.map[[2]], "bart2.5.tif" ,format = 'GTiff', overwrite = T)
writeRaster(bart.map[[3]], "bart97.5.tif" ,format = 'GTiff', overwrite = T)
writeRaster(uncertinity, "uncertainitybart.tif" ,format = 'GTiff', overwrite = T)


###geographical null model using geographical distance
null_model2 <- convHull(pres_train, lonlat =TRUE)
null_e <- evaluate(null_model2, p= pres_test, a=backg_test)

##================================================================================
##
## 5. Ensemble model - binomial regression with hierarchical Bayesian framework
##
##================================================================================
#extract data
coordinates(obs.data)<-~x+y
p <- obs.data
presence = p@coords
presvals <- extract(rreds, presence, cellnumber=TRUE)
prevals_coords =cbind(presence,presvals)
##create the pseudo-absences##
set.seed(10)
background <- spsample(bw,n= 1000,"random", iter= 5000) 
plot(background)
absence = background@coords
absvals <- extract(preds, absence, cellnumber=TRUE)
absvals_coords =cbind(absence,absvals)
pb <- c(rep(1, nrow(prevals_coords)), rep(0, nrow(absvals_coords)))
pa <- data.frame(cbind(pb, rbind(prevals_coords, absvals_coords)))

## Extract environmental values and cell number for observations
pa$Presences <- pa$pb
pa$Trials <- c(1)
## omit rows with missing data
pa.cc=na.omit(pa)


## Normalized continuous covariates
pa.norm <- pa.cc
Mean <- vector()
Sd <- vector()
for (i in c(5:10)) {
  m <- mean(pa.cc[,i],na.rm=TRUE)
  s <- sd(pa.cc[,i],na.rm=TRUE)
  Mean <- c(Mean,m)
  Sd <- c(Sd,s)
  pa.norm[,i] <- (pa.cc[,i]-m)/s
}
## Data-frame with mean and sd for each variable
df.mean.sd <- as.data.frame(rbind(Mean,Sd))
names(df.mean.sd) <- names(pa.norm)[c(5:10)]

## Raster stack for predictions (with normalized covariates)
env <- covars

for (i in c(5:10)) {
  var.name <- names(pa.norm)[i] ## Variable name
  w <- which(names(env)==var.name) ## Position in the stack 
  m <- df.mean.sd[1,var.name] ## Mean
  s <- df.mean.sd[2,var.name] ## Sd
  orig <- values(subset(env,w)) ## Original values
  trans <- (orig-m)/s ## Transformed values
  env[[w]][] <- trans
}

## Select only grid cells with no NA
env.df.pred <- as.matrix(env)
w <- complete.cases(env.df.pred) ## Note: w will be used to obtain the cell identifier for predictions in iCAR model
env.df.pred.complete <- as.data.frame(env.df.pred[w,])

## Make a cluster for parallel MCMCs
nchains <- 2
ncores <- nchains ## One core for each MCMC chains
cores<-detectcores()
clust <- makeCluster(2)
registerDoParallel(clust)

## Starting values and random seed
seed <- 1234
set.seed(seed)
beta.start <- runif(nchains,-1,1)
gamma.start <- runif(nchains,-1,1)
Vrho.start <- runif(nchains,0,10)
seed.mcmc <- round(runif(nchains,0,1e6))
pa.norm$Trials <- c(1)
pa.norm$Presences <- pa.norm$pb

##=============================================================================
##
## 5.1. Binomial regression
##=============================================================================

mod.binomial <- foreach (i=1:nchains, .packages="hSDM") %dopar% {
  mod <- hSDM.binomial(presences=pa.norm$Presences,
                       trials=pa.norm$Trials,
                       suitability=~ RandomForest + BRT + BART,
                       data=pa.norm,
                       suitability.pred=env.df.pred.complete,
                       burnin=2000,
                       mcmc=2000, thin=5,
                       beta.start=beta.start[i],
                       mubeta=0, Vbeta=1.0E6,
                       seed=seed.mcmc[i], verbose=1,
                       save.p=0)
  return(mod)
}

## Extract list of MCMCs from output
binomial.env.mcmc <- mcmc.list(lapply(mod.binomial,"[[","mcmc"))

sink(file="Results/binomial_mcmc_summary.txt")
summary(binomial.env.mcmc)
sink()
## Outputs summary
bionomial.env.stat <- summary(binomial.env.mcmc)$statistics
sink(file="Results/binomial_mcmc_summary.txt")
summary(binomial.env.mcmc)
cat(rep("\n",3))
gelman.diag(binomial.env.mcmc)
sink()
## Deviance
deviance.bionomial.env <- bionomial.env.stat["Deviance","Mean"]

## Plot trace and posterior distributions
pdf("Results/bionomial_mcmc_trace.pdf")
plot(binomial.env.mcmc)
dev.off()

## Prediction on the landscape
prob.p.bi <- subset(covars,1) ## create a raster for predictions
values(prob.p.bi)[w] <- mod.binomial[[1]]$theta.pred ## assign predicted values
values(prob.p.bi)[!w] <- NA ## set NA where no environmental data
## Plot the predictions
pdf(file="Results/binomialpredictions.pdf")
plot(prob.p.bi)
plot(pa.norm[pa.norm$pb==0,],pch=".",col=grey(0.5),add=TRUE)
plot(pa.norm[pa.norm$pb>0,],pch=3,add=TRUE)
dev.off()

## Export the results as GeoTIFF
writeRaster(prob.p.bi,filename="Binomial_pred.tif",overwrite=TRUE)



##=============================================================================
##
## 5.2. Binomial with iCAR (spatial autocorrelation)  
##=============================================================================


## Landscape and neighbors
ncells <- ncell(preds)
neighbors.mat <- adjacent(preds, cells=c(1:ncells), directions=8, pairs=TRUE, sorted=TRUE)
n.neighbors <- as.data.frame(table(as.factor(neighbors.mat[,1])))[,2]
adj <- neighbors.mat[,2]
cells.pred <- which(w) ## Vector w indicates the cells with environmental information (without NA)

## binomial icar model
## hSDM model using Binomial icar for perfect detection
mod.binomial.icar <- foreach (i=1:nchains, .packages="hSDM") %dopar% {
  mod <- hSDM.binomial.iCAR(presences=pa.norm$Presences,
                            trials=pa.norm$Trials,
                            suitability=~ RandomForest + BRT + BART,
                            data=pa.norm,
                            ## Spatial structure
                            spatial.entity=pa.norm$cells,
                            n.neighbors=n.neighbors,
                            neighbors=adj,
                            suitability.pred=env.df.pred.complete,
                            spatial.entity.pred=cells.pred,
                            burnin=2000,
                            mcmc=2000, thin=5,
                            beta.start=beta.start[i],
                            Vrho.start=Vrho.start[i],
                            ## Priors
                            priorVrho="Uniform",
                            mubeta=0, Vbeta=1.0E6,
                            Vrho.max=10,
                            seed=seed.mcmc[i], verbose=1,
                            save.p=1) # save probability per pixel
  return(mod)
}

## Extract list of MCMCs from output
binomial.icar.mcmc <- mcmc.list(lapply(mod.binomial.icar,"[[","mcmc"))

## Outputs summary
bionomial.icar.stat <- summary(binomial.icar.mcmc)$statistics
sink(file="results/binomial.icar_mcmc_summary.txt")
summary(binomial.icar.mcmc)
cat(rep("\n",3))
gelman.diag(binomial.icar.mcmc)
sink()
## Deviance
deviance.bionomial.icar <- bionomial.icar.stat["Deviance","Mean"]

## Plot trace and posterior distributions
pdf("results/bionomial.icar_mcmc_trace.pdf")
plot(binomial.icar.mcmc)
dev.off()
## Spatial random effects
rho <- subset(covars,1) ## create a raster
values(rho) <- mod.binomial.icar[[1]]$rho.pred
pdf(file="results/binomial.iCAR_random_effects.pdf")
plot(rho)
dev.off()
## Prediction on the landscape
prob.p.b <- subset(covars,1) ## create a raster for predictions
values(prob.p.b)[w] <- mod.binomial.icar[[1]]$theta.pred ## assign predicted values
values(prob.p.b)[!w] <- NA ## set NA where no environmental data
## Plot the predictions
pdf(file="Results/binomial.icar_predictions.pdf")
plot(prob.p.b)
plot(pa.norm[pa.norm$pb==0,],pch=".",col=grey(0.5),add=TRUE)
plot(pa.norm[pa.norm$pb>0,],pch=3,add=TRUE)
dev.off()

## Export the results as GeoTIFF
writeRaster(prob.p.b,filename="Results/binomial_icar_pred.tif",overwrite=TRUE)

## Extract predicted probability of presence
pa$prob.p <- extract(prob.p.b,pa.obs.data)



##===============================================
##
## Model comparison based on deviance
##
##===============================================

## Null model
mod.ZIB.null <- foreach (i=1:nchains, .packages="hSDM") %dopar% {
  mod <- hSDM.ZIB(presences=pa.norm$Presences,
                  trials=pa.norm$Trials,
                  suitability=~1,
                  observability=~1,
                  data=pa.norm,
                  suitability.pred=env.df.pred.complete,
                  burnin=2000,
                  mcmc=2000, thin=5,
                  beta.start=beta.start[i],
                  gamma.start=gamma.start[i],
                  mubeta=0, Vbeta=1.0E6,
                  mugamma=0, Vgamma=1.0E6,
                  seed=seed.mcmc[i], verbose=1,
                  save.p=0)
  return(mod)
}
## Stop cluster
stopCluster(clust)

## Extract list of MCMCs from output
ZIB.null.mcmc <- mcmc.list(lapply(mod.ZIB.null,"[[","mcmc"))

## Deviance
ZIB.null.stat <- summary(ZIB.null.mcmc)$statistics
deviance.null <- ZIB.null.stat["Deviance","Mean"]


prob.p.null <- subset(covars,1) ## create a raster for predictions
values(prob.p.null)[w] <- mod.ZIB.null[[1]]$prob.p.pred ## assign predicted values
values(prob.p.null)[!w] <- NA 
plot(prob.p.null)

## Extract predicted probability of presence
pa$prob.p <- extract(prob.p.null,pa.obs.data)


## TSS as a function of the threshold
OA <- vector()
TSS <- vector()
thresh.seq <- seq(0,1,by=0.01)
for (i in 1:length(thresh.seq)) {
  Index <- Index.fun(pa$Suit,pa$prob.p,thresh.seq[i])
  OA[i] <- Index$OA
  TSS[i] <- Index$TSS
}

## Identifying the threshold maximizing TSS
maxTSS <- max(TSS,na.rm=TRUE)
maxTSS








