rastlist <- list.files(path = "covars_5k", pattern='.tif$',  full.names=TRUE)
covars <- stack(rastlist)
plot(covars)
res(covars)
tmax <- getData(name = "worldclim",var = "tmax", res = 2.5, path = "covars_5k")
tmax_m <- mean(tmax)
res(tmax)
e <- extent(covars)
tmax_m <- crop(tmax_m, e)
writeRaster(tmax_m, "covars_5k/tmax_m.tif",format = 'GTiff', overwrite = T)
tmin <- getData(name = "worldclim",var = "tmin", res = 2.5, path = "covars_5k")
tmin_m <- mean(tmin)
tmin_m <- crop(tmin_m, e)
writeRaster(tmin_m, "covars_5k/tmin_m.tif",format = 'GTiff', overwrite = T)
ppt <- getData(name = "worldclim",var = "prec", res = 2.5, path = "covars_5k")
ppt_m <- mean(ppt)
ppt_m <- crop(ppt_m, e)
writeRaster(ppt_m, "covars_5k/ppt_m.tif",format = 'GTiff', overwrite = T)
clim = stack(tmax_m,tmin_m,ppt_m)
names(clim) = c("tmax","tmin","ppt")

covars <- stack(clim, covars)
plot(covars)

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

## training and testing data ##
set.seed(0)
k<-4
group <- kfold(obs.data, 4)

pres_train <- obs.data[group != 1, ]
pres_test <- obs.data[group == 1, ]

##create the pseudo-absences##
set.seed(10)
background <- spsample(bw,n= 1000,"random", iter= 5000) 
plot(background)
background <- background@coords
group <- kfold(background, 4)
backg_train <- background[group != 1, ]
backg_test <- background[group == 1, ]

#extract data
train <- rbind(pres_train, backg_train)
pb_train <- c(rep(1, nrow(pres_train)), rep(0, nrow(backg_train)))
envtrain <- extract(covars, train)
envtrain <- data.frame( cbind(pa=pb_train, envtrain) )
envtrain <- na.omit(envtrain)

### check covariate correlation #
correlation.matrix<-cor(envtrain, method=c("spearman"))

hc <- hclust(as.dist(1-abs(correlation.matrix)),
             method="centroid")
par(mfrow=c(1,1))
plot (hc, sub="", xlab="Variables", hang=-1, ylab="1-
     abs(correlation)")
abline(0.3,0)
rect.hclust(hc,h=0.3)
##Variance Inflation Factor (VIF) 
#load the car library
library(car)
model <- lm(pa ~acc_50k+  bushmeat+pop_d+def_1k
            + proximity_raster_wdpa + all_mammals +tmin+ppt+tmax, data = envtrain)

#view the output of the regression model
summary(model)
#calculate the VIF for each predictor variable in the model
vif(model)
#create vector of VIF values
vif_values <- vif(model)

#create horizontal bar chart to display each VIF value
barplot(vif_values, main = "VIF Values", horiz = TRUE, col = "steelblue")

#add vertical line at 5
abline(v = 10, lwd = 3, lty = 2)
#gdp and tmax dropped as VIF > 10

covars <- dropLayer(covars, tmax)
plot(covars)