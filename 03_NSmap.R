## defining weights
Wmammala <- raster("covars_5k/all_mammals.tif")
Wpop <- raster("covars_5k/pop_density.tif")
## load uncertainty raster
uncertainty <-  raster("uncertainty2.tif")
e<- raster::extent(uncertainty)
Wpop <- crop(Wpop, e)
Wpop <- resample(Wpop, uncertainty)
Wmammala <- resample(Wmammala, uncertainty)
#Wpop_log <- log10(Wpop)
#plot(Wpop_log)
#pop_log[Wpop_log==-Inf] <- NA

## calculate necessity for additional surveillance (standardized from 0 to 1)

NS <- (uncertainty *Wpop*Wmammala)
NS <-climateStability::rescale0to1(NS)
NS<- (NS*20)/10
plot(NS)
writeRaster(NS, "NS.tif",format = 'GTiff', overwrite = T)

#loop following for 100 (4*100) iterations
idx = which.max(NS)
xy1 = xyFromCell(NS,idx)
colnames(xy1) <- c('x', 'y')
xy1 <- data.frame(xy1)
coordinates(xy1) <-~x+y
extract(NS, xy1@coords, buffer=100000, fun=mean)
extract(NS, xy1@coords)
b = circles(xy1, d=100000, lonlat=T)
m <- polygons(b)
b1 = circles(xy1, d=75000, lonlat=T)
m1 <- polygons(b1)
b2 = circles(xy1, d=50000, lonlat=T)
m2 <- polygons(b2)
b3 = circles(xy1, d=25000, lonlat=T)
m3 <- polygons(b3)
d1 <- gDifference(m, m1)
d2 <- gDifference(m1, m2)
d3 <- gDifference(m2, m3)


r <- rasterize(d1, NS, mask = TRUE)
r<- r*1
NS <-merge(r,NS)

r1 <- rasterize(d2, NS, mask = TRUE)
r1 <-r1*0.75
NS <-merge(r1, NS)

r2 <- rasterize(d2, NS, mask = TRUE)
r2 <-r2*0.50
NS <-merge(r2,NS)

r3 <- rasterize(d3, NS, mask = TRUE)
r3 <- r3*0.25
NS <-merge(r3,NS)

r4 <- rasterize(m3, NS, mask = TRUE)
r4 <- r4*0
NS <-merge(r4,NS)
xy1@coords
