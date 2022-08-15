
#library
library(dismo)
library(rgeos)
library(maptools)
library(maps)
library(mapdata) 
library(randomForest)
library(readxl)
library(tidyverse)

setwd("D:/phymaturus")
phym <- read_excel("phym.xlsx")

# read files from pc
filename <- file.choose()
chi.adm1 <- readRDS(filename)
filename2 <- file.choose()
arg.adm1 <- readRDS(filename2)
chi.arg <- rbind(chi.adm1, arg.adm1)

bpath <- "D:\\phymaturus\\wc0.5\\bio_43" # define path to bio
# make a list of paths, only .bil files
raster_data<-list.files(path = bpath, pattern= "_43.bil")   
# read files as rasters
setwd(bpath) 
chi.all <- stack(raster_data)
nlayers(chi.all) # number of layers inthe rasterstack

# finally, set CRS
projection(chi.all) <- "+proj=longlat +datum=WGS84 +no_defs"

#crop layers into small frame
exx <- extent(-74,-70,-38.5,-36)
chi.all2 <- crop(chi.all, exx) # crop biovars
chi.arg2 <- crop(chi.arg,exx) # crop map 

#function rasterize
biolay <- raster::rasterize(chi.arg2, chi.all2,
                            fun='last', background=NA, mask=T)

names(biolay) <- gsub("_43","",names(biolay))

# variables in a 10,15 and 20% less precipitations scenario
# variables in a +1,1.5 and 2ºC scenario
names(biolay)
bio1f <- calc(biolay[["bio1"]],fun=function(x){x+15}) #Annual Mean Temperature
bio2f <- calc(biolay[["bio2"]],fun=function(x){x+15}) #Mean Diurnal Range (Mean of monthly (max temp - min temp))
bio8f <- calc(biolay[["bio8"]],fun=function(x){x+15}) #Mean Temperature of Wettest Quarter
bio10f <- calc(biolay[["bio10"]],fun=function(x){x+15})#Mean Temperature of Warmest Quarter 
bio11f <- calc(biolay[["bio11"]],fun=function(x){x+15})#Mean Temperature of Coldest Quarter
bio12f <- calc(biolay[["bio12"]],fun=function(x){x-(x*0.15)}) #Annual Precipitation
bio13f <- calc(biolay[["bio13"]],fun=function(x){x-(x*0.15)}) #Precipitation of Wettest Month
bio19f <- calc(biolay[["bio19"]],fun=function(x){x-(x*0.15)}) #Precipitation of Coldest Quarter

# this is the modifyed "future" brick (some variables (whit r) are let equal)
biolay.f <- biolay
biolay.f[["bio1"]] <- bio1f
biolay.f[["bio2"]] <- bio2f
biolay.f[["bio8"]] <- bio8f
biolay.f[["bio10"]] <- bio10f
biolay.f[["bio11"]] <- bio11f
biolay.f[["bio12"]] <- bio12f
biolay.f[["bio13"]] <- bio13f
biolay.f[["bio19"]] <- bio19f
names(biolay.f)

# loop for a series of calculations of a distribution model

library(geosphere)

# location of populations (data points)
phym.coor <- phym[,2:3]

set.seed(1024) # for repeteability
seedindex <- runif(15,1,1000)  # the first argument is the number of repetitions (15)

# empty lists and data frames to resume the loop output, 
# if the numner of repetitions change, here all "15's" should change accordingly
p.list <- vector("list",15)
p.tab <- data.frame(auc=rep(0,15),km2=rep(0,15))
f.list <- vector("list",15)
f.tab <- data.frame(auc=rep(0,15),km2=rep(0,15))

# Present -----------------------------------------------------------------

for(i in 1:15){
  set.seed(seedindex[i])
  group <- kfold(phym.coor, 3)
  pres_train <- phym.coor[group != 1, ]
  pres_test <- phym.coor[group == 1, ]
  
  # WARNING, USE JUST ONE!   background for training and testdata
  #for today data set use #1 -> biolay2
  #for future data set use #2 -> biolay.f
  
  backg <- randomPoints(biolay, n=100, extf = 0.95) #1
  #backg <- randomPoints(biolay.f, n=100, extf = 0.95) #2
  
  colnames(backg) = c("lon", "lat")
  group <- kfold(backg, 3)
  backg_train <- backg[group != 2, ]
  backg_test <- backg[group == 2, ]
  
  #------------------------------------------------------
  # PLOT FIRST RASTER TO CHECK BACKGROUND POINTS... 
  # -----------------------------------------------------
  # r = raster(biolay.f, 1)
  # plot(!is.na(r), col=c("white", "light grey"), legend=FALSE)
  # plot(ext, add=TRUE, col="red", lwd=2) # TO PLOT A SUBSQUARE AROUND THE SUBSETED SPACE
  # points(backg_train, pch="-", cex=0.5, col="yellow")
  # points(backg_test, pch="-", cex=0.5, col="black")
  # points(pres_train, pch="+", col="red")
  # points(pres_test, pch="+", col="blue")
  # ------------------------------------------------------
  
  # extract background data and ocurrence data from
  # layers and prepare data for model running
  
  pb_train <- c(rep(1, nrow(pres_train)), rep(0,nrow(backg_train))) 
  
  # warning, use only one... actual or future bkgr:
  envtrain.bk <- raster::extract(biolay, backg_train) # actual data set background
  #envtrain.bk <- extract(biolay.f, backg_train) # future data set background
  
  envtrain.pr <- raster::extract(biolay, pres_train) # actual climate presences
  envtrain <- rbind(envtrain.pr, envtrain.bk)
  envtrain <- data.frame(cbind(pa=pb_train, envtrain))
  head(envtrain)
  
  testpres <- data.frame( raster::extract(biolay, pres_test))
  colnames(testpres) <- names(biolay)
  
  # warning, use only one... actual or future bkgr:
  testbackg <- data.frame( raster::extract(biolay, backg_test)) #actual
  # testbackg <- data.frame( extract(biolay.f, backg_test)) #future
  
  
  colnames(testbackg) <- names(biolay)
  
  # random forest modelling 8vars
  library(randomForest)
  colnames(envtrain) <- c("pa", names(biolay))
  #first model uses pa as continuous for regression
  #model <- pa ~ bio1 + bio2 + bio3 + bio4 + bio5 + bio6 + bio7 + bio8 +
  #              bio9 + bio10 + bio11 + bio12 + bio13 + bio14 + bio15 + 
  #                bio16 + bio17 + bio18 + bio19
  
  #model <- pa ~ bio1 + bio2 + bio3 + bio4 + bio6 + bio8 + bio11 + bio12 + bio13 + bio15 + bio16 + bio18
  model <- pa ~ bio1 + bio2 + bio3 + bio4 + bio7 + bio8 + bio10 + bio11 + bio12 + bio13 + bio15 + bio19
  rf1 <- randomForest(model, data=envtrain)
  
  #model2 <- factor(pa) ~ bio2 + bio4 + bio5 + bio6 + bio12 + bio15 + bio18 + bio19
  #rf2 <- randomForest(model2, data=envtrain)
  
  erf1 <- evaluate(testpres, testbackg, rf1)
  #erf2 <- evaluate(testpres, testbackg, rf2)
  
  
  # WARNING  plot rf predictions (change model each time you need)
  pr <- predict(biolay, rf1)
  #pr2 <- predict(biolay.f, rf2)
  
  # extract the polygon from values above 0.85
  
  pr1 <- rasterToPolygons(pr,fun=function(x){x>0.85}) 
  pr1.a <- sum(areaPolygon(pr1))/1000000 # determine the area of the 0.85 polygon
  
  
  # WARNING change "f." of "p."  in the name of the 
  # objects depending on the desired model (future or present)
  p.list[[i]] <- pr
  p.tab[i,1] <- erf1@ auc
  p.tab[i,2] <- pr1.a
}

today <- brick(unlist(p.list)) #create a "raster brick" from all repetitions
today.mean <- mean(today) # find the mean value for each pixel
today.sd <- calc(today, fun = sd)

# mapa final today
ttmas <- today.mean+today.sd
ttmens <- today.mean-today.sd 

# Future ------------------------------------------------------------------

for(i in 1:15){
  set.seed(seedindex[i])
  group <- kfold(phym.coor, 3)
  pres_train <- phym.coor[group != 1, ]
  pres_test <- phym.coor[group == 1, ]
  
  # WARNING, USE JUST ONE!   background for training and testdata
  #for today data set use #1 -> biolay2
  #for future data set use #2 -> biolay.f
  
  #backg <- randomPoints(biolay, n=100, extf = 0.95) #1
  backg <- randomPoints(biolay.f, n=100, extf = 0.95, excludep = T) #2
  
  colnames(backg) = c("lon", "lat")
  group <- kfold(backg, 3)
  backg_train <- backg[group != 2, ]
  backg_test <- backg[group == 2, ]
  
  #------------------------------------------------------
  # PLOT FIRST RASTER TO CHECK BACKGROUND POINTS... 
  # -----------------------------------------------------
  # r = raster(biolay.f, 1)
  # plot(!is.na(r), col=c("white", "light grey"), legend=FALSE)
  # plot(ext, add=TRUE, col="red", lwd=2) # TO PLOT A SUBSQUARE AROUND THE SUBSETED SPACE
  # points(backg_train, pch="-", cex=0.5, col="yellow")
  # points(backg_test, pch="-", cex=0.5, col="black")
  # points(pres_train, pch="+", col="red")
  # points(pres_test, pch="+", col="blue")
  # ------------------------------------------------------
  
  # extract background data and ocurrence data from
  # layers and prepare data for model running
  
  pb_train <- c(rep(1, nrow(pres_train)), rep(0,nrow(backg_train))) 
  
  # warning, use only one... actual or future bkgr:
  #envtrain.bk <- raster::extract(biolay, backg_train) # actual data set background
  envtrain.bk <- raster::extract(biolay.f, backg_train) # future data set background
  
  envtrain.pr <- raster::extract(biolay, pres_train) # actual climate presences for train
  envtrain <- rbind(envtrain.pr, envtrain.bk)
  envtrain <- data.frame(cbind(pa = pb_train, envtrain))
  head(envtrain)
  
  testpres <- data.frame( raster::extract(biolay, pres_test)) #actual presences for test
  colnames(testpres) <- names(biolay)
  
  # warning, use only one... actual or future bkgr:
  # testbackg <- data.frame(raster::extract(biolay, backg_test)) #actual
  testbackg <- data.frame(raster::extract(biolay.f, backg_test)) #future
  
  
  colnames(testbackg) <- names(biolay)
  
  # random forest modelling 8vars
  #library(randomForest)
  colnames(envtrain) <- c("pa",names(biolay))
  #first model uses pa as continuous for regression
  #model2 <- pa ~ bio1 + bio2 + bio3 + bio4 + bio5 + bio6 + bio7 + bio8 +
  #  bio9 + bio10 + bio11 + bio12 + bio13 + bio14 + bio15 + 
  # bio16 + bio17 + bio18 + bio19
  
  #model2 <- pa ~ bio1 + bio2 + bio3 + bio4 + bio6 + bio8 + bio11 + bio12 + bio13 + bio15 + bio16 + bio18
  model2 <- pa ~ bio1 + bio2 + bio3 + bio4 + bio7 + bio8 + bio10 + bio11 + bio12 + bio13 + bio15 + bio19
  rf2 <- randomForest(model2, data=envtrain)
  
  #erf1 <- evaluate(testpres, testbackg, rf1)
  erf2 <- evaluate(testpres, testbackg, rf2)
  
  
  # WARNING  plot rf predictions (change model each time you need)
  #pr <- predict(biolay2, rf1)
  pr2 <- predict(biolay.f, rf2)
  
  # extract the polygon from values above 0.85
  
  
  #pr22 <- rasterToPolygons(pr2,fun=function(x){x>0.7}) 
  #pr2.a <- sum(areaPolygon(pr22))/1000000 # determine the area of the 0.85 polygon
  
  
  
  # WARNING change "f." of "p."  in the name of the 
  # objects depending on the desired model (future or present)
  f.list[[i]] <- pr2
  f.tab[i,1] <- erf2@ auc
  #f.tab[i,2] <- pr2.a
}

# future
future <- brick(unlist(f.list))
future.mean <- mean(future)
future.sd <- calc(future, fun = sd)

# mapa final future
ffmas <- future.mean+future.sd
ffmens <- future.mean-future.sd

rgb.palette <- colorRampPalette(c("snow2",
                                  "snow3","seagreen","orange","firebrick"),
                                space = "rgb")

# today and futureplot
tod.plo <- spplot(today.mean, 
                  col.regions=rgb.palette(200),
                  scales =list(draw=T),xlab="longitud",ylab="latitud",
                  main="Actual")

fut.plo <- spplot(today.mean, 
                  col.regions=rgb.palette(200),
                  scales =list(draw=T),xlab="longitud",ylab="latitud",
                  main="Future")

