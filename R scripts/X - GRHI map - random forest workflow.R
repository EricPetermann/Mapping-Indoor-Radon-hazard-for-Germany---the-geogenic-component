library(Smisc)
library(caret)
library(CAST)
library(party)
library(plyr)
library(GSIF)

#load GRP data
Data <- loadObject("GRP data.R")

#training indices of spatial folds
sb<-loadObject("spatialBlocks GRP FFS Tune.R")

# defining metrics to be evaluated
mySummary <- function (data,
                       lev = NULL,
                       model = NULL) {
  out <-
    c(
      Metrics::rmse(data$obs, data$pred),
      cor(data$obs, data$pred) ^ 2,
      Metrics::mape(data$obs, data$pred),
      Metrics::rmsle(data$obs, data$pred),
      Metrics::mae(data$obs, data$pred),
      Metrics::rrse(data$obs, data$pred),
      Metrics::percent_bias(data$obs, data$pred)
    )
  names(out) <- c("rmse", "R2","mape","rmsle","mae","rrse", "pb")
  out
}

#control parameters
fitControl <- trainControl(
  method = 'repeatedcv',            # repeated cross validation
  number = 10,                      # number of folds
  repeats=1,                        # number of repititions
  savePredictions = 'final',        # saves predictions for optimal tuning parameter
  summaryFunction=mySummary,        # results summary function
  index =sb                         # indices of training data over 10 folds, repeated 5 times
)


##########                          #########
###  Predictor selection                  ###  
#########                           #########

#define tuningGrid
rfGrid <-  expand.grid(mtry = 5)  # hold constant at mtry=5 for predictor selection

#variables to include
include<-which(!names(Data) %in% c("GRP"))

# start forward feature selection
FFS_GRP <-ffs(
  Data[,include],
  log(Data$GRP),
  metric = "rmse",
  maximize = FALSE,
  method = "cforest",
  tuneGrid = rfGrid,
  trControl = fitControl,
  controls = cforest_unbiased(ntree = 100, trace = TRUE)
) 


##########      #########
###       Tuning      ###  
#########       #########

#selected vars
vars_selected<-FFS_GRP$selectedvars
selected<-which(names(Data)%in%vars_selected)   # find column numbers with specific column names

# define hyperparameter grid
rfGrid<-expand.grid(mtry=seq(2,length(vars_selected),1))

Tune_GRP <-train(
  Data[,selected],
  log(Data$GRP),
  metric = "rmse",
  method = "cforest",
  tuneGrid = rfGrid,
  trControl = fitControl,
  maximize = FALSE,
  controls = cforest_unbiased(ntree = 500, trace = TRUE)
)


##########             #########
###       Performance       ###  
#########              #########

#cross-validation folds for 10 times repeated 10-fold spatial cross-validation
sb_pe <- loadObject("spatialBlocks GRP Performance.R")

#fit control
fitControl <- trainControl(
  method = 'repeatedcv',          # repeated cross validation
  number = 10,                    # 10 folds
  repeats=10,                     # 10 repitiona
  savePredictions = 'final',      # save final predictions
  summaryFunction=mySummary,      # performance metrics summary function
  index =sb_pe,                   # indices of training data 
  verboseIter = TRUE              # trace progress
)

# define optimal mtry
rf.tuned <- expand.grid(mtry=c(as.numeric(Tune_GRP$bestTune)))

Perf_RF_GRP <-train(
  Data[,selected],
  log(Data$GRP),
  metric = "rmse",
  method = "cforest",
  tuneGrid = rf.tuned,
  trControl = fitControl,
  maximize = FALSE,
  controls = cforest_unbiased(ntree=1000,trace = TRUE)
)

##########             #########
###       final model        ###  
#########              #########

#prediction grid
Grid<-loadObject("Predictor Grid GRP.R")

###unify levels for categorical variables for predictors attributed to observational data
Geo_levels <- union(levels(Data$GK250), levels(Grid$GK250))
Data$GK250=with(Data,factor(GK250,levels=Geo_levels))

include<-which(names(Data)%in%c(vars_selected,"GRP"))   # find column numbers with specific column names

Mod_GRP <- party::cforest(log(GRP)~.,data=Data[,include],
                             controls=cforest_unbiased(ntree = 1000, mtry = as.numeric(Tune_GRP$bestTune), 
                                                       trace = TRUE))

##########             #########
###       Mapping            ###  
#########              #########

# define geological data as factor
Grid$GK250=as.factor(Grid$GK250)

# find variable names for estimation grid
Grid.r<-raster::brick("Prediction Raster.tif")
Grid.r@data@names=names(Grid[3:40])

# prepare estimation grid
include<-which(names(Grid.r)%in%vars_selected)   # find column numbers with specific column names
Grid.GRP <- Grid.r[[include]]
plot(Grid.GRP)
Grid.df <- as.data.frame(Grid.GRP)
Grid.df$GK250<-as.factor(Grid.df$GK250)

#define extent
x=seq(280500,921500,1000)
y=seq(5235500,6106500,1000)
xy=expand.grid(x,y)
# define crs
UTM32=CRS("+init=epsg:25832")

#tiling
Grid.p=SpatialPointsDataFrame(xy[,1:2],xy,proj4string = UTM32,bbox=NULL) #create sp point object
tiles <- GSIF::getSpatialTiles(Grid.p, block.x=50000, overlap.percent = 5,return.SpatialPolygons = FALSE) #create 50 km*50 km tiles
tiles.pol <- GSIF::getSpatialTiles(Grid.p, block.x=50000,block.y=50000,overlap.percent = 5,
                                   return.SpatialPolygons = TRUE) #convert to spatial polygons
tile.pol <- SpatialPolygonsDataFrame(tiles.pol, tiles)            #convert to spatial polygons data frame

# here are the tiles...
plot(Grid.GRP[[7]], col=bpy.colors(20))
lines(tiles.pol, lwd=2)

#create folder for tiles; navigate
setwd("./tiles")

### function for mapping tiled data (including caterorical predictors)
dir="."
fun_RF <- function(i, tiles, dir="."){
  out.tif <- paste0("Tile_", i, ".tif")
  if(!file.exists(out.tif)){                      # check if tiles data already exists
    Grid_crop <- crop(Grid.GRP,extent(tiles.pol[i])) # crop estimation grid to extent of tile i
    Grid_tile_df <- as.data.frame(Grid_crop)         # convert to data.frame 
    NAs<-sum(is.na(Grid_tile_df))                 # sum NA values in tile  
    if(NAs!=length(unlist(Grid_tile_df))){        # check if tile contains not only NA;  if so continue with analysis, else next i
      Grid_tile_df$GK250 <- as.factor(Grid_tile_df$GK250)   #define factors
      Geo_levels <- union(levels(Data$GK250), levels(Grid_tile_df$GK250))  ##unify factor levels with union of levels
      Grid_tile_df$GK250 <- with(Grid_tile_df,factor(GK250,levels=as.factor(sort(as.numeric(Geo_levels))))) #rebuild factor levels with unified levels
      levels(Grid_tile_df$GK250) <- mapvalues(levels(Grid_tile_df$GK250),
                                           from = levels(Grid_tile_df$GK250),
                                           to = levels(Data$GK250))   #substitute factor values back to those used for model building
      Grid_tile_df$lnGRP <- predict(Mod_GRP,newdata=Grid_tile_df,type="response")[,1]  #make prediction
      ext=extent(Grid_crop)                     # define grid extent and resolution for tile
      x=seq(ext[1]+500,ext[2]-500,1000)         # set coordinates of cell midpoint (+/- 500 metres)
      y=seq(ext[3]+500,ext[4]-500,1000)         # set coordinates of cell midpoint (+/- 500 metres)
      xy=expand.grid(x,y)                       # create grid
      GRP.p<-SpatialPointsDataFrame(xy[,1:2],Grid_tile_df["lnGRP"],proj4string = UTM32,bbox=NULL) #convert dataframe to sp object
      gridded(GRP.p) <- TRUE                   # convert to sptial pixel dataframe
      GRP.rf<-flip(raster(GRP.p),direction=2)  # flipping raster required! 
      writeRaster(GRP.rf,out.tif,overwrite=TRUE) # save raster tile
    }
  }
}

# apply a function over  a raster (multicore usage possible)
x0 <- mclapply(1:length(tiles.pol), FUN=fun_RF, tiles=tiles)

### assembly tiles to a map
# find all tiles
f <-list.files(pattern = glob2rx("Tile_*.tif"))  # "*" is the wildcard

# make a list  
r <- lapply(f, FUN=raster) 
# as you have the arguments as a list call 'merge' with 'do.call'  
x <- do.call("merge",r) 
