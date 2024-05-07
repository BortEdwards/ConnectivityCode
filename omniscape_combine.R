# Post-run Omniscape processing:
# 1. - rescale, and stack outputs from different Omniscape window size/scale analyses
# 2. - bin stacked layers into quantiles by geographic area/pixel number
# 3. - trial code to normalize, reclassify and stack different window size/scale analyses (can be appropriate if not using the normalized Omniscape output)

########################
# Packages
library(dplyr)
library(terra)
library(raster)
########################

######################## 1. stack outputs from differnt Omniscape window size/scale analyses (mostly approproate for already normalized Ominscape outputs)

outpth <- "C:/Users/C836876692/Desktop/Omniscape Test/Omniscape CO-wide/Processed Output"
wopt <- list(gdal=c("COMPRESS=LZW", "TFW=YES", "BIGTIFF=YES"))

# Mean ntiles of Raw Current ------------------------------------------------
r2.5 <- "C:/Users/C836876692/OneDrive - Colostate/Omniscape/Omniscape CO-wide/Omniscape Output Rasters/AllCO2.5miles37blocks_normalized.tif"
r5 <- "C:/Users/C836876692/OneDrive - Colostate/Omniscape/Omniscape CO-wide/Omniscape Output Rasters/AllCO5miles75blocks_normalized.tif"
r10 <- "C:/Users/C836876692/OneDrive - Colostate/Omniscape/Omniscape CO-wide/Omniscape Output Rasters/AllCO10miles151blocks_normalized.tif"
r25 <- "C:/Users/C836876692/OneDrive - Colostate/Omniscape/Omniscape CO-wide/Omniscape Output Rasters/AllCO25miles375blocks_normalized.tif"
r50 <- "C:/Users/C836876692/OneDrive - Colostate/Omniscape/Omniscape CO-wide/Omniscape Output Rasters/AllCO50miles751blocks_normalized.tif"

#make raster objects
ras2.5 <- rast(r2.5)
ras5 <- rast(r5)
ras10 <- rast(r10)
ras25 <- rast(r25)
ras50 <- rast(r50)

#extract metadata based on 2.5 raster
testarea <- ext(ras2.5)
allrows <- nrow(ras2.5)
allcols <- ncol(ras2.5)

v2.5 <- values(ras2.5, mat=FALSE)
v5 <- values(ras5, mat=FALSE)
v10 <- values(ras10, mat=FALSE)
v25 <- values(ras25, mat=FALSE)
v50 <- values(ras50, mat=FALSE)

#calculate percentiles for each
ptile2.5 <- ntile(v2.5, 100)
ptile5 <- ntile(v5, 100)
ptile10 <- ntile(v10, 100)
ptile25 <- ntile(v25, 100)
ptile50 <- ntile(v50, 100)

# Create a list of raster objects containing percentile values
plist <- list(rast(nrows=allrows, ncols=allcols, crs=crs(ras2.5), extent=testarea, vals=ptile2.5),
              rast(nrows=allrows, ncols=allcols, crs=crs(ras2.5), extent=testarea, vals=ptile5),
              rast(nrows=allrows, ncols=allcols, crs=crs(ras2.5), extent=testarea, vals=ptile10),
              rast(nrows=allrows, ncols=allcols, crs=crs(ras2.5), extent=testarea, vals=ptile25),
              rast(nrows=allrows, ncols=allcols, crs=crs(ras2.5), extent=testarea, vals=ptile50))
names(plist) <- c("radius2.5", "radius5", "radius10", "radius25", "radius50")

#stack
ptstack <- rast(plist)

# Calculate the mean of the stacked percentile raster and write to a new file
combRast <- app(ptstack, fun=mean, filename=file.path(outpth, "CompiledNormalizedOutputs.tif"), wopt=wopt)

######################## 2. - bin composite layer into quantiles with equal pixel number

# Extract values from the raster
raster_values <- values(combRast)

# Calculate terciles
quantiles <- quantile(raster_values, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)

### alternate quantiles
## 5 bins (quintiles)
## quantiles <- quantile(raster_values, probs = seq(0, 1, by = 0.25), na.rm = TRUE)

## 10 bins (deciles)
## quantiles <- quantile(raster_values, probs = seq(0, 1, by = 0.1), na.rm = TRUE)

# Reclassify raster into quantiles
reclassified_raster <- classify(combRast, quantiles, include.lowest = TRUE)

# Plot the reclassified raster
plot(reclassified_raster)

# Save the reclassified raster as a .tif file
output_file <- "CompiledNormalizedOutputsTerciles.tif"
writeRaster(reclassified_raster, filename = output_file, overwrite = TRUE)


######################## 3. - trial code to normalize and rescale Omniscape output layers based on pixel values (can be appropriate if not using the already normalized Omniscape output)
# load files ---------------------------------------
n2.5 <- "C:/Users/C836876692/Desktop/Omniscape Test/Omniscape CO-wide/AllCO2.5miles37blocks/Output/AllCO2.5miles37blocks.tif"
n5 <- "C:/Users/C836876692/Desktop/Omniscape Test/Omniscape CO-wide/AllCO5miles75blocks/Output/AllCO5miles75blocks.tif"
n10 <- "C:/Users/C836876692/Desktop/Omniscape Test/Omniscape CO-wide/AllCO10miles151blocks/Output/AllCO10miles151blocks.tif"
n25 <- "C:/Users/C836876692/Desktop/Omniscape Test/Omniscape CO-wide/AllCO25miles375blocks/Output/AllCO25miles375blocks.tif"
n50 <- "C:/Users/C836876692/Desktop/Omniscape Test/Omniscape CO-wide/AllCO50miles751blocks/Output/AllCO50miles751blocks.tif"

#make raster objects
nras2.5 <- rast(n2.5)
nras5 <- rast(n5)
nras10 <- rast(n10)
nras25 <- rast(n25)
nras50 <- rast(n50)

#list
raslist <- list(nras2.5, nras5, nras10, nras25, nras50)
outlist <- list()

#Function to rescale and normalize raster value (setting min and max)
normrescale <- function(inras){
  setMinMax(inras, force=TRUE)
  inrange <- minmax(inras)
  orange1 <- 0.6 - inrange[1,1] #calculate range from min to 0.6
  orange2 <- inrange[2,1] - 1.4 #calculate range from 1.4 to max
  outras <- ifel(inras < 0.6, ((inras - inrange[1,1])*9)/orange1, #condition1: if < 0.6 rescale linearly 1-10
                 ifel(inras > 1.4, (((inras - 1.4)*10)/orange2) + 10, 10)) #condition2: if > 1.4 rescale linearly 10-20, remaining (0.6-1.4) = 10
  return(outras)
}

#Loop through each raster, apply the normalization and rescaling, and append to the output list
for(ras in raslist){
  rscale <- normrescale(ras)
  outlist <- append(outlist, rscale)
}

#assign names
names(outlist) <- c("radius2.5", "radius5", "radius10", "radius25", "radius50")

#Stack the rescaled and normalized current rasters
normout <- app(outlist, fun=mean, 
               filename=file.path(outpth, "CompiledNormalizedRescaledfromRawOutput.tif"), 
               wopt=wopt)
