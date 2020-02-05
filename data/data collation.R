library(tidyverse)

# load PDO/SST/salmon catch data from reversing-PDO paper

dat <- read.csv("data/salmon.and.covariate.data.csv")

# these PDO/SST data are winter (NDJFM) values, 
# and catches are lagged to ocean entry, long-transformed, and normalized

dat <- dat %>%
  select(-PDO2, -SST2) %>%
  pivot_wider(names_from=species, 
         values_from=catch)

# add npgo

# download.file("http://www.oces.us/npgo/npgo.php", "~npgo")
npgo <- read.table("data/npgo", skip=10, nrows=828, fill=T, col.names = c("Year", "month", "value"))

npgo$win.yr <- ifelse(npgo$month %in% 11:12, npgo$Year+1, npgo$Year)

npgo <- npgo %>%
  filter(month %in% c(11,12,1:3))

npgo1 <- tapply(npgo$value, npgo$win.yr, mean)

# and drop 2019 b/c we don't have complete observations
npgo1 <- npgo1[1:length(npgo1)-1]

npgo3 <- zoo::rollmean(npgo1, 3, fill=NA)

# adding these values for now, though I wonder if we shouldn't just use the pdo in order to keep things simple!

dat$NPGO1 <- npgo1[match(dat$Year, names(npgo1))]
dat$NPGO3 <- npgo3[match(dat$Year, names(npgo3))]

# now load more climate data
clim.dat <- read.csv("data/updated goa climate data.csv")

# drop the versions of SST here as we don't need them
clim.dat <- clim.dat %>%
  select(-c(10:13))

# and join them together 
names(dat)[1] <- "year"
dat <- left_join(dat, clim.dat)

# now add some more population time series

# first, chum salmon catch from 2018 PRSB paper
other.catch <- read.csv("data/total.goa.salmon.catch.csv")




###############
# breaking off from above to punch out a climate data set for nonstationary DFA

climate <- read.csv("data/updated goa climate data.csv")

head(climate)

climate$upwelling <- apply(climate[,6:8], 1, mean)

climate <- climate %>%
  filter(year <= 2012) %>%
  select(1:5, 14:15)

# add a time series of full-GOA winter SST

library(ncdf4)
library(maps)
library(mapdata)
library(fields)
library(chron)

# script for calculating GOA sst anomalies wrt 1951-1980
download.file("https://coastwatch.pfeg.noaa.gov/erddap/griddap/nceiErsstv5.nc?sst[(1900-01-01):1:(2019-12-01T00:00:00Z)][(0.0):1:(0.0)][(54):1:(62)][(200):1:(226)]", "data/temp.nc")
# load and process SST data
nc <- nc_open("data/temp.nc")

# extract dates

ncvar_get(nc, "time")   # seconds since 1-1-1970
raw <- ncvar_get(nc, "time")
h <- raw/(24*60*60)
d <- dates(h, origin = c(1,1,1970))

# extract study area
# 54-62 deg. N, 200-226 deg. E
x <- ncvar_get(nc, "longitude")
y <- ncvar_get(nc, "latitude")

SST <- ncvar_get(nc, "sst", verbose = F)

# Change data from a 3-D array to a matrix of monthly data by grid point:
# First, reverse order of dimensions ("transpose" array)
SST <- aperm(SST, 3:1)  

# Change to matrix with column for each grid point, rows for monthly means
SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))  

# Keep track of corresponding latitudes and longitudes of each column:
lat <- rep(y, length(x))   
lon <- rep(x, each = length(y))   
dimnames(SST) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

# need to drop Bristol Bay cells
BB <- c("N58E200", "N58E202", "N56E200")
SST[,BB] <- NA

month.sst <- rowMeans(SST, na.rm = T)

# get winter-only means
yr <- as.numeric(as.character(years(d)))
m <- months(d)

# set up winter year
win.yr <- yr
win.yr[m %in% c("Nov", "Dec")] <- win.yr[m %in% c("Nov", "Dec")]+1

winter.month.sst <- month.sst[m %in% c("Nov", "Dec", "Jan", "Feb", "Mar")]
win.yr <- win.yr[m %in% c("Nov", "Dec", "Jan", "Feb", "Mar")]

win.sst <- tapply(winter.month.sst, win.yr, mean)
win.sst <- win.sst[names(win.sst) %in% 1950:2012]
plot(1950:2012, win.sst, type="l")

climate$winter.sst <- win.sst

# need to add freshwater!

# load data from ecology paper
dat <- read.csv("data/GOA environmental data.csv")

head(dat)
tail(dat)

climate$freshwater <- dat$Freshwater

# for grins, compare SST time series!
plot(dat$SST, climate$winter.sst)
# similar, but different somehow!
# best to conitnue using the one I just downloaded...might be that the difference has
# something to do with ERSSTv5 vs ERSSTv4...also maybe data processing differences

# and save this example
write.csv(climate, "data/example 2 - goa climate variables.csv", row.names = F)
