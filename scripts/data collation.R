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
