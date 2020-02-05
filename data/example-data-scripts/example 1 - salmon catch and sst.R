# create a climate-biology data pair with nonstationary relationship linking the two:
# GOA total salmon catch for four spp. (estimated with a DFA trend)
# and winter SST (three-yr running mean)

library(tidyverse)
library(MARSS)

salmon <- read.csv("data/total.goa.salmon.catch.csv", row.names = 1)

# these data are annual total catch by species, log-transformed and lagged to ocean entry

# drop Chinook - too small a fishery
# and restrict to 1965-2012 - period of mature fisheries, including only two putative states (no heat wave years!)

salmon <- salmon %>%