# create a climate-biology data pair with nonstationary relationship linking the two:
# GOA total salmon catch for four spp. (estimated with a DFA trend)
# and winter SST (three-yr running mean)

library(tidyverse)
library(MARSS)

salmon <- read.csv("data/total.goa.salmon.catch.csv", row.names = 1)

# these data are annual total catch by species, log-transformed and lagged to ocean entry

# drop Chinook - too small a fishery
# and restrict to 1965-2012 - period of mature fisheries, including only two putative states (no heat wave years!)

salmon <- salmon[rownames(salmon) %in% 1965:2012,] %>%
  select(-Chinook)

# transpose and run DFA 
salmon <- as.matrix(t(salmon))

# set up forms of R matrices - these are the four candidate error structures
levels.R = c("diagonal and equal",
             "diagonal and unequal",
             "equalvarcov",
             "unconstrained")

# I'm using the default convergence criteria!

# make an object to save output
model.data = data.frame()

# fit models & store results
for(R in levels.R) {
  for(m in 1) {  # allowing only one shared trend 
    dfa.model = list(A="zero", R=R, m=m)
    kemz = MARSS(salmon, model=dfa.model,
                 form="dfa", z.score=TRUE)
    model.data = rbind(model.data,
                       data.frame(R=R,
                                  m=m,
                                  logLik=kemz$logLik,
                                  K=kemz$num.params,
                                  AICc=kemz$AICc,
                                  stringsAsFactors=FALSE))
    assign(paste("kemz", m, R, sep="."), kemz)
  } # end m loop
} # end R loop

# calculate delta-AICc scores, sort in descending order, and compare
model.data$dAICc <- model.data$AICc-min(model.data$AICc)
model.data <- model.data %>%
  arrange(dAICc)
model.data

# diagonal and equal the best...which I kind of expected!

# now run the best model, save the trend estimates and CIs, and plot loadings/trend
model.list = list(A="zero", m=1, R="diagonal and equal")
mod <- MARSS(salmon, model=model.list, z.score=TRUE, form="dfa")

# get CI and plot loadings...
modCI <- MARSSparamCIs(mod)

plot.CI <- data.frame(names=rownames(salmon), mean=modCI$par$Z, upCI=modCI$par.upCI$Z,
                      lowCI=modCI$par.lowCI$Z)

plot.CI$names <- reorder(plot.CI$names, plot.CI$mean)

dodge <- position_dodge(width=0.9)

ggplot(plot.CI, aes(x=names, y=mean)) + 
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymax=upCI, ymin=lowCI), position=dodge, width=0.5) +
  ylab("Loading") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1,  size=12), axis.title.x = element_blank()) +
  geom_hline(yintercept = 0) + 
  ggtitle("DFA loadings")

# now that's consitency!

# and the shared trend!
trend.CI <- broom::tidy(mod, type="states")
trend.CI$year <- 1965:2012

ggplot(trend.CI, aes(x=year, y=estimate)) + 
  geom_line() +
  geom_ribbon(aes(ymax=conf.high, ymin=conf.low), alpha=0.2) +
  ylab("Shared trend") + 
  theme(axis.title.x = element_blank()) +
  geom_hline(yintercept = 0) + 
  ggtitle("DFA trend")

# now load SST
SST <- read.csv("data/salmon.and.covariate.data.csv")

# and combine!
trend.CI <- trend.CI %>%
  select(year, estimate)

names(SST)[1] <- "year"

SST <- SST %>%
  filter(species=="Coho") %>%
  select(year, SST3) 

dat <- left_join(trend.CI, SST)

names(dat)[2:3] <- c("catch", "SST")

# check - looks familiar!
ggplot(dat, aes(SST, catch, label=year)) +
         geom_text()

# and save 
write.csv(dat, "data/example 1 - salmon catch and sst.csv", row.names = F)
