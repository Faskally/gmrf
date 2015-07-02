

# set proxy values
httr::set_config(httr::use_proxy(url="192.168.41.8", port=80))

# load packages
devtools::install_github("faskally/gmrf")
library(gmrf)
#pkg <- devtools::as.package("C:/work/repos/Faskally/gmrf")
#devtools::load_all(pkg)

example(smooth.construct.gmrf.smooth.spec)

# -----------------------------------------------
#
# lets try a spatial effect 
#
# -----------------------------------------------

## fit a regional smoother first
hma <- CLdata::hma
# drop shetland and orkney
hma <- hma[!hma $ HAName %in% c("Shetlands", "Orkneys", "Inner Hebrides", "Outer Hebrides"),]

plot(hma)
Q1 <- getQpoly(hma)

nb <- spdep::poly2nb(hma, queen = FALSE)
Q2 <- getQnb(nb)

nbmat <- spdep::nb2mat(nb, style = "B", zero.policy = TRUE)
Q3 <- getRegionalGMRF(nbmat)


## fit a regional smoother first
hma <- CLdata::hma
# drop shetland and orkney
hma <- hma[!hma $ HAName %in% c("Shetlands", "Orkneys"),]

# calculate adjadency matrix
hmaadj <- spdep::poly2nb(hma, queen = FALSE)
hmaadj <- spdep::nb2mat(hmaadj, style = "B", zero.policy = TRUE)

# add connections for inner and outer hebs
hmaadj[hma $ HAName == "Inner Hebrides", hma $ HAName == "Outer Hebrides"] <- 1
hmaadj[hma $ HAName == "Outer Hebrides", hma $ HAName == "Inner Hebrides"] <- 1
hmaadj[hma $ HAName == "Inner Hebrides", c(21, 42, 43)] <- 1
hmaadj[c(21, 42, 43), hma $ HAName == "Inner Hebrides"] <- 1

Q <- getRegionalGMRF(hmaadj)
colnames(Q) <- rownames(Q) <- hma $ HACode

par(mfrow = c(1,2))
# simulate spatial process
set.seed(334543)
x <- simQ(Q)
breaks <- seq(min(x)-0.001, max(x)+0.001, length = 11)
xcolgrp <- as.numeric(cut(x, breaks))
xcols <- heat.colors(10)[xcolgrp]
plot(hma, col = xcols)

# simulate observations
n <- 100
id <- sample(1:length(x), n, replace = TRUE)
y <- x[id] + rnorm(100) * 0.5
HA <- hma $ HACode[id]

g1 <- gam(y ~ s(HA, bs = "gmrf", xt = list(penalty = Q)), method="REML")
summary(g1)
pred <- predict(g1, newdata = list(HA = hma $ HACode))
pcolgrp <- as.numeric(cut(pred, breaks))
pcols <- heat.colors(10)[pcolgrp]
plot(hma, col = pcols)


# -----------------------------------------------
#
# lets try a spatial effect with discontinuities
#
# -----------------------------------------------

pkg <- devtools::as.package("C:/work/repos/Faskally/gmrf")
devtools::load_all(pkg)
library(sp)
library(magrittr)

 
# -------------------------------------
# get data
# -------------------------------------

# load catchments
ctm <- rgdal::readOGR("P:/vector/nationalgrid/Catchment/SEPA","Baseline_confluence_nested_catchments")
load("C:/work/repos/Faskally/gmrf/devel/ctm_data.rData")
ctm @ data <- ctm_data
hma <- CLdata::hma


# Select a few
hanames <- c("Esk Group", "Deveron Group")
ctm <- ctm[ctm $ HAName %in% hanames,]
# remove "605" so we have a singleton catchment
ctm <- ctm[rownames(ctm@data) != "605",]
hma <- hma[hma $ HAName %in% hanames,]

# -------------------------------------
# get GMRF
# -------------------------------------

nb <- spdep::poly2nb(ctm, queen = FALSE)
Q <- getQnb(nb)

Q %>% replace(. == 0, ".") %>% as.table(.)


# -------------------------------------
# Get constraint martrix and a factor for grouping
# -------------------------------------

constraint <- getCnb(Q)
ctmgrp <- getFactorsnb(Q)

# visualise groupings
plot(ctm, col = ctmgrp)


# -------------------------------------
# simulate some data
# -------------------------------------

# design matrix for catchment group means
X <- model.matrix(~ ctmgrp - 1)
cmu <- c(1, 2, 3) * 2

# spatial effect (in effect catchment within catchment group)
x <- simQ(Q)
# check groups sum to zero
tapply(x, ctmgrp, mean)

# simulate observations, 1 per region in this case 
y <- x + c(X %*% cmu) + rnorm(length(x))*0.5

# get region ids for observations
# note region ID should not be character
# it can either be numeric, or a factor.
cid <- factor(rownames(ctm @ data))

# collect data in a list
dat <- data.frame(y = y, cid = cid, ctmgrp)

# collect smoother details in a list
xt <- list(penalty = Q, constraint = constraint)

# plot simulation
breaks <- seq(min(y)-0.001, max(y)+0.001, length = 11)
plot(ctm, col = heat.colors(length(breaks)-1)[as.numeric(cut(y, breaks))])

# -------------------------------------
# fit a model
# -------------------------------------

# use REML this time
g1 <- gam(y ~ -1 + ctmgrp + s(cid, bs = "gmrf", xt = xt), method="REML", data = dat)

summary(g1)
# check groups sum to ctm group means
tapply(fitted(g1), dat$ctmgrp, mean)

# plot fitted values
plot(ctm, col = heat.colors(length(breaks)-1)[as.numeric(cut(fitted(g1), breaks))])



# -----------------------------------------------
#
# lets try a reduced rank spatial effect with discontinuities
#
# -----------------------------------------------

# run all of the previous code to get data and GMRF model

# what does




# use REML this time
g1 <- gam(y ~ -1 + ctmgrp + s(cid, k = 20, bs = "gmrf", xt = xt), method="REML", data = dat)

summary(g1)
# check groups sum to ctm group means
tapply(fitted(g1), dat$ctmgrp, mean)

# plot fitted values
plot(ctm, col = heat.colors(length(breaks)-1)[as.numeric(cut(fitted(g1), breaks))])

