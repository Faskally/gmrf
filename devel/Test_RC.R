

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
 

## try one for disconnected catchments
# load catchments
ctm <- rgdal::readOGR("P:/vector/nationalgrid/Catchment/SEPA","Baseline_confluence_nested_catchments")
load("C:/work/repos/Faskally/gmrf/devel/ctm_data.rData")
ctm @ data <- ctm_data
hma <- CLdata::hma
library(sp)


# Select a few
#hanames <- c("Esk Group", "Deveron Group", "Dee (Aberdeenshire)", "Don (Aberdeenshire)", "Spey","Ythan Group")
hanames <- c("Esk Group", "Deveron Group")
ctm <- ctm[ctm $ HAName %in% hanames,]
hma <- hma[hma $ HAName %in% hanames,]

# remove "605" so we have a singleton catchment
ctm <- ctm[rownames(ctm@data) != "605",]

plot(ctm, col = grey(0.8))
plot(hma, lwd = 2, add = TRUE)
text(coordinates(ctm), labels = rownames(ctm@data), cex = 0.5)

# get neighbours
nb <- spdep::poly2nb(ctm, queen = FALSE)
# look at the adjacency matrice
library(magrittr)
Dnb(nb) %>% replace(. == 0, ".") %>% as.table(.)


Q <- getQnb(nb)
Q %>% replace(. == 0, ".") %>% as.table(.)

# identify singletons
singles <- which(diag(Q) == 0)

# remove singletons from spatial matrix
Qsub <- Q[-singles, -singles]

# the null space of Q tells you how many groups
## the vectors of the null space correspond to the individual regions
nullQ <- MASS::Null(Qsub)
rownames(nullQ) <- rownames(ctm@data)[-singles]
Matrix(nullQ, sparse = TRUE)
null.space.dim <- ncol(nullQ)

ctm $ grp <- 0
ctm @ data $ grp[-singles] <- rowSums((nullQ != 0) * rep(1:null.space.dim, each = nrow(nullQ)))
plot(ctm, col = ctm$grp+1)

# create codes for groups of catchments
ctm $ ctmgrp <- NA
ctm @ data $ ctmgrp[singles] <- seq_along(singles)
ctm @ data $ ctmgrp[-singles] <- rowSums((nullQ != 0) * rep(1:rank, each = nrow(nullQ))) + length(singles)

# so the idea would be to constrain each group to sum to zero...
X <- model.matrix(~ factor(ctmgrp) - 1, ctm @ data)
cmu <- c(1, 2, 3) * 2

x <- simQ(Q)
tapply(x, ctm $ ctmgrp, mean) #  looks like groups are zero

y <- x + c(X %*% cmu) + rnorm(length(x))*0.5

breaks <- seq(min(y)-0.001, max(y)+0.001, length = 11)
xcolgrp <- as.numeric(cut(y, breaks))
xcols <- heat.colors(10)[xcolgrp]
plot(ctm, col = xcols)

cid <- rownames(ctm @ data)
cid <- 1:length(x)
rownames(Q) <- colnames(Q) <- cid
constraint <- matrix(0, null.space.dim + length(singles), nrow(Q))
constraint[1:null.space.dim,-singles] <- t(nullQ) %>% replace(. != 0, 1)
constraint[-(1:null.space.dim), singles] <- diag(length(singles))

ctmgrp <- ctm $ ctmgrp

g1 <- gam(y ~ -1 + factor(ctmgrp) + s(cid, bs = "gmrf", xt = list(penalty = Q, constraint = constraint)), method="REML")
summary(g1)
tapply(fitted(g1), ctm $ ctmgrp, mean) #  looks like groups are zero

# plot fitted values
xcolgrp <- as.numeric(cut(fitted(g1), breaks))
xcols <- heat.colors(10)[xcolgrp]
plot(ctm, col = xcols)

# extract just catchment spatial effect


##mgcv:::plot.mgcv.smooth
## but what happens if we drop off number
