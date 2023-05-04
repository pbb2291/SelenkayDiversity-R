# Shannon H Model - Possible Prediction Models
library(tidyverse)
library(mgcv)
library(arrow)
library(lme4)
library(gratia)
library(gridExtra)

# # # 1) Load XY Data
# Note: XY data prepared in "BirdDataPreperation.R"

# Datadir for location of in/out vars
datadir = '/n/home02/pbb/scripts/SelenkayDiversity/data/'

figd = '/n/home02/pbb/scripts/SelenkayDiversity/figs/mango'

tabled = '/n/home02/pbb/scripts/SelenkayDiversity/tables/mango/GAMs'

# toggle radius of inquiry
# radius = 130
# radius = 80
radius = 50

# RDF files keep relevant variable types (factors, characters, etc.)
XY = readRDS(paste0(datadir, 'in/XY_', radius,'mRadius.rds'))
XY_scale = readRDS(paste0(datadir, 'in/XY_scaled_', radius,'mRadius.rds'))

mod130.5.0 = gam(shannonH ~
                   s(sd_gapsize),
                 data=XY_scale,
                 family=gaussian,
                 select=TRUE,
                 method="REML")

summary(mod130.5.0)
mod130.5.0$aic
p130.5.0 = draw(mod130.5.0, residuals = TRUE)
p130.5.0
appraise(mod130.5.0)

# # # 

mod130.5.1 = gam(shannonH ~
                   s(sd_gapsize, by=Soil_f) +
                   s(X,Y, k=5),
                 data=XY_scale,
                 family=gaussian,
                 select=TRUE,
                 method="REML")

summary(mod130.5.1)
p130.5.1 = draw(mod130.5.1, residuals = TRUE)
p130.5.1
