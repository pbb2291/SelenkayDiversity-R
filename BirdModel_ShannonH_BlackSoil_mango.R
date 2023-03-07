
# library("devtools"); install_github("lme4/lme4",dependencies=TRUE)
library(devtools)
# devtools::install_github("sjPlot/devel")
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

# toggle radius of inquiry
radius = 130 

# RDF files keep relevant variable types (factors, characters, etc.)
XY = readRDS(paste0(datadir, 'in/XY_', radius,'mRadius.rds'))
XY_scale = readRDS(paste0(datadir, 'in/XY_scaled_', radius,'mRadius.rds'))

# Split by soils
XY_scale_black = XY_scale %>% filter(Soil=='Black')

# # # black
# These vars are the top 5 spearman corr
mod130.1.0.black = gam(shannonH ~ s(mean_gapsize, k=3) +
                       s(sd_gapsize, k=3) +
                       s(mean_stdpeakh, k=3) +
                       s(meanH_plot, k=3) +
                       s(Cover1p5m_plot, k=3) +
                       s(X,Y, k=3),
                       data=XY_scale_black,
                       family=gaussian,
                       select=TRUE,
                       method="REML")

summary(mod130.1.0.black)
p130.1.0.black = draw(mod130.1.0.black, residuals = TRUE)
p130.1.0.black

ggsave(plot = p130.1.0.black,
       filename=paste0(figd,
                       "/GAMs/GAM-130-1-0-black_PartialEffects_", radius ,"mRadius_mango.png"),
       width = 14, height = 8, units = "in", device='png', dpi=300)

a130.1.0.black = appraise(mod130.1.0.black)
a130.1.0.black

ggsave(plot = a130.1.0.black,
       filename=paste0(figd,
                       "/GAMs/GAM-130-1-0-black_Appraise_", radius ,"mRadius_mango.png"),
       width = 12, height = 8, units = "in", device='png', dpi=300)

gam.check(mod130.1.0.black, rep=500)
lines(c(3.0, 3.8), c(3.0, 3.8))
# # # This is great - looks like the sd_gapsize is doing the most work here

# Try a simpler model, just sd_gapsize and XY
mod130.1.1.black = gam(shannonH ~ s(sd_gapsize, k=5) +
                         s(X,Y, k=5),
                     data=XY_scale_black,
                     family=gaussian,
                     select=TRUE,
                     method="REML")

summary(mod130.1.1.black)
p130.1.1.black = draw(mod130.1.1.black, residuals = TRUE)
p130.1.1.black

ggsave(plot = p130.1.1.black,
       filename=paste0(figd,
                       "/GAMs/GAM-130-1-1-black_PartialEffects_", radius ,"mRadius_mango.png"),
       width = 14, height = 8, units = "in", device='png', dpi=300)

a130.1.1.black = appraise(mod130.1.1.black)
a130.1.1.black

ggsave(plot = a130.1.1.black,
       filename=paste0(figd,
                       "/GAMs/GAM-130-1-1-black_Appraise_", radius ,"mRadius_mango.png"),
       width = 12, height = 8, units = "in", device='png', dpi=300)

gam.check(mod130.1.1.black, rep=500)
lines(c(3.0, 3.8), c(3.0, 3.8))
# Ok, these don't work at all
# that said, they seem to peak around the 50 m radius
# (but remain high for 130 m radius)
