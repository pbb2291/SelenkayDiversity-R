# Initial Abundance Model
# For EA Presentation
# 11/14/22 - PB

library(tidyverse)
library(mgcv)
library(arrow)
library(gratia)
# # # 1) Load XY Data
# Note: XY data prepared in "BirdDataPreperation.R"

# Datadir for location of in/out vars
datadir = '/n/home02/pbb/scripts/SelenkayDiversity/data/'

# RDF files keep relevant variable types (factors, characters, etc.)
XY = readRDS(paste0(datadir, 'in/XY.rds'))
XY_scale = readRDS(paste0(datadir, 'in/XY_scaled.rds'))

# # # Examine Y variable

# Plot a histogram
hist(XY_scale$Abundance)

# Draw QQPlot for normality
qqnorm(XY_scale$Abundance)
qqline(XY_scale$Abundance)


# # # Fit Models
mod1 = gam(Abundance ~ s(mean.sdH., by=Soil_f) + Soil_f + Transect_f,
           data=XY_scale,
           family=poisson,
           select=TRUE,
           method="REML")

summary(mod1)
draw(mod1, residuals = TRUE)
appraise(mod1)
gam.check(mod1, rep=500)
lines(c(100, 220), c(100, 220))

# # # Fit Models

mod2 = gam(Abundance ~ s(herbh_plot) + s(cvpeakh_plot) ,
           data=XY_scale,
           family=poisson,
           select=TRUE,
           method="REML")

summary(mod2)


mod3 = gam(Abundance ~ sd.coverG.:Soil_f + s(Soil_f, bs="re"),
           data=XY_scale,
           family=poisson,
           select=TRUE,
           method="REML")

summary(mod3)

AIC(mod1, mod2, mod3)
draw(mod3, residuals = TRUE)
appraise(mod3)
gam.check(mod2, rep=500)
lines(c(3.0, 3.8), c(3.0, 3.8))
# Abundance MDI Top 6 X Vars:
#   
#   herbh_plot
# cv.FHD.
# cvpeakh_plot
# cover_grass_maxH
# mean.cvH._grass
# sd.cvH._shrub
# cv.cscore.
# ptoh_plot
# sd.cvH._grass
# max.ptoh.
# 
# Abundance Permutation Top 6 X Vars:
#   
#   ptoh_plot
# mean.sdH._shrub
# VDR_plot
# meanH_plot
# mean.cvH._shrub
# X75thPerc_plot
# iqr.FHD.
# max.herbh.
# max.VDRpeak.
# iqr.nlayers.


ggplot(XY_scale, aes(x = mean.sdH., y = Abundance,
                     group = Soil_f, colour = Transect_f)) +
  geom_point() +
  facet_wrap(~ Soil_f)


ggplot(XY_scale, aes(x = sd.coverG., y = Abundance,
                     group = Soil_f, colour = Transect_f)) +
  geom_point() +
  facet_wrap(~ Soil_f)


ggplot(XY_scale, aes(x = mean.nlayers., y = Abundance,
                     group = Soil_f, colour = Transect_f)) +
  geom_point()


ggplot(XY_scale, aes(x = sd.sdH._grass, y = Abundance,
                     group = Soil_f, colour = Transect_f)) +
  geom_point()

# # # Models by Soil type

XY_scale_black = XY_scale %>% filter(Soil == 'Black')

mod1_black = gam(Abundance ~ cv.coverG.:Transect_f + s(cv.coverG.),
                 data=XY_scale_black,
                 family=scat(link="identity"),
                 select=TRUE,
                 method="REML")

summary(mod1_black)



XY_scale_red = XY_scale %>% filter(Soil == 'Red')

mod1_red = gam(Abundance ~ cv.coverG.:Transect_f + s(cv.coverG.),
                 data=XY_scale_red,
                 family=scat(link="identity"),
                 select=TRUE,
                 method="REML")

summary(mod1_red)
# + s(sd.maxH.) + 
# + s(sd.maxH.) + 
  # s(mean.nlayers.) + s(mean.coverG.)
