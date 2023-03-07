# Initial Shannon H Model
# For EA Presentation
# 11/11/22 - PB

# Useful links:
# Random Effects in mgcv: https://fromthebottomoftheheap.net/2021/02/02/random-effects-in-gams/
# Intro to GAMs video- https://www.youtube.com/watch?v=sgw4cu8hrZM&t=161s

# Adding spatial autocorrelation with s(x,y) term:
# https://stats.stackexchange.com/questions/35510/why-does-including-latitude-and-longitude-in-a-gam-account-for-spatial-autocorre?newreg=af53a8f989cb4490a9b45f9498ab2a70
# https://onlinelibrary.wiley.com/doi/pdfdirect/10.1111/j.1466-8238.2006.00279.x
# https://noamross.github.io/gams-in-r-course/chapter3

# library("devtools"); install_github("lme4/lme4",dependencies=TRUE)
library(devtools)
# devtools::install_github("sjPlot/devel")
library(tidyverse)
library(mgcv)
library(arrow)
library(lme4)
library(gratia)

# # # 1) Load XY Data
# Note: XY data prepared in "BirdDataPreperation.R"

# Datadir for location of in/out vars
datadir = '/n/home02/pbb/scripts/SelenkayDiversity/data/'

# RDF files keep relevant variable types (factors, characters, etc.)
XY = readRDS(paste0(datadir, 'in/XY.rds'))
XY_scale = readRDS(paste0(datadir, 'in/XY_scaled.rds'))

# # # Examine Y variable

# Draw QQPlot for normality
qqnorm(XY_scale$shannonH)
qqline(XY_scale$shannonH)

# Plot a histogram
hist(XY_scale$shannonH)

# Black and red soils look like they have diff. distributions...
hist(XY_scale$shannonH[XY_scale$Soil=='Black'])
hist(XY_scale$shannonH[XY_scale$Soil=='Red'])

# # # Fit models - from simplest to most complex
mod_cvH_1.0 = gam(shannonH ~ s(mean.cvH._woody) + Transect_f,
             data=XY_scale,
             family=gaussian,
             select=TRUE,
             method="REML")

# Ok! Actually, looks quite good
draw(mod_cvH_1.0, residuals = TRUE)

summary(mod_cvH_1.0)
appraise(mod_cvH_1.0)
gam.check(mod_cvH_1.0, rep=500)
lines(c(3.0, 3.8), c(3.0, 3.8))

# Add another var
mod_cover0p5m_1.0 = gam(shannonH ~ s(Cover0p5m_plot) + Transect_f,
                            data=XY_scale,
                            family=gaussian,
                            select=TRUE,
                            method="REML")

# Ok! Actually, looks quite good
draw(mod_cover0p5m_1.0, residuals = TRUE)
summary(mod_cover0p5m_1.0)
appraise(mod_cover0p5m_1.0)
gam.check(mod_cover0p5m_1.0, rep=500)
lines(c(3.0, 3.8), c(3.0, 3.8))

#gam(shannonH ~ s(mean.cvH._woody) +
#      te(Cover0p5m_plot, as.numeric(Transect_f)) + Transect_f,

mod_cvH_cover0p5m_1.0 = gam(shannonH ~ s(mean.cvH._woody) +
               s(Cover0p5m_plot) + Transect_f,
             data=XY_scale,
             family=gaussian,
             select=TRUE,
             method="REML")


# Ok! Actually, looks quite good
draw(mod_cvH_cover0p5m_1.0, residuals = TRUE)
summary(mod_cvH_cover0p5m_1.0)
appraise(mod_cvH_cover0p5m_1.0)
gam.check(mod_cvH_cover0p5m_1.0, rep=500)
lines(c(3.0, 3.8), c(3.0, 3.8))

# Compare to a model built from the corr variables
mod_corr_2.0 = gam(shannonH ~ s(sd.coverG.) +
                     s(sd.sdH._grass) + s(X50thPerc_plot) + s(mean.cvH.) + Transect_f,
                   data=XY_scale,
                   family=gaussian,
                   select=TRUE,
                   method="REML")

summary(mod_corr_2.0)

# # # 2/8/23
# Split by soils
XY_scale_red = XY_scale %>% filter(Soil=='Red')
XY_scale_black = XY_scale %>% filter(Soil=='Black')


mod_corr_2.0_black = gam(shannonH ~ s(sd.coverG.) +
                     s(X50thPerc_plot) + Transect_f,
                   data=XY_scale_black,
                   family=gaussian,
                   select=TRUE,
                   method="REML")

summary(mod_corr_2.0_black)

mod_corr_1.0_black = gam(shannonH ~ s(mean.cvH._woody) +
                           s(Cover0p5m_plot) + Transect_f,
                         data=XY_scale_black,
                         family=gaussian,
                         select=TRUE,
                         method="REML")

draw(mod_corr_1.0_black, residuals = TRUE)
summary(mod_corr_1.0_black)
appraise(mod_corr_1.0_black)
gam.check(mod_corr_1.0_black, rep=500)
lines(c(3.0, 3.8), c(3.0, 3.8))


# # # RED 
mod_corr_2.0_red = gam(shannonH ~ s(sd.coverG.) +
                     s(sd.sdH._grass) + Transect_f,
                   data=XY_scale_red,
                   family=gaussian,
                   select=TRUE,
                   method="REML")

summary(mod_corr_2.0_red)
# doesn't work at all for red



# try using the RF chosen variabels for red
mod_1.0_red = gam(shannonH ~ s(sd.stdpeakh.) +
                    s(sd.sdH._shrub) + Transect_f,
                  data=XY_scale_red,
                  family=gaussian,
                  select=TRUE,
                  method="REML")

summary(mod_1.0_red)
# nope

mod_1.1_red = gam(shannonH ~ s(mean.cvH._woody) + Transect_f,
                  data=XY_scale_red,
                  family=gaussian,
                  select=TRUE,
                  method="REML")

summary(mod_1.1_red)
# barely works at all 

mod_corr_1.2_red = gam(shannonH ~ s(ptoh_plot) +
                         s(sd.VDR.) + Transect_f,
                       data=XY_scale_red,
                       family=gaussian,
                       select=TRUE,
                       method="REML")

summary(mod_corr_1.2_red)
# barely works... 

# so, model for black soils works great
# red soils is a no go right now! 2/8/23


