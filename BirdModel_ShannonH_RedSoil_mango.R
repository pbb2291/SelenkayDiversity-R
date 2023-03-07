
library(devtools)
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
XY_scale_red = XY_scale %>% filter(Soil=='Red')

# # # RED
# These vars are the top 5 spearman corr
mod130.1.0.red = gam(shannonH ~ s(gapsize_plot, k=3) +
                       s(maxpeakh_plot, k=3) +
                       s(ptoh_plot, k=3) +
                       s(VDRpeak_plot, k=3) +
                       s(X,Y, k=3),
                       data=XY_scale_red,
                       family=gaussian,
                       select=TRUE,
                       method="REML")

summary(mod130.1.0.red)
p130.1.0.red = draw(mod130.1.0.red, residuals = TRUE)
p130.1.0.red

ggsave(plot = p130.1.0.red,
       filename=paste0(figd,
                       "/GAMs/GAM-130-1-0-Red_PartialEffects_", radius ,"mRadius_mango.png"),
       width = 14, height = 8, units = "in", device='png', dpi=300)

a130.1.0.red = appraise(mod130.1.0.red)
a130.1.0.red

ggsave(plot = a130.1.0.red,
       filename=paste0(figd,
                       "/GAMs/GAM-130-1-0-Red_Appraise_", radius ,"mRadius_mango.png"),
       width = 12, height = 8, units = "in", device='png', dpi=300)

gam.check(mod130.1.0.red, rep=500)
lines(c(3.0, 3.8), c(3.0, 3.8))
# # # There is significance here, but it looks like the X,Y term
# is explaining most of the variation

# # # Simpler version
mod130.1.1.red = gam(shannonH ~ 
                       s(maxpeakh_plot, k=5) +
                       s(VDRpeak_plot, k=5) +
                       s(X,Y, k=5),
                     data=XY_scale_red,
                     family=gaussian,
                     select=TRUE,
                     method="REML")

summary(mod130.1.1.red)
p130.1.1.red = draw(mod130.1.1.red, residuals = TRUE)
p130.1.1.red

ggsave(plot = p130.1.1.red,
       filename=paste0(figd,
                       "/GAMs/GAM-130-1-1-Red_PartialEffects_", radius ,"mRadius_mango.png"),
       width = 14, height = 8, units = "in", device='png', dpi=300)

a130.1.1.red = appraise(mod130.1.1.red)
a130.1.1.red

ggsave(plot = a130.1.1.red,
       filename=paste0(figd,
                       "/GAMs/GAM-130-1-1-Red_Appraise_", radius ,"mRadius_mango.png"),
       width = 12, height = 8, units = "in", device='png', dpi=300)

gam.check(mod130.1.1.red, rep=500)
lines(c(3.0, 3.8), c(3.0, 3.8))
mod130.1.1.red$aic

# Using leave-one-out cross validation to get the RMSE
# see: https://stats.stackexchange.com/questions/35602/gam-cross-validation-to-test-prediction-error
library(caret)
# 
# XY_scale_red_EN = XY_scale_red %>% rename("Easting" = "X",
#                                           "Northing" = "Y")
# 
# b <- train(shannonH ~ 
#              maxpeakh_plot +
#              VDRpeak_plot + Easting:Northing, 
#            data = XY_scale_red_EN,
#            method = "gam",
#            trControl = trainControl(method = "LOOCV", number = 1, repeats = 1),
#            tuneGrid = data.frame(method = "GCV.Cp", select = FALSE)
# )

# # #



# Now try top 5 MI corr
mod130.2.0.red = gam(shannonH ~ s(mean_gapsize, k=3) +
                       s(sd_gapsize, k=3) +
                       s(mean_stdpeakh, k=3) +
                       s(meanH_plot, k=3) +
                       s(Cover1p5m_plot, k=3) +
                       s(X,Y, k=3),
                     data=XY_scale_red,
                     family=gaussian,
                     select=TRUE,
                     method="REML")

summary(mod130.2.0.red)
p130.2.0.red = draw(mod130.2.0.red, residuals = TRUE)
p130.2.0.red

ggsave(plot = p130.2.0.red,
       filename=paste0(figd,
                       "/GAMs/GAM-130-2-0-Red_PartialEffects_", radius ,"mRadius_mango.png"),
       width = 14, height = 8, units = "in", device='png', dpi=300)

a130.2.0.red = appraise(mod130.2.0.red)
a130.2.0.red

ggsave(plot = a130.2.0.red,
       filename=paste0(figd,
                       "/GAMs/GAM-130-2-0-Red_Appraise_", radius ,"mRadius_mango.png"),
       width = 12, height = 8, units = "in", device='png', dpi=300)

gam.check(mod130.2.0.red, rep=500)
lines(c(3.0, 3.8), c(3.0, 3.8))
# Ok, these don't work at all
# that said, they seem to peak around the 50 m radius
# (but remain high for 130 m radius)
