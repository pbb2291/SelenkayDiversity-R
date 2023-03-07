
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
radius = 50

# RDF files keep relevant variable types (factors, characters, etc.)
XY = readRDS(paste0(datadir, 'in/XY_', radius,'mRadius.rds'))
XY_scale = readRDS(paste0(datadir, 'in/XY_scaled_', radius,'mRadius.rds'))


# # # All Soils
# These vars are the top 5 spearman corr
mod50.1.0 = gam(shannonH ~ s(mean_gapsize, k=3) +
                         s(sd_gapsize, k=3) +
                         s(mean_stdpeakh, k=3) +
                         s(meanH_plot, k=3) +
                         s(Cover1p5m_plot, k=3) +
                         s(X,Y, k=3),
                       data=XY_scale,
                       family=gaussian,
                       select=TRUE,
                       method="REML")

summary(mod50.1.0)
p50.1.0 = draw(mod50.1.0, residuals = TRUE)
p50.1.0

ggsave(plot = p50.1.0,
       filename=paste0(figd,
                       "/GAMs/GAM-50-1-0_PartialEffects_", radius ,"mRadius_mango.png"),
       width = 14, height = 8, units = "in", device='png', dpi=300)

a50.1.0 = appraise(mod50.1.0)
a50.1.0

ggsave(plot = a50.1.0,
       filename=paste0(figd,
                       "/GAMs/GAM-50-1-0_Appraise_", radius ,"mRadius_mango.png"),
       width = 12, height = 8, units = "in", device='png', dpi=300)

gam.check(mod50.1.0, rep=500)
lines(c(3.0, 3.8), c(3.0, 3.8))

# # # Simpler model

# # # All Soils
# These vars are the top 5 spearman corr
mod50.2.0 = gam(shannonH ~
                  s(sd_gapsize, k=3) +
                  s(meanH_plot, k=3) +
                  s(X,Y, k=3),
                data=XY_scale,
                family=gaussian,
                select=TRUE,
                method="REML")

summary(mod50.2.0)
p50.2.0 = draw(mod50.2.0, residuals = TRUE)
p50.2.0

ggsave(plot = p50.2.0,
       filename=paste0(figd,
                       "/GAMs/GAM-50-2-0_PartialEffects_", radius ,"mRadius_mango.png"),
       width = 14, height = 8, units = "in", device='png', dpi=300)

a50.2.0 = appraise(mod50.2.0)
a50.2.0

ggsave(plot = a50.2.0,
       filename=paste0(figd,
                       "/GAMs/GAM-50-2-0_Appraise_", radius ,"mRadius_mango.png"),
       width = 12, height = 8, units = "in", device='png', dpi=300)

gam.check(mod50.2.0, rep=500)
lines(c(3.0, 3.8), c(3.0, 3.8))



# # # 

mod50.3.0 = gam(shannonH ~
                  s(sd_FHD, k=3) +
                  s(X25thPerc_plot, k=3) +
                  s(X50thPerc_plot, k=3) +
                  s(X100thPerc_plot, k=3) +
                  s(cscore_plot, k=3) +
                  s(X,Y, k=3),
                data=XY_scale,
                family=gaussian,
                select=TRUE,
                method="REML")

summary(mod50.3.0)
p50.3.0 = draw(mod50.3.0, residuals = TRUE)
p50.3.0

ggsave(plot = p50.3.0,
       filename=paste0(figd,
                       "/GAMs/GAM-50-3-0_PartialEffects_", radius ,"mRadius_mango.png"),
       width = 14, height = 8, units = "in", device='png', dpi=300)

a50.3.0 = appraise(mod50.3.0)
a50.3.0

ggsave(plot = a50.3.0,
       filename=paste0(figd,
                       "/GAMs/GAM-50-3-0_Appraise_", radius ,"mRadius_mango.png"),
       width = 12, height = 8, units = "in", device='png', dpi=300)

gam.check(mod50.3.0, rep=500)
lines(c(3.0, 3.8), c(3.0, 3.8))
