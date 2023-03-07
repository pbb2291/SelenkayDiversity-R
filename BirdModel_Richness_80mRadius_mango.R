# library("devtools"); install_github("lme4/lme4",dependencies=TRUE)
library(devtools)
# devtools::install_github("sjPlot/devel")
library(tidyverse)
library(mgcv)
library(arrow)
library(lme4)
library(gratia)
library(gridExtra)
library(gtsummary)
library(broom)

# # # 1) Load XY Data
# Note: XY data prepared in "BirdDataPreperation.R"

# Datadir for location of in/out vars
datadir = '/n/home02/pbb/scripts/SelenkayDiversity/data/'

figd = '/n/home02/pbb/scripts/SelenkayDiversity/figs/mango'

tabled = '/n/home02/pbb/scripts/SelenkayDiversity/tables/mango/GAMs'

# toggle radius of inquiry
radius = 80

# RDF files keep relevant variable types (factors, characters, etc.)
XY = readRDS(paste0(datadir, 'in/XY_', radius,'mRadius.rds'))
XY_scale = readRDS(paste0(datadir, 'in/XY_scaled_', radius,'mRadius.rds'))

XY_scale = XY_scale %>% mutate(logrich = log(Richness))
# # # All Soils
# 
# red
# horzcover_grass
# sd_sdH_vegtypegrass
# sd_cvH_vegtypegrass
# sd_FHD
# 25thPerc_plot
# 50th
# 75th
# meanH_plot
# 
# black
# sd_nlayers
# cv_nlayers
# meanH_plot
# stdH_plot
# Cover1p5m_plot

# both
# mean_PAI_G
# sd_PAI_G
# herbh_plot
# mean_maxHgrasslayer
# mean_CD_Ggrasslayer
# cv_CD_Ggrasslayer

# this is the only one thatcomes up as iportant 
# but it's really... not very informative
# s(sd_nlayers, k=5) +

# NOTE: log richness works a lot better
# follows Bae et al. 2019 - radar gam paper

rich80.1.0 = gam(logrich ~ s(meanH_plot,
                              by=Soil_f) +
                   s(sd_nlayers,
                     by=Soil_f) +
                   s(horzcover_grass, by=Soil_f, k=5) +
                   s(X,Y, k=5),
                       data=XY_scale,
                       family=gaussian,
                       select=TRUE,
                       method="REML")

summary(rich80.1.0)
p80.1.0 = draw(rich80.1.0, residuals = TRUE)
p80.1.0

ggsave(plot = p80.1.0,
       filename=paste0(figd,
                       "/GAMs/GAM-80-1-0-LogRichness_PartialEffects_", radius ,"mRadius_mango.png"),
       width = 14, height = 8, units = "in", device='png', dpi=300)

a80.1.0 = appraise(rich80.1.0)
a80.1.0

ggsave(plot = a80.1.0,
       filename=paste0(figd,
                       "/GAMs/GAM-80-1-0-LogRichness_Appraise_", radius ,"mRadius_mango.png"),
       width = 12, height = 8, units = "in", device='png', dpi=300)

gam.check(rich80.1.0, rep=500)
lines(c(3.0, 3.8), c(3.0, 3.8))

# Same some important results
vartbl_rich80.1.0 = tidy(rich80.1.0)
write_csv(vartbl_rich80.1.0,
          paste0(tabled, '/VarTable_LogRich80-1-0_', radius, 'mRadius.csv'))

# use sink to save summary outputs in a txt file
sink(paste0(tabled, '/ModelSummary_LogRich80-1-0_', radius, 'mRadius.txt'))
summary(rich80.1.0)
print(paste0("AIC = ", rich80.1.0$aic))
sink()



# library(stargazer)
# 
# latex = stargazer(rich80.1.0,summary=TRUE)



# # # By Soils Models

# RED
# horzcover_grass
# sd_sdH_vegtypegrass
# sd_cvH_vegtypegrass
# sd_FHD
# 25thPerc_plot
# 50th
# 75th
# meanH_plot

# Split by soils
XY_scale_red = XY_scale %>% filter(Soil=='Red')

# s(X75thPerc_plot) +
#   s(sd_sdH_vegtype_grass) +

rich80.1.0.red = gam(Richness ~ s(X75thPerc_plot) +
                   s(X,Y, k=5),
                 data=XY_scale_red,
                 family=poisson,
                 select=TRUE,
                 method="REML")

summary(rich80.1.0.red)
p80.1.0.red = draw(rich80.1.0.red, residuals = TRUE)
p80.1.0.red

ggsave(plot = p80.1.0,
       filename=paste0(figd,
                       "/GAMs/GAM-80-1-0-Red-Richness_PartialEffects_",
                       radius ,"mRadius_mango.png"),
       width = 14, height = 8, units = "in", device='png', dpi=300)

a80.1.0.red = appraise(rich80.1.0.red)
a80.1.0.red

ggsave(plot = a80.1.0.red,
       filename=paste0(figd,
                       "/GAMs/GAM-80-1-0-Red-Richness_Appraise_",
                       radius ,"mRadius_mango.png"),
       width = 12, height = 8, units = "in", device='png', dpi=300)

gam.check(rich80.1.0.red, rep=500)
lines(c(3.0, 3.8), c(3.0, 3.8))

library(gtsummary)
