# Initial Richness Model
# For EA Presentation
# 11/14/22 - PB

# Useful links:
# Random Effects in mgcv: https://fromthebottomoftheheap.net/2021/02/02/random-effects-in-gams/
# Intro to GAMs video- https://www.youtube.com/watch?v=sgw4cu8hrZM&t=161s

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

# Draw QQPlot for normality
qqnorm(XY_scale$Richness)
qqline(XY_scale$Richness)
# honestly, it doesn't really fit it. Long Tail

# Interesting- large values of H do not seem to fit into the normal dist.

# Plot a histogram
hist(XY_scale$Richness)

# Black and red soils look like they have diff. distributions...
hist(XY_scale$Richness[XY_scale$Soil=='Black'])
hist(XY_scale$Richness[XY_scale$Soil=='Red'])


# Boxplots of y var by transect
ggplot(XY_scale, aes(x = sd.coverG., y = Richness,
                     group = Transect_f, colour = Transect_f)) +
  geom_boxplot() +
  facet_wrap(~ Soil_f)

# # # Make plots of important variables
# using MI corr vars
# sd_sdHgrasslayer        0.650131
# sd_maxHgrasslayer       0.630723
# sd_herbh                0.560271
# max_maxHgrasslayer      0.552707
# sd_cvH_vegtype_shrub    0.548980
# mean_cscore             0.526452
# cv_maxHgrasslayer       0.523084
# mean_PAI_AboveG         0.518821
# max_VDRpeak             0.508966

# unscaled (just regular XY)
columns = c('sd_sdHgrasslayer',
            'sd_maxHgrasslayer',
            'sd_herbh',
            'max_maxHgrasslayer',
            'sd_cvH_vegtype_shrub',
            'mean_cscore',
            'cv_maxHgrasslayer',
            'mean_PAI_AboveG',
            'max_VDRpeak',
            'sd_CD_AboveG',
            'X', 'Y')
p <- list()
i = 0
for (c in columns) {
  # if (i == 0) {
  #   i = i+1
  #   }
  i = i+1
  
  XYplot = XY %>% select(c,
                         Richness,
                         Soil_f,
                         Treatment_f)
  
  p[[i]] = ggplot(data=XYplot) +
    geom_point(aes_string(x=c,
                          y='Richness',
                          shape='Treatment_f',
                          colour='Soil_f')) + 
    scale_colour_manual(values = c("grey30", "red"))
  
}
# https://stackoverflow.com/questions/9315611/grid-of-multiple-ggplot2-plots-which-have-been-made-in-a-for-loop
# https://cran.r-project.org/web/packages/egg/vignettes/Ecosystem.html
gridplot = do.call(grid.arrange, p)

ggsave(plot = gridplot,
       filename=paste0("/n/home02/pbb/scripts/SelenkayDiversity/figs/Scatter_LidarxRichness_MIcorr_mango.png"),
       width = 12, height = 8, device='png', dpi=300)


# # # 2) Initial Models

# Ran RF in Python to get a list of top vars for Shannon H
# Note that random effects were not included
# in which case, it may be better to run RF on black and red soils separately
# to come up with vars for random effects models.
# But, let's see how it goes first.
# #
# shannonH_vegan MDI Top 6 X Vars:
#   
#   herbh_plot
# max.ptoh.
# mean.cvH._grass
# cv.FHD.
# cover_ground_maxH
# sd.cvH._grass
# mean.VDR.
# cv.VDR.
# sd.VDR.
# sd.ptoh.
# 
# shannonH_vegan Permutation Top 6 X Vars:
#   
# sd.cvH._grass
# mean.VDR.
# max.ptoh.
# sd.cvH._woody
# iqr.VDR.
# mean.cvH.
# mean.sdH._grass
# sd.sdH._grass
# sd.coverG.
# sd.maxH.

# The model below is about as good as it gets...
mod1_rich = gam(Richness ~ 
             s(sd.coverG.) +
             s(sd.sdH._grass),
           data=XY_scale,
           family=gaussian,
           select=TRUE,
           method="REML")

summary(mod1_rich)

mod2 = gam(Richness ~ s(sd.coverG.) +
             s(sd.cvH._woody) +
             s(sd.sdH._grass),
           data=XY_scale,
           family=gaussian,
           select=TRUE,
           method="REML")

summary(mod2)

AIC(mod1, mod2)

ggplot(XY_scale, aes(x = herbh_plot, y = Richness,
                     group = Soil_f, colour = Soil_f)) +
  geom_point()+
  facet_grid(~Soil_f)

ggplot(XY_scale, aes(x = sd.sdH._grass, y = Richness,
                     group = Soil_f, colour = Transect_f)) +
  geom_point()

ggplot(XY_scale, aes(x = sd.cvH._woody, y = Richness,
                     group = Soil_f, colour = Transect_f)) +
  geom_point()

ggplot(XY_scale, aes(x = sd.maxH., y = Richness,
                     group = Soil_f, colour = Transect_f)) +
  geom_point() +
  facet_grid(~Soil_f)



# random effects
mod1 = gam(Richness ~ 
             sd.sdH._grass:Soil_f + s(sd.sdH._grass),
           data=XY_scale,
           family=gaussian,
           select=TRUE,
           method="REML")

summary(mod1)
