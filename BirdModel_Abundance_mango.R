# Initial Abundance Model
# For EA Presentation
# 11/14/22 - PB

library(tidyverse)
library(mgcv)
library(arrow)
library(gratia)
library(gridExtra)
# # # 1) Load XY Data
# Note: XY data prepared in "BirdDataPreperation.R"

# Datadir for location of in/out vars
datadir = '/n/home02/pbb/scripts/SelenkayDiversity/data/'

figd = '/n/home02/pbb/scripts/SelenkayDiversity/figs/mango'

# toggle radius of inquiry
radius = 30 

# RDF files keep relevant variable types (factors, characters, etc.)
XY = readRDS(paste0(datadir, 'in/XY_', radius,'mRadius.rds'))
XY_scale = readRDS(paste0(datadir, 'in/XY_scaled_', radius,'mRadius.rds'))


# # # Examine Y variable

# Plot a histogram
hist(XY_scale$Abundance)

# Draw QQPlot for normality
qqnorm(XY_scale$Abundance)
qqline(XY_scale$Abundance)


# # # Top vars - corr with MI index
# horzcover_grass         0.584746
# FHD_plot                0.566022
# mean_sdH                0.560303
# sd_cvH_vegtype_grass    0.523120
# mean_nlayers            0.480641
# sd_PAI_AboveG           0.475477
# mean_maxH               0.474805
# mean_cscore             0.472496
# ptoh_plot               0.467395
# sd_FHD                  0.465876

# make a bunch of plots with variables to explore
# micorr
columns = c('horzcover_grass',
        'FHD_plot',
        'mean_sdH',
        'sd_cvH_vegtype_grass',
        'mean_nlayers',
        'sd_PAI_AboveG',
        'mean_maxH',
        'mean_cscore',
        'ptoh_plot',
        'sd_FHD',
        'X', 'Y')
p <- list()
i = 0
for (c in columns) {
  # if (i == 0) {
  #   i = i+1
  #   }
  i = i+1
  
  XYplot = XY %>% select(c,
                         Abundance,
                         Soil_f,
                         Treatment_f)
  
  p[[i]] = ggplot(data=XYplot) +
    geom_point(aes_string(x=c,
                          y='Abundance',
                          shape='Treatment_f',
                          colour='Soil_f')) + 
    scale_colour_manual(values = c("grey30", "red"))
  
}
 # https://stackoverflow.com/questions/9315611/grid-of-multiple-ggplot2-plots-which-have-been-made-in-a-for-loop
# https://cran.r-project.org/web/packages/egg/vignettes/Ecosystem.html
gridplot = do.call(grid.arrange, p)

gridplot

ggsave(plot = gridplot,
       filename=paste0("/n/home02/pbb/scripts/SelenkayDiversity/figs/Scatter_LidarxAbundance_MIcorr_mango.png"),
       width = 12, height = 8, device='png', dpi=300)

# ggsave(paste0("/n/home02/pbb/scripts/SelenkayDiversity/figs/Scatter_LidarxAbundance_mango.tif"),
#        width = 12, height = 8, device='tiff', dpi=300)


# # # Fit Models
mod_1.0 = gam(Abundance ~ s(mean_sdH, by=Soil_f) + Soil_f + Transect_f,
           data=XY_scale,
           family=poisson,
           select=TRUE,
           method="REML")

summary(mod1)
draw(mod1, residuals = TRUE)
appraise(mod1)
gam.check(mod1, rep=500)
lines(c(100, 220), c(100, 220))


mod_1.1 = gam(Abundance ~ s(sd_cvH_vegtype_grass, by=Soil_f, k=5) + s(X,Y, k=5),
            data=XY_scale,
            family=poisson,
            select=TRUE,
            method="REML")

summary(mod_1.1)
draw(mod_1.1, residuals = TRUE)
appraise(mod_1.1)
gam.check(mod_1.1, rep=500)
lines(c(100, 280), c(100, 280))





# # # Fit Models

mod_1.0 = gam(Abundance ~ s(herbh_plot) + s(cvpeakh_plot) ,
           data=XY_scale,
           family=poisson,
           select=TRUE,
           method="REML")

summary(mod2)


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
