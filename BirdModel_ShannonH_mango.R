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

# using gratia
# https://gavinsimpson.github.io/gratia/reference/draw.gam.html

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
radius = 30 

# RDF files keep relevant variable types (factors, characters, etc.)
XY = readRDS(paste0(datadir, 'in/XY_', radius,'mRadius.rds'))
XY_scale = readRDS(paste0(datadir, 'in/XY_scaled_', radius,'mRadius.rds'))

# # # Examine Y variable

# Draw QQPlot for normality
qqnorm(XY_scale$shannonH)
qqline(XY_scale$shannonH)

# Plot a histogram
hist(XY_scale$shannonH)

# Black and red soils look like they have diff. distributions...
hist(XY_scale$shannonH[XY_scale$Soil=='Black'])
hist(XY_scale$shannonH[XY_scale$Soil=='Red'])

# # # Top vars - corr with MI index
# sd_CD_AboveG            0.711085
# sd_sdH_vegtype_grass    0.673687
# cv_CD_AboveG            0.667515
# sd_CD_G                 0.643478
# X50thPerc_plot          0.635744
# sd_nlayers              0.628162
# mean_cvH                0.619774
# Cover0p5m_plot          0.612713
# mean_CD_AboveG          0.607847
# cv_gapsize              0.601971
# make a bunch of plots with variables to explore
# mean_cvH_vegtype_woody added below
# columns = c('sd_CD_AboveG',
#             'sd_sdH_vegtype_grass',
#             'cv_CD_AboveG',
#             'sd_CD_G',
#             'sd_nlayers',
#             'mean_cvH',
#             'Cover0p5m_plot',
#             'mean_CD_AboveG',
#             'cv_gapsize',
#             "X50thPerc_plot",
#             'mean_cvH_vegtype_woody')

# shannonH - spearman corr
# X25thPerc_plot 	0.515102
# sd_FHD 	0.473229
# Cover1p5m_plot 	0.469964
# sd_VDRpeak 	0.468427
# horzcover_woody 	0.460936
# sd_cscore 	0.457575
# Cover0p5m_plot 	0.452869
# X50thPerc_plot 	0.450564
# mean_nlayers 	0.449988
# mean_VDRpeak 	0.449796

# shannonH spearman corr vars
columns = c('X25thPerc_plot',
            'sd_FHD',
            'Cover1p5m_plot',
            'sd_VDRpeak',
            'horzcover_woody',
            'sd_cscore',
            'Cover0p5m_plot',
            'X50thPerc_plot',
            'mean_nlayers',
            'mean_VDRpeak',
            'X', 'Y')
p <- list()
i = 0
for (c in columns) {
  # if (i == 0) {
  #   i = i+1
  #   }
  i = i+1
  
  XYplot = XY_scale %>% select(c,
                         shannonH,
                         Soil_f,
                         Treatment_f)
  
  p[[i]] = ggplot(data=XYplot) +
    geom_point(aes_string(x=c,
                          y='shannonH',
                          shape='Treatment_f',
                          colour='Soil_f')) + 
    scale_colour_manual(values = c("grey30", "red"))
  
}
# https://stackoverflow.com/questions/9315611/grid-of-multiple-ggplot2-plots-which-have-been-made-in-a-for-loop
# https://cran.r-project.org/web/packages/egg/vignettes/Ecosystem.html
gridplot = do.call(grid.arrange, p)

ggsave(plot = gridplot,
       filename=paste0(figd,
                       "/GAMs/Scatter_LidarScaledxShannon_",
                       radius,
                       "mRadius_mango_SpearmanVars.png"),
       width = 10, height = 8, units = "in", device='png', dpi=300)


# unscaled (just regular XY)
# MIcorr vars
# columns = c('sd_CD_AboveG',
#             'sd_sdH_vegtype_grass',
#             'cv_CD_AboveG',
#             'sd_CD_G',
#             'sd_nlayers',
#             'mean_cvH',
#             'horzcover_grass',
#             'Cover0p5m_plot',
#             'mean_CD_AboveG',
#             'cv_gapsize',
#             'mean_cvH_vegtype_woody',
#             'X',
#             'Y')

# Spearman vars
# columns = c('sd_FHD',
#             'Cover1p5m_plot',
#             'sd_VDRpeak',
#             'horzcover_woody',
#             'sd_cscore',
#             'Cover0p5m_plot',
#             'mean_nlayers',
#             'mean_VDRpeak',
#             'X', 'Y')

# pearson corr
# These are mainly linear vars
# spearman seems to grab more of the non-linear
columns = c('sd_CD_AboveG',
            'sd_CD_G',
            'sd_cvpeakh',
            'sd_VDRpeak',
            'sd_FHD',
            'sd_sdHgrasslayer',
            'sd_herbh',
            'sd_maxHgrasslayer',
            "FHD_plot", 
            'cv_cscore', 'X', 'Y')

p <- list()
i = 0
for (c in columns) {
  
  i = i+1
  
  XYplot = XY %>% select(c,
                         shannonH,
                         Soil_f,
                         Treatment_f)
  
  p[[i]] = ggplot(data=XYplot) +
    geom_point(aes_string(x=c,
                          y='shannonH',
                          shape='Treatment_f',
                          colour='Soil_f')) + 
    scale_colour_manual(values = c("grey30", "red"))
  
}
# https://stackoverflow.com/questions/9315611/grid-of-multiple-ggplot2-plots-which-have-been-made-in-a-for-loop
# https://cran.r-project.org/web/packages/egg/vignettes/Ecosystem.html
gridplot = do.call(grid.arrange, p)

ggsave(plot = gridplot,
       filename=paste0(figd,
                       "/GAMs/Scatter_LidarxShannon_",
                       radius,
                       "mRadius_mango_pearsonVars.png"),
       width = 10, height = 8, units = "in", device='png', dpi=300)


# 'sd_FHD',
# 'Cover1p5m_plot',
# 'sd_VDRpeak',
# 'horzcover_woody',
# 'sd_cscore',
mod_1.0 = gam(shannonH ~ s(X25thPerc_plot, k=5) +
                s(horzcover_woody, k=5) +
                s(X,Y, k=10) + 
                Transect_f,
              data=XY_scale,
              family=gaussian,
              select=TRUE,
              method="REML")
summary(mod_1.0)

mod_1.1 = gam(shannonH ~ 
               s(horzcover_woody, k=5, by=Soil_f) +
               s(mean_cvH_vegtype_woody, k=5, by=Soil_f) +
               s(X,Y, k=10),
             data=XY_scale,
             family=gaussian,
             select=TRUE,
             method="REML")
summary(mod_1.1)
p1.1 = draw(mod_1.1, residuals = TRUE)
p1.1

mod_1.1 = gam(shannonH ~ 
                s(horzcover_woody, k=5, by=Soil_f) +
                s(mean_cvH_vegtype_woody, k=5, by=Soil_f) +
                s(X,Y, k=10),
              data=XY_scale,
              family=gaussian,
              select=TRUE,
              method="REML")
summary(mod_1.1)
p1.1 = draw(mod_1.1, residuals = TRUE)
p1.1 

# # # Maybe the top model:


mod_2.1 = gam(shannonH ~ s(sd_CD_AboveG, k=5, by=Soil_f) +
                s(mean_cvH_vegtype_woody, k=5, by=Soil_f) +
                s(X, Y, k=10) + 
                Transect_f,
              data=XY_scale,
              family=gaussian,
              select=TRUE,
              method="REML")


p2.1 = draw(mod_2.1, residuals = TRUE)
p2.1 

ggsave(plot = p2.1,
       filename=paste0(figd,
                       "/GAMs/GAM2p1_PartialEffects_", radius ,"mRadius_mango.png"),
       width = 14, height = 8, units = "in", device='png', dpi=300)

summary(mod_2.1)
a2.1 = appraise(mod_2.1)
a2.1

ggsave(plot = a2.1,
       filename=paste0(figd,
                       "/GAMs/GAM2p1_Appraise_", radius ,"mRadius_mango.png"),
       width = 12, height = 8, units = "in", device='png', dpi=300)

gam.check(mod_2.1, rep=500)
lines(c(3.0, 3.8), c(3.0, 3.8))

# # other top models

mod_1.0 = gam(shannonH ~ s(mean_cvH_vegtype_woody, k=5, by=Soil_f) +
                              s(Cover0p25m_plot, k=5, by=Soil_f) +
                            s(X,Y),
                            data=XY_scale,
                            family=gaussian,
                            select=TRUE,
                            method="REML")


# Ok! Actually, looks quite good
draw(mod_1.0, residuals = TRUE)
summary(mod_1.0)
appraise(mod_1.0)
gam.check(mod_1.0, rep=500)
lines(c(3.0, 3.8), c(3.0, 3.8))



mod_2.0 = gam(shannonH ~ s(sd_CD_AboveG,k=5, by=Soil_f) +
                     s(mean_cvH_vegtype_woody, k=15, by=Soil_f) +
                     s(X,Y, k=5),
                   data=XY_scale,
                   family=gaussian,
                   select=TRUE,
                   method="REML")

summary(mod_2.0)
draw(mod_2.0, residuals = TRUE)
summary(mod_2.0)
appraise(mod_2.0)
gam.check(mod_2.0, rep=500)
lines(c(3.0, 3.8), c(3.0, 3.8))




mod_3.0 = gam(shannonH ~ s(mean_cvH_vegtype_woody, k=3, by=Soil_f) +
                s(Cover0p25m_plot, k=3, by=Soil_f) +
                s(sd_CD_AboveG,k=3, by=Soil_f) +
                s(X,Y, k=5),
              data=XY_scale,
              family=gaussian,
              select=TRUE,
              method="REML")

draw(mod_3.0, residuals = TRUE)
summary(mod_3.0)
appraise(mod_3.0)
gam.check(mod_3.0, rep=500)
lines(c(3.0, 3.8), c(3.0, 3.8))



mod_3.1 = gam(shannonH ~ s(mean_cvH_vegtype_woody, k=3, by=Soil_f) +
                s(Cover0p25m_plot, k=3, by=Soil_f) +
                s(sd_CD_AboveG,k=3, by=Soil_f),
              data=XY_scale,
              family=gaussian,
              select=TRUE,
              method="REML")

draw(mod_3.1, residuals = TRUE)
summary(mod_3.1)
appraise(mod_3.1)
gam.check(mod_3.1, rep=500)
lines(c(3.0, 3.8), c(3.0, 3.8))

# gam.vcomp(mod_1.0)

# # # Fit models - from simplest to most complex
mod_cvH_1.0 = gam(shannonH ~ s(mean_cvH_vegtype_woody) + Transect_f,
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

# included spatial autocorr in a different way - 2/13
mod_cvH_1.1 = gam(shannonH ~ s(mean_cvH_vegtype_woody) + s(X,Y),
                  data=XY_scale,
                  family=gaussian,
                  select=TRUE,
                  method="REML")

# Ok! Actually, looks quite good
draw(mod_cvH_1.1, residuals = TRUE)
summary(mod_cvH_1.1)
appraise(mod_cvH_1.1)
gam.check(mod_cvH_1.1, rep=500)
lines(c(3.0, 3.8), c(3.0, 3.8))

# Add another var
mod_cover0p5m_1.0 = gam(shannonH ~ s(Cover0p5m_plot) + Transect_f,
                            data=XY_scale,
                            family=gaussian,
                            select=TRUE,
                            method="REML")


draw(mod_cover0p5m_1.0, residuals = TRUE)
summary(mod_cover0p5m_1.0)
appraise(mod_cover0p5m_1.0)
gam.check(mod_cover0p5m_1.0, rep=500)
lines(c(3.0, 3.8), c(3.0, 3.8))

#gam(shannonH ~ s(mean.cvH._woody) +
#      te(Cover0p5m_plot, as.numeric(Transect_f)) + Transect_f,

mod_cvH_cover0p5m_1.0 = gam(shannonH ~ s(mean_cvH_vegtype_woody) +
               s(Cover0p25m_plot) +  Transect_f,
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


# with new spatial autocorr
mod_cvH_cover0p5m_1.1 = gam(shannonH ~ s(mean_cvH_vegtype_woody) +
                              s(sd_CD_AboveG) + s(X,Y),
                            data=XY_scale,
                            family=gaussian,
                            select=TRUE,
                            method="REML")


# Ok! Actually, looks quite good
draw(mod_cvH_cover0p5m_1.1, residuals = TRUE)
summary(mod_cvH_cover0p5m_1.1)
appraise(mod_cvH_cover0p5m_1.1)
gam.check(mod_cvH_cover0p5m_1.1, rep=500)
lines(c(3.0, 3.8), c(3.0, 3.8))

# Compare to a model built from the corr variables
mod_corr_2.0 = gam(shannonH ~ s(sd_CD_AboveG) + s(X50thPerc_plot) + Transect_f,
                   data=XY_scale,
                   family=gaussian,
                   select=TRUE,
                   method="REML")

summary(mod_corr_2.0)

mod_corr_2.1 = gam(shannonH ~ s(sd_CD_AboveG,k=5) +
                              s(mean_cvH_vegtype_woody, k=5) +
                              s(X,Y),
                   data=XY_scale,
                   family=gaussian,
                   select=TRUE,
                   method="REML")

summary(mod_corr_2.1)
draw(mod_corr_2.1, residuals = TRUE)
summary(mod_corr_2.1)
appraise(mod_corr_2.1)
gam.check(mod_corr_2.1, rep=500)
lines(c(3.0, 3.8), c(3.0, 3.8))


# # # 2/8/23
# Split by soils
XY_scale_red = XY_scale %>% filter(Soil=='Red')
XY_scale_black = XY_scale %>% filter(Soil=='Black')

# # # BLACK SOIL best models 2/9 
mod_corr_1.0_black = gam(shannonH ~ s(mean_cvH_vegtype_woody) +
                           s(Cover0p25m_plot) + Transect_f,
                         data=XY_scale_black,
                         family=gaussian,
                         select=TRUE,
                         method="REML")

draw(mod_corr_1.0_black, residuals = TRUE)
summary(mod_corr_1.0_black)
appraise(mod_corr_1.0_black)
gam.check(mod_corr_1.0_black, rep=500)
lines(c(3.0, 3.8), c(3.0, 3.8))

mod_corr_1.0_black = gam(shannonH ~ s(mean_cvH_vegtype_woody) +
                           s(sd_CD_AboveG) + Transect_f,
                         data=XY_scale_black,
                         family=gaussian,
                         select=TRUE,
                         method="REML")

draw(mod_corr_1.0_black, residuals = TRUE)
summary(mod_corr_1.0_black)
appraise(mod_corr_1.0_black)
gam.check(mod_corr_1.0_black, rep=500)
lines(c(3.0, 3.8), c(3.0, 3.8))


# # Spatial autocorre added
mod_corr_1.0_black = gam(shannonH ~ s(mean_cvH_vegtype_woody, k=5) +
                           s(sd_CD_AboveG, k=5) + s(X, Y, k=5),
                         data=XY_scale_black,
                         family=gaussian,
                         select=TRUE,
                         method="REML")

draw(mod_corr_1.0_black, residuals = TRUE)
summary(mod_corr_1.0_black)
appraise(mod_corr_1.0_black)
gam.check(mod_corr_1.0_black, rep=500)
lines(c(3.0, 3.8), c(3.0, 3.8))

mod_corr_1.1_black = gam(shannonH ~ s(mean_cvH_vegtype_woody, k=3) +
                           s(sd_CD_AboveG, k=5) + 
                           s(X, Y, k=5) + 
                           Transect_f,
                         data=XY_scale_black,
                         family=gaussian,
                         select=TRUE,
                         method="REML")

draw(mod_corr_1.1_black, residuals = TRUE)
summary(mod_corr_1.1_black)
appraise(mod_corr_1.1_black)
gam.check(mod_corr_1.1_black, rep=500)
lines(c(3.0, 3.8), c(3.0, 3.8))



mod_corr_1.0_black = gam(shannonH ~ s(mean_cvH_vegtype_woody) +
                           s(sd_CD_AboveG) + Transect_f,
                         data=XY_scale_black,
                         family=gaussian,
                         select=TRUE,
                         method="REML")

draw(mod_corr_1.0_black, residuals = TRUE)
summary(mod_corr_1.0_black)
appraise(mod_corr_1.0_black)
gam.check(mod_corr_1.0_black, rep=500)
lines(c(3.0, 3.8), c(3.0, 3.8))


mod_corr_1.0_black = gam(shannonH ~ s(mean_cvH_vegtype_woody) +
                           s(cv_CD_AboveG) + Transect_f,
                         data=XY_scale_black,
                         family=gaussian,
                         select=TRUE,
                         method="REML")

draw(mod_corr_1.0_black, residuals = TRUE)
summary(mod_corr_1.0_black)
appraise(mod_corr_1.0_black)
gam.check(mod_corr_1.0_black, rep=500)
lines(c(3.0, 3.8), c(3.0, 3.8))
### 
### 


mod_corr_2.0_black = gam(shannonH ~ s(sd_sdH_vegtype_grass) +
                     s(X25thPerc_plot) + Transect_f,
                   data=XY_scale_black,
                   family=gaussian,
                   select=TRUE,
                   method="REML")

summary(mod_corr_2.0_black)
draw(mod_corr_2.0_black, residuals=TRUE)

# Spearan black soil correlation
# shannonH
# mean_stdpeakh    0.765385
# Cover0p5m_plot   0.762308
# mean_gapsize     0.760000
# sd_CD_AboveG     0.751538
# mean_CD_AboveG   0.751538
# Cover0p25m_plot  0.746923
# Cover1p5m_plot   0.745385
# X25thPerc_plot   0.739231
# mean_nlayers     0.734615
# horzcover_woody  0.733846

mod_spearman_1.0_black = gam(shannonH ~ s(mean_stdpeakh, k=5) +
                           s(Cover0p5m_plot, k=5) + s(mean_gapsize, k=5) +
                             s(sd_CD_AboveG, k=5) + s(X,Y, k=5),
                         data=XY_scale_black,
                         family=gaussian,
                         select=TRUE,
                         method="REML")

draw(mod_spearman_1.0_black , residuals = TRUE)
summary(mod_spearman_1.0_black )
appraise(mod_spearman_1.0_black )
gam.check(mod_spearman_1.0_black , rep=500)
lines(c(3.0, 3.8), c(3.0, 3.8))

mod_spearman_1.1_black = gam(shannonH ~ 
                               s(sd_CD_AboveG, k=5) + s(X,Y, k=5),
                             data=XY_scale_black,
                             family=gaussian,
                             select=TRUE,
                             method="REML")

draw(mod_spearman_1.1_black , residuals = TRUE)
summary(mod_spearman_1.1_black )
appraise(mod_spearman_1.1_black )
gam.check(mod_spearman_1.1_black , rep=500)
lines(c(3.0, 3.8), c(3.0, 3.8))


mod_3.0_black = gam(shannonH ~ s(X25thPerc_plot) + s(sd_VDR) + Transect_f,
                         data=XY_scale_black,
                         family=gaussian,
                         select=TRUE,
                         method="REML")

summary(mod_3.0_black)
draw(mod_3.0_black, residuals = TRUE)
appraise(mod_3.0_black)
gam.check(mod_3.0_black, rep=500)
lines(c(3.0, 3.8), c(3.0, 3.8))

ggplot() + 
  geom_point(data=XY_scale,
             aes(x=X25thPerc_plot,
                 y=shannonH))

ggplot() + 
  geom_point(data=XY_scale,
             aes(x=sd_VDR,
                 y=shannonH))

ggplot() + 
  geom_point(data=XY_scale,
             aes(x=cv_CD_AboveG,
                 y=shannonH))

ggplot() + 
  geom_point(data=XY_scale,
             aes(x=cv_CD_AboveG,
                 y=shannonH))


ggplot() + 
  geom_point(data=XY_scale,
             aes(x=Cover0p05m_plot,
                 y=shannonH))

ggplot() + 
  geom_point(data=XY_scale,
             aes(x=Cover0p25m_plot,
                 y=shannonH))

ggplot() + 
  geom_point(data=XY_scale,
             aes(x=sd_maxpeakh,
                 y=shannonH))



# # # RED 
mod_corr_2.0_red = gam(shannonH ~ s(sd_CD_G) +
                     s(sd_sdH_vegtype_grass) + Transect_f,
                   data=XY_scale_red,
                   family=gaussian,
                   select=TRUE,
                   method="REML")

summary(mod_corr_2.0_red)
# doesn't work at all for red



# try using the RF chosen variabels for red
mod_1.0_red = gam(shannonH ~ s(sd_stdpeakh) +
                    s(sd_sdH_vegtype_shrub) + Transect_f,
                  data=XY_scale_red,
                  family=gaussian,
                  select=TRUE,
                  method="REML")

summary(mod_1.0_red)
# nope

mod_1.1_red = gam(shannonH ~ s(sd_PAI_AboveG) + Transect_f,
                  data=XY_scale_red,
                  family=gaussian,
                  select=TRUE,
                  method="REML")

summary(mod_1.1_red)
# barely works at all 

mod_corr_1.2_red = gam(shannonH ~ s(cv_CD_Ggrasslayer) +
                         s(cv_maxHgrasslayer) + Transect_f,
                       data=XY_scale_red,
                       family=gaussian,
                       select=TRUE,
                       method="REML")

summary(mod_corr_1.2_red)
# barely works... 

ggplot() + geom_point(data=XY_scale_red, aes(x=cv_maxHgrasslayer, y=shannonH))

# so, model for black soils works great
# red soils is a no go right now! 2/8/23

mod_5.0 =  gam(shannonH ~ s(Cover0p25m_plot) +
                 s(cvpeakh_plot) + Transect_f,
               data=XY_scale,
               family=gaussian,
               select=TRUE,
               method="REML")

summary(mod_5.0)


mod_5.0 =  gam(shannonH ~ s(cv_CD_AboveG) + Transect_f,
               data=XY_scale,
               family=gaussian,
               select=TRUE,
               method="REML")

summary(mod_5.0)
