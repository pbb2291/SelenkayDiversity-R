# Shannon H Model - Prediction Model Run on Site 01
# Followed Bottom of the heap code on exptrapolation to a t
# https://fromthebottomoftheheap.net/2020/06/03/extrapolating-with-gams/
# PB 4/20/23
library(tidyverse)
library(mgcv)
library(arrow)
library(lme4)
library(gratia)
library(gridExtra)
library(sf)

# # # 1) Load XY Data
# Note: XY data prepared in "BirdDataPreperation.R"

# Datadir for location of in/out vars
datadir = '/n/home02/pbb/scripts/SelenkayDiversity/data/'

figd = '/n/home02/pbb/scripts/SelenkayDiversity/figs/PredictionsSite01'

tabled = '/n/home02/pbb/scripts/SelenkayDiversity/tables/PredictionsSite01'

# Read in the 260 m pixel shapefile
# with sd_gapsize in it
shape = read_sf('/n/davies_lab/Users/pbb/SelenkayDiversity/data/PredictionSite01/grid_metrics/ModelMetrics.shp')

# toggle radius of inquiry
radius = 130
# radius = 80
# radius = 50

# RDF files keep relevant variable types (factors, characters, etc.)
XY = readRDS(paste0(datadir, 'in/XY_', radius,'mRadius.rds'))
XY_scale = readRDS(paste0(datadir, 'in/XY_scaled_', radius,'mRadius.rds'))

# use the knots argument to cover the range of the prediction data (shannonH)
# See https://fromthebottomoftheheap.net/2020/06/03/extrapolating-with-gams/ 
# for details
knots <- list(x = c(2.75, 3.29, 3.49, 3.8))

# First, train the model:
mod130.5.0 = gam(shannonH ~
                   s(sd_gapsize, bs='bs'),
                 data=XY,
                 family=gaussian,
                 select=FALSE,
                 method="REML",
                 knots=knots)

summary(mod130.5.0)
mod130.5.0$aic
p130.5.0 = draw(mod130.5.0, residuals = TRUE)
p130.5.0
appraise(mod130.5.0)

# # #

# library('rgl')
# plot3d(mod130.5.0)

# create the prediction
predictions = as_tibble(predict.gam(mod130.5.0,
                          shape,
                          type='response',
                          se.fit = TRUE))  %>%
  rename(fit_bs_default = fit, se_bs_default = se.fit)

# Merge into predictions
shape = shape %>% mutate(shannonH_fit = predictions$fit_bs_default)

# Save out as a new shapefile
st_write(shape,
         dsn = paste0(datadir, "/out/PredictionSite01/Shapefile/SelenkaySite01_ShannonPredictionGrid_WGS84UTM37S.shp"),
         layer="Site01Predictions",
         driver="ESRI Shapefile")

# # # PLOT

crit <- qnorm((1 - 0.95)/ 2, lower.tail = FALSE)

# bs_default <- basis(s(sd_gapsize, bs='bs'), knots=knots, data=XY)
# draw(bs_default)
new_data_bs_eg <- bind_cols(shape, predictions) %>%
  pivot_longer(fit_bs_default:se_bs_default, names_sep = '_',
               names_to = c('variable', 'spline', 'penalty')) %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  mutate(upr_ci = fit + (crit * se), lwr_ci = fit - (crit * se))


XYplot = ggplot() +
  geom_ribbon(data = new_data_bs_eg,
              mapping = aes(ymin = lwr_ci, ymax = upr_ci, x = sd_gapsize),
              inherit.aes = FALSE, alpha = 0.2) +
  scale_fill_manual(values = c("#C0D9D9")) + 
  geom_line(data = new_data_bs_eg,
            aes(x=sd_gapsize, y=fit, group = "Prediction")) +
  geom_point(data= XY, aes(x = sd_gapsize,
                           y = shannonH,
                           shape = Treatment,
                           colour = Soil), size=2.3) +
  scale_colour_manual(values = c("grey30", "#D55E00")) + 
  labs(shape='Protected Status',
       colour='Soil') +
  theme(axis.text.y = element_text(colour = "black", size = 12), 
        axis.text.x = element_text(colour = "black", size = 12), 
        legend.text = element_text(size = 12, colour ="black"), 
        legend.position = "right",
        title = element_text(face = "bold", size = 12, colour = "black"),
        legend.title = element_text(size = 12, colour = "black", face = "bold"), 
        legend.key=element_blank())+
  xlab('Variation (SD) in Canopy Gap Size [m]') +
  ylab('Shannon Diversity')

XYplot

ggsave(plot = XYplot,
       filename=paste0(figd, "/PredictionSite01_GAM-130-5-0_XYplot_", radius ,"mRadius_mango.png"),
       width = 9, height = 6, units = "in", device='png', dpi=300)

