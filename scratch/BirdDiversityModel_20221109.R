
# install.packages("mgcv")

library(tidyverse)
library(mgcv)
library(arrow)

# NOTE: Maybe try making the Shannon and Spp Richness again with vegan in R,
# just in case there was something weird with how you manually calculated it 
# in python.

# NOTE: Perform variable selection elsewhere, maybe with random forests, 
# and see if that does better for you

# Helpful sites:
# https://www.rdocumentation.org/packages/dplyr/versions/1.0.10
# https://www.rstudio.com/wp-content/uploads/2015/02/data-wrangling-cheatsheet.pdf

# Datadir for location of in/out vars
datadir = '/n/home02/pbb/scripts/SelenkayDiversity/data/'

# Open up X and Y vars
y = read_csv(file=paste0(datadir, 'in/BirdSurvey/y_vars/SpotSummaryStats_Birds_final.csv'))
X_raw = read_csv(file=paste0(datadir, 'in/Lidar/x_vars/LidarMetrics_30m_20221109.csv'))

# # # Data Prep.

# Regularize Vars
X_scale = X_raw %>%
  select_if(is.numeric) %>%
  scale() %>%
  data.frame()
# NOTE: be careful when scaling something like nlayers, which is actually an interger

# Make an XY df
XY_scale = X_scale %>%
  bind_cols(y)

# rename Shannon H col
XY_scale = XY_scale %>% rename("ShannonH" = "Shannon H")

# Write this scaled df to disk
write.csv(XY_scale, '~/scripts/SelenkayDiversity/data/in/XY_scaled_20221110.csv')

# Note: you have to manually select your columns
# can print them out using select() and "ends_with" "starts_with" or "contains"
select(XY_scale, contains("sd")) %>% colnames()

# # # Fit your first GAM!
# select only wholeplot vars to begin with
mod1 = gam(Abundance ~ s(mean.maxH.) + s(mean.sdH.) + s(sd.maxH.) + 
             s(mean.nlayers.) + s(mean.coverG.),
           data=XY_scale,
           family=poisson,
           select=TRUE,
           method="REML")

summary(mod1)

mod2 = gam(ShannonH ~ s(mean.maxH.) + s(mean.sdH.) + s(sd.maxH.) + 
             s(mean.nlayers.) + s(mean.coverG.),
           data=XY_scale,
           family=gaussian,
           select=TRUE,
           method="REML")

summary(mod2)

mod3 = gam(Richness ~ s(mean.maxH.) + s(mean.sdH.) + s(sd.maxH.) + 
             s(mean.nlayers.) + s(mean.coverG.),
           data=XY_scale,
           family=poisson,
           select=TRUE,
           method="REML")

summary(mod3)


mod4 = gam(ShannonH ~ s(nlayers_plot) +
             s(stdH_plot) +
             s(Cover0p5m_plot) +
             s(X50thPerc_plot)+
             s(X98thPerc_plot),
           data=XY_scale,
           family=gaussian,
           select=TRUE,
           method="REML")

summary(mod4)


mod5 = gam(ShannonH ~ 
             s(sd.coverG.) +
             s(sd.sdH._grass) +
             s(sd.sdH._shrub),
           data=XY_scale,
           family=gaussian,
           select=TRUE,
           method="REML")

summary(mod5)
# Note: std of nlayers and herbh both did well here.
# but sd of coverG is better

mod5 = gam(ShannonH ~ 
             s(sd.coverG.) +
             s(sd.sdH._grass) +
             s(sd.sdH._shrub),
           data=XY_scale,
           family=gaussian,
           select=TRUE,
           method="REML")

summary(mod5)

# Top 10 vars according to mdi index and RF
# herbh_plot          0.153698
# mean.cvH.           0.036883
# sd.cvH._grass       0.036696
# stdH_plot           0.036152
# cv.FHD.             0.031806
# sd.cvH._woody       0.028803
# cover_shrub_maxH    0.027524
# sd.cvH._shrub       0.025965
# mean.cvH._grass     0.025004
# sd.sdH._shrub       0.024486

# Not bad, but not great either!
mod6 = gam(ShannonH ~ 
             s(herbh_plot) +
             s(sd.cvH._woody) +
             s(cv.FHD.) +
             s(cover_shrub_maxH),
           data=XY_scale,
           family=gaussian,
           select=TRUE,
           method="REML")

summary(mod6)
# s(mean.cvH.)
# s(stdH_plot)

# Another run of RF showed these vars
# sd.sdH._grass       0.061094
# cover_shrub_maxH    0.057563
# sd.cvH._woody       0.033542
# cv.VDR.             0.033278
# sd.sdH._shrub       0.031265
# sd.coverG.          0.024531

mod7 = gam(ShannonH ~ 
             s(sd.sdH._grass) +
             s(sd.cvH._woody) +
             s(sd.coverG.),
           data=XY_scale,
           family=gaussian,
           select=TRUE,
           method="REML")

summary(mod7)


# Compare with AIC
AIC(mod6, mod7)

# OK! So abundance is great, richness is not great but doable, and diversity is 
# totally unrelated. Awesome!

# mod1 = gam(Richness ~ s(mean.maxH.) + s(mean.sdH.) + s(sd.maxH.) + 
#              s(mean.nlayers.) + s(mean.coverG.) +
#              s(cover_ground_maxH) + s(cover_grass_maxH) +
#              s(cover_woody_maxH) +
#              s(mean.ptoh.) + s(mean.VDR.) + s(mean.FHD.) +
#              s(mean.meanpeakh.),
#            data=XY_scale,
#            family=poisson,
#            select=TRUE,
#            method="REML")

# Train test split? 

# # # Model Fitting and Eval.

# Whole Plot Run

# Mean Vars Run

