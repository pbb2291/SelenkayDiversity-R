# Bird Models Shannon
# Scraps
# PB 11/15/2022



# # # 2) Initial Models

# Ran RF in Python to get a list of top vars for Shannon H
# Note that random effects were not included
# in which case, it may be better to run RF on black and red soils separately
# to come up with vars for random effects models.
# But, let's see how it goes first.
# #
# shannonH MDI Top 6 X Vars:
#   
#   cv.VDR.
# max.cvpeakh.
# max.ptoh.
# cv.cscore.
# sd.ptoh.
# max.VDRpeak.
# mean.cvH._woody
# sd.cvH.
# mean.cvH._shrub
# mean.VDR.
# 
# 
# shannonH Permutation Top 6 X Vars:
#   
# sd.cvH.
# sd.sdH._shrub
# cv.maxpeakh.
# mean.cvH._woody
# sd.ptoh.
# max.VDRpeak.
# FHD_plot
# max.cscore.
# mean.herbh.
# sd.FHD.

mod1 = gam(shannonH ~ s(cv.maxpeakh., by=Soil_f) +
             s(mean.cvH._woody) +
             Transect_f,
           data=XY_scale,
           family=gaussian,
           select=TRUE,
           method="REML")

summary(mod1)

mod2 = gam(shannonH ~ s(cv.maxpeakh., by=Soil_f) +
             s(mean.cvH._woody, by=Soil_f) + 
             s(Transect_f, bs="re"),
           data=XY_scale,
           family=gaussian,
           select=TRUE,
           method="REML")

summary(mod2)

mod2 = gam(shannonH ~ s(cv.maxpeakh., by=Soil_f) +
             s(mean.cvH._woody, by=Soil_f) + 
             cv.maxpeakh.:Transect_f,
           data=XY_scale,
           family=gaussian,
           select=TRUE,
           method="REML")

summary(mod2)

mod3 = gam(shannonH ~ sd.herbh.:Transect_f +
             s(mean.cvH._woody, by=Soil_f),
           data=XY_scale,
           family=gaussian,
           select=TRUE,
           method="REML")

summary(mod3)

# Compare models with AIC
AIC(mod2, mod3)

mod4 = gam(shannonH ~ s(mean.cvH._woody, by=Soil_f) +
             Cover0p5m_plot:Transect_f,
           data=XY_scale,
           family=gaussian,
           select=TRUE,
           method="REML")

summary(mod4)



mod5 = gam(shannonH ~ mean.cvH._woody:Transect_f +
             s(sd.herbh., by=Soil_f),
           data=XY_scale,
           family=gaussian,
           select=TRUE,
           method="REML")

summary(mod5)

AIC(mod3, mod5)


mod6 = gam(shannonH ~ sd.herbh.:Transect_f +
             s(mean.cvH._woody, by=Soil_f) + 
             s(sd.herbh., by=Soil_f),
           data=XY_scale,
           family=gaussian,
           select=TRUE,
           method="REML")

summary(mod6)


mod7 = gam(shannonH ~ s(mean.cvH._woody, by=Soil_f) +
             Cover0p05m_plot:Transect_f,
           data=XY_scale,
           family=gaussian,
           select=TRUE,
           method="REML")

summary(mod7)

mod8 = gam(shannonH ~ s(mean.cvH._woody, by=Soil_f) +
             mean.coverG.:Transect_f,
           data=XY_scale,
           family=gaussian,
           select=TRUE,
           method="REML")

summary(mod8)

AIC(mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8)


# try a tensor for an interaction effect
mod4.1 = gam(shannonH ~ s(mean.cvH._woody) +
               te(Cover0p5m_plot, as.numeric(Transect_f)) + Transect_f,
             data=XY_scale,
             family=gaussian,
             select=TRUE,
             method="REML")

summary(mod4.1)

gam.check(mod4.1, rep=500)
lines(c(3.0, 3.8), c(3.0, 3.8))

library(gratia)
appraise(mod4.1)


# test out which of these vars has a stronger relationship
mod4.2 = gam(shannonH ~ s(mean.cvH._woody) + Transect_f,
             data=XY_scale,
             family=gaussian,
             select=TRUE,
             method="REML")

summary(mod4.2)


mod4.3 = gam(shannonH ~ s(Cover0p5m_plot) + Transect_f,
             data=XY_scale,
             family=gaussian,
             select=TRUE,
             method="REML")

summary(mod4.3)


mod4.4 = gam(shannonH ~ s(sd.maxH.) + Transect_f,
             data=XY_scale,
             family=gaussian,
             select=TRUE,
             method="REML")

summary(mod4.4)

AIC(mod4.2, mod4.3, mod4.1, mod4.4)
# Combined model with tensor is best, 
# but let's try it with just 1 var


# Try this
mod4.te1 = gam(shannonH ~ s(mean.cvH._woody, by=Soil_f),
               data=XY_scale,
               family=gaussian,
               select=TRUE,
               method="REML")

summary(mod4.te1)

plot(mod4.te1, pages=1)

# try a tensor for an interaction effect
mod4.5 = gam(shannonH ~ s(Cover0p5m_plot) +
               s(Cover0p5m_plot, by=Soil_f) + Soil_f,
             data=XY_scale,
             family=gaussian,
             select=TRUE,
             method="REML")

summary(mod4.5)

gam.check(mod4.5)

# # # MOD 4 is the winner! 
# shannonH ~ shannonH ~ s(mean.cvH._woody, by=Soil_f) +
#                       Cover0p5m_plot:Transect_f


# ggplot(XY_scale, aes(x = Cover0p5m_plot, y = shannonH,
#                      group = Soil_f, colour = Transect_f)) +
#   geom_point() +
#   facet_wrap(~ Soil_f) + 
#   labs(colour = "Transect Number") + 
#   xlab("Canopy Cover above 0.5m") +
#   ylab("Shannon H")

ggplot(XY_scale, aes(x = Cover0p5m_plot, y = shannonH,
                     group = Soil_f, colour = Soil_f)) +
  geom_point(size=2.5) +
  labs(colour = "Soil Type") + 
  xlab("Canopy Cover above 0.5m") +
  ylab("Shannon H") + scale_color_hue(direction = -1)


ggplot(XY_scale, aes(x = mean.cvH._woody, y = shannonH,
                     group = Soil_f, colour = Soil_f)) +
  geom_point(size=2.5) +
  labs(colour = "Soil Type") + 
  xlab("Mean Vertical Variation in Woody Vegetation") +
  ylab("Shannon H") + scale_color_hue(direction = -1)


# Boxplots of y var by transect
ggplot(XY_scale, aes(x = Cover0p5m_plot, y = shannonH,
                     group = Transect_f, colour = Transect_f)) +
  geom_boxplot(varwidth=FALSE, alpha=0.6, width=.15, position = "dodge2") +
  facet_wrap(~ Soil_f, scales = "free")


# # # Plot mod 2 variables by soil type
ggplot(XY_scale, aes(x = cv.maxpeakh., y = shannonH,
                     group = Soil_f, colour = Transect_f)) +
  geom_point() +
  facet_wrap(~ Soil_f, scales="free")

ggplot(XY_scale, aes(x = mean.cvH._woody, y = shannonH,
                     group = Soil_f, colour = Transect_f)) +
  geom_point() +
  facet_wrap(~ Soil_f)

ggplot(XY_scale, aes(x = sd.coverG., y = shannonH,
                     group = Soil_f, colour = Transect_f)) +
  geom_point() +
  facet_wrap(~ Soil_f, scales="free")

ggplot(XY_scale, aes(x = sd.maxH., y = shannonH,
                     group = Soil_f, colour = Transect_f)) +
  geom_point() +
  facet_wrap(~ Soil_f)

ggplot(XY_scale, aes(x = Cover0p5m_plot, y = shannonH,
                     group = Soil_f, colour = Transect_f)) +
  geom_point() +
  facet_wrap(~ Soil_f)

# Try adding these in for red soil model
ggplot(XY_scale, aes(x = sd.VDRpeak., y = shannonH,
                     group = Soil_f, colour = Transect_f)) +
  geom_point() +
  facet_wrap(~ Soil_f)

ggplot(XY_scale, aes(x = sd.herbh., y = shannonH,
                     group = Soil_f, colour = Transect_f)) +
  geom_point() +
  facet_wrap(~ Soil_f)

ggplot(XY_scale, aes(x = Cover0p5m_plot, y = shannonH,
                     group = Soil_f, colour = Transect_f)) +
  geom_point() +
  facet_wrap(~ Soil_f)

# try these for red soil 
ggplot(XY_scale, aes(x = cv.maxH., y = shannonH,
                     group = Soil_f, colour = Transect_f)) +
  geom_point() +
  facet_wrap(~ Soil_f)

ggplot(XY_scale, aes(x = cover_shrub_maxH, y = shannonH,
                     group = Soil_f, colour = Transect_f)) +
  geom_point() +
  facet_wrap(~ Soil_f, scales="free")

ggplot(XY_scale, aes(x = mean.ptoh., y = shannonH,
                     group = Soil_f, colour = Transect_f)) +
  geom_point() +
  facet_wrap(~ Soil_f, scales="free")


# Boxplots of y var by transect
ggplot(XY_scale, aes(x = sd.coverG., y = shannonH,
                     group = Transect_f, colour = Transect_f)) +
  geom_boxplot(varwidth=FALSE, alpha=0.6, width=.15, position = "dodge2") +
  facet_wrap(~ Soil_f, scales = "free")

# # # # 

# # Older models you tried
# 
# 
# mod1_lmer = lmer(shannonH ~ (0 + cv.maxpeakh.|Transect_f) +
#                    (0 + mean.cvH._woody|Transect_f) + (1|Transect_f), XY_scale)
# summary(mod1_lmer)
# tab_model(mod1_lmer)

# mod1 = gam(shannonH ~
#              s(sd.coverG.)+ 
#              s(Soil_f, bs="re") +
#              s(Transect_f, bs="re"),
#            data=XY_scale,
#            family=gaussian,
#            select=TRUE,
#            method="REML")
# 
# summary(mod1)
# 
# mod2 = gam(shannonH ~ 
#              s(sd.maxH.),
#            data=XY_scale,
#            family=gaussian,
#            select=TRUE,
#            method="REML")
# 
# summary(mod2)
# 
# mod2 = gam(shannonH ~ 
#              s(sd.maxH.) + s(sd.cvH._woody),
#            data=XY_scale,
#            family=gaussian,
#            select=TRUE,
#            method="REML")
# 
# summary(mod2)
# 
# 
# mod3 = gam(shannonH ~ 
#              s(sd.cvH._woody) +
#              s(sd.coverG.),
#            data=XY_scale,
#            family=gaussian,
#            select=TRUE,
#            method="REML")

# summary(mod3)












# # # 3) Soils and Random Effects!

# Add in Soil as a factor
colnames(XY_scale)

# Soil
mod2_soil = gam(shannonH ~ sd.coverG.:Soil_f + s(sd.coverG.),
                data=XY_scale,
                family=gaussian,
                select=TRUE,
                method="REML")

summary(mod2_soil)

# Beta coefficients are totally different 
# whether you use sd.coverG., herbh_plot, or sd.maxH.
# This pretty much confirms - we need diff. models for black and red soils!

# # # Black Soil Model
# shannonH MDI Top 6 X Vars:
#   
# sd.coverG.
# sd.nlayers.
# iqr.FHD.
# cv.cscore.
# iqr.coverG.
# max.ptoh.
# VDRpeak_plot
# cv.maxH.
# mean.cvH._woody
# max.cvpeakh.
# 
# shannonH Permutation Top 6 X Vars:
#   
# sd.coverG.
# mean.sdH._shrub
# propMultiLayer
# sd.cvH._shrub
# max.meanpeakh.
# cv.herbh.
# sd.maxH.
# mean.VDRpeak.
# iqr.cscore.
# cover_tree_maxH

# select(XY_scale, contains("cover")) %>% colnames()

XY_scale_black = XY_scale %>% filter(Soil == 'Black')

mod1_black = gam(shannonH ~ s(sd.coverG.),
                 data=XY_scale_black,
                 family=scat(link="identity"),
                 select=TRUE,
                 method="REML")

summary(mod1_black)
# Note sd.FHD. and iqr.FHD. also seem to be important


# Add a random effect for Transect
mod1_black_re = gam(shannonH ~ sd.coverG.:Transect_f,
                    data=XY_scale_black,
                    family=gaussian,
                    select=TRUE,
                    method="REML")

summary(mod1_black_re)

AIC(mod1_black, mod1_black_re)
# Standard model without re is better! interesting

ggplot(XY_scale_black, aes(x = sd.coverG.,
                           y = shannonH,
                           colour = Transect_f)) +
  geom_point()


ggplot(XY_scale_black, aes(x = sd.FHD.,
                           y = shannonH,
                           colour = Transect_f)) +
  geom_point()


# # # Red Soil Model

# Vars for red soil: 
# shannonH MDI Top 6 X Vars:
#   
#   sd.coverG.
# max.ptoh.
# cover_shrub_maxH
# cv.FHD.
# sd.sdH._shrub
# mean.sdH._woody
# sd.VDR.
# mean.sdH._grass
# ptoh_plot
# cv.VDR.
# 
# 
# shannonH Permutation Top 6 X Vars:
#   
# cv.FHD.
# cv.VDR.
# cv.meanpeakh.
# sd.coverG.
# iqr.coverG.
# cv.maxH.
# mean.sdH._woody
# cv.nlayers.
# sd.stdpeakh.
# cvpeakh_plot

XY_scale_red = XY_scale %>% filter(Soil == 'Red')

# Soil
mod1_red = gam(shannonH ~ cover_shrub_maxH,
               data=XY_scale_red,
               family=gaussian,
               select=TRUE,
               method="REML")

summary(mod1_red)

mod1_red_re = gam(shannonH ~ cover_shrub_maxH:Transect_f +
                    s(Transect_f, bs="re"),
                  data=XY_scale_red,
                  family=gaussian,
                  select=TRUE,
                  method="REML")

summary(mod1_red_re)


residuals(mod1_red, type='working')

ggplot(XY_scale_red, aes(x = sd.coverG.,
                         y = shannonH,
                         colour = Transect_f)) +
  geom_point()
