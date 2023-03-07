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

# # # Plots of Important variables

# Shannon H
p1 = ggplot(XY, aes(x = Cover0p5m_plot, y = shannonH,
                          group = Soil_f, colour = Soil_f)) +
  geom_point(size=2.5) +
  labs(colour = "Soil Type") + 
  xlab("Canopy Cover above 0.5m") +
  ylab("Shannon H") + scale_color_hue(direction = -1)

p1

ggsave(plot=p1,
       filename="/n/home02/pbb/scripts/SelenkayDiversity/figs/Scatter_shannonh_cover0p5m.png",
       device="png",
       height=4, width=6)

# Shannon H
p2 = ggplot(XY, aes(x = `mean(cvH)_woody`, y = shannonH,
                    group = Soil_f, colour = Soil_f)) +
  geom_point(size=2.5) +
  labs(colour = "Soil Type") + 
  xlab("Mean Vertical Variation in Woody Vegetation") +
  ylab("Shannon H") + scale_color_hue(direction = -1)

p2

ggsave(plot=p2,
       filename="/n/home02/pbb/scripts/SelenkayDiversity/figs/Scatter_shannonh_meancvHwoody.png",
       device="png",
       height=4, width=6)


# # # Shannon H Model

# Add another var
mod_cvHwoody_1.0 = gam(shannonH ~ s(mean.cvH._woody) +
                         s(mean.cvH._woody, by=Soil_f) +
                         Transect_f + Soil_f,
                       data=XY_scale,
                       family=gaussian,
                       select=TRUE,
                       method="REML")

# Ok! Actually, looks quite good
draw(mod_cvHwoody_1.0, residuals = TRUE)
summary(mod_cvHwoody_1.0)
appraise(mod_cvHwoody_1.0)
gam.check(mod_cvHwoody_1.0, rep=500)
lines(c(3.0, 3.8), c(3.0, 3.8))


# Add another var
mod_cvHwoody_2.0 = gam(shannonH ~ s(mean.cvH._woody) + 
                         s(Cover0p5m_plot, by=Soil_f) +
                         s(mean.cvH._woody) +
                         Transect_f + Soil_f,
                       data=XY_scale,
                       family=gaussian,
                       select=TRUE,
                       method="REML")

# Ok! Actually, looks quite good
draw(mod_cvHwoody_2.0, residuals = TRUE)
summary(mod_cvHwoody_2.0)
appraise(mod_cvHwoody_2.0)
gam.check(mod_cvHwoody_2.0, rep=500)
lines(c(3.0, 3.8), c(3.0, 3.8))

## # Calc RMSE
# sqrt(mean((data$actual - data$predicted)^2))

pred = predict(mod_cvHwoody_2.0, data=XY_scale)
RMSE = sqrt(mean((pred - XY_scale$shannonH)^2))
rRMSE = RMSE/mean(XY_scale$shannonH)

# Now that you've selected your variables, 
# Split a training and test dataset

#make this example reproducible
set.seed(1)

# Use 70% of dataset as training set and remaining 30% as testing set
sample <- sample(c(TRUE, FALSE),
                 nrow(XY_scale),
                 replace=TRUE,
                 prob=c(0.7,0.3))

# From: https://www.statology.org/train-test-split-r/
train  <- XY_scale[sample, ]
test   <- XY_scale[!sample, ]

nrow(train)
nrow(test)

# Make a training model
# NOTE: had to leave out the transect term
# should be a more accurate view of the prediction error in any case
mod_train = gam(shannonH ~ s(mean.cvH._woody) + 
                  s(Cover0p5m_plot, by=Soil_f) +
                  s(mean.cvH._woody) + Soil_f,
                data=train,
                family=gaussian,
                select=TRUE,
                method="REML")

# Ok! Actually, looks quite good
draw(mod_train, residuals = TRUE)
summary(mod_train)
appraise(mod_train)
gam.check(mod_train, rep=500)
lines(c(3.0, 3.8), c(3.0, 3.8))

# Use test data to get RMSE
pred = predict(mod_train, newdata=test)
RMSE = sqrt(mean((pred - test$shannonH)^2))
rRMSE = RMSE/mean(test$shannonH)


RMSE_list = c()

# Run that like 1 million times if you want to get an idea of the dist.
for (x in 1:50) {
  
  #make this example reproducible
  set.seed(x)
  
  # Use 80% of dataset as training set and remaining 20% as testing set
  sample <- sample(c(TRUE, FALSE),
                   nrow(XY_scale),
                   replace=TRUE,
                   prob=c(0.8,0.2))
  
  # From: https://www.statology.org/train-test-split-r/
  train  <- XY_scale[sample, ]
  test   <- XY_scale[!sample, ]
  
  mod_train = gam(shannonH ~ s(mean.cvH._woody) + 
                    s(Cover0p5m_plot, by=Soil_f) +
                    s(mean.cvH._woody) + Soil_f,
                  data=train,
                  family=gaussian,
                  select=TRUE,
                  method="REML")
  
  # Use test data to get RMSE
  pred = predict(mod_train, newdata=test)
  RMSE = sqrt(mean((pred - test$shannonH)^2))
  rRMSE = RMSE/mean(test$shannonH)
  
  # Append
  RMSE_list = c(RMSE_list, RMSE)

}

# Print out values 
mean(RMSE_list)
max(RMSE_list)
min(RMSE_list)
sd(RMSE_list)

# # # PLOTS

# # # Plot shannon H inside and outside real quick
p1 = ggplot(data=XY_scale, mapping = aes(x = Treatment_f,
                                         y = shannonH,
                                         group = Treatment_f,
                                         colour = Soil_f,
)) + 
  geom_boxplot() +
  facet_wrap(~Soil_f) +
  labs(colour = "Soil Type") + 
  xlab("") +
  ylab("Shannon H") + scale_color_hue(direction = -1)

p1
ggsave(plot=p1,
       filename="/n/home02/pbb/scripts/SelenkayDiversity/figs/Boxplot_ShannonH_InOut_RedBlack.png",
       device="png",
       height=6, width=4)


ggplot(data=XY_scale, mapping = aes(x = Treatment_f,
                                    y = shannonH,
                                    group = Treatment_f)) + 
  geom_boxplot() +
  xlab("") +
  ylab("Shannon H") + scale_color_hue(direction = -1)



p = ggplot(data=XY_scale, mapping = aes(x = Treatment_f,
                                        y = Richness,
                                        group = Treatment_f,
                                        colour = Soil_f)) + 
  geom_boxplot() +
  xlab("") +
  labs(colour = "Soil Type") + 
  facet_wrap(~Soil_f) +
  ylab("Number of Bird Species") + scale_color_hue(direction = -1)

p

ggsave(plot=p,
       filename="/n/home02/pbb/scripts/SelenkayDiversity/figs/Boxplot_Richness_InOut_RedBlack.png",
       device="png",
       height=6, width=4)

p = ggplot(data=XY_scale, mapping = aes(x = Treatment_f,
                                        y = Abundance,
                                        group = Treatment_f,
                                        colour = Soil_f)) + 
  geom_boxplot() +
  xlab("") +
  labs(colour = "Soil Type") + 
  facet_wrap(~Soil_f) +
  ylab("Number of Birds Observed") + scale_color_hue(direction = -1)

p

ggsave(plot=p,
       filename="/n/home02/pbb/scripts/SelenkayDiversity/figs/Boxplot_Abundance_InOut_RedBlack.png",
       device="png",
       height=6, width=4)


p = ggplot(data=XY_scale, mapping = aes(x = Treatment_f,
                                        y = Cover0p5m_plot,
                                        group = Treatment_f,
                                        colour = Soil_f)) + 
  geom_boxplot() +
  xlab("") +
  labs(colour = "Soil Type") + 
  facet_wrap(~Soil_f) +
  ylab("Canopy Cover above 0.5 m") + scale_color_hue(direction = -1)

p

ggsave(plot=p,
       filename="/n/home02/pbb/scripts/SelenkayDiversity/figs/Boxplot_Cover0p5m_InOut_RedBlack.png",
       device="png",
       height=6, width=4)


# # # # 

p1 = ggplot(XY_scale, aes(x = Cover0p5m_plot, y = shannonH,
                          group = Treatment_f, colour = Treatment_f)) +
  geom_point(size=2.5) +
  labs(colour = "Inside/Outside Selenkay") + 
  xlab("Canopy Cover above 0.5m") +
  ylab("Shannon H") + scale_color_hue(direction = -1)

p1


p1 = ggplot(XY_scale, aes(x = mean.cvH._woody, y = shannonH,
                          group = Treatment_f, colour = Treatment_f)) +
  geom_point(size=2.5) +
  labs(colour = "Inside/Outside Selenkay") + 
  xlab("Mean Vertical Variation in Woody Vegetation") +
  ylab("Shannon H") + scale_color_hue(direction = -1)

p1
