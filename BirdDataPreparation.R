# Data Preperation - Bird Model
# 11/11/2022 - PB

library(tidyverse)
library(vegan)

# # # 1) Prepare XY Data

# Datadir for location of in/out vars
datadir = '/n/home02/pbb/scripts/halo-metadata-server/Selenkay/data/'

figd = '/n/home02/pbb/scripts/SelenkayDiversity/figs/mango'

# Adding labels for different data processing runs
# 'initial'
# mango - 2/8/23
label = 'mango'
# note - features were added and renamed in this run!
# Mainly, the grasslayer metrics of cover and percheights were added
# and all feature names with '()' in them were replaced with '_'

# Radius of lidar data inquiry for modelling
radius = 30

# Open up X and Y vars

# Also, open up spotsppdf to make shannonH values with vegan
spotsppdf = read_csv(file='~/scripts/SelenkayDiversity/data/in/BirdSurvey/SelenkayBirds_SpotSppDF.csv')
# spotsppdf = read_csv(file=paste0(datadir, 'in/SelenkayBirds_SpotSppDF.csv'))
# spotsppdf = read_csv(file=paste0(datadir, 'out/SpotSppDF_Birds_final.csv'))


# NOTE: discontinued using this file SpotSummaryStats_Birds_final.csv - 11/14/22
# as the Y-vars appear to be wrong for some sites (1A-E have very high abundances) 
# these were made in python in google collab - redid 11/14,
# and doubled checked the X vars below
# X vars look OK
# That's why I'm remaking the y vars here, in vegan
# y = read_csv(file=paste0(datadir, 'in/BirdSurvey/y_vars/SpotSummaryStats_Birds_final.csv'))

# Note: a quick way to look for duplicates of certain taxa/spp:
# spotsppdf %>% select(contains('shrike')) %>% colnames()

# Remove the nan and index columns
spotsppdf = subset(spotsppdf, select = -c(nan, ...1))

# make a numeric version (without Spot column)
spotsppdf_numeric = subset(spotsppdf, select= -c(Spot))

# Added 2/9 - filter out all spp with only 1 observation
spotsppdf_numeric = spotsppdf_numeric[spotsppdf_numeric %>%
                                        colSums() > 1]

# make y vars - abundances, shannonH, 
shannonH = diversity(spotsppdf_numeric)
abun = spotsppdf_numeric %>% rowSums()
rich = spotsppdf_numeric %>% apply(1, function(x) sum(x>0))
# Evenness index (Pielou’s evenness from vegan paper)
evenness = shannonH/log(rich)
# inverse simpson
invsimp <- diversity(spotsppdf_numeric, "invsimp")

# make a y df with Spot added back in
y = data.frame(shannonH = shannonH) %>%
  mutate(Spot = spotsppdf$Spot,
         Abundance = abun,
         Richness = rich,
         Evenness = evenness,
         Simpson = invsimp)

head(y)

# Add Spot col back into the spotsppdf for saving
spotsppdf_curated = spotsppdf_numeric %>%
  mutate(Spot = spotsppdf$Spot)

# save csv files
write.csv(y, '~/scripts/SelenkayDiversity/data/in/BirdSurvey_y.csv')
write.csv(spotsppdf_curated,
          '~/scripts/SelenkayDiversity/data/in/SpotSppDF_curated.csv')
write.csv(spotsppdf_numeric,
          '~/scripts/SelenkayDiversity/data/in/SpotSppDFnumeric_curated.csv')


# # # Load and clean X vars
X_raw = read_csv(file=paste0(datadir, 'out/', label,
                             '/LidarMetrics_', radius,
                             'm_', label,'.csv'))

# Rename X_raw Site col to Spot
X_raw = X_raw %>% rename("Spot" = "Site")

# Also, change 3a to 3A for merging below
X_raw["Spot"][X_raw["Spot"] == "3a"] <- "3A"

# Also, load in the x,y centroid locations in EPSG32736
centroids = read_csv(file='~/scripts/SelenkayDiversity/data/in/SelenkaySpotPolygons_XYCentroids_EPSG32737.csv')

# Change 3a to 3A for merging below
centroids["Spot"][centroids["Spot"] == "3a"] <- "3A"

# # # Join X and y 

# NOTE: the raw X DF is slightly offset from the y
# need to join them on Spot ID
#c(y['Spot'], X_raw['Spot'])

# Merge on spot Col
XY = inner_join(y, X_raw, by='Spot')

# Also, merge with easting northing coordinate data
XY = inner_join(XY, centroids, by='Spot')


# # # Make a scaled XY dataframe

# Merge with centroids
X_raw = inner_join(X_raw, centroids, by='Spot')

# Scale the Numeric Variables
X_scale = X_raw %>%
  select_if(is.numeric) %>%
  scale() %>%
  data.frame() %>%
  bind_cols(select(X_raw, c('Spot', 'Soil', 'Treatment')))

# Merge on spot Col, again, with y data
XY_scale = inner_join(y, X_scale)

# Make Treatment and Soil into Factors
# Also, make nlayers_plot into a factor (since it's an integer)
# Also, add a Transect Column as a factor 
# for using gsub, see: 
# https://stackoverflow.com/questions/14543627/extracting-numbers-from-vectors-of-strings
XY <- XY %>%
  mutate(Treatment_f = factor(Treatment),
         Soil_f = factor(Soil),
         nlayers_plot_f = factor(nlayers_plot)) %>%
  mutate(Transect_f = factor(as.numeric(gsub("\\D", "", Spot))))

XY_scale <- XY_scale %>%
  mutate(Treatment_f = factor(Treatment),
         Soil_f = factor(Soil),
         nlayers_plot_f = factor(nlayers_plot)) %>%
  mutate(Transect_f = factor(as.numeric(gsub("\\D", "", Spot))))

# Write Both DFs to disk as rdf and csv files

# RDF files keep relevant variable types (factors, characters, etc.)
write_rds(XY, paste0('~/scripts/SelenkayDiversity/data/in/XY_', radius, 'mRadius.rds'))
write_rds(XY_scale, paste0('~/scripts/SelenkayDiversity/data/in/XY_scaled_', radius, 'mRadius.rds'))

# Csv files for reference and human readability
write.csv(XY, paste0('~/scripts/SelenkayDiversity/data/in/XY_scaled_', radius, 'mRadius.csv'))
write.csv(XY_scale, paste0('~/scripts/SelenkayDiversity/data/in/XY_scaled_', radius, 'mRadius.csv'))

# Quick way to check the classes of each col
sapply(XY, class)
sapply(XY_scale, class)

# Make some plots
ggplot(XY, aes(x = sd_cvH_vegtype_woody, y = shannonH,
                     group = Soil_f, colour = Treatment, shape=Treatment)) +
  geom_point() + 
  scale_colour_manual(values = c("#FFC20A", "#0C7BDC")) +
  facet_wrap(~ Soil_f)

ggplot(XY, aes(x = sd_CD_AboveG, y = shannonH,
                     group = Soil_f, colour = Treatment, shape=Treatment)) +
  geom_point() + 
  scale_colour_manual(values = c("#FFC20A", "#0C7BDC")) +
  facet_wrap(~ Soil_f)

ggplot(XY, aes(x = sd_CD_AboveG,
               y = Simpson,
               group = Soil_f,
               colour = Treatment, shape=Treatment)) +
  geom_point(size=3) + 
  scale_colour_manual(values = c("#FFC20A", "#0C7BDC")) +
  facet_wrap(~ Soil_f)



# Make some plots
ggplot(XY, aes(x = X, y = shannonH,
               group = Soil_f,
               colour = Treatment,
               shape=Treatment)) +
geom_point(size=3) + 
  scale_colour_manual(values = c("#FFC20A", "#0C7BDC")) +
  facet_wrap(~ Soil_f)

ggplot(XY, aes(x = Y,
               y = shannonH,
               group = Soil_f,
               colour = Treatment,
               shape=Treatment)) +
  geom_point(size=3) + 
  scale_colour_manual(values = c("#FFC20A", "#0C7BDC")) +
  facet_wrap(~ Soil_f)

ggplot(XY, aes(x = X,
                     y = Simpson,
                     group = Soil_f,
                     colour = Treatment,
                     shape=Treatment)) +
  geom_point(size=3) + 
  scale_colour_manual(values = c("#FFC20A", "#0C7BDC")) +
  facet_wrap(~ Soil_f)

ggplot(XY, aes(x = Y,
                     y = Simpson,
                     group = Soil_f,
                     colour = Treatment,
                     shape=Treatment)) +
  geom_point(size=3) + 
  scale_colour_manual(values = c("#FFC20A", "#0C7BDC")) +
  facet_wrap(~ Soil_f)


# # # Acoustic data stuff below - 2/15 commented out for now
# # # # ACOUSTIC data
# ACIbyHour = read_csv(file=paste0(datadir,
#                                  'in/Acoustics/SelenkaySurvey1_ACIbyHour.csv'))
# 
# # Note: you can split it by dawn and dusk
# # dawn in bird survey data was 6-9am (6-9 hours)
# # dusk was 4-7pm (so 16-19 hours)
# 
# # Get the mean ACI by site
# ACIbySite = ACIbyHour %>%
#             filter( ( (TIME >=16) & (TIME <=23) ) ) %>%
#             select(c('Site', 'ACI')) %>%
#             group_by(Site) %>%
#             summarise(ACI_mean = mean(ACI),
#                       ACI_max = max(ACI),
#                       ACI_med = median(ACI)) %>%
#             rename("Spot"="Site")
# 
# # Merge with y vars on spot Col
# yACI = merge(y, ACIbySite, by='Spot')
# 
# # Now, plot em! 
# ggplot(yACI, aes(x = ACI_mean, y = shannonH, 
#                  colour = Spot)) +
#   geom_point()
# 
# ggplot(yACI, aes(x = ACI_mean, y = Richness, 
#                  colour = Spot)) +
#   geom_point()
# 
# ggplot(yACI, aes(x = ACI_mean, y = Abundance, 
#                  colour = Spot)) +
#   geom_point()
# 
# 
# # Now, plot max
# ggplot(yACI, aes(x = ACI_max, y = shannonH, 
#                  colour = Spot)) +
#   geom_point()
# 
# ggplot(yACI, aes(x = ACI_max, y = Richness, 
#                  colour = Spot)) +
#   geom_point()
# 
# ggplot(yACI, aes(x = ACI_max, y = Abundance, 
#                  colour = Spot)) +
#   geom_point()
# 
# 
# # median
# ggplot(yACI, aes(x = ACI_med, y = shannonH, 
#                  colour = Spot)) +
#   geom_point()
# 
# ggplot(yACI, aes(x = ACI_med, y = Richness, 
#                  colour = Spot)) +
#   geom_point()
# 
# ggplot(yACI, aes(x = ACI_med, y = Abundance, 
#                  colour = Spot)) +
#   geom_point()
# 
# 
# # Merge acoustics with lidar data
# XY_acoustic = merge(yACI, X_scale)
# 
# ggplot(XY_acoustic, aes(x = Cover0p5m_plot, y = ACI_mean, 
#                  colour = Spot)) +
#   geom_point()
# 
# ggplot(XY_acoustic, aes(x = mean.cvH._woody, y = ACI_mean, 
#                         colour = Spot)) +
#   geom_point()
# 
# ggplot(XY_acoustic, aes(x = mean.cvH._shrub, y = ACI_mean, 
#                         colour = Spot)) +
#   geom_point()
# 
# ggplot(XY_acoustic, aes(x = mean.cvH._grass, y = ACI_mean, 
#                         colour = Spot)) +
#   geom_point()


