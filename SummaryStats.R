# Summary Stats of Sites 
# PB 2/28/23

# Datadir for location of in/out vars
datadir = '/n/home02/pbb/scripts/SelenkayDiversity/data'

figd = '/n/home02/pbb/scripts/SelenkayDiversity/figs'

tabled = '/n/home02/pbb/scripts/SelenkayDiversity/tables'

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
XY = readRDS(paste0(datadir, '/in/XY_', radius,'mRadius.rds'))

# filter to only y vars, and soil, treatment, spot, etc.
y = XY %>% select( c('Abundance', 'Richness', 'shannonH', 'Evenness',
                     'Soil', 'Treatment'))

# make groupbys
y_soil = y %>% group_by('Soil') %>% summarise(avg=mean(Abundance, na.rm=TRUE))
y_treat = y %>% group_by('Treatment')
# y_soil_treat = y %>% group_by(c('Soil', 'Treatment'))



y_soil %>% apply(mean('Abundance'))

# summarize mean, std, etc. for different soil types and treatments

# output a table


# # # Some quick summary Stats
# 
# # total number of observations 
# # (after filtering for observations with just 1 value)
# df_spp %>% colSums() %>% data.frame() %>% colSums()
# 
# # shape of df (2nd number gives total number of species)
# df_spp %>% size_sum()
# 
# # # # Groupby summary stats
# 
# # Read XY df created in BirdDataPreperation.R
# XY = readRDS(paste0(datadir, 'in/XY.rds'))
# 
# # Merge the dfs on the spot column
# df_spp_spot = df_spp_spot %>%
#   inner_join(select(XY, c('Spot', 'Soil', 'Treatment')))
# 
# # # # RED SOILS
# df_spp_red = df_spp_spot %>%
#   filter(Soil=='Red')
# 
# # Group by soil, then other stuff
# df_spp_red = df_spp_red %>%
#   select( -c('Spot', 'Soil', 'Treatment'))
# 
# # take out any spp with no observations in them
# df_spp_red = df_spp_red[df_spp_red %>%
#                           colSums() > 0]
# 
# 
# # (after filtering for observations with just 1 value)
# df_spp_red %>%
#   colSums() %>%
#   data.frame() %>%
#   colSums()
# 
# # mean and sd of abundances 
# mean = df_spp_red %>% rowSums() %>% var()
# sd = df_spp_red %>% rowSums() %>% sd()

