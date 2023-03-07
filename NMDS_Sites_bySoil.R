# Initial NDMS Analysis of Bird Diversity Data
# Run on the 50 sites (5 per transect)
# And split by soil type (so, seperate NMDS for black and red soils)
# PB 
# 02/07/2023

# Nice article on picking colorblind friendly colors:
# https://davidmathlogic.com/colorblind/#%23FFC20A-%230C7BDC
# I picked orange and blue: "#FFC20A", "#0C7BDC"

# Nice article on NMS and stress values
# http://strata.uga.edu/8370/lecturenotes/multidimensionalScaling.html

library(vegan)
library(arrow)
library(tidyverse)

# dev.off()

datadir = '/n/home02/pbb/scripts/SelenkayDiversity/data/'

figd = '/n/home02/pbb/scripts/SelenkayDiversity/figs/mango'

# Load Data - 2/15 update
# with spot col
df_spp_spot = read.csv('~/scripts/SelenkayDiversity/data/in/SpotSppDF_curated.csv')
# only numeric
df_spp = read.csv('~/scripts/SelenkayDiversity/data/in/SpotSppDFnumeric_curated.csv')

# drop the X (row names) column
df_spp_spot = select(df_spp_spot, -c(X))
df_spp = select(df_spp, -c(X))

# # #  Discontinued 2/15 - now using inputs from BirdDataPreperation.R
# # Summary Stats of Bird Data from python, google collab (Rows = Site)
# df_birds = read.csv(paste0(datadir, "in/BirdSurvey/DFBirds_SummaryStatsbySpot.csv"))
# 
# # Bird Counts by Site (Rows) and Species (Cols)
# df_spot_spp = read.csv(paste0(datadir, "in/BirdSurvey/SelenkayBirds_SpotSppDF.csv"))
# 
# # Added 2/9 - filter out all spp with only 1 observation
# df_spot_spp = df_spot_spp[df_spot_spp %>% select(-c('Spot')) %>% colSums() > 1]
# # # 

# Read XY df created in BirdDataPreperation.R
XY = readRDS(paste0(datadir, 'in/XY.rds'))


# Add environmental variabels to the spot spp df
# by merging on the spot column
df_spot_spp = df_spp_spot %>%
  inner_join(select(XY, c('Spot', 'Soil', 'Treatment',
                     "Cover0p25m_plot",
                     'mean_cvH_vegtype_woody','mean_cvH_vegtype_grass',
                     "sd_sdH_vegtype_woody", "sd_sdH_vegtype_grass",
                     "sd_maxH", 'mean_maxH',
                     "cv_CD_AboveG",  "mean_cvH_vegtype_woody")))

# # # BLACK SOILS
df_spot_spp_black = df_spot_spp %>%
  filter(Soil=='Black')

# Drop the X (transect), spot label, and other xtra columns
df_spp_black = df_spot_spp_black %>%
  mutate(Spot=NULL, Treatment=NULL, Soil=NULL,
         "Cover0p25m_plot"=NULL,
         'mean_cvH_vegtype_woody'=NULL,'mean_cvH_vegtype_grass'=NULL,
         "sd_sdH_vegtype_woody"=NULL, "sd_sdH_vegtype_grass"=NULL,
         "sd_maxH"=NULL, 'mean_maxH'=NULL,
         "cv_CD_AboveG"=NULL,  "mean_cvH_vegtype_woody"=NULL)

# take out any spp with no observations in them
df_spp_black = df_spp_black[df_spp_black %>%
                              colSums() > 0]

# # # calculate basic summary stats

# total number of observations 
# (after filtering for observations with just 1 value)
df_spp_black %>% colSums() %>% data.frame() %>% colSums()

# shape of df after removing spp with 0 observations
# (2nd number gives total number of species)
df_spp_black %>%
  size_sum()

# make a proportions df
df_spp_p_black = decostand(df_spp_black, method="total")

# make a distance matrix
df_spp_distmat_black = vegdist(df_spp_p_black, method="bray")

# Creating easy to view matrix and writing .csv
df_spp_distmat_black <- 
  as.matrix(df_spp_distmat_black, labels = T)

write.csv(df_spp_distmat_black,
          paste0(datadir, "/out/bird_spp_distmat_bySpot_blacksoil.csv"))

# Running NMDS in vegan (metaMDS)
df_spp_NMS_black <-
  metaMDS(df_spp_distmat_black,
          distance = "bray",
          k = 3,
          maxit = 999, 
          trymax = 500,
          wascores = TRUE)

# Shepards test/goodness of fit
goodness(df_spp_NMS_black) # Produces a results of test statistics for goodness of fit for each point

stressplot(df_spp_NMS_black) # Produces a Shepards diagram

# stree value to report
df_spp_NMS_black$stress

# Plotting points in ordination space
plot(df_spp_NMS_black, "sites")   # Produces distance 
orditorp(df_spp_NMS_black, "sites")   # Gives points labels

# Instead, try a ggplot with labels

#extract NMDS scores (x and y coordinates)
df_scores_black = as.data.frame(scores(df_spp_NMS_black, display="site"))


# # # Add in Env Vars

# Add a spot col to the df_scores df and
# Merge the dfs on the spot column
df_scores_black = df_scores_black %>%
  mutate(Spot=df_spot_spp_black$Spot) %>%
  merge(select(df_spot_spp_black, c('Spot', 'Soil', 'Treatment')))


# # # PLOT
# Note, you can switch between the axes (NMDS1, 2, and 3) below
xx = ggplot(df_scores_black, aes(x = NMDS1, y = NMDS2, label= Spot)) + 
  geom_point(size = 4, aes(shape = Treatment, colour = Treatment)) + 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 12, colour = "black"), 
        title = element_text(face = "bold", size = 12, colour = "black"),
        legend.title = element_text(size = 12, colour = "black", face = "bold"), 
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", y = "NMDS2", shape = "Treatment") + 
  scale_colour_manual(values = c("#FFC20A", "#0C7BDC")) + 
  ggtitle('Black Soil') + coord_fixed(ratio = 1) + 
  stat_ellipse(data=df_scores_black,
               level = 0.95,
               aes(colour=Treatment))

# NOTE: took labels out 2/22
# geom_text(nudge_x = 0.015,nudge_y = 0.015) +

xx

ggsave(paste0(figd, "/NMDS/NMDS_Dims1v2_3Dims_bySpot_BlackSoils_nolabels.png"),
       xx,
       dpi=300,
       width = 6, height = 4,
       units = "in")



# # # PLOT
# Note, you can switch between the axes (NMDS1, 2, and 3) below
xx = ggplot(df_scores_black, aes(x = NMDS2, y = NMDS3, label= Spot)) + 
  geom_point(size = 4, aes(shape = Treatment, colour = Treatment)) + 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 12, colour = "black"), 
        title = element_text(face = "bold", size = 12, colour = "black"),
        legend.title = element_text(size = 12, colour = "black", face = "bold"), 
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS2", y = "NMDS3", shape = "Treatment") + 
  scale_colour_manual(values = c("#FFC20A", "#0C7BDC")) + 
  ggtitle('Black Soil') + coord_fixed(ratio = 1) + 
  stat_ellipse(data=df_scores_black,
               level = 0.95,
               aes(colour=Treatment))

# NOTE: took labels out 2/22
# geom_text(nudge_x = 0.015,nudge_y = 0.015) +

xx

ggsave(paste0(figd, "/NMDS/NMDS_Dims2v3_3Dims_bySpot_BlackSoils_nolabels.png"),
       xx,
       dpi=300,
       width = 6, height = 4,
       units = "in")


# # # PLOT
# Note, you can switch between the axes (NMDS1, 2, and 3) below
xx = ggplot(df_scores_black, aes(x = NMDS1, y = NMDS3, label= Spot)) + 
  geom_point(size = 4, aes(shape = Treatment, colour = Treatment)) + 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 12, colour = "black"), 
        title = element_text(face = "bold", size = 12, colour = "black"),
        legend.title = element_text(size = 12, colour = "black", face = "bold"), 
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", y = "NMDS3", shape = "Treatment") + 
  scale_colour_manual(values = c("#FFC20A", "#0C7BDC")) + 
  ggtitle('Black Soil') + coord_fixed(ratio = 1) + 
  stat_ellipse(data=df_scores_black,
               level = 0.95,
               aes(colour=Treatment))

# NOTE: took labels out 2/22
# geom_text(nudge_x = 0.015,nudge_y = 0.015) +

xx

ggsave(paste0(figd, "/NMDS/NMDS_Dims1v3_3Dims_bySpot_BlackSoils_nolabels.png"),
       xx,
       dpi=300,
       width = 6, height = 4,
       units = "in")

# Add in lidar Variables of interest with envfit
# https://jkzorz.github.io/2020/04/04/NMDS-extras.html
env = df_spot_spp %>%
  filter(Soil=='Black') %>%
  select("Cover0p25m_plot",
         'mean_cvH_vegtype_woody','mean_cvH_vegtype_grass',
         "sd_sdH_vegtype_woody", "sd_sdH_vegtype_grass",
         "sd_maxH", 'mean_maxH',
         "cv_CD_AboveG",  "mean_cvH_vegtype_woody")

# , 'Treatment'
en = envfit(df_spp_NMS_black, env, permutations = 999, na.rm = TRUE)
# make coordinates for plotting
en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en)

xx = ggplot() + 
  geom_point(data=df_scores_black,
             size = 4,
             mapping=aes(x = NMDS1, y = NMDS2,
                         shape = Treatment, colour = Treatment)) + 
  geom_text(data = df_scores_black,
            aes(x = NMDS1, y = NMDS2, label = Spot),
            colour = "black", 
            fontface = "bold",
            check_overlap = TRUE,
            size=2.2,
            nudge_x=0.02,
            nudge_y=0.015) +
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 12, colour = "black"), 
        title = element_text(face = "bold", size = 12, colour = "black"),
        legend.title = element_text(size = 12, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", y = "NMDS2", shape = "Treatment") + 
  scale_colour_manual(values = c("#FFC20A", "#0C7BDC")) + 
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = en_coord_cont, size =1, alpha = 0.5, colour = "grey30") + 
  geom_text(data = en_coord_cont,
            aes(x = NMDS1, y = NMDS2),
            colour = "black", 
            fontface = "bold",
            label = row.names(en_coord_cont),
            check_overlap = TRUE,
            size=2.2,
            nudge_y = 0.01,
            nudge_x = -0.03) +
  xlim(-0.4, 0.3) + 
  ylim(-0.25, 0.25) + 
  ggtitle('Black Soil') + coord_fixed(ratio = 1)

xx

ggsave(paste0(figd, "/NMDS/NMDSwithLidar_Dims1v2_2Dims_bySpot_BlackSoils.png"),
       xx,
       width = 6, height = 4, units = "in", device='png', dpi=300)



# # # RED SOILS
df_spot_spp_red = df_spot_spp %>%
  filter(Soil=='Red')

# Drop the X (transect) and spot label column
df_spp_red = df_spot_spp_red %>%
  mutate(Spot=NULL, Treatment=NULL, Soil=NULL,
         "Cover0p25m_plot"=NULL,
         'mean_cvH_vegtype_woody'=NULL,'mean_cvH_vegtype_grass'=NULL,
         "sd_sdH_vegtype_woody"=NULL, "sd_sdH_vegtype_grass"=NULL,
         "sd_maxH"=NULL, 'mean_maxH'=NULL,
         "cv_CD_AboveG"=NULL,  "mean_cvH_vegtype_woody"=NULL)

# take out any spp with no observations in them
df_spp_red = df_spp_red[df_spp_red %>%
                              colSums() > 0]

# # # calculate basic summary stats

# total number of observations 
# (after filtering for observations with just 1 value)
df_spp_red %>% colSums() %>% data.frame() %>% colSums()

# shape of df after removing spp with 0 observations
# (2nd number gives total number of species)
df_spp_red %>%
  size_sum()

# make a proportions df
df_spp_p_red = decostand(df_spp_red, method="total")

# make a distance matrix
df_spp_distmat_red = vegdist(df_spp_p_red, method="bray")

# Creating easy to view matrix and writing .csv
df_spp_distmat_red <- 
  as.matrix(df_spp_distmat_red, labels = T)

write.csv(df_spp_distmat_red, paste0(datadir, "/out/bird_spp_distmat_bySpot_redsoil.csv"))

# Running NMDS in vegan (metaMDS)
df_spp_NMS_red <-
  metaMDS(df_spp_distmat_red,
          distance = "bray",
          k = 3,
          maxit = 999, 
          trymax = 500,
          wascores = TRUE)

# Shepards test/goodness of fit
goodness(df_spp_NMS_red) # Produces a results of test statistics for goodness of fit for each point

stressplot(df_spp_NMS_red) # Produces a Shepards diagram

# Stress metric
df_spp_NMS_red$stress

# Instead, try a ggplot with labels

#extract NMDS scores (x and y coordinates)
df_scores_red = as.data.frame(scores(df_spp_NMS_red, display="site"))

# Add a spot col to the df_scores df and
# Merge the dfs on the spot column
df_scores_red = df_scores_red %>%
  mutate(Spot=df_spot_spp_red$Spot) %>%
  inner_join(select(df_spot_spp_red, c('Spot', 'Soil', 'Treatment')))

# # # PLOT
# Note, you can switch between the axes (NMDS1, 2, and 3) below
xx = ggplot(df_scores_red, aes(x = NMDS1, y = NMDS2, label= Spot)) + 
  geom_point(size = 4, aes(shape = Treatment, colour = Treatment)) + 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 12, colour = "black"), 
        title = element_text(face = "bold", size = 12, colour = "black"),
        legend.title = element_text(size = 12, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", y = "NMDS2", shape = "Treatment") + 
  scale_colour_manual(values = c("#FFC20A", "#0C7BDC")) + 
  ggtitle('Red Soil') + coord_fixed(ratio = 1)  + 
  stat_ellipse(data=df_scores_red,
               level = 0.95,
               aes(colour=Treatment))

# Removed text labels
# geom_text(nudge_x = 0.015,nudge_y = 0.015) +

xx

ggsave(paste0(figd, "/NMDS/NMDS_Dims1v2_3Dims_bySpot_RedSoils_nolabels.png"),
       xx,
       width = 6, height = 4, units = "in", device='png', dpi=300)


# # # PLOT
# Note, you can switch between the axes (NMDS1, 2, and 3) below
xx = ggplot(df_scores_red, aes(x = NMDS2, y = NMDS3, label= Spot)) + 
  geom_point(size = 4, aes(shape = Treatment, colour = Treatment)) + 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 12, colour = "black"), 
        title = element_text(face = "bold", size = 12, colour = "black"),
        legend.title = element_text(size = 12, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS2", y = "NMDS3", shape = "Treatment") + 
  scale_colour_manual(values = c("#FFC20A", "#0C7BDC")) + 
  ggtitle('Red Soil') + coord_fixed(ratio = 1)  + 
  stat_ellipse(data=df_scores_red,
               level = 0.95,
               aes(colour=Treatment))

# Removed text labels
# geom_text(nudge_x = 0.015,nudge_y = 0.015) +

xx

ggsave(paste0(figd, "/NMDS/NMDS_Dims2v3_3Dims_bySpot_RedSoils_nolabels.png"),
       xx,
       width = 6, height = 4, units = "in", device='png', dpi=300)


# # # PLOT
# Note, you can switch between the axes (NMDS1, 2, and 3) below
xx = ggplot(df_scores_red, aes(x = NMDS1, y = NMDS3, label= Spot)) + 
  geom_point(size = 4, aes(shape = Treatment, colour = Treatment)) + 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 12, colour = "black"), 
        title = element_text(face = "bold", size = 12, colour = "black"),
        legend.title = element_text(size = 12, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", y = "NMDS3", shape = "Treatment") + 
  scale_colour_manual(values = c("#FFC20A", "#0C7BDC")) + 
  ggtitle('Red Soil') + coord_fixed(ratio = 1)  + 
  stat_ellipse(data=df_scores_red,
               level = 0.95,
               aes(colour=Treatment))

# Removed text labels
# geom_text(nudge_x = 0.015,nudge_y = 0.015) +

xx

ggsave(paste0(figd, "/NMDS/NMDS_Dims1v3_3Dims_bySpot_RedSoils_nolabels.png"),
       xx,
       width = 6, height = 4, units = "in", device='png', dpi=300)

# # # Add in Env Vars


# Add in lidar Variables of interest with envfit
# https://jkzorz.github.io/2020/04/04/NMDS-extras.html
env_red = df_spot_spp %>%
  filter(Soil=='Red') %>%
  select("Cover0p25m_plot",
         'mean_cvH_vegtype_woody','mean_cvH_vegtype_grass',
         "sd_sdH_vegtype_woody", "sd_sdH_vegtype_grass",
         "sd_maxH", 'mean_maxH',
         "cv_CD_AboveG",  "mean_cvH_vegtype_woody")

# , 'Treatment'
en_red = envfit(df_spp_NMS_red,
                env_red,
                permutations = 999,
                na.rm = TRUE)

# make coordinates for plotting
# en_coord_cont_red = as.data.frame(scores(en_red, "vectors")) * ordiArrowMul(en_red)
# 
# xx = ggplot() + 
#   geom_point(data=df_scores_red,
#              size = 4,
#              mapping=aes(x = NMDS1,
#                          y = NMDS2,
#                          shape = Treatment,
#                          colour = Treatment)) + 
#   geom_text(data = df_scores_red,
#             aes(x = NMDS1,
#                 y = NMDS2,
#                 label = Spot),
#             colour = "black", 
#             fontface = "bold",
#             check_overlap = TRUE,
#             size=2.2,
#             nudge_x=0.02,
#             nudge_y=0.015) +
#   theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
#         axis.text.x = element_text(colour = "black", face = "bold", size = 11), 
#         legend.text = element_text(size = 11, face ="bold", colour ="black"), 
#         legend.position = "right", axis.title.y = element_text(face = "bold", size = 12), 
#         axis.title.x = element_text(face = "bold", size = 12, colour = "black"), 
#         legend.title = element_text(size = 12, colour = "black", face = "bold"), 
#         panel.background = element_blank(),
#         panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
#         legend.key=element_blank()) + 
#   labs(x = "NMDS1", y = "NMDS2", shape = "Treatment") + 
#   scale_colour_manual(values = c("#FFC20A", "#0C7BDC")) + 
#   geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
#                data = en_coord_cont, size =1, alpha = 0.5, colour = "grey30") + 
#   geom_text(data = en_coord_cont_red,
#             aes(x = NMDS1, y = NMDS2),
#             colour = "black", 
#             fontface = "bold",
#             label = row.names(en_coord_cont_red),
#             check_overlap = TRUE,
#             size=2.2,
#             nudge_y = 0.01,
#             nudge_x = -0.03) +
#   xlim(-0.4, 0.3) + 
#   ylim(-0.25, 0.25) + 
#   ggtitle('Red Soil') + coord_fixed(ratio = 1)
# 
# xx
# 
# ggsave(paste0(figd, "/NMDS/NMDSwithLidar_Dims1v2_2Dims_bySpot_RedSoils.png"),
#        width = 6, height = 6, device='png', dpi=300)


# # # Plot some explanatory vars

# Make a red soil Df
XY_red = XY %>% filter(Soil=="Red")

ggplot(data = XY_red) + geom_point(aes(y=shannonH,
                                       x=horzcover_grass,
                                       shape = Treatment,
                                       colour = Treatment))

ggplot() + geom_point(aes(y=XY_red$'shannonH',
                          x=XY_red$'sd_cvH_vegtype_woody',
                          shape = XY_red$'Treatment',
                          colour = XY_red$'Treatment'))

ggplot() + geom_point(aes(y=XY_red$'shannonH',
                          x=XY_red$'mean_cvH_vegtype_woody',
                          shape = XY_red$'Treatment',
                          colour = XY_red$'Treatment'))

# # BLACK

# Make a black soil Df
XY_black = XY %>% filter(Soil=="Black")

ggplot() + geom_point(aes(y=XY_black$'shannonH',
                          x=XY_black$'sd_herbh',
                          shape = XY_black$'Treatment',
                          colour = XY_black$'Treatment'))


ggplot() + geom_point(aes(y=XY_black$'shannonH',
                          x=XY_black$'mean_cvH_vegtype_woody',
                          shape = XY_black$'Treatment',
                          colour = XY_black$'Treatment'))

# # # 
# # # 
# # # 
# OUTSIDE and INSIDE 
# Basic summary stats
# added 2/28/23

# # # BLACK SOILS
df_spot_spp_in = df_spot_spp %>%
  filter(Treatment=='Inside')

# Drop the X (transect), spot label, and other xtra columns
df_spp_in = df_spot_spp_in %>%
  mutate(Spot=NULL, Treatment=NULL, Soil=NULL,
         "Cover0p25m_plot"=NULL,
         'mean_cvH_vegtype_woody'=NULL,'mean_cvH_vegtype_grass'=NULL,
         "sd_sdH_vegtype_woody"=NULL, "sd_sdH_vegtype_grass"=NULL,
         "sd_maxH"=NULL, 'mean_maxH'=NULL,
         "cv_CD_AboveG"=NULL,  "mean_cvH_vegtype_woody"=NULL)

# take out any spp with no observations in them
df_spp_in = df_spp_in[df_spp_in %>%
                              colSums() > 0]

# # # calculate basic summary stats

# total number of observations 
# (after filtering for observations with just 1 value)
df_spp_in %>% colSums() %>% data.frame() %>% colSums()

# shape of df after removing spp with 0 observations
# (2nd number gives total number of species)
df_spp_in %>%
  size_sum()

