# NDMS Analysis of Bird Diversity Data
# Run on the 10 Transects
# PB 
# initial run: 09/13/2022
# revised: 02/07/2023

# TBD - add the acoustic metrics in here, might be cool

# Most of this code is lifted straight from these tutorials:
# Running NDMS: https://rpubs.com/CPEL/NMDS
# Plotting: https://jkzorz.github.io/2019/06/06/NMDS.html
# Could do a follow-up indicator spp analysis following:
# https://jkzorz.github.io/2019/07/02/Indicator-species-analysis.html
# https://cran.r-project.org/web/packages/vegan/vignettes/intro-vegan.pdf

# Other useful refs:
# https://jkzorz.github.io/2020/04/04/NMDS-extras.html
# https://dplyr.tidyverse.org/reference/mutate.html
# https://github.com/rstudio/cheatsheets/blob/main/data-transformation.pdf

library(vegan)
library(arrow)
library(tidyverse)
# install.packages(plot_ly)
# library(plot_ly)

# https://cran.r-project.org/web/packages/vegan/vegan.pdf
# https://rdrr.io/cran/vegan/f/inst/doc/intro-vegan.pdf
datadir = '/n/home02/pbb/scripts/SelenkayDiversity/data/'

# Full bird dataframe
# df = read_parquet(paste0(datadir, "in/BirdSurveyData_All.pqt"))

# Summary Stats of Bird Data from python, google collab (Rows = Site)
df_birds = read.csv(paste0(datadir, "in/BirdSurvey/DFBirds_SummaryStats.csv"))
                       
# Bird Counts by Site (Rows) and Species (Cols)
df_site_spp = read.csv(paste0(datadir, "in/BirdSurvey/SelenkayBirds_SiteSppDF.csv"))

# Bird Counts by Site (Rows) and Species (Cols) Proportions from python (just double checked they're the same as the decostand func)
# df_p_py = read.csv(paste0(datadir, "in/SelenkayBirds_SiteSppProportionsDF.csv"))

# Drop the X (transect) column
df_spp = df_site_spp %>%
  mutate(X=NULL)

# make a proportions df
df_spp_p = decostand(df_spp, method="total")

# make a distance matrix
df_spp_distmat = vegdist(df_spp_p, method="bray")

# Creating easy to view matrix and writing .csv
df_spp_distmat <- 
  as.matrix(df_spp_distmat, labels = T)

write.csv(df_spp_distmat, paste0(datadir, "/out/bird_spp_distmat.csv"))

# Running NMDS in vegan (metaMDS)
df_spp_NMS <-
  metaMDS(df_spp_distmat,
          distance = "bray",
          k = 3,
          maxit = 999, 
          trymax = 500,
          wascores = TRUE)

# Shepards test/goodness of fit
goodness(df_spp_NMS) # Produces a results of test statistics for goodness of fit for each point

stressplot(df_spp_NMS) # Produces a Shepards diagram

# Plotting points in ordination space
plot(df_spp_NMS, "sites")   # Produces distance 
orditorp(df_spp_NMS, "sites")   # Gives points labels

# Instead, try a ggplot with labels

#extract NMDS scores (x and y coordinates)
df_scores = as.data.frame(scores(df_spp_NMS, display="site"))

# Add group cols to df
df_scores = df_scores %>%
  mutate(Treatment = df_birds$Treatment) %>%
  mutate(Soil = df_birds$Soil) %>%
  mutate(TransectNum = df_birds$X + 1)

# # # PLOT
# Note, you can switch between the axes (NMDS1, 2, and 3) below
xx = ggplot(df_scores, aes(x = NMDS1, y = NMDS2, label= TransectNum)) + 
  geom_point(size = 4, aes(shape = Treatment, colour = Soil)) + 
  geom_text(nudge_x = 0.015,nudge_y = 0.015,) +
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "Soil", y = "NMDS2", shape = "Treatment")  + 
  scale_colour_manual(values = c("#000000", "#D55E00")) 

xx

ggsave(paste0("/n/home02/pbb/scripts/SelenkayDiversity/figs/SelenkayBirdDiv_NMDS_Dims1v2_k=3Dims.png"))


# Now try in ggplot
# Note, you can switch between the axes (NMDS1, 2, and 3) below
xx = ggplot(df_scores, aes(x = NMDS1, y = NMDS3, label= TransectNum)) + 
  geom_point(size = 4, aes(shape = Treatment, colour = Soil)) + 
  geom_text(nudge_x = 0.015,nudge_y = 0.015,) +
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "Soil", y = "NMDS3", shape = "Treatment")  + 
  scale_colour_manual(values = c("#000000", "#D55E00")) 

xx

ggsave(paste0("/n/home02/pbb/scripts/SelenkayDiversity/figs/SelenkayBirdDiv_NMDS_Dims1v3_k=3Dims.png"))


# Now try in ggplot
# Note, you can switch between the axes (NMDS1, 2, and 3) below
xx = ggplot(df_scores, aes(x = NMDS2, y = NMDS3, label= TransectNum)) + 
  geom_point(size = 4, aes(shape = Treatment, colour = Soil)) + 
  geom_text(nudge_x = 0.015,nudge_y = 0.015,) +
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS2", colour = "Soil", y = "NMDS3", shape = "Treatment")  + 
  scale_colour_manual(values = c("#000000", "#D55E00")) 

xx

ggsave(paste0("/n/home02/pbb/scripts/SelenkayDiversity/figs/SelenkayBirdDiv_NMDS_Dims2v3_k=3Dims.png"))


## Cool!
# 2/7/23 - After a re-run with cleaned bird data,
# The major axes of variation (1 and 2) appear to be along soil type gradient
# The treatment effect (inside vs out) is definately there, but it's secondary.
# This probably means you need 2 seperate models. 
# Need to run an NMDS on all 50 sites, that might help. 

