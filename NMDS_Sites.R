#  NDMS Analysis of Bird Diversity Data
# Run on the 50 sites (5 per transect)
# PB 
# 02/07/2023

# Vegan Cheatsheet
# https://raw.githubusercontent.com/rstudio/cheatsheets/main/vegan.pdf

# Vegan paper
# https://mirror.linux.duke.edu/cran/web/packages/vegan/vignettes/diversity-vegan.pdf
# https://cran.ism.ac.jp/web/packages/vegan/vegan.pdf

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
# Vegan
# https://cran.r-project.org/web/packages/vegan/vegan.pdf
# https://rdrr.io/cran/vegan/f/inst/doc/intro-vegan.pdf

library(vegan)
library(arrow)
library(tidyverse)

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

# # # Some quick summary Stats

# total number of observations 
# (after filtering for observations with just 1 value)
df_spp %>% colSums() %>% data.frame() %>% colSums()

# shape of df (2nd number gives total number of species)
df_spp %>% size_sum()


# # # 2/15- discontinued doing spot spp df processing here
# everything is done in BirdDataPreperation.R

# # Summary Stats of Bird Data from python, google collab (Rows = Site)
# df_birds = read.csv(paste0(datadir, "in/BirdSurvey/DFBirds_SummaryStatsbySpot.csv"))
# 
# # Bird Counts by Site (Rows) and Species (Cols)
# df_site_spp = read.csv(paste0(datadir, "in/BirdSurvey/SelenkayBirds_SpotSppDF.csv"))
# 
# # Added 2/9 - filter out all spp with only 1 observation
# df_site_spp = df_site_spp[df_site_spp %>% select(-c('Spot')) %>% colSums() > 1]
# 
# # Drop the X (transect) and spot label column
# # Also, remove this pesky 'nan' column
# df_spp = df_site_spp %>%
#   mutate(X=NULL, Spot=NULL) %>% 
#   subset(select=-c(nan, ...1))
# # #

# make a proportions df
df_spp_p = decostand(df_spp, method="total")

# make a distance matrix
df_spp_distmat = vegdist(df_spp_p, method="bray")

# Creating easy to view matrix and writing .csv
df_spp_distmat <- 
  as.matrix(df_spp_distmat,
            labels = T)

write.csv(df_spp_distmat,
          paste0(datadir, "/out/bird_spp_distmat_bySpot.csv"))

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

# Stress metric
df_spp_NMS$stress

# Plotting points in ordination space
plot(df_spp_NMS, "sites")   # Produces distance 
orditorp(df_spp_NMS, "sites")   # Gives points labels

# Instead, try a ggplot with labels

#extract NMDS scores (x and y coordinates)
df_scores = as.data.frame(scores(df_spp_NMS, display="site"))


# # # Add in Env Vars

# Read XY df created in BirdDataPreperation.R
XY = readRDS(paste0(datadir, 'in/XY.rds'))

# Add a spot col to the df_scores df and
# Merge the dfs on the spot column
df_scores = df_scores %>%
  mutate(Spot=df_spp_spot$Spot) %>%
  inner_join(select(XY, c('Spot', 'Soil', 'Treatment')))


# # # PLOT
# Note, you can switch between the axes (NMDS1, 2, and 3) below
xx = ggplot(df_scores, aes(x = NMDS1, y = NMDS2, label= Spot)) + 
  geom_point(size = 4, aes(shape = Treatment, colour = Soil)) +
  labs(x = "NMDS1", colour = "Soil", y = "NMDS2", shape = "Treatment")  + 
  scale_colour_manual(values = c("grey30", "#D55E00")) +
  coord_fixed(ratio = 1) + 
  stat_ellipse(data=df_scores,
               level = 0.95,
               aes(colour=Soil)) + 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 12, colour = "black"), 
        title = element_text(face = "bold", size = 12, colour = "black"),
        legend.title = element_text(size = 12, colour = "black", face = "bold"), 
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) 
# took out labels
# geom_text(nudge_x = 0.02,nudge_y = 0.015) +
  # + 
  #   theme(axis.text.y = element_text(colour = "black", size = 11, face = "bold"), 
  #         axis.text.x = element_text(colour = "black", face = "bold", size = 11), 
  #         legend.text = element_text(size = 12, face ="bold", colour ="black"), 
  #         legend.position = "right", axis.title.y = element_text(face = "bold", size = 12), 
  #         axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
  #         legend.title = element_text(size = 14, colour = "black", face = "bold"), 
  #         panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
  #         legend.key=element_blank()) + 
xx

ggsave(filename=paste0(figd, "/NMDS/SelenkayBirdDiv_NMDS_Dims1v2_k=3Dims_bySpot_nolabels.png"),
       plot=xx,
       width = 8,
       height = 6,
       device='png',
       dpi=300)


# # # # EXTRA plots if running a k=3 axis NMDS
xx = ggplot(df_scores, aes(x = NMDS2, y = NMDS3, label= Spot)) + 
  geom_point(size = 4, aes(shape = Treatment, colour = Soil)) + 
  theme(axis.text.y = element_text(colour = "black", size = 11, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 11), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 12), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS2", colour = "Soil", y = "NMDS3", shape = "Treatment")  + 
  scale_colour_manual(values = c("grey30", "#D55E00")) +
  coord_fixed(ratio = 1) + 
  stat_ellipse(data=df_scores,
               level = 0.95,
               aes(colour=Soil))

xx

ggsave(filename=paste0(figd, "/NMDS/SelenkayBirdDiv_NMDS_Dims2v3_k=3Dims_bySpot_nolabels.png"),
       plot=xx,
       width = 8,
       height = 6,
       device='png',
       dpi=300)


# # # # EXTRA plots if running a k=3 axis NMDS
xx = ggplot(df_scores, aes(x = NMDS1, y = NMDS3, label= Spot)) + 
  geom_point(size = 4, aes(shape = Treatment, colour = Soil)) + 
  theme(axis.text.y = element_text(colour = "black", size = 11, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 11), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 12), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "Soil", y = "NMDS3", shape = "Treatment")  + 
  scale_colour_manual(values = c("grey30", "#D55E00")) +
  coord_fixed(ratio = 1) + 
  stat_ellipse(data=df_scores,
               level = 0.95,
               aes(colour=Soil))

xx

ggsave(filename=paste0(figd, "/NMDS/SelenkayBirdDiv_NMDS_Dims1v3_k=3Dims_bySpot_nolabels.png"),
       plot=xx,
       width = 8,
       height = 6,
       device='png',
       dpi=300)

# Interesting - at the spot level, we don't see much variation by Treatment
# mainly by soil
# Try splitting by soil type (seperate files)


# # # Run anosim to find differences between soil types and inside/outside

# Add site, treatment and spot into df_site_spp to align them correctly
df_site_spp_env = df_spp_spot %>%
  inner_join(select(XY, c('Soil', 'Treatment', 'Spot')),
             by='Spot')



library(wesanderson)
# if you want a wes anderson color palette
# https://github.com/karthik/wesanderson
pal <- wes_palette("Zissou1",
                   300,
                   type = "continuous")


# Make and save a spp heat map plot

# Way to save R plots: 
# http://www.sthda.com/english/wiki/creating-and-saving-graphs-r-base-graphs
# 1. Open png file
png(paste0(figd, "/NMDS/SpotSppHeatMap.png"),
    height = 14, width = 16, units = 'in', res=300)

# Note: Used inner_join above to keep the row order from the first df
# merge will use the order of the 2nd df 
# which would throw things off, since spot 10 comes before
# plot community data table
# https://search.r-project.org/CRAN/refmans/vegan/html/vegemite.html
tabasco(select(df_site_spp_env, -c('Soil', 'Treatment', 'Spot')),
              labCol = df_site_spp_env$Spot,
              col=terrain.colors(n=50))

# 3. Close the file
dev.off()


# Note - have to save the above plot manually...
# ggsave(filename=paste0(figd, "/NMDS/SpotSppHeatMap.png"),
#        width = 12,
#        height = 8,
#        device='png',
#        dpi=300)

# not sure how to interpret anosim... 
# https://sites.google.com/site/mb3gustame/hypothesis-tests/anosim
# https://www.researchgate.net/post/Can_anyone_help_me_in_understanding_and_clearly_interpreting_ANOSIM_Analysis_of_Similarityand_SIMPER_Similarity_percentage_analysisresults
# I guess, an R value close to 1 means they're more similar
# close to 0, they are dissimilar
# the pvalue tells you that it's significant... 
# https://rdrr.io/rforge/vegan/man/anosim.html
df_spp_distmat = vegdist(df_spp_p, method="bray")
df_anosim_soil = anosim(x=df_spp_distmat,
                        grouping=df_site_spp_env$Soil)

summary(df_anosim_soil)
df_anosim_soil
plot(df_anosim_soil)


df_anosim_treatment = anosim(x=df_spp_distmat,
                        grouping=df_site_spp_env$Treatment) 
df_anosim_treatment
summary(df_anosim_treatment)
plot(df_anosim_treatment)

# # # Aggregating for beta diversity
# sum abundances of each spp by soil and treatment

df_agg_treatment = aggregate(df_spp,
                             list(df_site_spp_env$Treatment),
                             FUN=sum)

df_agg_soil = aggregate(df_spp,
                        list(df_site_spp_env$Soil),
                        FUN=sum)

# # # Compute beta diversity
# https://rdrr.io/rforge/vegan/man/betadiver.html
# wtf is it? idk
# https://methodsblog.com/2015/05/27/beta_diversity/
beta_soil = betadiver(df_agg_soil)
beta_treatment = betadiver(df_agg_treatment)

# NOTE: maybe you should actually be using betadispers...
# https://rdrr.io/cran/vegan/man/betadisper.html

betadisp_treatment = betadisper(df_spp_distmat,
                                df_site_spp_env$Treatment)
summary(betadisp_treatment)
betadisp_treatment
plot(betadisp_treatment)
anova(betadisp_treatment)
tukeyHSD_treatment = TukeyHSD(betadisp_treatment)
plot(tukeyHSD_treatment)


betadisp_soil = betadisper(df_spp_distmat,
                                df_site_spp_env$Soil)
summary(betadisp_soil)
betadisp_soil
plot(betadisp_soil)
anova(betadisp_soil)
tukeyHSD_soil = TukeyHSD(betadisp_soil)
plot(tukeyHSD_soil)

# # # Model Beta diversities

df_beta_dist = betadiver(df_spp, "w")

quantile(df_beta_dist)


# Example from the vegan paper - Oksansen 2022
betamod_soil = with(select(df_site_spp_env, c(Soil, Treatment)),
           betadisper(df_beta_dist, Soil))

betamod_soil
boxplot(betamod_soil)
summary(betamod_soil)
