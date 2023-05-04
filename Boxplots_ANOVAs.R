# Boxplots and ANOVAs for paper
# 2/22/2023 - PB

library(tidyverse)
library(vegan)
library(car)
library(broom)
# broom is useful for tidying up tables for output
# https://stackoverflow.com/questions/61106658/one-way-anova-for-loop-how-do-i-iterate-through-multiple-columns-of-a-dataframe
# https://stackoverflow.com/questions/67014880/2-way-anova-loop-export-to-csv-or-excel-table

# dev.off()

datadir = '/n/home02/pbb/scripts/SelenkayDiversity/data/'

figd = '/n/home02/pbb/scripts/SelenkayDiversity/figs/mango'
tabled = '/n/home02/pbb/scripts/SelenkayDiversity/tables/mango/ANOVA'

# toggle radius of inquiry
radius = 130 

# RDF files keep relevant variable types (factors, characters, etc.)
XY = readRDS(paste0(datadir, 'in/XY_', radius,'mRadius.rds'))
XY_scale = readRDS(paste0(datadir, 'in/XY_scaled_', radius,'mRadius.rds'))

# # #  Boxplots of y vars
# and Important metrics


ggplot(XY, aes(x = Treatment,
                     y = Evenness,
                     fill=Treatment)) +
  geom_boxplot() + 
  facet_wrap(~ Soil_f) +
  scale_fill_manual(values = c("#FFC20A", "#0C7BDC")) +  
  theme(axis.text.y = element_text(colour = "black", size = 11), 
        axis.text.x = element_text(colour = "black", size = 11),
        strip.text.x = element_text(colour = c("black", "red"),
                                    face = "bold", 
                                    size = 11),
        legend.text = element_text(size = 12, face = "bold", 
                                   colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold",
                                                               size = 12), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"))

ggsave(paste0(figd, "/boxplots/Boxplots_Evenness_SplitbySoil.png"),
       width = 6, height = 4, device='png', dpi=300)


ggplot(XY, aes(x = Treatment,
                     y = Simpson,
                     fill=Treatment)) +
  geom_boxplot() + 
  facet_wrap(~ Soil_f) +
  scale_fill_manual(values = c("#FFC20A", "#0C7BDC")) +  
  theme(axis.text.y = element_text(colour = "black", size = 11), 
        axis.text.x = element_text(colour = "black", size = 11),
        strip.text.x = element_text(colour = c("black", "red"),
                                    face = "bold", 
                                    size = 11),
        legend.text = element_text(size = 12, face = "bold", 
                                   colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold",
                                                               size = 12), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"))

ggsave(paste0(figd, "/boxplots/Boxplots_SimpsonInv_SplitbySoil.png"),
       width = 6, height = 4, device='png', dpi=300)


ggplot(XY, aes(x = Treatment,
                     y = Abundance,
                     fill=Treatment)) +
  geom_boxplot() + 
  facet_wrap(~ Soil_f) +
  scale_fill_manual(values = c("#FFC20A", "#0C7BDC")) +  
  theme(axis.text.y = element_text(colour = "black", size = 11), 
        axis.text.x = element_text(colour = "black", size = 11),
        strip.text.x = element_text(colour = c("black", "red"),
                                    face = "bold", 
                                    size = 11),
        legend.text = element_text(size = 12, face = "bold", 
                                   colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold",
                                                               size = 12), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"))

ggsave(paste0(figd, "/boxplots/Boxplots_Abundance_SplitbySoil.png"),
       width = 6, height = 4, device='png', dpi=300)

#

ggplot(XY, aes(x = Treatment,
                     y = Richness,
                     fill=Treatment)) +
  geom_boxplot() + 
  facet_wrap(~ Soil_f) +
  scale_fill_manual(values = c("#FFC20A", "#0C7BDC")) +  
  theme(axis.text.y = element_text(colour = "black", size = 11), 
        axis.text.x = element_text(colour = "black", size = 11),
        strip.text.x = element_text(colour = c("black", "red"),
                                    face = "bold", 
                                    size = 11),
        legend.text = element_text(size = 12, face = "bold", 
                                   colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold",
                                                               size = 12), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"))

ggsave(paste0(figd, "/boxplots/Boxplots_Richness_SplitbySoil.png"),
       width = 6, height = 4, device='png', dpi=300)

#

ggplot(XY, aes(x = Treatment,
                     y = shannonH,
                     fill=Treatment)) +
  geom_boxplot() + 
  facet_wrap(~ Soil_f) +
  scale_fill_manual(values = c("#FFC20A", "#0C7BDC")) +  
  theme(axis.text.y = element_text(colour = "black", size = 11), 
        axis.text.x = element_text(colour = "black", size = 11),
        strip.text.x = element_text(colour = c("black", "red"),
                                    face = "bold", 
                                    size = 11),
        legend.text = element_text(size = 12, face = "bold", 
                                   colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold",
                                                               size = 12), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"))

ggsave(paste0(figd, "/boxplots/Boxplots_Shannon_SplitbySoil.png"),
       width = 6, height = 4, device='png', dpi=300)

# # # ADDED PLOTS NOT SPLIT BY SOIL

ggplot(XY, aes(x = Treatment,
               y = Richness,
               fill=Treatment)) +
  geom_boxplot() + 
  scale_fill_manual(values = c("#FFC20A", "#0C7BDC")) +  
  theme(axis.text.y = element_text(colour = "black", size = 11), 
        axis.text.x = element_text(colour = "black", size = 11),
        strip.text.x = element_text(colour = c("black", "red"),
                                    face = "bold", 
                                    size = 11),
        legend.text = element_text(size = 12, face = "bold", 
                                   colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold",
                                                               size = 12), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"))

ggsave(paste0(figd, "/boxplots/Boxplots_Richness_AllSoils.png"),
       width = 6, height = 4, device='png', dpi=300)

# Quick anova
formula = as.formula(paste0('Richness', '~ Treatment_f'))
aov_res <- aov(formula,
               data = XY)
aov_tbl = broom::tidy(aov_res)
aov_tbl

#

ggplot(XY, aes(x = Treatment,
               y = shannonH,
               fill=Treatment)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FFC20A", "#0C7BDC")) +  
  theme(axis.text.y = element_text(colour = "black", size = 11), 
        axis.text.x = element_text(colour = "black", size = 11),
        strip.text.x = element_text(colour = c("black", "red"),
                                    face = "bold", 
                                    size = 11),
        legend.text = element_text(size = 12, face = "bold", 
                                   colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold",
                                                               size = 12), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"))

ggsave(paste0(figd, "/boxplots/Boxplots_Shannon_AllSoils.png"),
       width = 6, height = 4, device='png', dpi=300)


# Quick anova
formula = as.formula(paste0('shannonH', '~ Treatment_f'))
aov_res <- aov(formula,
               data = XY)
aov_tbl = broom::tidy(aov_res)
aov_tbl

# 
ggplot(XY, aes(x = Treatment,
               y = Abundance,
               fill=Treatment)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FFC20A", "#0C7BDC")) +  
  theme(axis.text.y = element_text(colour = "black", size = 11), 
        axis.text.x = element_text(colour = "black", size = 11),
        strip.text.x = element_text(colour = c("black", "red"),
                                    face = "bold", 
                                    size = 11),
        legend.text = element_text(size = 12, face = "bold", 
                                   colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold",
                                                               size = 12), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"))

ggsave(paste0(figd, "/boxplots/Boxplots_Abundance_AllSoils.png"),
       width = 6, height = 4, device='png', dpi=300)


# Quick anova
formula = as.formula(paste0('Abundance', '~ Treatment_f'))
aov_res <- aov(formula,
               data = XY)
aov_tbl = broom::tidy(aov_res)
aov_tbl
# # # Now, also make boxplots of some environmental vars

# These are the top correlated variables with shannonH
# Plus a few simple ones (maxH, sdH, etc.)
# columns = c('sd_CD_AboveG',
#             'sd_sdH_vegtype_grass',
#             'cv_CD_AboveG',
#             'mean_CD_AboveG',
#             'sd_CD_G',
#             'mean_nlayers',
#             'sd_nlayers',
#             'mean_cvH',
#             'horzcover_grass',
#             'Cover0p5m_plot',
#             'mean_CD_AboveG',
#             'cv_gapsize',
#             'mean_maxH',
#             'sd_maxH')
columns = c('Abundance', 
            'Richness',
            'Evenness',
            'shannonH',
            'Simpson',
            'mean_gapsize',
            'sd_gapsize',
            'cv_gapsize',
            'gapsize_plot',
            'mean_maxH',
            'sd_maxH',
            'cv_maxH',
            'mean_CD_AboveG',
            'sd_CD_AboveG',
            'Cover0p05m_plot',
            'Cover1p5m_plot',
            'mean_nlayers',
            'nlayers_plot',
            'sd_nlayers',
            'cv_nlayers',
            'horzcover_grass',
            'horzcover_woody',
            'X', 'Y')

for (c in columns) {
  
  p1 = ggplot(XY, aes_string(x = 'Treatment',
                             y = c,
                             fill='Treatment')) +
    geom_boxplot() + 
    facet_wrap(~ Soil_f) +
    scale_fill_manual(values = c("#FFC20A", "#0C7BDC")) +  
    theme(axis.text.y = element_text(colour = "black", size = 11), 
          axis.text.x = element_text(colour = "black", size = 11),
          strip.text.x = element_text(colour = c("black", "red"),
                                      face = "bold", 
                                      size = 11),
          legend.text = element_text(size = 12, face = "bold", 
                                     colour ="black"), 
          legend.position = "right", axis.title.y = element_text(face = "bold",
                                                                 size = 12), 
          axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
          legend.title = element_text(size = 14, colour = "black", face = "bold"))
  
  ggsave(paste0(figd, "/boxplots/Boxplots_", c, "_", radius,"mRadius.png"),
         p1,
         width = 6, height = 4,
         device='png', dpi=300)
  
}


p1

# # # ANOVAs


# Loop for multiple vars and saving outputs
columns = c('Abundance', 
            'Richness',
            'Evenness',
            'shannonH',
            'Simpson',
            'mean_gapsize',
            'sd_gapsize',
            'cv_gapsize',
            'gapsize_plot',
            'mean_maxH',
            'sd_maxH',
            'cv_maxH',
            'mean_CD_AboveG',
            'sd_CD_AboveG',
            'Cover0p05m_plot',
            'Cover1p5m_plot',
            'mean_nlayers',
            'nlayers_plot',
            'sd_nlayers',
            'cv_nlayers',
            'horzcover_grass',
            'horzcover_woody')

#  Run an  ANOVA test
# with the variables nested
# https://www.statology.org/nested-anova-in-r/

for (c in columns) {
  
  #https://stackoverflow.com/questions/26889240/looping-many-one-sided-anova-in-r
  formula = as.formula(paste0(c, '~ Soil_f / Treatment_f'))
  
  # Check for equal varainces
  # https://statsandr.com/blog/anova-in-r/
  lt = leveneTest(formula,
                  data = XY)
  
  lt_tbl = broom::tidy(lt)
  
  write_csv(lt_tbl,
            paste0(tabled, '/Nested/LevenesTest_', c, '_', radius, 'mRadius.csv'))
  
  # # Compute the analysis of variance
  aov_res <- aov(formula,
                 data = XY)
  
  aov_tbl = broom::tidy(aov_res)
  
  write_csv(aov_tbl,
            paste0(tabled, '/Nested/Anova_', c, '_', radius, 'mRadius.csv'))
  
  # Tukey's Honest Sig diff test
  thsd = TukeyHSD(aov_res)
  thsd_tbl = broom::tidy(thsd)
  write_csv(thsd_tbl,
            paste0(tabled, '/Nested/TukeysHSD_', c, '_', radius, 'mRadius.csv'))
  
  
}

# Unnested anova 
# (using SoilTreatment & 4 groups )
# Add a SoilTreatment group for Anova testing
XY = XY %>% mutate(SoilTreatment = factor(paste(Soil, Treatment, sep="")))

for (c in columns) {
  
  #https://stackoverflow.com/questions/26889240/looping-many-one-sided-anova-in-r
  formula = as.formula(paste0(c, '~ SoilTreatment'))
  
  # Check for equal varainces
  # https://statsandr.com/blog/anova-in-r/
  lt = leveneTest(formula,
                  data = XY)
  
  lt_tbl = broom::tidy(lt)
  
  write_csv(lt_tbl,
            paste0(tabled, '/SoilTreatment/LevenesTest_', c, '_', radius, 'mRadius.csv'))
  
  # # Compute the analysis of variance
  aov_res <- aov(formula,
                 data = XY)
  
  aov_tbl = broom::tidy(aov_res)
  
  write_csv(aov_tbl,
            paste0(tabled, '/SoilTreatment/Anova_', c, '_', radius, 'mRadius.csv'))
  
  # Tukey's Honest Sig diff test
  thsd = TukeyHSD(aov_res)
  thsd_tbl = broom::tidy(thsd)
  write_csv(thsd_tbl,
            paste0(tabled, '/SoilTreatment/TukeysHSD_', c, '_', radius, 'mRadius.csv'))
  
}

# # # TESTING

# 
# # Basic anovaworkflow for 1 var
# NOT nested approach 
# # Add a SoilTreatment group for Anova testing
# XY = XY %>% mutate(SoilTreatment = factor(paste(Soil, Treatment, sep="")))
# 
# # Check for equal varainces
# # https://statsandr.com/blog/anova-in-r/
# lt = leveneTest(Evenness ~ SoilTreatment,
#            data = XY)
# lt
# 
# # # Compute the analysis of variance
# aov_even <- aov(Evenness ~ SoilTreatment,
#                  data = XY)
# 
# summary(aov_even)
# 
# # can use broom here to make a nice table
# aov_even_tbl = broom::tidy(aov_even)
# aov_even_tbl
# write_csv(aov_even_tbl,
#           paste0(tabled, '/Anova_', 'Evenness', '.csv'))
# 
# # Tukey's Honest Sig diff test
# thsd = TukeyHSD(aov_even)
# thsd
# thsd_tbl = broom::tidy(thsd)
# thsd_tbl


# # Split into red and black dfs for anovas
# XY_black = XY %>%
#   filter(Soil=='Black')
# 
# XY_red = XY %>%
#   filter(Soil=='Red')
# 
# # Check for normality
# # Shapiro-Wilk or Kolmogorov-Smirnov test.
# 
# # Check for equal varainces
# # https://statsandr.com/blog/anova-in-r/
# lt = leveneTest(Evenness ~ Treatment_f,
#            data = XY_black)
# lt
# # if lt < 0.05, then groups don't have equal variances
# # and you should use welch's ttest or K-W anova
# # otherwise, good to go
# 
# # Compute the analysis of variance
# aov_black <- aov(Evenness ~ Treatment_f,
#                data = XY_black)
# 
# summary(aov_black)
# 
# # Kruskal wallis anova
# kruskal.test(variable ~ group, data = dat)