# Looking for indicator spp in black and red soils, inside and out
# PB 
# 02/16/2023

library(vegan)
library(indicspecies)
library(tidyverse)

# # #  Load Data and set Directories- 2/15 update


figd = '/n/home02/pbb/scripts/SelenkayDiversity/figs/mango'

# with spot col
df_spp_spot = read.csv('~/scripts/SelenkayDiversity/data/in/SpotSppDF_curated.csv')

# only numeric
df_spp = read.csv('~/scripts/SelenkayDiversity/data/in/SpotSppDFnumeric_curated.csv')

# drop the X (row names) column
df_spp_spot = select(df_spp_spot, -c(X))
df_spp = select(df_spp, -c(X))

# # # Run indicator spp analysis

# # Get env variables (Spot, Soil, Treatment)

# Load in XY df
datadir = '/n/home02/pbb/scripts/SelenkayDiversity/data/'
XY = readRDS(paste0(datadir, 'in/XY.rds'))

# Add site, treatment and spot into df_site_spp to align them correctly
df_spp_env = df_spp_spot %>%
  inner_join(select(XY, c('Soil', 'Treatment', 'Spot')),
             by='Spot')

# Add a new df that combines multiple treatments
# so "black inside", "red outside", etc.
df_spp_env = df_spp_env %>% mutate(SoilTreatment = factor(paste(Soil, Treatment, sep="")))

# make a proportions df
df_spp_p = decostand(df_spp, method="total")

# use multipatt function to identify different spp. inside/outside
# https://esajournals.onlinelibrary.wiley.com/doi/10.1890/08-1823.1

inv = multipatt(df_spp_p,
                cluster = df_spp_env$SoilTreatment,
                control = how(nperm=9999))

summary(inv, indvalcomp = TRUE,  alpha=0.05)

# # # Also do a Rank Abundance Curve

# Aggregate by soil treatment
df_SoilTreatment = aggregate(df_spp,
                             list(df_spp_env$SoilTreatment),
                             FUN=sum)

RAcurve = radfit(df_SoilTreatment %>% select(-c(Group.1)),
                 model='Null')

plot(RAcurve)

# make a better plot
library(wesanderson)

pal = wes_palette(n=4, name="Moonrise2")
  
p = ggplot() + geom_line(aes(x=order(RAcurve$`1`$y, decreasing=TRUE),
                           y=RAcurve$`1`$y,
                           colour='BlackInside')) + 
    geom_line(aes(x=order(RAcurve$`2`$y, decreasing=TRUE),
                    y=RAcurve$`2`$y,
                    colour='BlackOutside')) + 
    geom_line(aes(x=order(RAcurve$`3`$y, decreasing=TRUE),
                  y=RAcurve$`3`$y,
                  colour='RedInside')) + 
    geom_line(aes(x=order(RAcurve$`4`$y, decreasing=TRUE),
                  y=RAcurve$`4`$y,
                  colour='RedOutside')) + 
    scale_y_continuous(trans='log10') +
    annotation_logticks(sides="l") +
    labs(title="Rank Abundance Curves",
         x ="Species Rank", y = "Abundance",
         colour='') +
    scale_colour_manual(values=c('#000033',
      '#0000CC',
      '#FF0033',
      '#993300'
      ))

p

ggsave(filename=paste0(figd, "/RankAbundanceCurves.png"),
         plot=p,
         width = 6,
         height = 5,
         device='png',
         dpi=300)
  

# # # Scatter Plot of Shannon and spp richness

# make a numeric version (without Spot column)
spotsppdf_numeric = subset(df_spp_spot, select= -c(Spot))

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

# # # Scatter plots of Richness, shannon, even, and abundance
ggplot() + geom_point(aes(x=abun,
                          y=evenness,
                          colour=df_spp_env$SoilTreatment)) +
  scale_colour_manual(values=c('#000033',
                               '#0000CC',
                               '#FF0033',
                               '#993300'
  ))

ggplot() + geom_boxplot(aes(x=abun,
                            colour=df_spp_env$Soil)) +
  scale_colour_manual(values=c('#000000',
                               '#FF0033'
  )) + 
  facet_grid(~df_spp_env$Treatment)


