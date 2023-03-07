# Final 130model

# library("devtools"); install_github("lme4/lme4",dependencies=TRUE)
library(devtools)
# devtools::install_github("sjPlot/devel")
library(tidyverse)
library(mgcv)
library(arrow)
library(lme4)
library(gratia)
library(gridExtra)

# # # 1) Load XY Data
# Note: XY data prepared in "BirdDataPreperation.R"

# Datadir for location of in/out vars
datadir = '/n/home02/pbb/scripts/SelenkayDiversity/data/'

figd = '/n/home02/pbb/scripts/SelenkayDiversity/figs/mango'

tabled = '/n/home02/pbb/scripts/SelenkayDiversity/tables/mango/GAMs'

# toggle radius of inquiry
radius = 130

# RDF files keep relevant variable types (factors, characters, etc.)
XY = readRDS(paste0(datadir, 'in/XY_', radius,'mRadius.rds'))
XY_scale = readRDS(paste0(datadir, 'in/XY_scaled_', radius,'mRadius.rds'))



mod130.5.0 = gam(shannonH ~
                   s(sd_gapsize, by=Soil_f) +
                   s(VDRpeak_plot, by=Soil_f) +
                   s(X,Y, k=5),
                 data=XY_scale,
                 family=gaussian,
                 select=TRUE,
                 method="REML")

summary(mod130.5.0)
p130.5.0 = draw(mod130.5.0, residuals = TRUE)
p130.5.0


ggsave(plot = p130.5.0,
       filename=paste0(figd,
                       "/GAMs/MixedGAM-130-5-0_PartialEffects_", radius ,"mRadius_mango.png"),
       width = 14, height = 8, units = "in", device='png', dpi=300)

a130.5.0 = appraise(mod130.5.0)
a130.5.0

ggsave(plot = a130.5.0,
       filename=paste0(figd,
                       "/GAMs/MixedGAM-130-5-0_Appraise_", radius ,"mRadius_mango.png"),
       width = 12, height = 8, units = "in", device='png', dpi=300)

gam.check(mod130.5.0, rep=500)
lines(c(3.0, 3.8), c(3.0, 3.8))

# compute RMSE
mod130.5.0_RMSE = sqrt(mean(mod130.5.0$residuals**2))

# Same some important results
vartbl_mod130.5.0 = tidy(mod130.5.0)
write_csv(vartbl_mod130.5.0,
          paste0(tabled, '/VarTable_Shannon130-5-0-Mixed_', radius, 'mRadius.csv'))

# use sink to save summary outputs in a txt file

sink(paste0(tabled, '/ModelSummary_Shannon130-5-0-Mixed_', radius, 'mRadius.txt'))
summary(mod130.5.0)
print(paste0("AIC = ", mod130.5.0$aic))
print(paste0("RMSE = ", mod130.5.0_RMSE))
sink()

# # # Plots of fInal model (mod130-5-0)
# draw the basis functions (this is just looking at it with beta splines instead)
bs_default <- basis(s(sd_gapsize, bs='bs'), data=XY_scale)
draw(bs_default)

# draw the basis functions
# from https://fromthebottomoftheheap.net/2020/06/03/extrapolating-with-gams/
tp_default <- basis(s(sd_gapsize), data=XY_scale)
draw(tp_default)

# get prediction values with SE 
p = as_tibble(predict(mod130.5.0, se.fit = TRUE))

# critical value for 95% CIs
crit = qnorm((1 - 0.95) / 2,
             lower.tail = FALSE)

p_CI = data.frame(p) %>% mutate(upr_ci = fit + (crit * se.fit),
                                lwr_ci = fit - (crit * se.fit),
                                sd_gapsize = XY$sd_gapsize)

# Plot x vs y
# built off of plotting code in https://fromthebottomoftheheap.net/2020/06/03/extrapolating-with-gams/
XYplot =  ggplot() +
  geom_ribbon(data = p_CI,
              mapping = aes(ymin = lwr_ci,
                            ymax = upr_ci,
                            x = sd_gapsize,
                            fill='95% CI'),
              inherit.aes = FALSE, alpha = 0.75) +
  scale_fill_manual(values = c("#C0D9D9")) + 
  geom_point(data = XY,
             mapping = aes(x = sd_gapsize,
                           y = shannonH,
                           shape = Treatment,
                           colour = Soil),
             size=3) +
  scale_colour_manual(values = c("grey30", "#D55E00")) + 
  geom_line(data = p_CI,
            aes(y = fit,
                x = sd_gapsize))  + 
  theme(axis.text.y = element_text(colour = "black", size = 12), 
        axis.text.x = element_text(colour = "black", size = 12), 
        legend.text = element_text(size = 12, colour ="black"), 
        legend.position = "right",
        title = element_text(face = "bold", size = 12, colour = "black"),
        legend.title = element_text(size = 12, colour = "black", face = "bold"), 
        legend.key=element_blank()) +
  xlab('Variation (SD) in Canopy Gap Size [m]') +
  ylab('Shannon Diversity') +
  labs(fill='Model Fit',
       shape='Protected Status')
# panel.background = element_blank(),
# panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
# axis.title.x = element_text(face = "bold", size = 12, colour = "black"), 
# axis.title.y = element_text(face = "bold", size = 12, colour = "black"), 
XYplot

ggsave(plot = XYplot,
       filename=paste0(figd,
                       "/GAMs/GAM-130-5-0-Mixed_XYplot_", radius ,"mRadius_mango.png"),
       width = 9, height = 6, units = "in", device='png', dpi=300)

# # # also make a 1-to-1 plot
p1to1 =  ggplot() +
  geom_point(mapping = aes(x = p_CI$fit,
                           y = XY$shannonH,
                           shape = XY$Treatment,
                           colour = XY$Soil),
             size=3) +
  scale_colour_manual(values = c("grey30", "#D55E00")) + 
  geom_abline()  + 
  theme(axis.text.y = element_text(colour = "black", size = 12), 
        axis.text.x = element_text(colour = "black", size = 12), 
        legend.text = element_text(size = 12, colour ="black"), 
        legend.position = "right",
        title = element_text(face = "bold", size = 12, colour = "black"),
        legend.title = element_text(size = 12, colour = "black", face = "bold"), 
        legend.key=element_blank()) +
  xlab('Fitted') +
  ylab('Observed') +
  labs(shape='Protected Status',
       colour='Soil') + 
  coord_fixed(xlim=c(3, 3.8),
              ylim=c(3, 3.8))

# coord_fixed(ratio = 1)
# panel.background = element_blank(),
# panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
# axis.title.x = element_text(face = "bold", size = 12, colour = "black"), 
# axis.title.y = element_text(face = "bold", size = 12, colour = "black"), 
p1to1

ggsave(plot = p1to1,
       filename=paste0(figd,
                       "/GAMs/GAM-130-5-0-Mixed_1to1Plot_", radius ,"mRadius_mango.png"),
       width = 9, height = 6, units = "in", device='png', dpi=300)
