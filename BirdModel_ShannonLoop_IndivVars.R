# # # Model Loop for Shannon Index
# PB 5/4/23 - updated to include individual variable models as well

library(tidyverse)
library(mgcv)
library(gratia)
library(gridExtra)
library(broom)
library(ggplot2)
library(qpcR) # for aikike weights

# Datadir for location of in/out vars
datadir = '/n/home02/pbb/scripts/SelenkayDiversity/data/'

# Note: figures disabled in this run
# figd = '/n/home02/pbb/scripts/SelenkayDiversity/figs/mango/GAMs/ShannonLoop'

tabled = '/n/home02/pbb/scripts/SelenkayDiversity/tables/mango/GAMs/ShannonLoop_IndivVars'

for (radius in c(10, 20, 30, 50, 80, 130)){
  
  # # # 1) Load XY Data
  # Note: XY data prepared in "BirdDataPreperation.R"
  # RDF files keep relevant variable types (factors, characters, etc.)
  XY = readRDS(paste0(datadir, 'in/XY_', radius,'mRadius.rds'))
  XY_scale = readRDS(paste0(datadir, 'in/XY_scaled_', radius,'mRadius.rds'))
  
  # make a logrichness variable
  # XY_scale = XY_scale %>% mutate(logrich = log(Richness))
  
  # # # Grab top 5 important variables 
  # (calculated with spearman corr in python)
  # Note these are hard-coded
  if (radius == 10){
    
    # top 3 black, then top 3 red
    m = c("sd_CD_AboveG",
          "sd_herbh",
          "sd_PAI_AboveG",
          "VDRpeak_plot",
          "cv_CD_AboveG",
          "cv_CD_AboveGgrasslayer")
    
  } else if (radius == 20) {
    
    # top 3 black, then top 3 red
    m = c("Cover1p5m_plot",
          "X25thPerc_plot",
          "Cover0p5m_plot",
          "cv_FHD",
          "cv_cvpeakh",
          "cv_CD_AboveG")
    
  } else if (radius == 30) {
    
    # top 3 black, then top 3 red
    m = c("mean_stdpeakh",
          "Cover0p5m_plot",
          "mean_gapsize",
          "sd_ptoh",
          "cv_FHD",
          "cv_CD_G")
    
  } else if (radius == 50) {
    
    # top 3 black, then top 3 red
    m = c("meanH_plot",
          "Cover1p5m_plot",
          "sd_CD_AboveG",
          "cvpeakh_plot",
          "cv_cscore",
          "cv_FHD")
  
  } else if (radius == 80) {
    
    # top 3 black, then top 3 red
    m = c("Cover1p5m_plot",
          "mean_gapsize",
          "mean_stdpeakh",
          "meanpeakh_plot",
          "ptoh_plot",
          "stdpeakh_plot")
    
  } else if (radius == 130) {
    
    # top 3 black, then top 3 red
    m = c("mean_gapsize",
          "sd_gapsize",
          "mean_stdpeakh",
          "VDRpeak_plot",
          "gapsize_plot",
          "maxpeakh_plot")
    
    # Top 6 vars - just doesn't work as well!
    # m = c("X100thPerc_plot",
    #       "cscore_plot",
    #       "X98thPerc_plot",
    #       "maxpeakh_plot",
    #       "gapsize_plot",
    #       "stdpeakh_plot")
    
  }
  
  # Initialize List of aic values for weight calculation
  aic = c()
  # also a list of labels for the models
  modlabels = c()
  # Also a list of radii so that you can stick it all in a df after
  radlist = c()
  
  # for each individual variable at this radius
  for (v in m) {
    
    # # Fit a fixed effects model
    formula = as.formula(paste0('shannonH ~ s(', v,') + s(X,Y, k=5)'))
    
    # using scaled t dist - long tailed t here
    # note - select = false for individual model fitting
    mod = gam(formula,
              data=XY_scale,
              family=gaussian,
              select=FALSE,
              method="REML")
    
    # compute RMSE
    mod_RMSE = sqrt(mean(mod$residuals**2))
    
    # Same some important results
    vartbl_mod = tidy(mod)
    write_csv(vartbl_mod,
              paste0(tabled, '/', radius, 
                     'mRadius/ModVarTable_Shannon_',
                     v, '_',
                     radius, 'mRadius.csv'))
    
    # use sink to save summary outputs in a txt file
    sink(paste0(tabled,'/', radius, 
                'mRadius/ModSummary_Shannon_',
                v, '_',
                radius, 'mRadius.txt'))
    print(summary(mod))
    print(paste0("AIC = ", mod$aic))
    print(paste0("RMSE = ", mod_RMSE))
    sink()
    
    # add aic to list, as well as variable and radius to label the model
    aic = aic %>% append(mod$aic)
    modlabels = modlabels %>% append(v)
    radlist = radlist %>% append(radius)
    
  }
  
  
  # # FIT the fixed effects model
  # https://stackoverflow.com/questions/26889240/looping-many-one-sided-anova-in-r
  formula = as.formula(paste0('shannonH ~ s(',
                              m[1], ', k=5) + s(',
                              m[2], ', k=5) + s(',
                              m[3], ', k=5) + s(',
                              m[4], ', k=5) + s(',
                              m[5], ', k=5) + s(',
                              m[6], ', k=5) + s(X,Y, k=5)'))
  
  mod = gam(formula,
            data=XY_scale,
            family=gaussian,
            select=TRUE,
            method="REML")
  
  # compute RMSE
  mod_RMSE = sqrt(mean(mod$residuals**2))
  
  # Same some important results
  vartbl_mod = tidy(mod)
  write_csv(vartbl_mod,
            paste0(tabled, '/', radius, 
                   'mRadius/ModVarTable_Shannon_',
                   radius, 'mRadius.csv'))
  
  # use sink to save summary outputs in a txt file
  sink(paste0(tabled,'/', radius, 
              'mRadius/ModSummary_Shannon_',
              radius, 'mRadius.txt'))
  print(summary(mod))
  print(paste0("AIC = ", mod$aic))
  print(paste0("RMSE = ", mod_RMSE))
  sink()
  
  # add aic to list, as well as variable and radius to label the model
  aic = aic %>% append(mod$aic)
  modlabels = modlabels %>% append('Combined')
  radlist = radlist %>% append(radius)
  
  # Now that all the models have been fit (should be 6 of them)
  # Calculate aic weights
  aicweights = akaike.weights(aic)
  
  # add model number to df
  aicweights = data.frame(aicweights) %>% 
    mutate(model_label = modlabels,
           radius = radius)
  
  write_csv(aicweights,
            paste0(tabled,'/',
                   radius,
                   'mRadius/AllModels_', radius,
                   '_AkaikeWeights.csv'))
  
  # # # PLOTS
  # get prediction values
  # fitted = data.frame(fitted = predict(mod))
  # 
  # # # # Save out a 1-to-1 plot
  # p1to1 =  ggplot() +
  #   geom_point(mapping = aes(x = fitted$fitted,
  #                            y = XY$shannonH,
  #                            shape = XY$Treatment,
  #                            colour = XY$Soil),
  #              size=3) +
  #   scale_colour_manual(values = c("grey30", "#D55E00")) + 
  #   geom_abline()  + 
  #   theme(axis.text.y = element_text(colour = "black", size = 12), 
  #         axis.text.x = element_text(colour = "black", size = 12), 
  #         legend.text = element_text(size = 12, colour ="black"), 
  #         legend.position = "right",
  #         title = element_text(face = "bold", size = 12, colour = "black"),
  #         legend.title = element_text(size = 12, colour = "black", face = "bold"), 
  #         legend.key=element_blank()) +
  #   xlab('Fitted') +
  #   ylab('Observed') +
  #   labs(shape='Protected Status',
  #        colour='Soil') + 
  #   coord_fixed(xlim=c(3, 3.8),
  #               ylim=c(3, 3.8))
  # 
  # # coord_fixed(ratio = 1)
  # 
  # # p1to1
  # 
  # ggsave(plot = p1to1,
  #        filename=paste0(figd,'/', radius, 
  #                        "mRadius/Shannon_1to1Plot_",
  #                        radius,"mRadius_mango.png"),
  #        width = 9, height = 6, units = "in", device='png', dpi=300)
  # 
  # 
  # # # # Also, fit a random effects model 
  # # this time, using the top 2 vars from each soil type
  # # with a RE of soil type
  # mixformula = as.formula(paste0('shannonH ~ s(',
  #                             m[1], ', by=Soil_f, k=5) + s(',
  #                             m[2], ', by=Soil_f, k=5) + s(',
  #                             m[4], ', by=Soil_f, k=5) + s(',
  #                             m[5], ', by=Soil_f, k=5) + s(X,Y, k=5)'))
  # 
  # mixmod = gam(mixformula,
  #           data=XY_scale,
  #           family=gaussian,
  #           select=TRUE,
  #           method="REML")
  # 
  # # compute RMSE
  # mixmod_RMSE = sqrt(mean(mixmod$residuals**2))
  # 
  # # Same some important results
  # vartbl_mixmod = tidy(mixmod)
  # write_csv(vartbl_mixmod,
  #           paste0(tabled, '/', radius, 
  #                  'mRadius/MixModVarTable_Shannon_',
  #                  radius, 'mRadius.csv'))
  # 
  # # use sink to save summary outputs in a txt file
  # sink(paste0(tabled,'/', radius, 
  #             'mRadius/MixModSummary_Shannon_',
  #             radius, 'mRadius.txt'))
  # print(summary(mixmod))
  # print(paste0("AIC = ", mixmod$aic))
  # print(paste0("RMSE = ", mixmod_RMSE))
  # sink()
  # 
  # 
  # # # # PLOTS
  # # get prediction values
  # mixfitted = data.frame(fitted = predict(mixmod))
  # 
  # # # # Save out a 1-to-1 plot
  # mixp1to1 =  ggplot() +
  #   geom_point(mapping = aes(x = mixfitted$fitted,
  #                            y = XY$shannonH,
  #                            shape = XY$Treatment,
  #                            colour = XY$Soil),
  #              size=3) +
  #   scale_colour_manual(values = c("grey30", "#D55E00")) + 
  #   geom_abline()  + 
  #   theme(axis.text.y = element_text(colour = "black", size = 12), 
  #         axis.text.x = element_text(colour = "black", size = 12), 
  #         legend.text = element_text(size = 12, colour ="black"), 
  #         legend.position = "right",
  #         title = element_text(face = "bold", size = 12, colour = "black"),
  #         legend.title = element_text(size = 12, colour = "black", face = "bold"), 
  #         legend.key=element_blank()) +
  #   xlab('Fitted') +
  #   ylab('Observed') +
  #   labs(shape='Protected Status',
  #        colour='Soil') + 
  #   coord_fixed(xlim=c(3, 3.8),
  #               ylim=c(3, 3.8))
  # 
  # # coord_fixed(ratio = 1)
  # 
  # # p1to1
  # 
  # ggsave(plot = mixp1to1,
  #        filename=paste0(figd,'/', radius, 
  #                        "mRadius/MixMod_Shannon_1to1Plot_",
  #                        radius,"mRadius_mango.png"),
  #        width = 9, height = 6, units = "in", device='png', dpi=300)
  
  
}
