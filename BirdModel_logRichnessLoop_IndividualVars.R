# # # Model Loop for logrichness
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

# Note figures disabled in this run
# figd = '/n/home02/pbb/scripts/SelenkayDiversity/figs/mango/GAMs/logRichnessLoop/LogRichnessLoop_IndivVars'

tabled = '/n/home02/pbb/scripts/SelenkayDiversity/tables/mango/GAMs/logRichnessLoop_IndivVars'

for (radius in c(10, 20, 30, 50, 80, 130)){
  
  # # # 1) Load XY Data
  # Note: XY data prepared in "BirdDataPreperation.R"
  # RDF files keep relevant variable types (factors, characters, etc.)
  XY = readRDS(paste0(datadir, 'in/XY_', radius,'mRadius.rds'))
  XY_scale = readRDS(paste0(datadir, 'in/XY_scaled_', radius,'mRadius.rds'))
  
  # make a logrichness variable
  XY_scale <- XY_scale %>% mutate(logrich = log(Richness))
  XY <- XY %>% mutate(logrich = log(Richness))
  
  # # # Grab top 5 important variables 
  # (calculated with spearman corr in python)
  # Note these are hard-coded
  if (radius == 10){
    
    # top 3 black, then top 3 red
    m = c("Cover0p5m_plot",
          "Cover0p25m_plot",
          "sd_nlayers",
          "mean_cvH_tree",
          "VDRpeak_plot",
          "sd_cvH_vegtype_tree"
          )
    
  } else if (radius == 20) {
    
    # top 3 black, then top 3 red
    m = c("sd_nlayers",
          "cv_nlayers",
          "horzcover_tree",
          "horzcover_grass",
          "mean_PAI_G",
          "mean_CD_Ggrasslayer")
    
  } else if (radius == 30) {
    
    # top 3 black, then top 3 red
    m = c("sd_nlayers",
          "cv_nlayers",
          "mean_stdpeakh",
          "horzcover_grass",
          "mean_PAI_G",
          "cv_cvpeakh")
    
  } else if (radius == 50) {
    
    # top 3 black, then top 3 red
    m = c("sd_nlayers",
          "meanH_plot",
          "cv_nlayers",
          "horzcover_grass",
          "cv_CD_AboveG",
          "cvpeakh_plot")
  
  } else if (radius == 80) {
    
    # top 3 black, then top 3 red
    m = c("meanH_plot",
          "stdH_plot",
          "Cover1p5m_plot",
          "horzcover_grass",
          "cv_VDRpeak",
          'cv_cvpeakh')
    
  } else if (radius == 130) {
    
    # top 3 black, then top 3 red
    m = c("X98thPerc_plot",
          "sd_gapsize",
          "stdH_plot",
          "horzcover_grass",
          "mean_cvH_tree",
          "mean_cvH_vegtype_woody"
    )
    
  }
  
  # Initialize List of aic values for weight calculation
  aic = c()
  # also a list of labels for the models
  modlabels = c()
  # Also a list of radiuses so that you can stick it all in a df after
  radlist = c()
  
  # for each individual variable at this radius
  for (v in m) {
    
    # # Fit a fixed effects model
    formula = as.formula(paste0('logrich ~ s(', v,') + s(X,Y, k=5)'))
    
    # using scaled t dist - long tailed t here
    # note - select = false for individual model fitting
    mod = gam(formula,
              data=XY_scale,
              family=scat,
              select=FALSE,
              method="REML")
    
    # compute RMSE
    mod_RMSE = sqrt(mean(mod$residuals**2))
    
    # Same some important results
    vartbl_mod = tidy(mod)
    write_csv(vartbl_mod,
              paste0(tabled, '/', radius, 
                     'mRadius/ModVarTable_logRich_',
                     v, '_',
                     radius, 'mRadius.csv'))
    
    # use sink to save summary outputs in a txt file
    sink(paste0(tabled,'/', radius, 
                'mRadius/ModSummary_logRich_',
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
  
  # Fit a combined model 
  # and use variable selection (similar to ridge regression)
  formula = as.formula(paste0('logrich ~ s(',
                              m[1], ', k=5) + s(',
                              m[2], ', k=5) + s(',
                              m[3], ', k=5) + s(',
                              m[4], ', k=5) + s(',
                              m[5], ', k=5) + s(',
                              m[6], ', k=5) + s(X,Y, k=5)'))
  
  # using scaled t dist - long tailed t here
  mod = gam(formula,
            data=XY_scale,
            family=scat,
            select=TRUE,
            method="REML")
  
  # compute RMSE
  mod_RMSE = sqrt(mean(mod$residuals**2))
  
  # Same some important results
  vartbl_mod = tidy(mod)
  write_csv(vartbl_mod,
            paste0(tabled, '/', radius, 
                   'mRadius/ModVarTable_logRich_CombinedModel_',
                   radius, 'mRadius.csv'))
  
  # use sink to save summary outputs in a txt file
  sink(paste0(tabled,'/', radius, 
              'mRadius/ModSummary_logRich_CombinedModel_',
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
  # NOTE: really weird stuff happening
  # only getting 35 values back for XY$Richness...
  # Model will only predict 35 values... 
  # no idea why
  # fitted = data.frame(fitted = predict(mod))
  # 
  # dfplot = fitted %>% mutate(logrich = log(XY$Richness),
  #                            Treatment = XY$Treatment_f,
  #                            Soil = XY$Soil_f)
  # 
  # # # # Save out a 1-to-1 plot
  # p1to1 =  ggplot() +
  #   geom_point(mapping = aes(x = fitted,
  #                            y = XY_scale$logrich,
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
  # p1to1
  # 
  # ggsave(plot = p1to1,
  #        filename=paste0(figd,'/', radius, 
  #                        "mRadius/logrich_1to1Plot_",
  #                        radius,"mRadius_mango.png"),
  #        width = 9, height = 6, units = "in", device='png', dpi=300)
  # 
  # 
  
  # # # TBD ISSUE HERE - Mixed model doesn't work because of the strange 35 vars thing
  # figure out another day 3/3/23
  # # # Also, fit a random effects model
  # this time, using the top 2 vars from each soil type
  # with a RE of soil type
  # mixformula = as.formula(paste0('logrich ~ s(',
  #                             m[1], ', by=Soil_f, k=5) + s(',
  #                             m[2], ', by=Soil_f, k=5) + s(',
  #                             m[4], ', by=Soil_f, k=5) + s(',
  #                             m[5], ', by=Soil_f, k=5) + s(X,Y, k=5)'))
  # 
  # mixmod = gam(mixformula,
  #           data=XY_scale,
  #           family=scat,
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
  #                  'mRadius/MixModVarTable_logrich_',
  #                  radius, 'mRadius.csv'))
  # 
  # # use sink to save summary outputs in a txt file
  # sink(paste0(tabled,'/', radius,
  #             'mRadius/MixModSummary_logrich_',
  #             radius, 'mRadius.txt'))
  # print(summary(mixmod))
  # print(paste0("AIC = ", mixmod$aic))
  # print(paste0("RMSE = ", mixmod_RMSE))
  # sink()

  # 
  # # # # PLOTS
  # # get prediction values
  # mixfitted = data.frame(fitted = predict(mixmod))
  # 
  # # # # Save out a 1-to-1 plot
  # mixp1to1 =  ggplot() +
  #   geom_point(mapping = aes(x = mixfitted$fitted,
  #                            y = XY$logrichH,
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
  #                        "mRadius/MixMod_logrich_1to1Plot_",
  #                        radius,"mRadius_mango.png"),
  #        width = 9, height = 6, units = "in", device='png', dpi=300)
  
  
}
