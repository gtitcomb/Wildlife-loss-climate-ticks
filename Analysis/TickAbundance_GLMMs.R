### Interacting effects of wildlife loss and climate on ticks and tick-borne disease

## Georgia Titcomb1,2, Brian Allan2.3, Tyler Ainsworth1, Lauren Henson7, Tyler Hedlund3, Robert M. Pringle2,5, Todd M. Palmer2,6, Laban Njoroge4, Michael G. Campana7, Robert Fleischer7, John Naisikie Mantas2, Hillary S. Young1,2

## 1 Department of Ecology, Evolution and Marine Biology, University of California, Santa Barbara
## 2 Mpala Research Centre, Box 555, Nanyuki, Kenya 
## 3 Department of Entomology, University of Illinois at Urbana-Champaign
## 4 Mammal Section, National Museums of Kenya, Nairobi, Kenya
## 5 Department of Ecology and Evolutionary Biology, Princeton University, Princeton, NJ 08544
## 6 Department of Biology, University of Florida, Gainesville, FL 32611
## 7 Center for Conservation Genomics, Smithsonian Conservation Biology Institute, National Zoological Park, Washington, DC 20008 

## Corresponding author: georgiatitcomb@gmail.com
## Running Head: Defaunation and climate effects on tick-borne disease
## Keywords: ticks, tick-borne disease, defaunation, climate, rainfall, exclosure, Coxiella burnetii, Rickettsia

############################################################################################

# The following code was used to generate generalized linear mixed models of tick abundance

######## Setup

rm(list = ls())

# Load packages
require(lme4)
require(glmmADMB) 
require(MuMIn)
require(bbmle)
require(coefplot2)

# Create internal functions

#plots of data
plots <- function(){
  plot(Total~Treatment)
  plot(Total~Rain)
  plot(Total~Treatment*Rain)
  plot(Total~Replicate)
  plot(Total~Period)
  hist(Total)
}

#model diagnostics
Validation <- function(bestmodel, tickdata){
  Model.Residuals <- residuals(bestmodel, type = 'pearson')
  Model.Fitted <- fitted(bestmodel)
  scatter.smooth(Model.Fitted ~ Total)
  scatter.smooth(Model.Residuals~Model.Fitted)
  hist(Model.Residuals)
  qqnorm(Model.Residuals)
  Validate <- cbind(tickdata, Model.Residuals)
  boxplot(Model.Residuals~Treatment)
  boxplot(Model.Residuals~Rain)
  boxplot(Model.Residuals~(Rain*Treatment))
  boxplot(Model.Residuals~Period)
  boxplot(Model.Residuals~Replicate)
}

### Set working directory to folder containing datafiles ###

########################################################
####### All ticks for CONT vs. LMH, all months #########


# Setup
  Allticks.LC <- read.csv("ALL_LC.csv", header = T, sep = ",")
  Allticks.LC <- transform(Allticks.LC, Rain = scale(Rain), Period = factor(Period), Replicate = factor(Replicate))
  attach(Allticks.LC)  

  plots()

  # Create GLMM
      ALL.LC <- glmmadmb(formula = Total~ Treatment + Rain + Treatment*Rain + (1|Replicate) + (1|Period), data = Allticks.LC, family = "nbinom")

  # Compare all models
      ALL.LC.dredged <- dredge(ALL.LC)
      ALL.LC.dredged
      # Only Model 8 and has deltaAICc < 2
    
      ALL.LC.models <- get.models(ALL.LC.dredged, subset = T)
    
  # Summary of the top model (8) 
    summary(ALL.LC.models$`8`)
   
  # Diagnostic plots for this model
    Validation(ALL.LC.models$`8`, Allticks.LC)
    

# Visualization
    coefplot2(ALL.LC.models$`8`)
    
detach(Allticks.LC)
    
   
############################################################  
####### Ticks for all plots, a subset of months #########
  
  Allticks.ALL <- read.csv("ALL_ALL.csv", header = T, sep = ",")
  Allticks.ALL <- transform(Allticks.ALL, Rain=scale(Rain), Period = factor(Period), Replicate = factor(Replicate))
  attach(Allticks.ALL)
  
  # Investigate the data/distribution
  plots()
  
# Create GLMM
  
       ALL.ALL.m <- glmmadmb(formula = Total~ Treatment + Rain + Treatment*Rain + (1|Replicate) + (1|Period), data = Allticks.ALL, family = "nbinom")
      
  # Compare all models
  
  ALL.ALL.dredged <- dredge(ALL.ALL.m)
  ALL.ALL.dredged
  ALL.ALL.models <- get.models(ALL.ALL.dredged, subset = T)
  
  
  # Summary of the top model (8) 
  summary(ALL.ALL.models$`8`)
  
  
  # Diagnostic plots for this model
  
  Validation(ALL.ALL.models$`8`, Allticks.ALL)
  
 
  detach(Allticks.ALL)


  
############################################################################################################# 
  
# The following code repeats these analyses for each tick species to determine species-specific differences 
  
  ####### RHPV for all months, CONT and LMH only #########
  
  RHPV.LC <- read.csv("RHPV_LC.csv", header = T, sep = ",")
  RHPV.LC <- transform(RHPV.LC, Rain = scale(Rain), Period = factor(Period), Replicate = factor(Replicate))
  attach(RHPV.LC)
  
  # Investigate the data/distribution
  plots()
  
  # Create GLMM. We used zero inflated models here, as nearly 40% of the values are 0. The nonzero data still follows the negative binomial distribution.
  RHPV.LC.m <- glmmadmb(formula = Total~ Treatment + Rain + Treatment*Rain + (1|Replicate) + (1|Period), data = RHPV.LC, family = "nbinom", zeroInflation = T)
  
  # Compare models
  RHPV.LC.dredged <- dredge(RHPV.LC.m)
  RHPV.LC.dredged
  RHPV.LC.models <- get.models(RHPV.LC.dredged, subset = T)

  # Summary of the top model (8) 
  summary(RHPV.LC.models$`8`)
  
  # Diagnostic plots for this model
  Validation(RHPV.LC.models$`8`, RHPV.LC)
  
 
  # Visualization
  coefplot2(RHPV.LC.models$`8`)
  
  detach(RHPV.LC)
  
 
  
  ####### RHPR for all months, CONT and LMH only #########
  
  RHPR.LC <- read.csv("RHPR_LC.csv", header = T, sep = ",")
  RHPR.LC <- transform(RHPR.LC, Rain = scale(Rain), Period = factor(Period), Replicate = factor(Replicate))
  attach(RHPR.LC)
  
  # Investigate the data/distribution
  plots()
  
  # Create GLMM
  RHPR.LC.m <- glmmadmb(formula = Total ~ Treatment + Rain + Treatment*Rain + (1|Replicate) + (1|Period), data = RHPR.LC, family = "nbinom1", zeroInflation = T)

  # Compare models
  RHPR.LC.dredged <- dredge(RHPR.LC.m)
  RHPR.LC.dredged
  RHPR.LC.models <- get.models(RHPR.LC.dredged, subset = T)

  # Summary of the top models (3 and 4) 
  summary(RHPR.LC.models$`3`)
  summary(RHPR.LC.models$`4`)
  
  # Diagnostic plots for this model
  Validation(RHPR.LC.models$`3`, RHPR.LC)
  Validation(RHPR.LC.models$`4`, RHPR.LC)
  
  
  # Visualization
  coefplot2(RHPR.LC.models$`3`)
 
  detach(RHPR.LC)
  
  
####### RHPU for all months, CONT and LMH only #########
  
  RHPU.LC <- read.csv("RHPU_LC.csv", header = T, sep = ",")
  RHPU.LC <- transform(RHPU.LC, Rain = scale(Rain), Period = factor(Period), Replicate = factor(Replicate))
  attach(RHPU.LC)
  
  # Investigate the data/distribution
  plots()
  
  # Create GLMM
  
 RHPU.LC.m <- glmmadmb(formula = Total ~ Treatment + Rain + Treatment*Rain + (1|Replicate) + (1|Period), data = RHPU.LC, family = "nbinom", zeroInflation = T)
 
  
  # Compare models: warning generated indicates that this model may be a poor fit
  
  RHPU.LC.dredged <- dredge(RHPU.LC.m)
  RHPU.LC.dredged
  RHPU.LC.models <- get.models(RHPU.LC.dredged, subset = T)
  
  
  # Summary of the top models (3 and 4) 
  summary(RHPU.LC.models$`3`)
  summary(RHPU.LC.models$`4`)
  
  # Diagnostic plots for this model
  Validation(RHPU.LC.models$`3`, RHPU.LC)
  Validation(RHPU.LC.models$`4`, RHPU.LC)

  
  # Visualization
  coefplot2(RHPU.LC.models$`4`)
  
  
  detach(RHPU.LC)
  
  
#################################################################################
  
  # The following code repeats these analyses for each tick species to determine species-specific differences 
  
####### RHPV for all treatments, subset of months #########
  
  RHPV.ALL <- read.csv("RHPV_ALL.csv", header = T, sep = ",")
  RHPV.ALL <- transform(RHPV.ALL, Rain = scale(Rain), Period = factor(Period), Replicate = factor(Replicate))
  attach(RHPV.ALL)
  
  # Investigate the data/distribution
  plots()
  
  # Create GLMM

 RHPV.ALL.m <- glmmadmb(formula = Total~ Treatment + Rain + Treatment*Rain + (1|Replicate) + (1|Period), data = RHPV.ALL, family = "nbinom", zeroInflation = T)
      
   # Compare Models
  
  RHPV.ALL.dredged <- dredge(RHPV.ALL.m)
  RHPV.ALL.dredged
  RHPV.ALL.models <- get.models(RHPV.ALL.dredged, subset = T)
  
  
  # Summary of the top model (8) 
  summary(RHPV.ALL.models$`8`)
  
  # Diagnostic plots for this model
  Validation(RHPV.ALL.models$`8`, RHPV.ALL)
  
 
  # Visualization
  coefplot2(RHPV.ALL.models$`8`)
  
  
  detach(RHPV.ALL)
  
  
  
####### RHPR for all treatments, subset of months #########
  
  RHPR.ALL <- read.csv("RHPR_ALL.csv", header = T, sep = ",")
  RHPR.ALL <- transform(RHPR.ALL, Rain = scale(Rain), Period = factor(Period), Replicate = factor(Replicate))
  attach(RHPR.ALL)
  
  # Investigate the data/distribution
  plots()
  
  # Create GLMM
  RHPR.ALL.m <- glmmadmb(formula = Total ~ Treatment + Rain + Treatment*Rain + (1|Replicate) + (1|Period), data = RHPR.ALL, family = "nbinom", zeroInflation = T)

 
  # Compare Models
  RHPR.ALL.dredged <- dredge(RHPR.ALL.m)
  RHPR.ALL.dredged
  RHPR.ALL.models <- get.models(RHPR.ALL.dredged, subset = T)
  
  
  # Summary of the top models (8) 
  summary(RHPR.ALL.models$`8`)    # This is odd because RHPR only increases in Meso
  summary(RHPR.ALL.models$`4`)
  summary(RHPR.ALL.models$`3`)
  
  # Diagnostic plots for this model
  Validation(RHPR.ALL.models$`8`, RHPR.ALL)
  
  # Visualization
  coefplot2(RHPR.ALL.models$`8`)
  
  
  detach(RHPR.ALL)
  
  
####### RHPU for all treatments, a subset of months #########
  
  RHPU.ALL <- read.csv("RHPU_ALL.csv", header = T, sep = ",")
  RHPU.ALL <- transform(RHPU.ALL, Rain = scale(Rain), Period = factor(Period), Replicate = factor(Replicate))
  attach(RHPU.ALL)
  
  # Investigate the data/distribution
  plots()
  
  # Create GLMM
  
  RHPU.ALL.m <- glmmadmb(formula = Total~ Treatment + Rain + Treatment*Rain + (1|Replicate) + (1|Period), data = RHPU.ALL, family = "nbinom")

  # Compare Models
  
  RHPU.ALL.dredged <- dredge(RHPU.ALL.m)
  RHPU.ALL.dredged
  RHPU.ALL.models <- get.models(RHPU.ALL.dredged, subset = T)
  
  
  # Summary of the top models (3 and 4) 
  summary(RHPU.ALL.models$`3`)
  summary(RHPU.ALL.models$`4`)

  
  # Diagnostic plots for this model
  Validation(RHPU.ALL.models$`3`, RHPU.ALL)
  Validation(RHPU.ALL.models$`4`, RHPU.ALL)
  

  
    detach(RHPU.ALL)


  
  
  
 
  