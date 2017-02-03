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

# The following code was used to generate generalized linear mixed models of tick infection with Coxiella burnetii and Rickettsia sp.

###### Load required packages

require(glmmADMB)
require(MuMIn)

### Set working directory to folder containing datafiles ###

# Setup
pathogens <- read.csv("Pathogens_CB_R.csv", header = T, sep = ",")
pathogens <- transform(pathogens, Species = factor(Species), Replicate = factor(Replicate), Rainfall = scale(Rainfall), Infection_CB = factor(Infection_CB), Infection_R = factor(Infection_R))

##### C. burnetii

# Create GLMM. Species, Treatment, Rainfall, and the interaction of Treatment and Rainfall are fixed effects; replicates at rainfall levels are random effects
CB.glmm <- glmmadmb(Infection_CB ~ Species + Treatment + Rainfall + Treatment*Rainfall + (1|Replicate), data = pathogens, family = "binomial", link = "logit")

# Compare models
dd.CB <- dredge(CB.glmm)
dd.CB
top.CB <- get.models(dd.CB, subset = T)
summary(top.CB$`1`)
summary(top.CB$`3`)
summary(top.CB$`2`)

# Plots
plot(Infection_CB ~ Species, data = pathogens)
plot(Infection_CB ~ Treatment, data = pathogens)
plot(Infection_CB ~ Rainfall, data = pathogens)
plot(Infection_CB ~ Replicate, data = pathogens)



##### Rickettsia sp.

# Create GLMM
R.glmm <- glmmadmb(Infection_R ~ Species + Treatment + Rainfall + Treatment*Rainfall + (1|Replicate), data = pathogens, family = "binomial", link = "logit")

# Compare models
dd.R <- dredge(R.glmm)
dd.R
top.R <- get.models(dd.R, subset = T)
summary(top.R$`4`)
summary(top.R$`3`)
summary(top.R$`2`)

# Plots
plot(Infection_R ~ Species, data = pathogens)
plot(Infection_R ~ Treatment, data = pathogens)
plot(Infection_R ~ Rainfall, data = pathogens)
plot(Infection_R ~ Replicate, data = pathogens)