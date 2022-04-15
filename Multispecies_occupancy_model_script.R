
# This script produces multispecies distribution map (or species richness map).
# It uses camera trap data but can also be applied to other dataset (e.g., transect survey).
# It accounts for imperfect detection and hence the map produced is the true distribution map.
# This modelling framework uses Bayesian approach.
# Ugyen Penjor, 2021 ugyenpenjor.bt@gmail.com 
# Link to relevant paper: https://doi.org/10.1016/j.ecolind.2020.107085
# Penjor et al., (2021) Vulnerability of mammal communities to the combined impacts of anthropic land-use and climate change in the Himalayan conservation landscape of Bhutan

rm(list=ls(all=TRUE))
ls()

# Load required packages
library(jagsUI)
library(wiqid)

# Load data
load("MSOM_species_detection.RData")
load("MSOM_covariates.RData")
ls()

# Abbreviations
# dset = distance to settlement
# agri = percentage of land under agriculture
# roa = density of road (highway + farm road)
# tem = mean annual temperature of a site
# prec = mean annual precipitation
# ps = precipitation seasonality
# covariates with 'S' indicate standardised (x - mean / sd)
# Ymat = detection history matrices for all observed species (45 in this case)
# nrep = number of sampling occasions
# nSobs = number of observed species
# nSites = number of sites


# Bundle data for JAGS
jagsData <- list(y=Ymat, nSites=nSites, nSobs=nSobs, nRep=nRep, nCamera=nCamera, nLandtype=nLandtype,
                      fore=forestS, set=setS, agri=agriS, roa=roaS, tem=temS, tem2=temS^2, prec=precS, prec2=precS^2, ps=psS, ps2=psS^2,
                      trail=trail, landtype=landtype, camera=camera, tau=1/(2.25^2))
str(jagsData)

# Write model in BUGS (Bayesian analysis using Gibbs Sampler) language
######## BUGS model code######## 

modelText <- "model {

for(k in 1:nSobs) { # loop through species

  # Likelihood
  for(i in 1:nSites) { # loop over sites
    
    # Ecological model
    logit(psi[i, k]) <- lpsi[k] +  
                        betalpsiFore[k] * fore[i] + betalpsiSet[k] * set[i] + betalpsiAgri[k] * agri[i] + 
                        betalpsiRoa[k] * roa[i] + betalpsiTem[k] * tem[i] + betalpsiTem2[k] * tem2[i] +
                        betalpsiPre[k] * prec[i] +betalpsiPre2[k] * prec2[i] + betalpsiPS[k] * ps[i] + betalpsiPS2[k] * ps2[i]
    z[i, k] ~ dbern(psi[i, k])
    
    # Observation model
    logit(p[i, k]) <- lp[k] + betalpLand[landtype[i]] + betalpCamera[camera[i]] + betalpTrail[k] * trail[i] 
    y[i, k] ~ dbin(z[i, k] * p[i, k], nRep[i])
  } # i loop ends here
  
  # Priors  (species level) 
  lpsi[k] ~ dnorm(mu.lpsi, tau.lpsi)
 
  # Correlated prior for psi and p 
  mu.eta[k] <- mu.lp + rho * sd.lp/sd.lpsi * (lpsi[k] - mu.lpsi)
  lp[k] ~ dnorm(mu.eta[k], tau.eta)
  
  # priors on regression parameters
  betalpsiFore[k] ~ dnorm(mu.betalpsiFore, tau.betalpsiFore)
  betalpsiSet[k] ~ dnorm(mu.betalpsiSet, tau.betalpsiSet)
  betalpsiAgri[k] ~ dnorm(mu.betalpsiAgri, tau.betalpsiAgri)
  betalpsiRoa[k] ~ dnorm(mu.betalpsiRoa, tau.betalpsiRoa)
  betalpsiTem[k] ~ dnorm(mu.betalpsiTem, tau.betalpsiTem)
  betalpsiTem2[k] ~ dnorm(mu.betalpsiTem2, tau.betalpsiTem2)
  betalpsiPre[k] ~ dnorm(mu.betalpsiPre, tau.betalpsiPre)
  betalpsiPre2[k] ~ dnorm(mu.betalpsiPre2, tau.betalpsiPre2)
  betalpsiPS[k] ~ dnorm(mu.betalpsiPS, tau.betalpsiPS)
  betalpsiPS2[k] ~ dnorm(mu.betalpsiPS2, tau.betalpsiPS2)
  betalpTrail[k] ~ dnorm(mu.betalpTrail, tau.betalpTrail)
} # k loop ends here

  # Hyperpriors (community level)
  # Occupancy
  lpsi.mean ~ dbeta(1, 1)
  mu.lpsi <- logit(lpsi.mean)  
  sd.lpsi ~ dt(0, tau, 4) T(0, )  # half-t priors with nu = 4; T(0, ) allows positive values only 
  tau.lpsi <- pow(sd.lpsi, -2)

  mu.betalpsiFore ~ dnorm(0, 0.1)
  sd.betalpsiFore ~ dt(0, tau, 4) T(0, )
  tau.betalpsiFore <- 1/(sd.betalpsiFore^2)

  mu.betalpsiSet ~ dnorm(0, 0.1)
  sd.betalpsiSet ~ dt(0, tau, 4) T(0, )
  tau.betalpsiSet <- 1/(sd.betalpsiSet^2)

  mu.betalpsiAgri ~ dnorm(0, 0.1)
  sd.betalpsiAgri ~ dt(0, tau, 4) T(0, )
  tau.betalpsiAgri <- 1/(sd.betalpsiAgri^2)

  mu.betalpsiRoa ~ dnorm(0, 0.1)
  sd.betalpsiRoa ~ dt(0, tau, 4) T(0, )
  tau.betalpsiRoa <- 1/(sd.betalpsiRoa^2)

  mu.betalpsiTem ~ dnorm(0, 0.1)
  sd.betalpsiTem ~ dt(0, tau, 4) T(0, )
  tau.betalpsiTem <- 1/(sd.betalpsiTem^2)

  mu.betalpsiTem2 ~ dnorm(0, 0.1)
  sd.betalpsiTem2 ~ dt(0, tau, 4) T(0, )
  tau.betalpsiTem2 <- 1/(sd.betalpsiTem2^2)

  mu.betalpsiPre ~ dnorm(0, 0.1)
  sd.betalpsiPre ~ dt(0, tau, 4) T(0, )
  tau.betalpsiPre <- 1/(sd.betalpsiPre^2)

  mu.betalpsiPre2 ~ dnorm(0, 0.1)
  sd.betalpsiPre2 ~ dt(0, tau, 4) T(0, )
  tau.betalpsiPre2 <- 1/(sd.betalpsiPre2^2)

  mu.betalpsiPS ~ dnorm(0, 0.1)
  sd.betalpsiPS ~ dt(0, tau, 4) T(0, )
  tau.betalpsiPS <- 1/(sd.betalpsiPS^2)

  mu.betalpsiPS2 ~ dnorm(0, 0.1)
  sd.betalpsiPS2 ~ dt(0, tau, 4) T(0, )
  tau.betalpsiPS2 <- 1/(sd.betalpsiPS2^2)

  # Detection
  lp.mean ~ dbeta(1, 1)
  mu.lp <- logit(lp.mean)
  sd.lp ~ dt(0, tau, 4) T(0, )
  tau.lp <- 1/(sd.lp^2)

  # For random land type effect
  betalpLand[1] <- 0
  for (l in 2:nLandtype) {      # loop over land types
    betalpLand[l] ~ dnorm(mu.betalpLand, tau.betalpLand)
  }
  mu.betalpLand ~ dnorm(0, 0.1)
  sd.betalpLand  ~ dt(0, tau, 4) T(0, )
  tau.betalpLand  <- 1/(sd.betalpLand^2)

  # For random camera effect
  betalpCamera[1] <- 0
  for (c in 2:nCamera) {      # loop over camera models
    betalpCamera[c] ~ dnorm(mu.betalpCamera, tau.betalpCamera)
  }
  mu.betalpCamera ~ dnorm(0, 0.1)
  sd.betalpCamera ~ dt(0, tau, 4) T(0, )
  tau.betalpCamera <- 1/(sd.betalpCamera^2)

  mu.betalpTrail ~ dnorm(0, 0.1)
  sd.betalpTrail ~ dt(0, tau, 4) T(0, )
  tau.betalpTrail <- 1/(sd.betalpTrail^2)

  rho ~ dunif(-1, 1)
  tau.eta <- tau.lp/(1 - rho^2)

  # Derived quantities 
    for(k in 1:nSobs) {
      Nocc.fs[k] <- sum(z[, k]) # No. of occupied sites for each species
    }
  
    for(i in 1:nSites) {
      Nsite[i] <- sum(z[i, ]) # No of occuring species at each site
    }
  
  # Derived parameters: Mean occupancy and detection per species
    for(k in 1:nSobs) {
      species.psi[k] <- ilogit(lpsi[k]) # inverse logit
      species.p[k] <- ilogit(lp[k])
    }
  
}"
writeLines(modelText, con="msom_model_code.jags")

################################################################################################################################################

# Initial values - to train the model to pick appropriate intial steps (Gibbs sampling)
inits <- function() list(z = matrix(1, nSites, nSobs))

# Parameters required
wanted <- c("mu.lpsi", "sd.lpsi", 
             "mu.betalpsiFore", "sd.betalpsiFore", "mu.betalpsiSet", "sd.betalpsiSet",
             "mu.betalpsiAgri", "sd.betalpsiAgri",
             "mu.betalpsiRoa", "sd.betalpsiRoa",
             "mu.betalpsiTem", "sd.betalpsiTem", "mu.betalpsiTem2", "sd.betalpsiTem2",
             "mu.betalpsiPre", "sd.betalpsiPre", "mu.betalpsiPre2", "sd.betalpsiPre2",
             "mu.betalpsiPS", "sd.betalpsiPS", "mu.betalpsiPS2", "sd.betalpsiPS2",
             "mu.lp", "sd.lp", 
             "mu.betalpCamera", "sd.betalpCamera", "mu.betalpTrail", "sd.betalpTrail", "mu.betalpLand", "sd.betalpLand",
             "lpsi", "lp",
             "betalpsiFore", "betalpsiSet", 
             "betalpsiAgri",
             "betalpsiRoa",
             "betalpsiTem", "betalpsiTem2",
             "betalpsiPre", "betalpsiPre2",
             "betalpsiPS", "betalpsiPS2",
             "betalpCamera", "betalpTrail", "betalpLand",
             "z")

# Call JAGS from R ( ~ 12 mins) # with 100K iterations, it takes days to run the model
# Requires a good number of cores (at least 12) to run this heavily parameterised model
out1 <- jags.basic(jagsData, inits, wanted, "msom_model_code.jags",
                        n.chains=3, n.iter=1500, n.adapt=400, n.burnin=500, n.thin=2, parallel=TRUE, DIC=FALSE)

# This is a heavily parameterised model, so takes time to run. 
# To save time, you can load the results (that I ran for you already!) and prepare species richness map.
load("MSOM_out1_results.RData")

# Convert the output to matrix 
all0 <- as.matrix(out1)
str(all0)
nms <- colnames(all0)

### SPECIES RICHNESS MAP (MULTISPECIES DISTRIBUTION MODELLING CORRECTED FOR IMPERFECT DETECTION)

library(raster)
library(rgdal)
library(sp)

# Standardise pixels (here I have used 1 km by 1 km resolution for the country size of ~39000 km2)
# A finer resolution will be better for small study area but it is computationally expensive because it works on each pixel.
FOREST <- (Bhutan$Sample_tif4.treeCover90m_EasternHimalayas_BNG_Band_1 - mean(forest)) / sd(forest)
DSET <- (Bhutan$Sample_tif4.EuclideanDist_SettlementGEODESIC_Band_1 - mean(dset)) / sd(dset)
ROAD <- (Bhutan$Sample_tif4.EuclideanDist_road90mGEODESIC_Band_1 - mean(roa)) / sd(roa)
AGRI_2 <- (Bhutan$Sample_tif2.agriPLAND2km_Band_1 - mean(agri)) / sd(agri)
TEMP <- (Bhutan$Sample_tif4.CHELSA_bio1a_90m_extbyMask_Band_1 - mean(tem)) / sd(tem)
TEMP2 <- TEMP^2
PRECIP <- (Bhutan$Sample_tif4.CHELSA_bio12a_90m_extbyMask_Band_1 - mean(prec)) / sd(prec)
PRECIP2 <- PRECIP^2
PS <- (Bhutan$Sample_tif4.CHELSA_bio10_15_PrecSea_90m_Band_1 - mean(ps)) / sd(ps)
PS2 <- PS^2

nkm2 <- nrow(Bhutan) # 1 km2 pixels across Bhutan
nsamp <- nrow(all0) # total sample size, too much?
select.samp <- sort(sample(1:nsamp, 100)) # only taking out 100 for faster run
nsamp2 <- length(select.samp)

# Get the coefficients and assign them to each pixel 
# Fine resolution (e.g., 250m x 250 or even 90m by 90m if your study area is small) mapping with take time!
INTERCEPT <- all0[select.samp, grep("^lpsi", nms)]
BETAFOREST <- all0[select.samp, grep("^betalpsiFore", nms)]
BETADSET <- all0[select.samp, grep("^betalpsiSet", nms)]
BETAAGRI <- all0[select.samp, grep("^betalpsiAgri", nms)]
BETAROAD <- all0[select.samp, grep("^betalpsiRoa", nms)]
BETATEMP <- all0[select.samp, grep("^betalpsiTem", nms)]
BETATEMP2 <- all0[select.samp, grep("^betalpsiTem2", nms)]
BETAPREC <- all0[select.samp, grep("^betalpsiPre", nms)]
BETAPREC2 <- all0[select.samp, grep("^betalpsiPre2", nms)]
BETAPS <- all0[select.samp, grep("^betalpsiPS", nms)]
BETAPS2 <- all0[select.samp, grep("^betalpsiPS2", nms)]

# Create an empty array to hold the predicted values
zBH <- array(NA, dim = c(nkm2, nSobs, nsamp2))

# Assign prediction to each pixel
for(i in 1:nkm2) {
  cat(paste("\nPixel", i, "\n"))
  for(u in 1:length(select.samp)) {
    psi <- plogis( INTERCEPT[u, ] + 
                     BETAFOREST[u, ] * FOREST[i] + 
                     BETADSET[u, ] * DSET[i] + 
                     BETAAGRI[u, ] * AGRI_2[i] + 
                     BETAROAD[u, ] * ROAD[i] + 
                     BETATEMP[u, ] * TEMP[i] + BETATEMP2[u, ] * TEMP2[i] +
                     BETAPREC[u, ] * PRECIP[i] + BETAPREC2[u, ] * PRECIP2[i] +
                     BETAPS[u, ] * PS[i] + BETAPS2[u, ] * PS2[i] )
    zBH[i,,u] <- rbinom(nSobs, 1, psi)
  }
}

# Compute posterior distribution of species richness by collapsing z array
SR <- apply(zBH, c(1, 3), sum)
pmSR <- apply(SR, 1, mean)
sdSR <- apply(SR, 1, sd)

head(pmSR); head(sdSR)

# Posterior mean - Species richness
myCol <- colorRampPalette(c('white', 'mediumpurple3', 'yellow', 'red'))

ppdm <- rasterFromXYZ(data.frame(x=Bhutan$X1x1_points.POINT_X, y=Bhutan$X1x1_points.POINT_Y, z=pmSR))
plot(ppdm, col= myCol(100), axes=F, box=F, main='Posterior mean of species richness', legend.shrink=0.35,
     legend.args=list(text="Species richness", side=4, line=-2, font=2))
# Adding scale, boundary, axes, and north arrow can be performed later with few lines of codes.

# Posterior SD - this accounts for prediction uncertainty
ppdsd <- rasterFromXYZ(data.frame(x=Bhutan$X1x1_points.POINT_X, y=Bhutan$X1x1_points.POINT_Y, z=sdSR))
plot(ppdsd, col= myCol(100), axes=F, box=F, main='Posterior SD of species richness', mar=c(1,1,0,1))

# Save output raster for post processing (can import in GIS software)
writeRaster(ppdm, "species_richness_plot.tif")

######################################################################################################
################################################ END #################################################
######################################################################################################
