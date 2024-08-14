# Load packages
library(Hmsc)
library(tidyverse)

# Load data
load('HMSCdat_full.RData')

# Assign the number of chains to run
nChains = 2

# Assign whether to run a test version or the full analysis
test.run = T

# Assign runtime parameters for monte-carlo markov chain sampling based on testing or running the full analysis
if (test.run){
  # With this option, the analysis runs fast but results are not reliable
  thin = 1
  samples = 5
  transient = 10
  verbose = 1
  npar = 2
} else {
  # With this option, the analysis evaluates slowly (and may nto run at all based on memory limitations) but it produces the full analysis
  thin = 100
  samples = 500
  transient = 25000
  verbose = 1
  npar = 2
}

# Initalize and run the HMSC model
hmsc_mods <- hmsc_dat %>%
  mutate(
  ref = map(Reference, ~.$Reference),
  bay = map(Bay, ~.$Bay),
  
  # Assign an index for each geographical location; some lat longs are duplicated across the dataset, so we need to assign a unique identifier for each location
  geo_ind = map(spat, ~group_indices(group_by(., Latitude, Longitude))),
  
  # Initialize the HMSC model
  hmsc_mod = pmap(list(geo_ind,bay,ref, predictors,spat,Year, Month,bio),function(i,j,u,v,w,x,y,z){
    # Assign the date as a decimal year
    date_YM <- as.matrix(x + y/12)
    
    # Convert spatial data to matrix
    spat <- as.matrix(cbind((w$Longitude), (w$Latitude)))
    
    # Assign unique row names to spatial and temporal data (required for HMSC)
    rownames(spat) <- i
    rownames(date_YM) <- date_YM
    
    # Set study design dataframe with location, year, decimal date, and sample ID
    studyDesign <- data.frame(Loc = as.factor(i),Time = as.factor(date_YM),Bay = as.factor(j),Sample = as.factor(u), Year = as.factor(x))
    
    # Assign row names of study design matrix
    rownames(studyDesign) <- u
    
    # Initialize random levels for the HMSC model
    rl1 <- HmscRandomLevel(units = studyDesign$Bay)
    rl2 <- HmscRandomLevel(sData = unique(date_YM))
    rl2 <- setPriors(rl2,nfMin = 1, nfMax = 8)
    rl3 <- HmscRandomLevel(sData = unique(spat), sMethod = "NNGP", nNeighbours = 20)
    rl3 <- setPriors(rl3,nfMin = 1, nfMax = 8)
    rl4 <- HmscRandomLevel(units = studyDesign$Sample)
    
    # Set spike-slab covariate selection parameters
    qq <- 0.1 # Prior probability required for inclusion set at 0.1
    nc <- ncol(v)
    ns <- ncol(z)
    
    
    XSelect <- list()
    for (k in 1:nc){
      covGroup <- k
      spGroup <- 1:ns
      q <- rep(qq,max(spGroup))
      XSelect[[k]] <- list(covGroup = covGroup,
                           spGroup = spGroup,
                           q = q)
    }
    
    # Initialize the HMSC model with a lognormal poisson
    m <- Hmsc(Y = z, 
              XData = v, 
              XFormula = ~., 
              XSelect = XSelect,
              studyDesign = studyDesign,
              ranLevels = list(Bay = rl1, 
                               Time = rl2, 
                               Loc  = rl3,
                               Sample = rl4),
              distr = "lognormal poisson")
    # Compute data parameters
    m_datpar <- computeDataParameters(m)
    return(list(mod = m, datpar = m_datpar))}))

# Isolate mode and data parameters
mod <- hmsc_mods$hmsc_mod[[1]]$mod
datpar <- hmsc_mods$hmsc_mod[[1]]$datpar

# Clean up workspace
rm(list = setdiff(ls(),list('mod','datpar', 'thin','samples', 'transient','verbose','npar')))

# Run the monte-carlo markov chain for each model
mod <- sampleMcmc(mod, 
                  thin = thin, 
                  samples = samples, 
                  transient = transient,
                  nChains = nChains,
                  nParallel = npar, 
                  verbose = verbose, 
                  dataParList = datpar,
                  updater = list(GammaEta=FALSE,
                                 Gamma2=FALSE),
                  useSocket = F)

# Save the model
save(mod, file = "HMSC_full.RData")
