# Load packages
library(Hmsc)
library(tidyverse)

# Load data
load('HMSCdat_bay.RData')

# Assign which bays to process; this option was used to avoid memory limitations
bay <- c('AP','CK','TB','CH','JX','IR','TQ')

# Assign the number of chains to run
nChains = 4

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
  # With this option, the analysis evaluates slowly (and may not run at all based on memory limitations), but produces the full analysis
  thin = 100
  samples = 500
  transient = thin*samples
  verbose =1
  npar = 4
}

# Initalize and run the HMSC model for each bay
hmsc_mods <- hmsc_dat %>%
  filter(Bay %in% bay) %>%
  mutate(
    ref = map(Reference, ~.$Reference),
    
    # Assign an index for each geographical location; some lat longs are duplicated across the dataset, so we need to assign a unique identifier for each location
    geo_ind = map(spat, ~group_indices(group_by(., Latitude, Longitude))),
    
    # Initialize the HMSC model for each bay
    hmsc_mod = pmap(list(geo_ind,ref, predictors,spat,Year, Month,bio),function(i,u,v,w,x,y,z){
      
      # Filter the biotic data for hyper-rare species (occurring in less than 0.5% of samples)
      z <- z %>%
        select(where(~sum(.x != 0)/length(.x) > 0.005))
      
      # Assign date as decimal year
      date_YM <- as.matrix(x + y/12)
      
      # Convert spatial data to matrix
      spat <- as.matrix(cbind((w$Longitude), (w$Latitude)))
      
      # Assign unique row names to spatial and temporal data (required for HMSC)
      rownames(spat) <- i
      rownames(date_YM) <- date_YM
      
      # Set study design dataframe with location, year, decimal date, and sample ID
      studyDesign <- data.frame(Loc = as.factor(i),Year = as.factor(as.matrix(x)),Time = as.factor(date_YM),Sample = as.factor(u))
      
      # Assign row names of study design matrix
      rownames(studyDesign) <- u
      
      # Initalize random levels for the HMSC model
      rl1 <- HmscRandomLevel(units = studyDesign$Year) # Year
      rl2 <- HmscRandomLevel(sData = unique(date_YM)) # Time
      rl2 <- setPriors(rl2,nfMin = 1, nfMax = 8) # Set priors for the time random level
      rl3 <- HmscRandomLevel(sData = unique(spat), sMethod = "NNGP", nNeighbours = 20) # Location
      rl3 <- setPriors(rl3,nfMin = 1, nfMax = 8) # Set priors for the location random level
      rl4 <- HmscRandomLevel(units = studyDesign$Sample) # Sample
      
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
                ranLevels = list(Year = rl1,
                                 Time = rl2, 
                                 Loc  = rl3,
                                 Sample = rl4),
                distr = "lognormal poisson")
      
      # Compute data parameters
      m_datpar <- computeDataParameters(m)
      return(list(mod = m, datpar = m_datpar))
    }),
    
    # Run the monte-carlo markov chain for each model
    hmsc_mod = walk2(hmsc_mod,Bay, function(dat,bay){
      # Garbage collect prior to running the MCMC
      gc()
      mod <- dat$mod
      datpar <- dat$datpar
      print(paste("Starting MCMC Sampling:",bay))
      
      # Run MCMC with given runtime parameters. Parallel processing was used for full analysis; analyses run using a Linux machine (useSocket = F); if running on Windows, change to useSocket = T
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
      
      # Save each model with MCMC samples as they finish running
      save(mod, file = paste0('../Results/Bays/Bay Data/mod',bay,'.RData'))
      gc()
    })
  )
