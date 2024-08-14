# Load packages
library(tidyverse)
library(SpatialEpi)

# Read and process data
hmsc_dat <- read_csv('../Data/samples.csv') %>%
  # Nest the biotic, sample, and spatial data into df's within tibble
  nest(sample_info = -c(Bay,starts_with("Bio_"),Longitude,Latitude),
       bio = starts_with("Bio_"),
       geo = c(Longitude, Latitude)) %>%
  # Convert spatial data to cartesian coordinates
  mutate(geo      = map(geo,~latlong2grid(.))) %>%
  mutate(
    # Prepare the data for HMSC
    all_dat = pmap(list(bio,sample_info,geo), function(bio, samp, spat){
      dat <- cbind(bio, samp, spat) %>%
        select(Reference,
               starts_with("Bio_"),
               StartDepth,
               vis,
               Temperature,
               pH,
               Salinity,
               DissolvedO2,
               starts_with("Bottom_"),
               starts_with("Veg_"),
               starts_with("Bycatch_"),
               starts_with("Shore_"),
               Latitude,
               Longitude,
               Month,
               Year)
    }),
    # Isolate predictors
    predictors = map(all_dat,function(dat){
      
      hab <- select(dat,StartDepth,
                    vis,
                    Temperature,
                    pH,
                    Salinity,
                    DissolvedO2,
                    starts_with("Bottom_"),
                    starts_with("Veg_"),
                    starts_with("Bycatch_"),
                    starts_with("Shore_"))
    }),
    # Isolate biotic data
    bio = map(all_dat, ~select(.,starts_with("Bio_"))),
    # Isolate spatial data
    spat = map(all_dat, ~select(., Latitude, Longitude)),
    # Isolate Month and Year
    Month = map(all_dat, ~select(., Month)),
    Year = map(all_dat, ~select(., Year)),
    # Isolate reference
    Reference = map(all_dat, ~select(., Reference))) %>%
  select(Reference,Bay,predictors,bio,spat,Month,Year) %>%
  unnest(cols = c(Reference,predictors, bio, spat, Month, Year)) %>%
  nest(predictors = c(StartDepth:DissolvedO2, 
                      starts_with("Bottom_"),
                      starts_with("Veg_"),
                      starts_with("Bycatch_"),
                      starts_with("Shore_")),
       bio = starts_with("Bio_"),
       spat = c(Latitude,Longitude),
       Month = Month,
       Year = Year,
       Bay = Bay,
       Reference = Reference)

save(hmsc_dat, file = "Processing/Full/Data/HMSCdat_full.RData")

