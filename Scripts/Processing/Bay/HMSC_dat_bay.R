# Load packages
library(tidyverse)
library(SpatialEpi)

# Read and process the data
hmsc_dat <- read_csv('../Data/samples.csv') %>%
  group_by(Bay) %>%
  # Nest the biotic, sample, and geo data for each bay
  nest(sample_info = -c(Bay,starts_with("Bio_"),Longitude,Latitude),
       bio = starts_with("Bio_"),
       geo = c(Longitude, Latitude)) %>%
  ungroup() %>%
  # Convert geo lat long to cartesian coordinates
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
               Year) %>%
        mutate(temp_poly = Temperature^2, 
               sal_poly = Salinity^2) %>%
        drop_na(StartDepth:last_col())
      return(dat)
    }),
    # Isolate predictors
    predictors = map(all_dat,function(dat){
      
      hab <- select(dat,StartDepth,
                    vis,
                    Temperature,
                    temp_poly,
                    pH,
                    Salinity,
                    sal_poly,
                    DissolvedO2,
                    starts_with("Bottom_"),
                    starts_with("Veg_"),
                    starts_with("Bycatch_"),
                    starts_with("Shore_")) %>%
              select(where(~sum(.x != 0)/length(.x) > 0.005))
      
      return(hab)
    }),
    # Isolate biotic data
    bio = map(all_dat, ~select(.,starts_with("Bio_"))),
    # Isolate spatial data
    spat = map(all_dat, ~select(., Latitude, Longitude)),
    # Isolate temporal data
    Month = map(all_dat, ~select(., Month)),
    Year = map(all_dat, ~select(., Year)),
    # Isolate reference data
    Reference = map(all_dat, ~select(., Reference))) %>%
  select(Reference,Bay,predictors,bio,spat,Month,Year)

save(hmsc_dat, file = "Processing/Bay/Data/HMSCdat_bay.RData")

