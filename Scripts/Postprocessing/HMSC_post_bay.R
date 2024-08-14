#### SETUP ####
# Load packages
library(Hmsc)
library(ggtern)
library(tidyverse)
library(corrplot)
library(colorspace)
library(patchwork)
library(tidytext)
library(RColorBrewer)
library(viridis)

# Assign which bays to process (used in for loops)
bays <- c('AP','CK','TB','CH','JX','IR','TQ')

#### POST-PROCESSING OF MCMC OUTPUT ####
# Evaluate model fit and compute variance partitioning for each bay
for(i in bays) {
  print(i)
  gc()
  
  # Load model
  load(paste0('../Results/Bays/Data/mod',i,'.RData'))
  
  # Compute predicted values
  preds <- computePredictedValues(mod)
  
  # Evaluate fit of predicted values with initial biotic input
  MF <- evaluateModelFit(hM = mod, predY = preds)
  
  # Assign labels to predictors
  preds <- data.frame(predictors = colnames(mod$XData)) %>%
    mutate(predtype = factor(case_when(str_detect(predictors,'Bottom_') ~ 'Substrate',
                                       str_detect(predictors, 'Shore_') & str_detect(predictors, 'Over') ~ 'Overhanging shoreline',
                                       str_detect(predictors, "Shore_") ~ 'Non-overhanging shoreline',
                                       str_detect(predictors,'Veg_')  ~ 'Vegetation',
                                       str_detect(predictors,'Bycatch_') ~ 'Bycatch',
                                       .default = 'Environmental'),
                             levels = c('Environmental','Substrate','Vegetation','Bycatch','Overhanging shoreline', 'Non-overhanging shoreline')))
  
  # Compute variance partitioning across predictors
  group <- as.numeric(preds$predtype)
  vp <- computeVariancePartitioning(mod, group = group, groupnames = levels(preds$predtype))
  
  # Save model, variance partitioning, and model fit objects for each bay
  save(mod,MF,vp, file = paste0('../Results/Bays/Data/mod',i,'.RData'))
  rm(list = setdiff(ls(),'bays'))
}

# Initalize variance partioning/model fit list
Varpar_MF <- list()

# Combine variance partitioning and model fit objects for each bay into one list object
for(i in bays){
  load(paste0('../Results/Bays/Data/mod',i,'.RData'))
  out <- list(vp,MF)
  Varpar_MF[[i]] <- out
}

# Initialize diagnostic and effective sample size lists
ess <- list()
diagnost <- list()

# Compute MCMC diagnostics and effective sample size for each bay
for (i in bays){
  print(i)
  load(paste0('../Results/Bays/Data/mod',i,'.RData'))
  mpost <- convertToCodaObject(mod)
  ess[[i]] <- effectiveSize(mpost$Beta)
  diagnost[[i]] <- gelman.diag(mpost$Beta,multivariate = F)
}


# Initialize posterior estimates and Omega correlation matrix
PostEst  <- tibble()
OmegaCor <- list()

# Load species list
bio_list <- read_csv('../Data/tbl_corp_ref_species_list.csv')

# Assign a lookup field to match species names to NODC codes
lookup <- paste0('Bio_',bio_list$NODCCODE)
names(lookup) <- bio_list$Scientificname

# Compute posterior estimates and Omega correlation matrix for each bay
for(i in bays){
  print(i)
  load(paste0('../Results/Bays/Data/mod',i,'.RData'))
  mod <- alignPosterior(mod)
  betaPost <- getPostEstimate(mod, 'Beta')
  alphaTimePost <- getPostEstimate(mod, 'Alpha',r = 2)
  alphaLocPost  <- getPostEstimate(mod, 'Alpha', r = 3)
  
  lambdaPost_time <- crossprod(getPostEstimate(mod, 'Lambda', r = 2)$mean)
  lambdaPost_loc <- crossprod(getPostEstimate(mod, 'Lambda', r = 3)$mean)
  lambdaPost_codist <- crossprod(getPostEstimate(mod, 'Lambda', r = 4)$mean)
  
  
  # If all alpha values are zero, then the lambda matrix is the sum of the codist and time matrices (i.e., there is no temporal structure, and this element is assigned to random-level site effects)
  lambda_all <- if(all(alphaTimePost$mean == 0)) lambdaPost_codist + lambdaPost_time else lambdaPost_codist
  
  # If all alpha values are zero, then the lambda matrix is the sum of the location lambda matrix and the lambda from the previous step (i.e., there is no spatial structure, and this element is assigned to random-level site effects)
  lambda_all <- if(all(alphaLocPost$mean == 0)) lambda_all +lambdaPost_loc else lambda_all
  
  # Convert covariance matrix to correlation matrix
  OmegaCor_mod <- cov2cor(lambda_all)
  
  # Assign row and column names to correlation matrix
  rownames(OmegaCor_mod) <- mod$spNames
  colnames(OmegaCor_mod) <- mod$spNames
  
  # Assign scientific name to omega correlation matrix
  OmegaCor[[i]] <- data.frame(OmegaCor_mod) %>%
    mutate(NODCCODE = str_extract(row.names(.),'Bio_(.+)',group = 1)) %>%
    left_join(bio_list) %>%
    select(Scientificname,starts_with('Bio_')) %>%
    rename(any_of(lookup)) %>%
    column_to_rownames(var = 'Scientificname') %>%
    as.matrix(.)
    
  # Save posterior estimates for each bay as a tibble  
  out <- tibble(Bay = i, beta = betaPost,alphaTime = alphaTimePost, alphaLoc = alphaLocPost)
  
  # Append output to initialized tibble
  PostEst <- bind_rows(PostEst,out)
}

# Nest posterior estimates by bay
PostEst_tbl <- group_by(PostEst,Bay) %>%
  nest(beta = beta, alphaTime = alphaTime, alphaLoc = alphaLoc)

# Save variance patitioning, post estimate tibble, and omega correlation list as RData files
save(Varpar_MF, file = '../Results/Bays/VarPar_MF.RData')
save(PostEst_tbl, OmegaCor, file = '../Results/Bays/Data//PostEst.RData')

# Clean up workspace
rm(list = ls())

### PREPARE VISUALIZATIONS ####
# Load data required for visualizations
load('../Results/Bays/Data/VarPar_MF.RData')
load('Processing/Bay/Data/HMSCdat_bay.RData')
load('../Results/Bays/Data/PostEst.RData')

# Calculate frequence of occurrence and abundance for each species
freqoccur <- hmsc_dat %>%
  mutate(freq = map(bio, ~summarize(., across(everything(),~sum(.x>0)/length(.x))) %>%
  pivot_longer(everything(),names_to = "Species",values_to = "freq"))) %>%
  select(Bay, freq)

abun <-  hmsc_dat %>%
  mutate(abun = map(bio, ~summarize(., across(everything(),~sum(.x, na.rm = FALSE)/length(.x))) %>%
                      pivot_longer(everything(),names_to = "Species",values_to = "abun"))) %>%
  select(Bay, abun)

# Prep variance partitioning and model fits for visualization
Varpar_MF_tbl <- tibble(Bay = names(Varpar_MF),
                    R2  = map(Varpar_MF,
                              ~.[[2]]$SR2),
                    varpar = map(Varpar_MF,
                                 ~.[[1]]$vals)) %>%
  left_join(freqoccur) %>%
  left_join(abun) %>%
  left_join(PostEst_tbl) %>%
  mutate(varpar = pmap(list(freq,abun,varpar, R2, alphaTime,alphaLoc), \(freq,abun,vp, R2,aT,aL) {
                    vp <- t(vp)
                    R2 <- data.frame(R2, row.names = rownames(vp))
                    bind_cols(R2 = R2,vp) %>%
                      rownames_to_column(var = "Species") %>%
                      mutate(env = Environmental,
                             ovshore = `Overhanging shoreline`,
                             shore = `Non-overhanging shoreline`,
                             bycatch = Bycatch,
                             veg = Vegetation,
                             sub = Substrate,
                             spa = `Random: Loc`, 
                             codist = `Random: Sample`, 
                             time = `Random: Time`,
                             year = `Random: Year`,
                             hab = env + ovshore + shore + bycatch + veg + sub,
                             # If all time alpha values are zero, then the dispersal element includes only the year effect, otherwise it is the sum of the year and time effects
                             disp = if(all(aT$alphaTime$mean == 0)) year else time + year,
                             # If all location alpha values are zero, then the dispersal element includes only the previous dispersal effects, otherwise it is the sum of the dispersal and location effects
                             disp = if(all(aL$alphaLoc$mean == 0)) disp else disp + spa,
                             # If all time alpha values are zero, then the interaction element includes the time effect, otherwise it is only the sample-level effect
                             inter = if(all(aT$alphaTime$mean == 0)) codist + time else codist,
                             # If all location alpha values are zero, then the interaction element includes the location effect, otherwise it is only the sample-level effect and possible time effect
                             inter = if(all(aL$alphaLoc$mean == 0)) inter + spa else inter
                             ) %>%
                      
                      left_join(freq) %>%
                      left_join(abun) %>%
                      select(Species,
                             R2,
                             env,
                             ovshore,
                             shore,
                             bycatch,
                             veg,
                             sub,
                             spa,
                             codist,
                             time,
                             year,
                             freq,
                             abun,
                             hab,
                             disp,
                             inter)
                    })) %>%
  
  select(Bay, varpar) %>%
  # Unnest the tibble so that each bay has its own set of rows with the variance partitioning output values as columns
  unnest(cols = c(varpar)) %>%
  select(Bay,Species,R2,env,ovshore,shore,bycatch,veg,sub,spa,codist,time,year, freq, abun,hab,disp,inter) %>%
  drop_na() %>%
  # Assign bay as factor and order levels north -> south, west -> east
  mutate(Bay = factor(Bay, levels = c('AP','CK','TB','CH','JX','IR','TQ')))

# Save the variance partition/model fit tibble as 
saveRDS(Varpar_MF_tbl, file = '../Results/Bays/Data/VarPar_MF_tbl.RDS')

### MODEL FIT TRENDS ####
# Make plot for frequency of occurrence vs model fit
VarPar_trends_freqR2 <- Varpar_MF_tbl %>%
  mutate(Bay = factor(Bay,
                      levels = c('AP','CK','TB','CH','JX','IR','TQ'),
                      labels = c('Apalachicola Bay',
                                 'Cedar Key',
                                 'Tampa Bay',
                                 'Charlotte Harbor',
                                 'Northeast Florida',
                                 'Northern Indian River',
                                 'Southern Indian River'))) %>%
  ggplot(aes(x = freq, y = R2, color = Bay, fill = Bay)) +
  geom_point() +
  geom_smooth(method = 'lm', alpha = 0.2) +
  scale_fill_viridis_d(option = 'plasma') +
  scale_color_viridis_d(option = 'plasma') +
  labs(x = 'Frequency of occurrence',
       y = 'Model fit') +
  theme_bw() +
  
  theme(legend.title = element_blank(),
        )
# Make plot for showing log abundance vs model fit
VarPar_trends_abunR2 <- Varpar_MF_tbl %>%
  mutate(Bay = factor(Bay,
                      levels = c('AP','CK','TB','CH','JX','IR','TQ'),
                      labels = c('Apalachicola Bay',
                                 'Cedar Key',
                                 'Tampa Bay',
                                 'Charlotte Harbor',
                                 'Northeast Florida',
                                 'Northern Indian River',
                                 'Southern Indian River'))) %>%
  ggplot(aes(x = log(abun), y = R2, color = Bay, fill = Bay)) +
  geom_point() +
  geom_smooth(method = 'lm', alpha = 0.2) +
  scale_fill_viridis_d(option = 'plasma') +
  scale_color_viridis_d(option = 'plasma') +
  labs(x = 'Log average abundance',
       y = 'Model fit') +
  theme_bw() +
  
  theme(legend.title = element_blank(),
  )

# Combine plots into two-panel patchwork
(VarPar_trends_R2 <- VarPar_trends_freqR2 + VarPar_trends_abunR2 + plot_layout(guides = 'collect'))

# Save plot
ggsave('../Results/Bays/Figures/freq_R2.svg', VarPar_trends_R2, width = 12, height = 6)


#### BARPLOT ####
# load species list
bio_list <- read_csv('../Data/tbl_corp_ref_species_list.csv')

# Compute variance partitioning averages for each variable type, scaled by total model fit
VarPar_avgs <- Varpar_MF_tbl %>%
  group_by(Bay) %>%
  mutate(across(env:inter, ~.x*R2)) %>%
  summarise(across(env:inter, mean))


# Prepare variance partitioning tibble for barplots
VarPar_bar_dat <- Varpar_MF_tbl %>% 
  pivot_longer(cols = c(env:year,hab:inter), names_to = 'varpart', values_to = 'varpct') %>%
  # Scale variance partitioning by model fit, assign each variable type to dispersal, interactions, or habitat
  mutate(varpct = varpct*R2,
         vartype = case_when(varpart %in% c('env','ovshore','shore','bycatch','veg','sub','bay') ~ 'hab',
                             varpart == 'time' ~ 'disp',
                             varpart == 'year' ~ 'disp',
                             varpart == 'spa'  ~ NA,
                             varpart == 'inter' ~ 'inter')) %>%
  # Drop any NA values (in this case, "spatial" variation that was looped in with codistribution due to lack of spatial structure)
  drop_na() %>%
  
  # Match NODC codes to species names, order factors for proper plotting
  mutate(NODCCODE = str_sub(Species,start = 5),
         vartype = factor(vartype, levels = c('hab','disp','inter'))) %>%
  left_join(bio_list) %>%
  arrange(desc(freq), desc(vartype)) %>%
  
  # Dynamically assign labels to variance partitioning components for legend, including the average variance partitioning values
  mutate(varpart_labs = case_when(varpart == 'inter' ~ paste0('Interactions (',round(100*min(VarPar_avgs$inter), 1),'-',round(100*max(VarPar_avgs$inter), 1),'%)'),
                                  varpart == 'time' ~ paste0('Temporal structure (',round(100*min(VarPar_avgs$time), 1),'-',round(100*max(VarPar_avgs$time), 1),'%)'),
                                  varpart == 'year' ~ paste0('Annual variability (',round(100*min(VarPar_avgs$year), 1),'-',round(100*max(VarPar_avgs$year), 1),'%)'),
                                  varpart == 'env' ~ paste0('Environment (',round(100*min(VarPar_avgs$env), 1),'-',round(100*max(VarPar_avgs$env), 1),'%)'),
                                  varpart == 'ovshore' ~ paste0('Overhanging shoreline (',round(100*min(VarPar_avgs$ovshore), 1),'-',round(100*max(VarPar_avgs$ovshore), 1),'%)'),
                                  varpart == 'shore' ~ paste0('Non-overhanging shoreline (',round(100*min(VarPar_avgs$shore), 1),'-',round(100*max(VarPar_avgs$shore), 1),'%)'),
                                  varpart == 'veg' ~ paste0('Vegetation (',round(100*min(VarPar_avgs$veg), 1),'-',round(100*max(VarPar_avgs$veg), 1),'%)'),
                                  varpart == 'sub' ~ paste0('Substrate (',round(100*min(VarPar_avgs$sub), 1),'-',round(100*max(VarPar_avgs$sub), 1),'%)'),
                                  varpart == 'bycatch' ~ paste0('Bycatch (',round(100*min(VarPar_avgs$bycatch), 1),'-',round(100*max(VarPar_avgs$bycatch), 1),'%)')),
         varpart_labs = factor(varpart_labs, levels = unique(varpart_labs))) %>%
  drop_na()

# Assign colors to variance partitioning components
n_hab <- sum(unique(VarPar_bar_dat$varpart) %in% c('env','ovshore','shore','bycatch','veg','sub'))
n_disp <- sum(unique(VarPar_bar_dat$varpart) %in% c('time','year'))
n_inter <- sum(unique(VarPar_bar_dat$varpart) %in% c('inter'))

color_palette <- c(rep("purple4", n_inter),
                   rep("darkblue", n_disp),
                   rep("darkgreen",n_hab))

fill_palette <- c('darkmagenta',
                  rev(brewer.pal(n = max(3,2*n_disp), 'Blues'))[1:n_disp],
                  rev(brewer.pal(n = 9, 'Greens'))[1:n_hab])

max_freq <- max(VarPar_bar_dat$freq)

# Create bar plots for each bay, plot in a patchwork framework
VarPar_bar <- VarPar_bar_dat %>%
                 mutate(Bay = factor(Bay, levels = c('AP','CK','TB','CH','JX','IR','TQ'),
                                     labels = c('Apalachicola Bay',
                                                'Cedar Key',
                                                'Tampa Bay',
                                                'Charlotte Harbor',
                                                'Northeast Florida',
                                                'Northern Indian River',
                                                'Southern Indian River'))) %>%
                 arrange(Bay) %>%
                 group_by(Bay) %>%
                 nest() %>%
                 mutate(barplot = map2(Bay,data, \(Bay,data){
                   bar <- data %>%
                     mutate(Species = factor(Scientificname, levels = unique(Scientificname))) %>%
                     arrange(desc(freq)) %>%
                     ggplot() +
                       geom_bar(aes(x = Species, y = varpct, color = varpart_labs, fill = varpart_labs), position = "stack", stat = "identity") +
                       scale_color_manual(values = color_palette, guide = guide_legend(title = expression(underline('Variance component')),nrow = 3 )) +
                       scale_fill_manual(values = fill_palette, guide = guide_legend(title = expression(underline('Variance component')),nrow = 3 )) +
                       scale_y_continuous(limits = c(0,NA),expand = c(0,0,0.1,0)) +
                       labs(title = Bay, y = '%Var') +
                       theme_minimal() +
                       theme(
                         axis.text.x = element_blank(),
                         panel.grid = element_blank(),
                         axis.title.x = element_blank(),
                         axis.title.y = element_text(size = 18),
                         axis.text.y = element_text(size = 16)
                       )
                     freq <- data %>%
                       mutate(Species = factor(Scientificname, levels = unique(Scientificname))) %>%
                       arrange(desc(freq)) %>%
                       ggplot() +
                        geom_point(aes(x = Species, y = 0, color = freq), size = 4) +
                        scale_color_viridis_c(limits = c(0,round(max_freq,1)), breaks = c(0,0.2,0.4,0.6,0.8), guide = guide_legend(reverse = T,title = expression(underline('Occurrence')),override.aes = list(size = 10))) +
                        theme_minimal() +
                        theme(axis.title.x = element_blank(),
                              axis.title.y = element_blank(),
                              axis.text.x = element_text(size = 16, angle = -45, hjust = 0, vjust = 1),
                              axis.text.y = element_blank(),
                              panel.grid = element_blank(),
                              plot.margin = margin(0,0,0,0))
                     bar / freq +
                     plot_layout(heights = c(10,1.5))
                   }))

(plot_all <- patchwork::wrap_plots(VarPar_bar$barplot, 
                                   ncol = 1,
                                   guides = 'collect',
                                   byrow = F) & 
    theme(legend.box = 'vertical',
          legend.position = 'bottom',
          legend.direction = 'horizontal',
          legend.key.size = unit(0.5, 'in'),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 16),
          title = element_text(size = 18),
          plot.margin = unit(c(0,0.75,0,0), 'in')))

# Save bar plot
ggsave(plot_all,file = '../Results/Bays/Figures/varparBar_all.svg', height = 5, width = 3, scale = 6)

#### TERNARY PLOTS ####

# Prepare ternary plots as faceted plots for each bay
VarPar_tern_all <- Varpar_MF_tbl %>% 
  mutate(Bay = factor(Bay,
                      levels = c('AP','CK','TB','CH','JX','IR','TQ'),
                      labels = c('Apalachicola Bay',
                                 'Cedar Key',
                                 'Tampa Bay',
                                 'Charlotte Harbor',
                                 'Northeast Florida',
                                 'Northern Indian River',
                                 'Southern Indian River'))) %>%
  ggtern(aes(x = hab, z = disp, y = inter)) +
  stat_density_tern(
    geom='polygon',
    aes(fill=after_stat(level)),
    linewidth = 0.001,
    breaks = seq(1,10,2)) +
  geom_point(aes(color = freq, size = R2), alpha = 0.8) +
  
  scale_T_continuous(limits=c(0,1.0),
                     breaks=seq(0,1,by=0.2),
                     labels=seq(0,1,by=0.2)) +
  scale_L_continuous(limits=c(0.0,1),
                     breaks=seq(0,1,by=0.2),
                     labels=seq(0,1,by=0.2)) +
  scale_R_continuous(limits=c(0.0,1.0),
                     breaks=seq(0,1,by=0.2),
                     labels=seq(0,1,by=0.2)) +
  labs(title = "",
       x = "H",
       xarrow = "Habitat",
       y = "I",
       yarrow = "Interactions",
       z = "D", 
       zarrow = "Dispersal") +
  
  theme_bw() +
  scale_fill_viridis_c(option = 'C', alpha = 0.2, guide = guide_legend(reverse = T, title = expression(underline('Kernel density')))) +
  scale_color_viridis_c(limits = c(0,0.8), breaks = c(0,0.2,0.4,0.6,0.8), guide = guide_legend(reverse = T,title = expression(underline('Occurrence')), order = 2, override.aes = list(size = 5))) +
  scale_size_area(breaks = seq(0,1,0.2), limits = c(0,1), guide = guide_legend(reverse = T, title = expression(underline(R^2)), order = 1)) +
  theme(panel.grid = element_line(color = "darkgrey", ),
        tern.axis.line = element_line(color = 'darkgrey'),
        axis.text = element_text(size = 8),
        strip.background = element_blank(),
        strip.text = element_text(size = 12,face = 'bold')) +
  facet_wrap(vars(Bay), ncol = 2, dir = 'v')

# Save ternary plots
ggsave('../Results/Bays/Figures/varpar_tern.svg',VarPar_tern_all, height = 9, width = 6)

#### HABITAT ASSOCIATIONS ####

# Load data
load('Postprocessing/Bay/Data/HMSCdat_bay.RData')
bio_list <- read_csv('../Data/tbl_corp_ref_species_list.csv')

# Prepare posterior estimates for visualization
modPreds <- hmsc_dat %>%
  select(Bay,predictors) %>%
  mutate(predictors = map(predictors,  ~c('Int',colnames(.))))

# Extract beta estimates for each bay
PostEst_Beta <- PostEst_tbl %>%
  select(Bay,beta) %>%
  left_join(modPreds) %>%
  mutate(beta = map2(beta,predictors, \(beta,preds){
    est <- data.frame(beta$beta$mean) %>%
      mutate(predictors = preds) %>%
      pivot_longer(cols = starts_with('Bio_'),names_to = 'spp', values_to = 'mean')
    support <- data.frame(beta$beta$support) %>%
      mutate(predictors = preds) %>%
      pivot_longer(cols = starts_with('Bio_'),names_to = 'spp',values_to = 'support')
    support_neg <- data.frame(beta$beta$supportNeg) %>%
      mutate(predictors = preds) %>%
      pivot_longer(cols = starts_with('Bio_'),names_to = 'spp',values_to = 'support_neg')
    left_join(est,left_join(support,support_neg)) %>%
      # Filter out all predictors with support less than 0.95 (positive or negative) and intercepts
      filter((support > 0.95 | support_neg > 0.95),predictors != 'Int') %>%
      # Assign negative values to support for negative coefficients
      mutate(support = ifelse(support > support_neg, support, -support_neg)) %>%
      select(-support_neg)
  })) %>%
  select(Bay,beta) %>%
  # Unnest bay-specific columns to produce a full table of beta estimates across bays
  unnest(cols = 'beta') %>%
  right_join(Varpar_MF_tbl, by = join_by(Bay, spp == Species)) %>%
  mutate(Bay = factor(Bay,
                      levels = c('AP','CK','TB','CH','JX','IR','TQ'),
                      labels = c('Apalachicola Bay',
                                 'Cedar Key',
                                 'Tampa Bay',
                                 'Charlotte Harbor',
                                 'Northeast Florida',
                                 'Northern Indian River',
                                 'Southern Indian River')),
         NODCCODE = str_sub(spp,start = 5)) %>%
  left_join(bio_list) %>%
  arrange(Bay,freq) %>%
  
  # Assign predictor names based on coded column name
  mutate(Commonname = factor(Commonname, levels = unique(Commonname)),
         predtype = factor(case_when(str_detect(predictors,'Bottom_') ~ 'Bottom Type',
                                     str_detect(predictors, 'Shore_') ~ 'Shoreline',
                                     str_detect(predictors,'Veg_')  ~ 'Vegetation',
                                     str_detect(predictors,'Bycatch_') ~ 'Bycatch',
                                     is.na(predictors) ~ '',
                                     .default = 'Environmental'),
                           levels = c('','Environmental','Bottom Type','Vegetation','Bycatch','Shoreline','Estuary')),
         shoretype = factor(case_when(str_detect(predictors, 'Shore') & str_detect(predictors, 'Over') ~ 'Overhanging',
                                      str_detect(predictors, 'Shore') ~ 'Non-Overhanging',
                                      .default = ""),
                            levels = c('','Overhanging','Non-Overhanging')),
         predictor_labels = case_when(predictors == 'Bottom_D' ~ 'Detritus',
                                      predictors == 'Bottom_E' ~ 'Peat',
                                      predictors == 'Bottom_H' ~ 'Shell',
                                      predictors == 'Bottom_M' ~ 'Mud',
                                      predictors == 'Bottom_O' ~ 'Oyster',
                                      predictors == 'Bottom_R' ~ 'Rocks',
                                      predictors == 'Bottom_S' ~ 'Sand',
                                      predictors == 'Bycatch_AD' ~ 'Drift algae',
                                      predictors == 'Bycatch_AM' ~ 'Mixed algae',
                                      predictors == 'Bycatch_HA' ~ 'Halodule',
                                      predictors == 'Bycatch_MU' ~ 'Mud',
                                      predictors == 'Bycatch_GM' ~ 'Mixed seagrass',
                                      predictors == 'Bycatch_ST' ~ 'Sticks and branches',
                                      predictors == 'Bycatch_OY' ~ 'Oysters',
                                      predictors == 'Bycatch_AG' ~ 'Filamentous algae',
                                      predictors == 'Bycatch_TH' ~ 'Thalassia',
                                      predictors == 'Bycatch_CT' ~ 'Ctenophores',
                                      predictors == 'Bycatch_UL' ~ 'Ulva',
                                      predictors == 'DissolvedO2' ~ 'Dissolved oxygen',
                                      predictors == 'pH' ~ 'pH',
                                      predictors == 'Salinity' ~ 'Salinity',
                                      predictors == 'sal_poly' ~ 'Salinity (poly)',
                                      predictors == 'StartDepth' ~ 'Depth',
                                      predictors == 'Temperature' ~ 'Temperature',
                                      predictors == 'temp_poly' ~ 'Temperature (poly)',
                                      predictors == 'vis' ~ 'Water clarity',
                                      predictors %in% c('Shore_APOver', 'Shore_AP') ~ 'Australian pine',
                                      predictors %in% c('Shore_BMOver','Shore_BM') ~ 'Black mangrove',
                                      predictors %in% c('Shore_BPOver', 'Shore_BP') ~ 'Brazilian pepper',
                                      predictors %in% c('Shore_BUOver', 'Shore_BU') ~ 'Bulrush',
                                      predictors %in% c('Shore_BWOver', 'Shore_BW') ~ 'Buttonwood',
                                      predictors %in% c('Shore_CS') ~ 'Cattails',
                                      predictors %in% c('Shore_DS') ~ 'Dead oysters',
                                      predictors %in% c('Shore_JUOver', 'Shore_JU') ~ 'Juncus',
                                      predictors %in% c('Shore_MAOver', 'Shore_MA') ~ 'Mixed mangroves',
                                      predictors %in% c('Shore_MGOver', 'Shore_MG') ~ 'Marsh grasses',
                                      predictors %in% c('Shore_NO') ~ 'Sand',
                                      predictors %in% c('Shore_OSOver') ~ 'Shrub/tree',
                                      predictors %in% c('Shore_OY') ~ 'Oysters',
                                      predictors %in% c('Shore_RMOver', 'Shore_RM') ~ 'Red mangrove',
                                      predictors %in% c('Shore_RO') ~ 'Rocks',
                                      predictors %in% c('Shore_RR') ~ 'Rip rap',
                                      predictors %in% c('Shore_SGOver', 'Shore_SG') ~ 'Sea grapes',
                                      predictors %in% c('Shore_SN') ~ 'Spartina',
                                      predictors %in% c('Shore_SW') ~ 'Seawall',
                                      predictors %in% c('Shore_TC') ~ 'Cypress',
                                      predictors %in% c('Shore_TDOver', 'Shore_TD') ~ 'Dead tree',
                                      predictors %in% c('Shore_TVOver', 'Shore_TV') ~ 'Terrestrial vegetation',
                                      predictors %in% c('Shore_WMOver', 'Shore_WM') ~ 'White mangrove',
                                      predictors %in% c('Shore_WR') ~ 'Wrack',
                                      predictors == 'Veg_CA' ~ 'Caulerpa',
                                      predictors == 'Veg_GM' ~ 'Mixed seagrass',
                                      predictors == 'Veg_HA' ~ 'Halodule',
                                      predictors == 'Veg_NO' ~ 'None',
                                      predictors == 'Veg_SY' ~ 'Syringodium',
                                      predictors == 'Veg_TH' ~ 'Thalassia',
                                      predictors == 'Veg_VA' ~ 'Valisneria',
                                      predictors == 'Veg_CA' ~ 'Caulerpa',
                                      predictors == 'Veg_RU' ~ 'Ruppia',
                                      .default = 'No associations')) %>%
  filter(predictor_labels != 'No associations') %>%
  arrange(desc(freq))

# Separate out frequency of occurrence data for plotting as seperate patch
PostEst_Beta_points <- distinct(PostEst_Beta,Bay,Scientificname, freq)

# Make plot for eastern Gulf of Mexico bays
eGOM_plot <- (PostEst_Beta_points %>%
                filter(Bay %in% c('Apalachicola Bay','Cedar Key','Tampa Bay','Charlotte Harbor')) %>%
               # Plot points for each species, colored by frequency of occurrence
               ggplot() +
               geom_point(aes(x = 0,y = reorder_within(Scientificname,freq,Bay), color = freq), size = 3.75) +
               facet_grid(rows = vars(Bay), scales = 'free', space = 'free', switch = 'y') +
               scale_color_viridis_c(limits = c(0,.8), breaks = c(0,0.2,0.4,0.6,0.8,1), guide = guide_legend(reverse = T,title = expression(underline('Occurrence')), override.aes = list(size = 10))) +
               theme_minimal() +
               theme(axis.title = element_blank(),
                     axis.text.y = element_text(size = 12),
                     axis.text.x = element_blank(),
                     strip.placement = 'outside',
                     strip.text.y = element_text(size = 32),
                     panel.grid = element_blank(),
                     plot.margin = margin(0,0,0,0)) +
               scale_y_reordered()
) +
  # Plot habitat associations as "heat plot" of positive and negative associations, faceted by type of predictor and bay
  (PostEst_Beta %>% 
     filter(Bay %in% c('Apalachicola Bay','Cedar Key','Tampa Bay','Charlotte Harbor')) %>%
     group_by(predictors) %>%
     mutate(nassoc = n()) %>%
     ungroup() %>%
     arrange(desc(nassoc)) %>%
     mutate(predictor_labels = factor(predictor_labels, levels = unique(predictor_labels))) %>%
     ggplot(aes(x = predictor_labels, y = reorder_within(Scientificname,freq,Bay) , fill = support)) +
     geom_tile() +
     scale_fill_gradientn(colors = brewer.pal(5, "RdBu"), guide = 'none') +
     scale_x_discrete(position = 'bottom') +
     facet_grid(rows = vars(Bay), cols = vars(shoretype,predtype), scales = 'free', space = 'free', switch = 'y') +
     theme_minimal() +
     theme(axis.title = element_blank(),
           axis.text.y = element_blank(),
           axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 32),
           strip.placement = 'outside',
           strip.text.x = element_text(size = 32),
           strip.text.y = element_blank(),
           panel.grid.major.x = element_blank(),
           panel.grid.major.y = element_line(color = 'gray50'),
           panel.spacing.x = unit(2, "lines"),
           plot.margin = margin(0, 2.25, 0, 0, "in"))) +
  scale_y_reordered() +
  plot_layout(widths = c(1,40), guides = 'collect') & theme(legend.direction = 'horizontal',
                                                            legend.position = 'bottom',
                                                            legend.title = element_text(size = 30),
                                                            legend.text = element_text(size = 20))

# Save plots as svg and pdf
ggsave('../Results/Bays/Figures/Abiotic_eGOM.svg',plot = eGOM_plot,width = 16,height = 16,units = 'in',scale = 2)
ggsave('../Results/Bays/Figures/Abiotic_eGOM.pdf',plot = eGOM_plot,width = 16,height = 16,units = 'in',scale = 2)

# Make plot for atlantic bays
atl_plot <- (PostEst_Beta_points %>%
               filter(Bay %in% c('Northeast Florida','Northern Indian River','Southern Indian River')) %>%
               # Plot points for each species, colored by frequency of occurrence
               ggplot() +
               geom_point(aes(x = 0,y = reorder_within(Scientificname,freq,Bay), color = freq), size = 3.75) +
               facet_grid(rows = vars(Bay), scales = 'free', space = 'free', switch = 'y') +
               scale_color_viridis_c(limits = c(0,.8), breaks = c(0,0.2,0.4,0.6,0.8,1), guide = guide_legend(reverse = T,title = expression(underline('Occurrence')), override.aes = list(size = 10))) +
               theme_minimal() +
               theme(axis.title = element_blank(),
                     axis.text.y = element_text(size = 12),
                     axis.text.x = element_blank(),
                     strip.placement = 'outside',
                     strip.text.y = element_text(size = 32),
                     panel.grid = element_blank(),
                     plot.margin = margin(0,0,0,0)) +
               scale_y_reordered()
               ) +
  # Plot habitat associations as "heat plot" of positive and negative associations, faceted by type of predictor and bay
(PostEst_Beta %>% 
  filter(Bay %in% c('Northeast Florida','Northern Indian River','Southern Indian River')) %>%
   group_by(predictors) %>%
   mutate(nassoc = n()) %>%
   ungroup() %>%
   arrange(desc(nassoc)) %>%
   mutate(predictor_labels = factor(predictor_labels, levels = unique(predictor_labels))) %>%
   ggplot(aes(x = predictor_labels, y = reorder_within(Scientificname,freq,Bay) , fill = support)) +
   geom_tile() +
   scale_fill_gradientn(colors = brewer.pal(5, "RdBu"), guide = 'none') +
   scale_x_discrete(position = 'bottom') +
   facet_grid(rows = vars(Bay), cols = vars(shoretype,predtype), scales = 'free', space = 'free', switch = 'y') +
   theme_minimal() +
   theme(axis.title = element_blank(),
         axis.text.y = element_blank(),
         axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 32),
         strip.placement = 'outside',
         strip.text.x = element_text(size = 32),
         strip.text.y = element_blank(),
         panel.grid.major.x = element_blank(),
         panel.grid.major.y = element_line(color = 'gray50'),
         panel.spacing.x = unit(2, "lines"),
         plot.margin = margin(0, 2.25, 0, 0, "in"))) +
  scale_y_reordered() +
  plot_layout(widths = c(1,40), guides = 'collect') & theme(legend.direction = 'horizontal',
                                                            legend.position = 'bottom',
                                                            legend.title = element_text(size = 30),
                                                            legend.text = element_text(size = 20))

# Save plots as svg and pdf
ggsave('../Results/Bays/Figures/Abiotic_atl.svg',plot = atl_plot,width = 16,height = 13,units = 'in',scale = 2)
ggsave('../Results/Bays/Figures/Abiotic_atl.pdf',plot = atl_plot,width = 16,height = 13,units = 'in',scale = 2)


# Make separate plots for each bays
ind_plots <- PostEst_Beta %>%
  group_by(Bay,predictors) %>%
  mutate(nassoc = n()) %>%
  ungroup() %>%
  group_by(Bay) %>%
  nest() %>%
  ungroup() %>%
  mutate(plots = map2(Bay,data, function(bay,dat) {
    # prep data 
    plot_dat <- dat %>%
      arrange(desc(nassoc)) %>%
      mutate(predictor_labels = factor(predictor_labels, levels = unique(predictor_labels))) %>%
      arrange((freq)) %>%
      mutate(Scientificname = factor(Scientificname, levels = unique(Scientificname)))
    # Make heat tile plot
   tile_plot <- ggplot(plot_dat, aes(x = predictor_labels, y = Scientificname , fill = support)) +
      geom_tile() +
      scale_fill_gradientn(colors = brewer.pal(5, "RdBu"), guide = 'none') +
      scale_x_discrete(position = 'bottom') +
      theme_minimal() +
      theme(axis.title = element_blank(),
            axis.text.y = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16),
            strip.placement = 'outside',
            strip.text.x = element_text(size = 16, face = 'bold'),
            panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_line(color = 'gray50'),
            panel.spacing.x = unit(1, "lines"),
            plot.margin = margin(0, 0, 0, 0, "in")) +
      facet_grid(cols = vars(shoretype,predtype), scales = 'free', space = 'free')
    # Make frequency plot
   freq_plot <- ggplot(plot_dat, aes(x = 0, y = Scientificname, color = freq)) +
      geom_point(size = 3.75) +
      scale_color_viridis_c(limits = c(0,.8), breaks = c(0,0.2,0.4,0.6,0.8,1), guide = guide_legend(reverse = T,title = expression(underline('Occurrence')), override.aes = list(size = 10)))+
      theme_minimal() +
      theme(axis.title = element_blank(),
            axis.text.y = element_text(size = 16),
            axis.text.x = element_blank(),
            panel.grid = element_blank(),
            plot.margin = margin(0,0,0,0))
   
    # Combine plots in patchwork
   full_plot <- freq_plot + tile_plot + 
     plot_layout(widths = c(1,100), guides = 'collect') & theme(legend.direction = 'horizontal',
                                                               legend.position = 'bottom',
                                                               legend.title = element_text(size = 20),
                                                               legend.text = element_text(size = 18))
  }))

# Save plots as svg
walk2(ind_plots$plots, ind_plots$Bay,~ggsave(paste0('../Results/Bays/Figures/Abiotic_',.y,'.svg'),plot = .x,width = 2.5,height = 1,units = 'in',scale = 12))

#### SPECIES ASSOCIATIONS ####
# Load data
Varpar_MF_tbl <- readRDS(file = '../Results/Bays/Data/VarPar_MF_tbl.RDS')

# Prepare interactions data for plotting
intxs <- Varpar_MF_tbl %>%
  mutate(inter = inter*R2) %>%
  select(Bay,inter, freq) %>%
  group_by(Bay) %>% 
  nest() %>%
  arrange(Bay) %>%
  .$data

# Load posterior estimates (which includes omega correlation list)
load('../Results/Bays/Data/PostEst.RData')

# Prepare and plot omega correlations
OmegaCor_plot <-map2(OmegaCor, intxs, function(x,y) {
  # add species names as row names to interactions vector
  y <- y %>% mutate(spp = row.names(x)) %>% column_to_rownames('spp')
  # Create diagonal matrix from interactions vector
  y <- diag(y$inter)
  
  # Assign dimension names for plotting
  dimnames(y) <- dimnames(x)
  
  # Order species by hierarchical clustering of pairwise correlation strength
  ord <- corrMatOrder(x, order = 'hclust')
  
  # Assign ordering to interactions matrix
  y <- y[ord,ord]
  y[y == 0] <- NA
  
  # Create correlation plots, stacking correlations and interactions
  fig <- {corrplot(y,diag = T,add = F,is.corr = F, method = 'color',cl.pos = 'b',col = COL1('YlGn',200),tl.col = 'black',cl.cex = 1)
          corrplot(x,diag = F, order = 'hclust',add = T,tl.pos = 'n',cl.cex = 1)
          recordPlot()
  }
  
}) %>%
  # Save plots as svg and pdf
  iwalk(~ggsave(filename = paste0('../Results/Bays/Figures/SpeciesAssociations_',.y,'.svg'), plot = replayPlot(.x),width = 7, height = 6)) %>%
  iwalk(~ggsave(filename = paste0('../Results/Bays/Figures/SpeciesAssociations_',.y,'.pdf'), plot = replayPlot(.x),width = 7, height = 6))
