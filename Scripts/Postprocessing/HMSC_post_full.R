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

# Load MCMC model
hmsc_mod <- readRDS('../Results/All/Model/HMSC_full.rds')

#### POST PROCESSING OF MCMC OUTPUT ####
# Compute predicted values
preds <- computePredictedValues(hmsc_mod)

# Evaluate model fit
MF <- evaluateModelFit(hM = hmsc_mod, predY = preds)

# Pre-treat data for variance partitioning
preds <- data.frame(predictors = colnames(hmsc_mod$XData)) %>%
  mutate(predtype = factor(case_when(str_detect(predictors,'Bottom_') ~ 'Substrate',
                                     str_detect(predictors, 'Shore_') & str_detect(predictors, 'Over') ~ 'Overhanging shoreline',
                                     str_detect(predictors, "Shore_") ~ 'Non-overhanging shoreline',
                                     str_detect(predictors,'Veg_')  ~ 'Vegetation',
                                     str_detect(predictors,'Bycatch_') ~ 'Bycatch',
                                     .default = 'Environmental'),
                           levels = c('Environmental','Substrate','Vegetation','Bycatch','Overhanging shoreline', 'Non-overhanging shoreline')))

group <- as.numeric(preds$predtype)

# Compute variance partitioning
vp <- computeVariancePartitioning(hmsc_mod, group = group, groupnames = levels(preds$predtype))

# Save variance partitioning results
Varpar_MF <- list(vp,MF)
save(Varpar_MF, file = 'VarPar_MF_ALL.RData')

# Get posterior estimates
betaPost        <- getPostEstimate(hmsc_mod, 'Beta')
alphaLocBayPost <- getPostEstimate(hmsc_mod, 'Alpha')
alphaTimePost   <- getPostEstimate(hmsc_mod, 'Alpha')
alphaLocPost    <- getPostEstimate(hmsc_mod, 'Alpha')
etaBayPost      <-  getPostEstimate(hmsc_mod, 'Eta')
lambdaBayPost   <-  getPostEstimate(hmsc_mod, 'Lambda')
PostEst         <- list(beta = betaPost,alphaTime = alphaTimePost, alphaLoc = alphaLocPost, alphaLocBay = alphaLocBayPost,  eta = etaBayPost, lambda = lambdaBayPost)

modPreds <- hmsc_mod$covNames

bays <- levels(hmsc_mod[["rL"]][["Bay"]][["pi"]])

# Save posterior estimates
save(PostEst,modPreds, bays, file = 'PostEst_ALL.RData')

# Get diagnostics and effective sample size
mpost <- convertToCodaObject(hmsc_mod, start = start)
essBeta <- effectiveSize(mpost$Beta)
diagnostBeta <- gelman.diag(mpost$Beta,multivariate = F)
save(mpost, essBeta, diagnostBeta, file = 'mpost.RData')

# Clean up workspace
rm(list = ls())

#### PREPARE VISUALIZATIONS ####
# Load data
load('../Results/All/Data/VarPar_MF_ALL.RData')
load('../Results/All/Data/PostEst_ALL.RData')
load('../Results/All/HMSCdat_full.RData')
load('../Results/All/Data/mpost.RData')

# Load species list and prepare for matching
bio_list <- read_csv('../Data/tbl_corp_ref_species_list.csv')
lookup <- paste0('Bio_',bio_list$NODCCODE)
names(lookup) <- bio_list$Scientificname

# Extract posterior estimates
alphaLocBay        <- PostEst$alphaLocBay
alphaLoc           <- PostEst$alphaLoc
alphaTime          <- PostEst$alphaTime
lambdaPost_spatbay <- PostEst$lambda$spatbay
lambdaPost_time    <- PostEst$lambda$time
lambdaPost_loc     <- PostEst$lambda$loc
lambdaPost_codist  <- PostEst$lambda$codist

# Assign lambda matrix based on whether temporal and/or spatial structure exists (alpha > 0)
lambda_all <- if(all(alphaTime$mean == 0)) rbind(lambdaPost_codist, lambdaPost_time) else lambdaPost_codist
lambda_all <- if(all(alphaLoc$mean == 0)) rbind(lambda_all,lambdaPost_loc) else lambda_all

# Convert covariance matrix to correlation matrix
OmegaCor_mod <- cov2cor(crossprod(lambda_all))
rownames(OmegaCor_mod) <- colnames(hmsc_dat$bio[[1]])
colnames(OmegaCor_mod) <- colnames(hmsc_dat$bio[[1]])

# Save correlation matrix
OmegaCor <- data.frame(OmegaCor_mod) %>%
  mutate(NODCCODE = str_extract(row.names(.),'Bio_(.+)',group = 1)) %>%
  left_join(bio_list) %>%
  select(Scientificname,starts_with('Bio_')) %>%
  rename(any_of(lookup)) %>%
  column_to_rownames(var = 'Scientificname') %>%
  as.matrix(.)
saveRDS(OmegaCor,'../Results/All/Data/OmegaCor.RDS')

# Calculate frequency of occurrence and abundance for each species
freqoccur <- hmsc_dat %>%
  mutate(freq = map(bio, ~summarize(., across(everything(),~sum(.x>0)/length(.x))) %>%
  pivot_longer(everything(),names_to = "Species",values_to = "freq"))) %>%
  select(freq) %>%
  unnest(cols = c(freq))

abun <-  hmsc_dat %>%
  mutate(abun = map(bio, ~summarize(., across(everything(),~sum(.x, na.rm = FALSE)/length(.x))) %>%
                      pivot_longer(everything(),names_to = "Species",values_to = "abun"))) %>%
  select(abun) %>%
  unnest(cols = c(abun))

# Prep variance partitioning and model fits for visualization
vp <- t(Varpar_MF[[1]]$vals)
MF <- Varpar_MF[[2]]

Varpar_MF_tbl <- bind_cols(
                    R2  = MF$SR2,
                    vp) %>%
  mutate(R2 = ifelse(is.na(R2), 0, R2)) %>%
  bind_cols(freqoccur) %>%
  left_join(abun) %>%
  mutate(env = Environmental,
         ovshore = `Overhanging shoreline`,
         shore = `Non-overhanging shoreline`,
         bycatch = Bycatch,
         veg = Vegetation,
         sub = Substrate,
         spa = `Random: Loc`, 
         codist = `Random: Sample`, 
         time = `Random: Time`,
         bay = `Random: Bay`,
         bayloc = `Random: BayLoc`,
         year = `Random: Year`
         ) %>%
  select(Species,R2,env,ovshore,shore,bycatch,veg,sub,spa,codist,time,year, bayloc, bay, freq, abun) %>%
  drop_na() %>%
  
  mutate(codist = case_when(all(PostEst$alphaLoc$mean == 0) ~ codist + spa,
                            .default = codist),
         bay = if(all(PostEst$alphaLocBay$mean == 0)) bay + bayloc else bay,
         hab = env + ovshore + shore + bycatch + veg + sub + bay,
         hab = if(all(PostEst$alphaLocBay$mean == 0)) hab + bayloc else hab,
         disp = if(all(PostEst$alphaLoc$mean == 0)) time + year else spa + time + year,
         disp = if(all(PostEst$alphaLocBay$mean == 0)) disp else disp + bayloc,
         inter = if(all(PostEst$alphaLoc$mean == 0)) spa + codist else codist)
# Save data for visualizations
saveRDS(Varpar_MF_tbl, '../Results/All/Data/VarPar_MF_tbl.RDS')

#### MODEL FIT TRENDS ####
# Plot model fit by frequency of occurrence and abundance
VarPar_trends_freqR2 <- Varpar_MF_tbl %>%
  ggplot(aes(x = freq, y = R2)) +
  geom_point(color = 'darkblue') +
  geom_smooth(method = 'lm', alpha = 0.2, fill = 'darkblue', color = 'darkblue') +
  labs(x = 'Frequency of occurrence',
       y = 'Model fit') +
  theme_bw() +
  
  theme(legend.title = element_blank(),
        )

VarPar_trends_abunR2 <- Varpar_MF_tbl %>%
  ggplot(aes(x = log(abun), y = R2)) +
  geom_point(color = 'darkorange') +
  geom_smooth(method = 'lm', alpha = 0.2, fill = 'darkorange', color = 'darkorange') +
  labs(x = 'Log average abundance',
       y = 'Model fit') +
  theme_bw() +
  
  theme(legend.title = element_blank(),
  )

(VarPar_trends_R2 <- VarPar_trends_freqR2 + VarPar_trends_abunR2)

# Save model fit trend plot
ggsave('../Results/All/Figures/freq_R2_all.svg', VarPar_trends_R2, width = 12, height = 6)
ggsave('../Results/All/Figures/freq_R2_all.pdf', VarPar_trends_R2, width = 12, height = 6)

#### BARPLOT ####
# Calculate average variance partitioning values
VarPar_avgs <- Varpar_MF_tbl %>%
  mutate(across(env:bay, ~.x*R2)) %>%
  summarise(across(env:bay, mean))
  

# Prepare data for variance partitioning bar plot
VarPar_bar_dat <- Varpar_MF_tbl %>% 
  pivot_longer(cols = env:bay, names_to = 'varpart', values_to = 'varpct') %>%
  # Scale variance paritioning by model fit, assign each variable type to dispersal, interactions, or habitat
  mutate(varpct = varpct*R2,
         vartype = case_when(varpart %in% c('env','ovshore','shore','bycatch','veg','sub','bay') ~ 'hab',
                             varpart == 'time' ~ 'disp',
                             varpart == 'year' ~ 'disp',
                             varpart == 'spa' & all(PostEst$alphaLoc$mean == 0) ~ NA,
                             varpart == 'spa' & any(PostEst$alphaLoc$mean != 0) ~ 'disp',
                             varpart == 'bayloc' & all(PostEst$alphaLocBay$mean == 0) ~ NA,
                             varpart == 'bayloc' & any(PostEst$alphaLocBay$mean != 0) ~ 'disp',
                             varpart == 'codist' ~ 'inter')) %>%
  # Drop any NA values
  drop_na() %>%
  
  # Match NODC codes to species names, order factors for proper plotting
  mutate(NODCCODE = str_sub(Species,start = 5),
         vartype = factor(vartype, levels = c('hab','disp','inter'))) %>%
  left_join(bio_list) %>%
  arrange(desc(freq), desc(vartype)) %>%
  
  # Dynamically assign labels to variance partitioning components for legend, including the average variance partitioning values
  mutate(Species = factor(Scientificname, levels = unique(Scientificname)),
         varpart_labs = case_when(varpart == 'codist' ~ paste0('Interactions (', round(100*VarPar_avgs$codist, 1),'%)'),
                             varpart == 'time' ~ paste0('Temporal structure (', round(100*VarPar_avgs$time, 1),'%)'),
                             varpart == 'year' ~ paste0('Annual variability (', round(100*VarPar_avgs$year, 1),'%)'),
                             varpart == 'bayloc' ~ paste0('Broad spatial structure (', round(100*VarPar_avgs$bayloc, 1),'%)'),
                             varpart == 'spa' ~ paste0('Fine spatial structure (', round(100*VarPar_avgs$spa, 1),'%)'),
                             varpart == 'env' ~ paste0('Environment (', round(100*VarPar_avgs$env, 1),'%)'),
                             varpart == 'ovshore' ~ paste0('Overhanging shoreline (', round(100*VarPar_avgs$ovshore, 1),'%)'),
                             varpart == 'shore' ~ paste0('Non-overhanging shoreline (', round(100*VarPar_avgs$shore, 1),'%)'),
                             varpart == 'veg' ~ paste0('Vegetation (', round(100*VarPar_avgs$veg, 1),'%)'),
                             varpart == 'sub' ~ paste0('Substrate (', round(100*VarPar_avgs$sub, 1),'%)'),
                             varpart == 'bycatch' ~ paste0('Bycatch (', round(100*VarPar_avgs$bycatch, 1),'%)'),
                             varpart == 'bay' ~ paste0('Estuary (', round(100*VarPar_avgs$bay, 1),'%)')),
         varpart_labs = factor(varpart_labs, levels = unique(varpart_labs)))

# Assign colors to variance partitioning components
n_hab <- sum(unique(VarPar_bar_dat$varpart) %in% c('env','ovshore','shore','bycatch','veg','sub','bay'))
n_disp <- sum(unique(VarPar_bar_dat$varpart) %in% c('spa','time','year','bayloc'))
n_inter <- sum(unique(VarPar_bar_dat$varpart) %in% c('codist'))

color_palette <- c(rep("purple4", n_inter),
                   rep("darkblue", n_disp),
                   rep("darkgreen",n_hab))

fill_palette <- c('darkmagenta',
                  rev(brewer.pal(n = max(3,2*n_disp), 'Blues'))[1:n_disp],
                  rev(brewer.pal(n = 9, 'Greens'))[1:n_hab])

# Create variance partitioning bar plot
VarPar_bar <- (VarPar_bar_dat %>%
  ggplot() +
  geom_bar(aes(x = Species, y = varpct, color = varpart_labs, fill = varpart_labs), position = "stack", stat = "identity") +
  scale_color_manual(values = color_palette, guide = guide_legend(title = expression(underline('Variance component')) )) +
  scale_fill_manual(values = fill_palette, guide = guide_legend(title = expression(underline('Variance component')) )) +
  scale_y_continuous(limits = c(0,NA),expand = c(0,0,0.1,0)) +
  labs(y = 'Proportion of variance explained') +
  theme_minimal() +
  theme(
    # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size = 16),
    ))

VarPar_freq <- VarPar_bar_dat %>%
     ggplot() +
  geom_point(aes(x = Species, y = 0, color = freq), size = 4) +
  scale_color_viridis_c(limits = c(0,NA), breaks = c(0,0.2,0.4,0.6), guide = guide_legend(reverse = T,title = expression(underline('Occurrence')),override.aes = list(size = 8))) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 16, angle = -45, hjust = 0, vjust = 1),
        axis.text.y = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(0,0,0,0),
        legend.position = 'left')

(VarPar_plot <- VarPar_bar / VarPar_freq  + 
    plot_layout(heights = c(20,1), guides = 'collect') &
    theme(legend.box = 'vertical',
          legend.position = 'bottom',
          legend.direction = 'horizontal',
          legend.key.size = unit(0.5, 'in'),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 16),
          title = element_text(size = 18),
          plot.margin = unit(c(0,0.75,0,0), 'in')))

# Save bar plot
ggsave('../Results/All/Figures/varpar_bar_all.svg',VarPar_plot, height = 1.25, width = 3.25, scale = 8)  
ggsave('../Results/All/Figures/varpar_bar_all.pdf',VarPar_plot, height = 1, width = 3.25, scale = 8)  

#### TERNARY PLOT ####
# Prepare ternary plot
VarPar_tern <- Varpar_MF_tbl %>% 
  mutate(hab = env + ovshore + shore + bycatch + veg + sub + bay,
         hab = if(all(PostEst$alphaLocBay$mean == 0)) hab + bayloc else hab,
         disp = if(all(PostEst$alphaLoc$mean == 0)) time + year else spa + time + year,
         disp = if(all(PostEst$alphaLocBay$mean == 0)) disp else disp + bayloc,
         inter = if(all(PostEst$alphaLoc$mean == 0)) spa + codist else codist) %>%
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
  scale_color_viridis_c(limits = c(0,NA), breaks = c(0,0.2,0.4,0.6), guide = guide_legend(reverse = T,title = expression(underline('Occurrence')), order = 2, override.aes = list(size = 5))) +
  scale_size_area(breaks = seq(0,1,0.2), limits = c(0,1), guide = guide_legend(reverse = T, title = expression(underline(R^2)), order = 1)) +
  theme(panel.grid = element_line(color = "darkgrey", ),
        tern.axis.line = element_line(color = 'darkgrey'),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 12))

# Save ternary plot
ggsave('../Results/All/Figures/varpar_tern_all.svg',VarPar_tern, height = 2, width = 3, scale = 3)
ggsave('../Results/All/Figures/varpar_tern_all.pdf',VarPar_tern, height = 2, width = 3, scale = 3)

#### HABITAT ASSOCIATIONS ####

# Load data
load('../Results/All/Data/PostEst_ALL.RData')
bio_list <- read_csv('../Data/tbl_corp_ref_species_list.csv')

# Extract beta estimates for bay-level associations
bay_supp           <- (PostEst$bay_assoc$supp > 0.9) * PostEst$bay_assoc$supp
bay_negsupp        <- -1*(PostEst$bay_assoc$neg_supp > 0.9) * PostEst$bay_assoc$neg_supp
bay_supp           <- bay_supp + bay_negsupp
colnames(bay_supp) <- colnames(PostEst$beta$mean)
bay_supp           <- data.frame(Bay = bays, bay_supp)


bay_tbl <- bay_supp %>%
  pivot_longer(cols = starts_with("Bio_"),names_to = "Species",values_to = 'support') %>%
  filter(support != 0) %>%
  mutate(predictors = Bay
  ) %>%
  select(-Bay) %>%
  right_join(Varpar_MF_tbl, by = join_by(Species)) %>%
  mutate(NODCCODE = str_sub(Species,start = 5)) %>%
  left_join(bio_list) %>%
  drop_na()

# Extract beta estimates for predictor associations
est <- data.frame(modPreds,PostEst$beta$mean) %>%
  pivot_longer(cols = starts_with("Bio_"),names_to = "Species",values_to = 'est')

support <- data.frame(modPreds,PostEst$beta$support) %>%
  pivot_longer(cols = starts_with("Bio_"),names_to = "Species",values_to = 'support')

supportNeg <- data.frame(modPreds,PostEst$beta$supportNeg) %>%
  pivot_longer(cols = starts_with("Bio_"),names_to = "Species",values_to = 'supportNeg')
  
# Prepare posterior estimates for habitat associations
betaPost_tbl <- left_join(est,left_join(support,supportNeg)) %>%
  filter(modPreds != '(Intercept)',(support > 0.75 | supportNeg > 0.75)) %>%
  mutate(support = ifelse(support > supportNeg, support, -supportNeg)) %>%
  select(-supportNeg) %>%
  right_join(Varpar_MF_tbl, by = join_by(Species)) %>%
  mutate(NODCCODE = str_sub(Species,start = 5)) %>%
  left_join(bio_list) %>%
  rename(predictors = modPreds) %>%
  bind_rows(bay_tbl) %>%
  arrange(freq) %>%

  #  Assign predictor names based on coded column name
  mutate(Scientificname = factor(Scientificname, levels = unique(Scientificname)),
         predtype = factor(case_when(str_detect(predictors,'Bottom_') ~ 'Bottom Type',
                                     str_detect(predictors, 'Shore_') & str_detect(predictors, 'Over') ~ 'Overhanging\nShoreline',
                                     str_detect(predictors, 'Shore_') ~ 'Non-overhanging\nShoreline',
                                     str_detect(predictors,'Veg_')  ~ 'Vegetation',
                                     str_detect(predictors,'Bycatch_') ~ 'Bycatch',
                                     predictors %in% c('AP','CH','CK','TB','IR','TQ','JX') ~ 'Estuary',
                                     is.na(predictors) ~ ' ',
                                     .default = 'Environmental'),
                           levels = c(' ','Environmental','Bottom Type','Vegetation','Bycatch','Overhanging\nShoreline',"Non-overhanging\nShoreline",'Estuary')),
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
                                      predictors == 'AP' ~ 'Apalachicola Bay',
                                      predictors == 'CK' ~ 'Cedar Key',
                                      predictors == 'TB' ~ 'Tampa Bay',
                                      predictors == 'CH' ~ 'Charlotte Harbor',
                                      predictors == 'JX' ~ 'Northeast Florida',
                                      predictors == 'IR' ~ 'Northern IRL',
                                      predictors == 'TQ' ~ 'Southern IRL',
                                      .default = 'No associations')) %>%
  arrange(desc(freq))

# Separate out frequency of occurrence for each species
PostEst_Beta_points <- distinct(betaPost_tbl,Scientificname, freq)

# Make plot
plot <- (PostEst_Beta_points %>%
  ggplot() +
  geom_point(aes(x = 0,y = Scientificname, color = freq), size = 8) +
  scale_color_viridis_c(limits = c(0,NA), breaks = c(0,0.2,0.4,0.6), guide = guide_legend(reverse = T,title = expression(underline('Occurrence')))) +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_blank(),
        strip.placement = 'outside',
        strip.text = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(0,0,0,0),
        legend.position = 'left')) +
   (betaPost_tbl %>%
      filter(Scientificname %in% PostEst_Beta_points$Scientificname) %>%
      group_by(predictors) %>%
      mutate(nassoc = n()) %>%
      ungroup() %>%
  
   
   ggplot(aes(reorder_within(x = predictor_labels, by = desc(nassoc), within = predtype), y = Scientificname, fill = support)) +

  geom_tile() +
    scale_fill_gradientn(colors = brewer.pal(5, "RdBu"), guide = 'none') +
  facet_grid(cols = vars(predtype), scales = 'free', space = 'free') +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(size = 20, angle = 45,hjust = 1, vjust = 1),
        axis.text.y = element_blank(),
        strip.placement = 'outside',
        strip.text = element_text(size = 22, face = 'bold',vjust = 0),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = 'gray50'),
        panel.spacing.x = unit(1.5, "lines"),
        plot.margin = margin(0, 0, 0, 0, "in"))) +
  plot_layout(widths = c(1,60),guides = 'collect') +
  scale_x_reordered() &
  theme(legend.box = 'vertical',
        legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.key.size = unit(0.5, 'in'),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        title = element_text(size = 18),
        plot.margin = unit(c(0,0.75,0,0), 'in'))

# Save plot
ggsave('../Results/All/Figures/Abiotic_all.svg',plot = plot,width = 4 ,height = 2.5, units = 'in',scale = 10)
ggsave('../Results/All/Figures/Abiotic_all.pdf',plot = plot,width = 4 ,height = 2.5, units = 'in',scale = 10)

#### SPECIES ASSOCIATIONS ####
# Load data
OmegaCor <- readRDS('../Results/All/Data/OmegaCor.RDS')
Varpar_MF_tbl <- readRDS('../Results/All/Data/VarPar_MF_tbl.RDS')

# Prepare interactions data for plotting
intxs <- Varpar_MF_tbl %>%
  mutate(hab = env + ovshore + shore + bycatch + veg + sub + bay,
         hab = if(all(PostEst$alphaLocBay$mean == 0)) hab + bayloc else hab,
         disp = if(all(PostEst$alphaLoc$mean == 0)) time + year else spa + time + year,
         disp = if(all(PostEst$alphaLocBay$mean == 0)) disp else disp + bayloc,
         inter = if(all(PostEst$alphaLoc$mean == 0)) spa + codist else codist,
         inter = inter*R2) %>%
  mutate(NODCCODE = str_sub(Species,start = 5)) %>%
  left_join(bio_list) %>%
  mutate(spp = Scientificname) %>%
  select(spp,inter)

intxs <- intxs %>% column_to_rownames('spp')
inter <- diag(intxs$inter)
  
dimnames(inter) <- dimnames(OmegaCor)

# Order species by hierarchical clustering of pairwise correlation strength
ord <- corrMatOrder(OmegaCor, order = 'hclust')

# Assign ordering to interactions matrix
inter <- inter[ord,ord]

# Create correlation plots, stacking correlations and interactions
fig_sp <- {corrplot(inter,diag = T,add = F,is.corr = F, method = 'color',cl.pos = 'b',col = COL1('YlGn',200),tl.col = 'black',tl.cex = 0.8, cl.cex = 1)
           corrplot(OmegaCor,diag = F, order = 'hclust',add = T,tl.pos = 'n',cl.cex = 1)
           recordPlot()
}

# Save plot data
save(fig_sp, file = '../Results/All/Figures/SpeciesAssociations.RData')

# Save plots 
svg('../Results/All/Figures/SpeciesAssociations_all.svg', width = 14, height = 13)
replayPlot(fig_sp)
dev.off()

pdf('../Results/All/Figures/SpeciesAssociations_all.pdf', width = 14, height = 13)
replayPlot(fig_sp)
dev.off()
