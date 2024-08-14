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
library(sf)
library(spData)
library(cowplot)


# Load geospatial data for Florida
fl <- st_read('Figures/MapFiles/FL_Shoreline_12k.shp') %>%
  mutate(state = 'FL')

# Load geospatial data for bays
bays <- st_read('Figures/MapFiles/FIM_ALL_BAYS_Sampling_Grids_20201103.shp') %>%
  filter(REGION %in% c('Apalachicola','Cedar_Key','Tampa_Bay','Charlotte_Harbor','Jacksonville','Indian_River','Tequesta'))

#### MAP ####

# Assign geospatial data to bay labels
bay_labels <- data.frame(bay = c('Apalachicola Bay\n1998',
                                 'Cedar Key\n1996',
                                 'Tampa Bay\n1989',
                                 'Charlotte Harbor\n1989',
                                 'Northeast Florida\n2001',
                                 'Northern Indian\nRiver Lagoon\n1990',
                                 'Southern Indian\nRiver Lagoon\n1997'),
                         x = c(-85,-83.75,-83.35,-83,-81,-79.75,-79.5),
                         y = c(29.25,29,27.65,26.65,30.75,28.5,27.5))

# Assign geospatial data to coast labels
coast_labels <- data.frame(coast = c('Gulf of Mexico','Atlantic Ocean'),
                            x = c(-86,-78),
                            y = c(28,28))

# Make plot for Florida and bays
fl_plot <- ggplot(data = fl) +
  geom_sf() +
  geom_sf(data = bays, fill = 'gray20',color = 'gray20') + 
  coord_sf(xlim = c(-88,-77)) +
  geom_text(data = bay_labels,aes(x = x, y = y, label = bay)) +
  geom_text(data = coast_labels,aes(x = x, y = y, label = coast),fontface = 'italic') +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(fill = NA,color = 'black'),
        axis.title = element_blank())

# Make plot for US states inset
us_plot <- us_states %>%
  mutate(highlight = ifelse(NAME == 'Florida',T,F)) %>%
  ggplot() +
  geom_sf(aes(fill = highlight)) +
  scale_fill_manual(values = c('white','black'),guide = F) +
  theme_void() +
  theme(panel.border = element_rect(fill = NA, color = 'black'),
        plot.background = element_rect(fill = 'white'))

# Make final map
(all_plot <- fl_plot + 
    inset_element(us_plot,bottom = 0.8,left = 0.755,right = 1,top = 1)) 

# Save map
ggsave(all_plot, filename = '../Results/Combo/Figures/FL_map.svg', width = 3, height = 2, scale = 3.5)

#### POLYGON TERNARY PLOT ####
# Load plot data
VarPar_MF_tbl_bay <- readRDS('../Results/Bays/Data/VarPar_MF_tbl.RDS')
VarPar_MF_tbl_reg <- readRDS('../Results/All/Data/VarPar_MF_tbl.RDS')
load('Processing/Bay/Data/HMSCdat_bay.RData')

# Combine bay and regional data
VarPar_MF_tbl <- VarPar_MF_tbl_reg %>%
  mutate(Bay = 'Regional') %>%
  bind_rows(VarPar_MF_tbl_bay) %>%
  mutate(Bay = factor(as.character(Bay), levels = c('AP','CK','TB','CH','JX','IR','TQ','Regional'),
                      labels = c('Apalachicola Bay',
                                 'Cedar Key',
                                 'Tampa Bay',
                                 'Charlotte Harbor',
                                 'Northeast Florida',
                                 'Northern Indian River',
                                 'Southern Indian River',
                                 'Across-estuary')))
# Assign color scheme
cols <- c(viridis_pal(option = 'C')(7),'black')

# Make ternary plot
(VarPar_tern_plot_combo <- VarPar_MF_tbl %>%
  ggtern(aes(x = hab, z = disp, y = inter, color = Bay, fill = Bay)) +
  stat_mean_ellipse(alpha = 0.25, geom = 'polygon') +
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
  theme_classic() +
  theme_arrowlong() +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  theme(
    axis.title = element_text(size = 20),
    tern.axis.text = element_text(size = 0),
    tern.axis.arrow.text = element_text(size = 18),
    tern.axis.arrow = element_line(linewidth = 1.75),
    legend.title = element_blank()))

# Save plot
ggsave('../Results/varpar_tern_combo.svg',VarPar_tern_plot_combo, height = 2, width = 3,scale = 2.25)

#### MODEL FIT TRENDS ####
# Make plot for frequency of occurrence vs model fit
VarPar_trends_freqR2 <- VarPar_MF_tbl %>%
  ggplot(aes(x = freq, y = R2, color = Bay, fill = Bay)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = 'lm', alpha = 0.2) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  labs(x = 'Frequency of occurrence',
       y = 'Model fit') +
  theme_bw() +
  theme(legend.title = element_blank(),
  )

# Make plot for log abundance vs model fit
VarPar_trends_abunR2 <- VarPar_MF_tbl %>%
  ggplot(aes(x = log(abun), y = R2, color = Bay, fill = Bay)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = 'lm', alpha = 0.2) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  labs(x = 'Log average abundance',
       y = NULL) +
  theme_bw() +
  theme(legend.title = element_blank(),
  )

# Combine plots
(VarPar_trends_R2 <- VarPar_trends_freqR2 + VarPar_trends_abunR2 + plot_layout(guides = 'collect'))

# Save plot
ggsave('../Results/Combo/Figures/freq_R2.svg', VarPar_trends_R2, width = 2, height = 1, scale = 3.5)


#### FACETED TERNARY PLOTS ####
# Make ternary plot faceted by bay/region
(VarPar_plot_all <- VarPar_MF_tbl %>% 
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
  theme(panel.grid = element_line(color = "darkgrey",linewidth = 0.2 ),
        tern.axis.line = element_line(color = 'darkgrey'),
        tern.axis.title = element_text(size = 12,vjust = 1),
        tern.axis.title.T = element_text(size = 12,vjust = 0),
        
        axis.text = element_text(size = 10),
        strip.background = element_blank(),
        strip.text = element_text(size = 12,face = 'bold')) +
  facet_wrap(vars(Bay), ncol = 2, dir = 'v'))

# Save plot
ggsave('../Results/Combo/Figures/varpar_tern.svg',VarPar_plot_all, height = 3, width = 2, scale = 2.6)

#### OMEGA CORRELATION VIOLIN PLOT ####
# Load data
load('../Results/Bays/Data/PostEst.RData')

# Extract lower tridiagonal, convert to dataframe for each bay
OmegaCor_bay <- map2_df(OmegaCor,names(OmegaCor), function(x,y) {
  data.frame(Bay = y, cors = x[lower.tri(x)])
  })

# Load regional data
OmegaCor_all <- readRDS('../Results/All/Data/OmegaCor.RDS') 

# Extract lower tridiagonal, convert to dataframe
OmegaCor_all <- data.frame(Bay = 'Regional', cors = OmegaCor_all[lower.tri(OmegaCor_all)])

# Assign color scheme
cols <- c(viridis_pal(option = 'C')(7),'black')

# Function to calculate mean and standard deviation
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

# Make plot
OmegaCor_plot <- bind_rows(OmegaCor_bay, OmegaCor_all) %>%
  mutate(Bay = factor(as.character(Bay), levels = c('AP','CK','TB','CH','JX','IR','TQ','Regional'),
                      labels = c('Apalachicola Bay',
                                 'Cedar Key',
                                 'Tampa Bay',
                                 'Charlotte Harbor',
                                 'Northeast Florida',
                                 'Northern Indian River',
                                 'Southern Indian River',
                                 'Across-estuary'))) %>%
  ggplot(aes(y = Bay, x = cors, fill = Bay)) +
  geom_violin(alpha = 0.5,draw_quantiles = T) +
  stat_summary(fun.data = data_summary, show.legend = F) +
  scale_fill_manual(values = cols) +
  scale_y_discrete(limits = rev,guide = 'none') +
  labs(x = 'Pairwise associations',
       y = NULL) +
  theme_classic() +
  theme(legend.title = element_blank())

# Save plot 
ggsave(filename = '../Results/Combo/Figures/OmegaCor_violin.svg',plot = OmegaCor_plot, height = 3, width = 2.75, scale = 1.5)
