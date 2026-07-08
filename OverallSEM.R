################################################
#############  Overall SEM #################
################################################

# The following script uses the data file 'BirdsTopDown_Overall.csv'
# uploaded at 10.5281/zenodo.21266745 to recreate the PSEM 
# across forest types and bird families from
# the manuscript 'Bird richness and abundance suppress 
# invasive insect pests in forests in the eastern United States' 


library(tidyverse)
library(ggplot2)
library(piecewiseSEM)
library(gridExtra)
library(lme4)

#### Load in Data #####
tree_insect_bird_merged <- read_csv('BirdsTopDown_Overall.csv')


### Linear Models #####
height_ <- lm(mean_cv_height ~  
                treediversity_specrich_total.x +  
                lat + 
                mean_stand_age + 
                tmean_breeding_season + 
                precip_breeding_season + 
                mean_sd_crown_class +
                log_num_pixels_forest,
              data=tree_insect_bird_merged)


bird_abundance = lm(log_bird_abundance ~ county_bird_richness +
                      treediversity_specrich_total.x +
                      mean_sd_crown_class +
                      log_num_pixels_forest + 
                      mean_cv_height,
                    data = tree_insect_bird_merged)

insect_richness = glm(Total_pests ~ tmean_breeding_season + 
                        precip_breeding_season  +
                        lat  +
                        log_num_pixels_forest + 
                        pop_density + 
                        prop_forest +
                        clumpiness +
                        mean_cv_height + 
                        mean_stand_age +
                        mean_sd_crown_class +
                        treediversity_specrich_total.x + 
                        log_bird_abundance,
                      family = 'poisson',
                      data = tree_insect_bird_merged)

topdown_height_abundance <- psem(height_, 
                                 bird_abundance,
                                 insect_richness)

summary(topdown_height_abundance, conserve = TRUE)


SEM_coeffecients = piecewiseSEM::coefs(topdown_height_abundance)

SEM_coeffecients_tibble = tibble(response =  SEM_coeffecients$Response,
                                 predictor = SEM_coeffecients$Predictor,
                                 estimate = SEM_coeffecients$Estimate, 
                                 std.error =SEM_coeffecients$Std.Error ,
                                 signficant = ifelse(SEM_coeffecients$P.Value<.05,1,0),
                                 df =SEM_coeffecients$DF,
                                 p.value =SEM_coeffecients$P.Value ,
                                 std.estimate = SEM_coeffecients$Std.Estimate)



bottomuptopdown = rev( c('Tree Species Richness','Mean Stand Age', 'Mean Standard Deviation \n Crown Class Code',
                         'Mean C.V. Height','Log Bird Abundance'))

overall_effect_size_plot <- SEM_coeffecients_tibble %>%
  filter(response == 'Total_pests') %>%
  filter(predictor %in% c('log_bird_abundance',
                          'treediversity_specrich_total.x',
                          'mean_stand_age',
                          'mean_sd_crown_class',
                          'mean_cv_height')) %>%
  mutate(is.negative = as.factor(ifelse(std.estimate < 0,1,0)),
         absolute_value =  abs(std.estimate)) %>%
  ggplot(aes(y = predictor, x = absolute_value, fill = is.negative)) +
  geom_col() +
  scale_x_continuous('Absolute Value \n Scale Standardized Estimate') +
  scale_fill_manual(name = NULL,
                    labels = c('Positive','Negative'),
                    values = c("#22A884FF","#414487FF")) +
  scale_y_discrete(name = 'Predictor of Pest Richness',
                   label = bottomuptopdown)  +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45,vjust = .65),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 14),
        legend.position = 'bottom')

overall_effect_size_plot

