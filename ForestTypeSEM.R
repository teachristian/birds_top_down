################################################
#############  Forest Type SEM #################
################################################

# The following script uses the data file 'BirdsTopDown_Overall.csv'
# uploaded at 10.5281/zenodo.21266745 to recreate the SEM 
# within forest types from the manuscript 'Bird richness and abundance suppress 
# invasive insect pests in forests in the eastern United States' 

library(piecewiseSEM) 
library(tidyverse)

source('MECLabRepo/scripts/birds_tree_diversity/multigroup.R')
tree_insect_bird_merged <- read_csv('BirdsTopDown_Overall.csv')


### Height ####
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

height_abundance <- psem(height_, bird_abundance,insect_richness)
summary(height_abundance, conserve = TRUE)

height_multigroup = multigroup2(height_abundance ,
                                standardize = 'scale',
                                group = 'broad_group')



getStandardEstimates <- function(group){
  forest_type = tree_insect_bird_merged %>%
    filter(broad_group == group)
  # variables needs to match how the variables are entered in the SEM
  pest_glm = glm(Total_pests ~ tmean_breeding_season + 
                   precip_breeding_season  +
                   lat + 
                   log_num_pixels_forest +
                   pop_density + 
                   prop_forest +
                   clumpiness + 
                   mean_cv_height + 
                   mean_stand_age +
                   mean_sd_crown_class  + 
                   treediversity_specrich_total.x  +
                   log_bird_abundance,
                 family = 'poisson',
                 data = forest_type)
  R2 <- cor(forest_type$Total_pests, predict(pest_glm, type = "response"))^2 # non-linear predictions
  
  sd.yhat <- sqrt(var(predict(pest_glm, type = "link"))/R2)
  variables_values = coef(pest_glm)[-1] # removing intercept term
  variables_names = names(variables_values) 
  standard_estimates = vector(mode = 'numeric',length = length(variables_names))
  multigroup_ = height_multigroup$group.coefs[group][[1]][13:24,]
  
  for(i in seq(1,length(variables_values))){
    variable_data = forest_type %>% pull(variables_names[i])
    
    # check if there is a standard estimate already calculatesd for that variables
    no_standard_estimates = variables_names[which(multigroup_[,8] == '-')]
    
    if(variables_names[i] %in% no_standard_estimates){
      standard_estimates[i] = variables_values[i] * sd(variable_data/sd.yhat)
    }
    else{
      standard_estimates[i] = multigroup_[which(multigroup_$Predictor == variables_names[i]),'Std.Estimate']
    }
    
    standard_estimate_table = tibble(predictor = variables_names,
                                     std.estimate = as.numeric(standard_estimates),
                                     estimate = as.numeric(multigroup_$Estimate),
                                     p.value = multigroup_$P.Value,
                                     std.error = multigroup_$Std.Error,
                                     significant = as.factor(ifelse(multigroup_$P.Value < .05, 1,0)),
                                     broad_group = rep(group, length(variables_names)))
    
    
    
  }
  return(standard_estimate_table)
  
}



broad_group_names = unique(tree_insect_bird_merged$broad_group)

heightdiff_stdest = map_df(broad_group_names, ~getStandardEstimates(.x))


descendinggroups =  tree_insect_bird_merged %>%
  group_by(broad_group) %>%
  dplyr::summarise(num_counties = n()) %>%
  arrange(desc(num_counties)) %>%
  pull(broad_group)


heightdiff_stdest <- heightdiff_stdest %>%
  mutate(broad_group = factor(broad_group, levels = descendinggroups ))


bottomuptopdown = rev( c('Tree Species Richness','Mean Stand Age',
                         'Mean St. Dev. \n Crown Class Code',
                         'Mean C.V. Height','Log Bird Abundance'))

heightdiff_effect_sizes_stdest_plot <- heightdiff_stdest %>%
  filter(predictor %in% c('log_bird_abundance',
                          'treediversity_specrich_total.x',
                          'mean_stand_age',
                          'mean_sd_crown_class',
                          'mean_cv_height')) %>%
  mutate(is.negative = as.factor(ifelse(std.estimate < 0,1,0)))%>%
  ggplot(aes(y = predictor, x = abs(std.estimate), fill = is.negative, alpha = as.numeric(significant))) +
  geom_col() +
  theme_bw() +
  facet_wrap(~broad_group) +
  scale_x_continuous('Absolute Value \n Scale Standardized Estimate') +
  scale_alpha(guide = 'none') +
  scale_fill_manual(name = NULL,
                    values = c("#22A884FF", "#414487FF" ),
                    label = c('Postive','Negative')) +
  scale_y_discrete(labels= bottomuptopdown,
                   name = 'Predictor of Pest Richness')  +
  theme_bw() + 
  labs(tag = "Fisher's C = 10.87
    \nP-value = 0.99
         \n24 degrees of freedom") +
  theme(legend.position = 'bottom',
        axis.text = element_text(size = 12,hjust = 0),
        axis.ticks.y = element_line(linewidth = 2,),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        plot.tag.position = c(.80,.3),
        plot.tag = element_text(hjust = 0))


heightdiff_effect_sizes_stdest_plot

ggsave('figures/birdstopdown/figure_6.png',
       plot =   heightdiff_effect_sizes_stdest_plot,
       dpi = 300)
