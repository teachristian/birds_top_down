################################################
#############  Bird Family SEM #################
################################################

# The following script uses the data file 'BirdsTopDown_BirdFamily.csv'
# uploaded at 10.5281/zenodo.21266745 to recreate the SEM 
# within bird family from the manuscript 'Bird richness and abundance suppress 
# invasive insect pests in forests in the eastern United States' 

library(piecewiseSEM)
library(tidyverse)

source('MECLabRepo/scripts/birds_tree_diversity/multigroup.R')
tree_insect_bird_merged_ <- read_csv('BirdsTopDown_BirdFamily.csv')


height_ <- lm(mean_cv_height ~  
                treediversity_specrich_total.x +  
                lat + 
                mean_stand_age + 
                tmean_breeding_season + 
                precip_breeding_season + 
                mean_sd_crown_class +
                log_num_pixels_forest,
              data=tree_insect_bird_merged_)


bird_abundance = lm(log_bird_abundance ~  richness +
                      treediversity_specrich_total.x +
                      mean_sd_crown_class +
                      log_num_pixels_forest + 
                      lat +
                      tmean_breeding_season + 
                      mean_cv_height,
                    data = tree_insect_bird_merged_)

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
                      data = tree_insect_bird_merged_)

height_abundance <- psem(height_, bird_abundance,insect_richness,
                         mean_cv_height %~~% clumpiness)
summary(height_abundance, conserve = TRUE)

height_multigroup = multigroup2(height_abundance,
                                group = 'Family_A')




getStandardEstimates_Family <- function(family){
  bird_family = tree_insect_bird_merged_%>%
    filter(Family_A == family)
  
  pest_glm = glm(Total_pests ~ tmean_breeding_season + 
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
                 data = bird_family)
  
  R2 <- cor(bird_family$Total_pests, predict(pest_glm, type = "response"))^2 # non-linear predictions
  
  sd.yhat <- sqrt(var(predict(pest_glm, type = "link"))/R2)
  variables_values = coef(pest_glm)[-1] # removing intercept term
  variables_names = names(variables_values)
  standard_estimates = vector(mode = 'numeric',length = length(variables_names))
  
  for(i in seq(1,length(variables_values))){
    print(i)
    variable_data = bird_family %>% pull(variables_names[i])
    if(variables_names[i] %in% c('mean_cv_height','mean_stand_age','mean_sd_crown_class','treediversity_specrich_total.x')){
      multigroup_ = height_multigroup$group.coefs[family][[1]][15:26,]
      standard_estimates[i] = multigroup_[which(multigroup_$Predictor == variables_names[i]),'Std.Estimate']
    }
    else{
      standard_estimates[i] = as.numeric(variables_values[i] * sd(variable_data/sd.yhat) )
    }
  }
  
  multigroup_height = height_multigroup$group.coefs[family][[1]][15:26,]
  standard_estimate_table = tibble(predictor = variables_names,
                                   std.estimate = standard_estimates,
                                   estimate = as.numeric(multigroup_height$Estimate),
                                   p.value = multigroup_height$P.Value,
                                   std.error = as.numeric(multigroup_height$Std.Error),
                                   significant = as.factor(ifelse(multigroup_height$P.Value < .05, 1,0)),
                                   family_name = rep(family, length(variables_names)))
  
  return(standard_estimate_table)
}


family_names = unique(tree_insect_bird_merged_$Family_A)

heightdiff_stdest_family = map_df(family_names, getStandardEstimates_Family)

family_desc = heightdiff_stdest_family %>%
  filter(predictor == 'log_bird_abundance') %>%
  arrange(abs(as.numeric(std.estimate))) %>%
  pull(family_name)

family_standardestimate_plot <- heightdiff_stdest_family %>%
  mutate(family_name  =factor(family_name,levels = family_desc)) %>%
  filter(predictor %in% c('log_bird_abundance')) %>%
  mutate(std.estimate = as.numeric(std.estimate)) %>%
  ggplot(aes(shape = significant, color = significant,alpha = significant)) +
  geom_vline(aes(xintercept = 0), linetype = 2, color = 'gray50') +
  geom_point(aes(y= family_name, x = std.estimate), size = 4) +
  scale_color_manual(name = NULL,
                     label = c('Signficant','Not Signficant'),
                     values = rev(c('#000004',"#b73779"))) +
  scale_alpha_discrete(name = NULL,
                       label = c('Signficant','Not Signficant'),
                       range= c(1,.25)) + 
  scale_shape_manual(name = NULL,
                     label = c('Signficant','Not Signficant'),
                     values = c(16,1)) + 
  scale_y_discrete(name = 'Family') +  
  scale_x_continuous(name = 'Observation-Empirical \n Standardized Coefficient',
                     limits = c(-.05,.05)) + 
  theme_bw() +
  labs(tag = "Fisher's C = 17.44
    \nP-value = 0.49
         \n18 degrees of freedom") +
  theme(legend.position = 'none', 
        plot.tag.position = c(.80,.3),
        plot.tag = element_text(hjust = 0), 
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15))


family_standardestimate_plot




