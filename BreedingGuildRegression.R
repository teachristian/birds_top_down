##### Guild Regression Model ####
#################################

library(tidyverse)
library(tidybayes)
library(brms)
library(ggridges)
library(modelr)
set.seed(12345)

#### Load in Data #####
tree_insect_bird_merged <- read_csv('BirdsTopDown_BirdGuild.csv')

# centering with 'scale()'
center_scale <- function(x) {
  return(x - mean(x,na.rm =TRUE)/ sd(x,na.rm =TRUE))
}

levels(tree_insect_bird_merged$habitat_A) = c('woodland','scrub','urban','grassland')




tree_insect_bird_merged <- tree_insect_bird_merged %>%
  mutate(mean_cv_height_centered = center_scale(mean_cv_height),
         bird_richness_centered = center_scale(bird_richness),
         mean_cv_diameter_centered = center_scale(mean_cv_diameter),
         tmean_breeding_season_centered = center_scale(tmean_breeding_season),
         pop_density_centered = center_scale(pop_density),
         prop_forest_centered = center_scale(prop_forest),
         clumpiness_centered = center_scale(clumpiness),
         precip_breeding_season_centered = center_scale(precip_breeding_season),
         tree_richness_centered = center_scale(treediversity_specrich_total.x),
         mean_stand_age_centered = center_scale(mean_stand_age),
         mean_biomass_centered = center_scale(mean_biomass),
         sd_crown_class_centered = center_scale(mean_sd_crown_class))


woodland_merged <-tree_insect_bird_merged %>% filter(habitat_A =='woodland')  %>%
  sample_frac(.80)
test_woodland <- tree_insect_bird_merged %>% filter(habitat_A =='woodland')  %>%
  anti_join(woodland_merged)


scrub_merged <- tree_insect_bird_merged %>% filter(habitat_A =='scrub')%>%
  sample_frac(.80)

test_scrub <- tree_insect_bird_merged %>% filter(habitat_A =='scrub')  %>%
  anti_join(scrub_merged)


urban_merged <- tree_insect_bird_merged %>% filter(habitat_A =='urban')%>%
  sample_frac(.80)
test_urban <- tree_insect_bird_merged %>% filter(habitat_A =='urban')  %>%
  anti_join(urban_merged)

grassland_merged <- tree_insect_bird_merged %>% filter(habitat_A =='grassland')%>%
  sample_frac(.8)

test_grassland <- tree_insect_bird_merged %>% filter(habitat_A =='grassland')  %>%
  anti_join(grassland_merged)



# This model includes the richness term which was the most important predictor
# of bird abundance in the PSEM. I think that the coefficient of this term will likely
# be positive at that very least zero. I'm thinking of a prior that 
# cuts off at 0, but is relatively uniform until coefficient is equal to 1. 
# This results in an estimated coefficient equivalent to the linear model
# 0.295 vs .3 in the bayesian model and the intercept is also equivalent
# as the linear model below. 


model_1 <-  brm(data = tree_insect_bird_merged,
                family = gaussian(),
                formula = bf(log_bird_abundance ~ 1 + bird_richness_centered),
                prior = c(prior(normal(0, 10), class = Intercept),
                          prior(normal(0, 1), class = b),
                          prior(cauchy(0, 5), class = sigma)),
                iter = 2000,
                warmup = 1000,
                chains = 4,
                cores = 4,
                fit = 'fit1_richness',
                refresh = 0) 

summary(lm(log_bird_abundance ~ 1 + bird_richness_centered, data = tree_insect_bird_merged))
model_1 <- add_criterion(model_1, c("loo", "waic"))


model_2 <-  brm(data = tree_insect_bird_merged,
                family = gaussian(),
                formula = bf(log_bird_abundance ~ 1 + bird_richness_centered +
                                sd_crown_class_centered),
                prior = c(prior(normal(0, 10), class = Intercept),
                          prior(normal(0, 1), class = b),
                          prior(cauchy(0, 5), class = sigma)),
                iter = 2000,
                warmup = 1000,
                chains = 4,
                cores = 4,
                fit = 'fit2_richness_crowncanopy',
                refresh = 0) 
summary(lm(log_bird_abundance ~ 1 + bird_richness_centered +  sd_crown_class_centered, data = tree_insect_bird_merged))

model_2 <- add_criterion(model_2, c("loo", "waic"))

print(model_2)
plot(model_2)

# This model includes the richness term along side the 
# tree diversity metrics starting with the metric that was
# most important to the bird abundance. We have the ability
# to see if which diversity metrics are important overall.


model_3b <-  brm(data = tree_insect_bird_merged,
                 family = gaussian(),
                 formula = bf(log_bird_abundance ~ 1 + bird_richness_centered +
                                 sd_crown_class_centered + 
                                mean_cv_height_centered +
                                mean_cv_diameter_centered  +
                                mean_stand_age_centered +
                                mean_biomass_centered),
                 prior = c(prior(normal(0, 10), class = Intercept),
                           prior(normal(0, 1), class = b),
                           prior(cauchy(0, 5), class = sigma)),
                 iter = 2000,
                 warmup = 1000,
                 chains = 4,
                 cores = 4,
                 fit = 'fit3b_richness_allstructuralmetrics',
                 refresh = 0) 
print(model_3b)
plot(model_3b)

# So the only structural metric that is important is the sd_canopy_cover more than any other 
# structural metric, which is evident by the coefficients of zero on all the other metrics. 
# so moving forward our investigation is limited to the  sd_canopy cover. 

model_4 <-  brm(data = tree_insect_bird_merged,
                family = gaussian(),
                formula = bf(log_bird_abundance ~ 1 + bird_richness_centered +
                                sd_crown_class_centered + habitat_A),
                prior = c(prior(normal(0, 10), class = Intercept),
                          prior(normal(0, 1), class = b),
                          prior(cauchy(0, 5), class = sigma)),
                iter = 2000,
                warmup = 1000,
                chains = 4,
                cores = 4,
                fit = 'fit4_richness_crowncanopy_habitat',
                refresh = 0) 
print(model_4)
plot(model_4)
summary(lm(log_bird_abundance ~ 1 + bird_richness_centered +  sd_crown_class_centered
           +  habitat_A, data = tree_insect_bird_merged))

model_4 <- add_criterion(model_4, c("loo", "waic"))

# This model looks into adding the habitat as an explanatory variables
# the results suggest that the differences in abundances are non-zero
# for the urban/grassland group than the succession-scrub group.
# What I want to know is if the coefficient of  sd_crown_class_centered or any of 
# the other habitat structural metrics is different for woodlands versus the other group. 

# One of the things we wanted to know is which characteristics are important,
# that data will tell us whether a metric is important component of the response. 
# like the other interaction effects intervals contain 0. 

model_5 <-  brm(data = tree_insect_bird_merged,
                family = gaussian(),
                formula = bf(log_bird_abundance ~ 1 + bird_richness_centered +
                                sd_crown_class_centered + habitat_A + habitat_A* sd_crown_class_centered),
                prior = c(prior(normal(0, 10), class = Intercept),
                          prior(normal(0, 1), class = b),
                          prior(cauchy(0, 5), class = sigma)),
                iter = 2000,
                warmup = 1000,
                chains = 4,
                cores = 4,
                fit = 'fit5_habitatinteraction',
                refresh = 0) 
print(model_5)
plot(model_5)
summary(lm(log_bird_abundance ~ 1 + bird_richness_centered +   sd_crown_class_centered+
             habitat_A +  habitat_A* sd_crown_class_centered, data = tree_insect_bird_merged))
model_5 <- add_criterion(model_5, c("loo", "waic"))


l <- loo_compare(model_0, model_1,model_2, model_4,model_5,
                 criterion = "loo")

#Model 4 that includes habitat but does not include an interaction effect. What
# I learned form the models that included all the habitat is that breeding habitat
# does have effect on the structural metric. What  I ended up deciding to do was
#to seperate out the habitat types because
# what we wanted is to compare the estimates on the specific parameter.


##### Final Models #####
woodland_1 <-  brm(data = woodland_merged,
                   family = gaussian(),
                   formula = bf(log_bird_abundance ~ 1 + bird_richness_centered +
                                   sd_crown_class_centered),
                   prior = c(prior(normal(0, 10), class = Intercept),
                             prior(normal(0, 1), class = b),
                             prior(cauchy(0, 5), class = sigma)),
                   iter = 2000,
                   warmup = 1000,
                   chains = 4,
                   cores = 4,
                   fit = 'woodland_1',
                   refresh = 0) 

print(woodland_1)



scrub_1 <-  brm(data = scrub_merged,
                family = gaussian(),
                formula = bf(log_bird_abundance ~ 1 + bird_richness_centered +
                                sd_crown_class_centered),
                prior = c(prior(normal(0, 10), class = Intercept),
                          prior(normal(0, 1), class = b),
                          prior(cauchy(0, 5), class = sigma)),
                iter = 2000,
                warmup = 1000,
                chains = 4,
                cores = 4,
                fit = 'scrub_1',
                refresh = 0) 

print(scrub_1)

urban_1 <- brm(data = urban_merged,
               family = gaussian(),
               formula = bf(log_bird_abundance ~ 1 + bird_richness_centered +
                               sd_crown_class_centered),
               prior = c(prior(normal(0, 10), class = Intercept),
                         prior(normal(0, 1), class = b),
                         prior(cauchy(0, 5), class = sigma)),
               iter = 2000,
               warmup = 1000,
               chains = 4,
               cores = 4,
               fit = 'urban_1',
               refresh = 0) 


print(urban_1)


grassland_1 <- brm(data = grassland_merged,
                   family = gaussian(),
                   formula = bf(log_bird_abundance ~ 1 + bird_richness_centered +
                                   sd_crown_class_centered),
                   prior = c(prior(normal(0, 10), class = Intercept),
                             prior(normal(0, 1), class = b),
                             prior(cauchy(0, 5), class = sigma)),
                   iter = 2000,
                   warmup = 1000,
                   chains = 4,
                   fit = 'grassland1', 
                   cores = 4,
                   refresh = 0) 


print(grassland_1)
