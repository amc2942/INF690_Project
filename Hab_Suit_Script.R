# TEAM ECOLOGY
# Andy Carlino, Rowan Christie, Johnathan Tenny




################################################################################
# SETUP

# Necessary packages
library(stats)
library(rstudioapi)
#library(tidyverse)
library(raster)
#library(plotly)
library(terra)
library(data.table)
library(animation)
library(ggplot2)

# Path to where data would be stored
#setwd("C:\Users\andre\OneDrive\Documents\GitHub\INF690_Project")

################################################################################
#################### PROCESS RASTER FILES ######################################
data_dir = 'D:/LANDIS-4FRI/Project-Arizona-4FRI-master_V6/Input_files/'

biomass_files = dir(paste0(data_dir,'outputs/biomass'),full.names = T)
biomass_files = data.table(Filename=biomass_files)

get_biomass_year = function(filename){
  noext = strsplit(filename,'.',fixed=T)[[1]][1] # Remove extension
  split = strsplit(noext,'-')[[1]] # Split on hyphens
  year = split[length(split)] # Get the last item, which should be the year
  return(year)
}

get_biomass_species = function(filename){
  noext = strsplit(filename,'.',fixed=T)[[1]][1] # Remove extension
  split = strsplit(noext,'-')[[1]] # Split on hyphens
  species = split[length(split)-1] # Get second to last item
  return(species)
}

biomass_files$Year = as.integer(sapply(biomass_files$Filename,get_biomass_year))
biomass_files$Species = sapply(biomass_files$Filename,get_biomass_species)
biomass_files = biomass_files[order(Year,Species),]

total_biomass = raster::stack(biomass_files[Species=='TotalBiomass',Filename])
ponderosa_biomass = raster::stack(biomass_files[Species=='PinuPond',Filename])



################## BIOMASS BY AGE #########################################
biomass_by_age_files = dir(paste0(data_dir,'outputs/spp-biomass-by-age'),full.names = T)
biomass_by_age_files = data.table(Filename=biomass_by_age_files)

get_biomass_by_age_age = function(filename){
  noext = strsplit(filename,'.',fixed=T)[[1]][1] # Remove extension
  split = strsplit(noext,'-')[[1]] # Split on hyphens
  species = split[length(split)-1] # Get second to last item
  return(species)
}

get_biomass_by_age_species = function(filename){
  noext = strsplit(basename(filename),'.',fixed=T)[[1]][1] # Remove extension
  split = strsplit(noext,'-')[[1]] # Split on hyphens
  species = split[length(split)-2] # Get third to last item
  return(species)
}


biomass_by_age_files$Year = as.integer(sapply(biomass_by_age_files$Filename,get_biomass_year))
biomass_by_age_files$Species = sapply(biomass_by_age_files$Filename,get_biomass_by_age_species)
biomass_by_age_files$Age = sapply(biomass_by_age_files$Filename,get_biomass_by_age_age)
biomass_by_age_files = biomass_by_age_files[order(Year,Age,Species),]

old_pondo = raster::stack(biomass_by_age_files[Species=="PinuPond" & Age == 'ageclass5',Filename])

############################# TIME SINCE DISTURBANCE ######################################
fire_files = dir(paste0(data_dir,'fire'),full.names = T,pattern = 'severity')
harvest_files = dir(paste0(data_dir,'harvest'),full.names=T,pattern = 'biomass-removed')
disturbance_files = data.table(FireFilename=fire_files,HarvestFilename=harvest_files)

disturbance_files$Year = as.integer(sapply(disturbance_files$FireFilename,get_biomass_year))
disturbance_files = disturbance_files[order(Year),]

# Get time since last fire in year 0
last_disturbance = rast(paste0(data_dir,'DFFS-output/TimeOfLastFire-1.img'))

first = T
for (year in disturbance_files$Year){
  # Read the fire severity and harvest intensity rasters for this year
  fire_rast = rast(disturbance_files[Year==year,FireFilename])
  harvest_rast = rast(disturbance_files[Year==year,HarvestFilename])

  # If there was any fire OR harvest, set "disturbed" to 1
  disturbed = fire_rast + harvest_rast
  disturbed = disturbed > 1

  # Convert to a raster that stores the year of disturbance
  disturbed = disturbed * year
  disturbed[disturbed==0] = NaN

  # Get years since last disturbance relative to this year
  last_disturbance = max(last_disturbance,disturbed,na.rm = T)

  if (first) {
    years_since_disturbance = year - last_disturbance
    first = F
  }
  else{
    years_since_disturbance = c(years_since_disturbance,year-last_disturbance)
  }
}




############################# RASTER PROCESSING ###################################
# Function to plot each layer in the stack
plot_raster <- function(layer) {
  plot(layer, main = "")
}

animate_raster_stack = function(raster_stack,output_gif_name=NULL,interval = .5){
  # Create animation
  if (!is.null(output_gif_name)){
    saveGIF({
      for (i in 1:dim(total_biomass)[3]) {
        plot_raster(raster_stack[[i]])
        title(main = paste("Year:", i))
        Sys.sleep(interval)  # Adjust sleep duration as needed
      }
    }, interval = interval, movie.name = output_gif_name)
  }
  else{
    for (i in 1:dim(total_biomass)[3]) {
      plot_raster(raster_stack[[i]])
      title(main = paste("Year:", i))
      Sys.sleep(interval)  # Adjust sleep duration as needed
    }
  }
}

#animate_raster_stack(total_biomass,'Total_biomass_test.gif')
#animate_raster_stack(years_since_disturbance,'Years_since_disturbance.gif')

# Normalize the data
normalize_between_2sd <- function(data) {
  # Scale to a value between 0 and 1, clip extreme values to either 0 or 1
  data_vals = values(data)

  mean_val <- mean(data_vals)
  sd_val <- sd(data_vals)

  normalized_data <- base::scale(data_vals, center = rep(mean_val,dim(data_vals)[2]), scale = rep(sd_val * 2,dim(data_vals)[2]))
  normalized_data = pmin(pmax(normalized_data, 0), 1)

  values(data) = normalized_data

  return(data)
}

total_biomass = normalize_between_2sd(total_biomass)
ponderosa_biomass = normalize_between_2sd(ponderosa_biomass)
years_since_disturbance = normalize_between_2sd(years_since_disturbance)
old_pondo = normalize_between_2sd(old_pondo)

# Summarizing data
#landis_data_mean <- calc(landis_data, mean, na.rm = TRUE)

# Consider neighboring cells in 3x3 window
#landis_data_spatial_mean <- focal(landis_data_mean, w = matrix(1, 3, 3), fun = mean, na.rm = TRUE)

################################################################################
################################################################################
################### WILDLIFE HABITAT SUITABILITY MODELS ########################
################################################################################
################################################################################



################################################################################
# Feature Selection, starting with Mexican Spotted Owl (mso)

##### DEFINE MODEL ######
# Creating a function that calculates habitat suitability based on the above features: 
# We would need to adjust variable names to match specific outputs
calculate_suitability_mso <- function(ponderosa_biomass, old_ponderosa_biomass, total_biomass, time_since_disturbance) {
  # Assign relative variable weights based on importance
  weight_ponderosa_mso <- 0.1 # Percent of biomass that is ponderosa
  #weight_canopy_mso <- 0.35
  weight_total_biomass_mso <- 0.35
  weight_biomass_age5 <- 0.35
  weight_time_since_dist_mso <- -0.2
  # Calculate suitability score as weighted sum of selected variables 
  suitability_score_mso <- ((weight_ponderosa_mso * ponderosa_biomass) +
                          (weight_biomass_age5 * old_ponderosa_biomass) +
                          (weight_total_biomass_mso * total_biomass) +
                          (weight_time_since_dist_mso * time_since_disturbance))
  return(suitability_score_mso)
}

##### APPLY MODEL ######
# Initialize a raster with matching properties
mso_suitability_raster_stack = total_biomass

# Fill values using suitiability function and input data
values(mso_suitability_raster_stack) = calculate_suitability_mso(values(ponderosa_biomass),
                                                                 values(old_pondo),
                                                                 values(total_biomass),
                                                                 values(years_since_disturbance))

# Write raster stack to disk
raster::writeRaster(mso_suitability_raster_stack,'mso_suitability_raster_stack.tif')

# Visualize time series of mso habitat
animate_raster_stack(mso_suitability_raster_stack,'mso_suitability.gif')

# Get trend
mso_suitability_simple_sum = apply(values(mso_suitability_raster_stack),FUN=sum,MARGIN=2)
plot(seq(1,45),mso_suitability_simple_sum[1:45],type = 'l')

# Make up a "Good habitat" threshold using 75th percentile of habitat quality in year 1
year1_suitability_values = values(mso_suitability_raster_stack[[1]])
thresh = quantile(year1_suitability_values[year1_suitability_values>0],na.rm = T,probs=c(.75))[[1]]
# Add up area greater than threshold in each year
mso_suitable_area = apply(values(mso_suitability_raster_stack>=thresh),FUN=sum,MARGIN=2)
mso_suitable_area = data.table(Year = seq(1:55), MSOSuitableArea=mso_suitable_area)

ggplot(mso_suitable_area[1:45],aes(Year,MSOSuitableArea)) +
  geom_line() +
  labs(title = 'Suitable Habitat for Mexican Spotted Owl',
        subtitle = 'Under Fast Implementation of 4FRI') +
  ylab('Hectares of suitable habitat') +
  xlab('Year of Simulation')

ggsave('mso_suitability_trend.png')

####################### MORE STATS AND MODELLING ######################

# Splitting into training and testing data for Cross-validation
# This is done temporally, training on 80% of observations and testing on 20%. If it makes more sense, we can also decide to adjust this to a spatial split, reducing the raster into individual grid cells and training on 80% and testing on 20% at any specific time step.
split_index_mso <- floor(0.8 * nlayers(mso_features))

training_data_mso <- mso_features[[1:split_index_mso]]
testing_data_mso <- mso_features[[(split_index_mso + 1):nlayers(mso_features)]]





################################################################################
# Model Fitting:

# Fitting a Logistic Regression model 
# Logistic regression is often used for habitat suitability, as it shows how well the selected variables predict the probability of a specific grid cell being suitable habitat. A higher probability indicates better habitat suitability based on our selected variables. 

response_mso <- landis_data["suitability_score_mso"]
train_with_response <- stack(training_data_mso, response_mso)

# Running the glm() to run a Logistic Regression
logistic_model_mso <- glm(suitability_score_mso ~ ., 
                          data = as.data.frame(train_with_response),
                          family = "binomial")
summary(logistic_model_mso)





################################################################################
# Evaluating the Model (Comparing training and testing data)

# Making predictions
test_with_response <- stack(testing_data_mso, response_mso)
actual_labels <- test_with_response$suitability_score_mso
predicted_probabilities <- predict(logistic_model_mso, 
                                   newdata = as.data.frame(test_with_response),
                                   type = "response")

# Convert probability to binary predictions (adjust threshold to appropriate level): 
predicted_labels <- ifelse(predicted_probabilities > 0.5, 1, 0)

# Create a confusion matrix to display model accuracy:
# Confusion matrices report true positive (tp), true negative (tn), false positive (fp) and false negative (fn) counts
conf_matrix_mso <- table(Actual = actual_labels, Predicted = predicted_labels)

# Calculate and report Accuracy and Precision:
# Accuracy is overall correctness of predictions: (tp + tn) / (tp + tn + fp + fn)
# Precision is ratio of correctly predicted positive observations to total predicted positives: tp / (tp + fp)
accuracy_mso <- sum(diag(conf_matrix_mso)) / sum(conf_matrix_mso)
precision_mso <- conf_matrix_mso[2, 2] / sum(conf_matrix_mso[, 2])

cat("Model Accuracy for Mexican Spotted Owl:", accuracy_mso, "\n")
cat("Model Precision for Mexican Spotted Owl:", precision_mso, "\n")





################################################################################
# Visualizing the model: 

# Here is an example of what the data could look like visually for one time step and one grid cell. These are just example data and modeled purely as an example, but could be something to include in our paper as proof of concept:
# Feel free to run this part of the code, it works and produces a 3D plot with plotly. I'll also include a PNG of the plot in our repository so that we can use it in our reports if we want. 
# Also, I'm only looking at three variables so my computer doesn't break

set.seed(690)
example_data <- expand.grid(
  canopy_cover_ex = runif(50, 0, 100),
  basal_area_ex = runif(50, 0, 30),
  time_since_dist_ex = runif(50, 0, 10)
)

weights_ex <- c(0.02, 0.02, -0.02)
log_odds <- ((weights_ex[1] * example_data$canopy_cover_ex) +
  (weights_ex[2] * example_data$basal_area_ex) +
  (weights_ex[3] * example_data$time_since_dist_ex)
)

example_data$Predicted_Suitability <- 1 / (1 + exp(-log_odds))

plot_ly(data = example_data, 
        x = ~canopy_cover_ex, 
        y = ~basal_area_ex, 
        z = ~time_since_dist_ex,
        color = ~Predicted_Suitability, 
        colors = colorRamp(c("white","yellow", "green", "blue")),
        type = "scatter3d", 
        mode = "markers",
        marker = list(size = 2, opacity = 0.5)
) %>% 
  layout(
    scene = list(
      xaxis = list(title = "Canopy Cover"),
      yaxis = list(title = "Basal Area"),
      zaxis = list(title = "Time since last Disturbance"),
      margin = list(l = 50, r = 50, b = 60, t = 60)),
    title = "Predicted Habitat Suitability for Mexican Spotted Owls"
)
