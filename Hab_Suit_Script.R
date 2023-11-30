# TEAM ECOLOGY
# Andy Carlino, Rowan Christie, Johnathan Tenny




################################################################################
# SETUP

# Necessary packages
library(stats)
library(rstudioapi)
library(tidyverse)
library(raster)
library(plotly)

# Path to where data would be stored
setwd("C:\Users\andre\OneDrive\Documents\GitHub\INF690_Project")

################################################################################





################################################################################
################################################################################
################### WILDLIFE HABITAT SUITABILITY MODELS ########################
################################################################################
################################################################################



# NOTE: None of this code actually runs right now, it assumes we are pulling the raster data from our directory, loading that into our environment, and using that to then run the model. This is really just to show proof of concept, on how we would hypothetically create a model, using the Mexican Spotted Owl as our example system. 



################################################################################
# Data Preprocessing (Loading in data from raster files)

# Reading LANDIS-II Outputs: Assumes .tif raster files
raster_list <- list.files(pattern = "*.tif", full.names = TRUE)
landis_data <- stack(raster_list)

# Summarizing data
landis_data_mean <- calc(landis_data, mean, na.rm = TRUE)

# Consider neighboring cells in 3x3 window
landis_data_spatial_mean <- focal(landis_data_mean, w = matrix(1, 3, 3), fun = mean, na.rm = TRUE)





################################################################################
# Feature Selection, starting with Mexican Spotted Owl (mso)
# P. Ponderosa dominated, Canopy cover, basal area, time since disturbance

# Names should be adjusted to match variable names of landis outputs
mso_features <- landis_data[c("ponderosa_distribution", "canopy_cover", "basal_area", "time_since_disturbance")]

# Creating a function that calculates habitat suitability based on the above features: 
# We would need to adjust variable names to match specific outputs
calculate_suitability_mso <- function(ponderosa_distribution, canopy_cover, basal_area, time_since_disturbance) {
  # Assign relative variable weights based on importance
  weight_ponderosa_mso <- 0.1
  weight_canopy_mso <- 0.35
  weight_basal_area_mso <- 0.35
  weight_time_since_dist_mso <- -0.2
  # Calculate suitability score as weighted sum of selected variables 
  suitability_score_mso <- ((weight_ponderosa * ponderosa_distribution) + 
                          (weight_canopy * canopy_cover) + 
                          (weight_basal_area * basal_area) + 
                          (weight_time_since_dist * time_since_disturbance))
  return(suitability_score_mso)
}

# And now, apply the function to calculate the suitability for each grid cell:
landis_data$suitability_score_mso <- calculate_suitability_mso(
  landis_data$ponderosa_distribution,
  landis_data$canopy_cover,
  landis_data$basal_area, 
  landis_data$time_since_disturbance
)

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
