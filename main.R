library(roahd)
library(DepthProc)
library(depthTools)
library(caret)
library(fda)
library(fda.usc)
library(roahd)
library(DepthProc)
library(depthTools)
library(caret)


## epi hipo function
epi_hipo_func <- function(data1, data2, index = 'epi') {
  
  all <- rbind(data1,data2)
  n = nrow(all)
  
  mix<- numeric(n)  
  miy<- numeric(n) 
  
  
  for(i in 1:n) {
    
    row <- all[i,]
    
    new_data1 <- rbind(data1, row)
    new_data2 <- rbind(data2, row)
    
    if (index == 'epi'){
      mi_data1 <- MEI(new_data1)
      mi_data2 <- MEI(new_data2)
      
      mix[i] <- mi_data1[length(mi_data1)]
      miy[i] <- mi_data2[length(mi_data2)]
    }
    
    if (index == 'hipo'){  
      mi_data1 <- MHI(new_data1)
      mi_data2 <- MHI(new_data2)
      
      mix[i] <- mi_data1[length(mi_data1)]
      miy[i] <- mi_data2[length(mi_data2)]
    }
    
  }
  
  return(list(mix, miy))
}

## plotting function
plot_function <- function(mix, miy, g, index = 'epi', plot_name) {
  xlab_text <- ifelse(index == 'epi', "MEI G1", "MHI G1")
  ylab_text <- ifelse(index == 'epi', "MEI G2", "MHI G2")
  
  p <- ggplot(appended, aes(x = mix, y = miy, shape = factor(g))) +
    geom_point(size = 2) +
    scale_shape_manual(values = c("1" = 1, "0" = 4)) +  # Set point shapes
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    xlab(xlab_text) +
    ylab(ylab_text) +
    labs(shape = "Group") +
    theme_minimal() +  # Use a minimal theme
    theme(legend.position = c(0.2, 0.7),  # Adjust the legend position
          legend.justification = c(1, 0))
  
  postscript(plot_name, 
             width = 5, height = 4, horizontal = FALSE)  # Open EPS device
  print(p)
  dev.off() 
}

##dd function
dep_func<-function(data1, data2){
  s <- rbind(data1, data2)
  n <- nrow(s)
  depx <- numeric(n)
  depy <- numeric(n)
  for (i in 1:n){
    depx[i] <- depth.FM(s[i,], data1)$dep 
    
    depy[i] <- depth.FM(s[i,], data2)$dep
  }
  return(list(depx, depy))
}

# Function to split dataset
split_dataset <- function(data, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  row_indices <- seq_len(nrow(data))
  validation_index <- createDataPartition(row_indices, p = 0.80, list = FALSE)
  validation_data <- data[-validation_index, ]
  dataset_data <- data[validation_index, ]
  return(list(dataset_data, validation_data))
}

# Function to apply classifier and append data
apply_classifier <- function(data1, data2, classifier_func, index = NULL) {
  if ("index" %in% names(formals(classifier_func))) {
    result_train <- classifier_func(data1[[1]], data2[[1]], index = index)
    result_test <- classifier_func(data1[[2]], data2[[2]], index = index)
  } else {
    result_train <- classifier_func(data1[[1]], data2[[1]])
    result_test <- classifier_func(data1[[2]], data2[[2]])
  }
  
  g_train <- c(rep(0, nrow(data1[[1]])), rep(1, nrow(data2[[1]])))
  g_test <- c(rep(0, nrow(data1[[2]])), rep(1, nrow(data2[[2]])))
  
  appended_train <- data.frame(cbind(result_train[[1]], result_train[[2]], g_train))
  appended_test <- data.frame(cbind(result_test[[1]], result_test[[2]], g_test))
  
  colnames(appended_train) <- c("x", "y", "g")
  colnames(appended_test) <- c("x", "y", "g")
  
  return(list(appended_train, appended_test))
}

# Function for model training
train_model <- function(method, x, y, metric, control) {
  set.seed(7)
  fit <- train(x, y, method = method, metric = metric, trControl = control)
  return(fit)
}


# Function to split result data
split_result_data <- function(result_data) {
  x_data <- result_data[, 1:2]
  y_factor <- factor(result_data[, 3])
  return(list(x_data, y_factor))
}

#####  Defining Synthetic Data #####
set.seed(45)

N = 200
P = 1e2

grid = seq( -1, 1, length.out = P )
grid2 = seq( -2, 2, length.out = P )

C1 = exp_cov_function( grid, alpha = 0.2, beta = 0.3 )
C2 = exp_cov_function( grid2, alpha = 0.2, beta = 0.3 )

centerline1 = (grid )^2-1
centerline2 = (grid )^2

m = generate_gauss_fdata( N, centerline1, C1 )
p = generate_gauss_fdata( N, centerline2, C2 )


####### Berkeley Growth Study Data #########
data(growth)
m = growth$hgtf
p = growth$hgtm
m = t(m)
p = t(p)

####### Tecator Data  #########
data("tecator", package = "fda.usc")
m=tecator$absorp.fdata$data[tecator$y$Fat >= 20,]
p=tecator$absorp.fdata$data[tecator$y$Fat <= 20,]

####### MCO Data #########
data(MCO)
m=MCO$intact$data[MCO$classintact == 1,]
p=MCO$intact$data[MCO$classintact == 2,]

postscript("/Users/clesmesr/Downloads/7_growth_sim.eps", width = 6, height = 5, horizontal = FALSE, onefile = FALSE)

matplot(m, type="l", lty = 1,  col="black", xlab = "age", ylab = "Height",  ylim = c(70, 200))


for (j in 1:ncol(p))
{lines(p[,j],col="grey",type="l")}


legend("topleft", legend = c("Girls − 0", 'Boys − 1'), 
       col = c("black", rep("grey", ncol(p))), lty = c(1, rep(1, ncol(p))))


# Close the EPS device
dev.off()


# Split datasets
seed_value <- 45
split_data1 <- split_dataset(m, seed = seed_value)
split_data2 <- split_dataset(p, seed = seed_value)

# EE-Classifier Hipo
epi_hipo_result <- apply_classifier(split_data1, split_data2, epi_hipo_func, index = 'hipo')
result_training_eeh <- epi_hipo_result[[1]]
result_validation_eeh <- epi_hipo_result[[2]]

# EE-Classifier Epi
epi_hipo_result <- apply_classifier(split_data1, split_data2, epi_hipo_func, index = 'epi')
result_training_eee <- epi_hipo_result[[1]]
result_validation_eee <- epi_hipo_result[[2]]

# DD-Classifier
dd_results <- apply_classifier(split_data1, split_data2, dep_func)
result_training_dd <- dd_results[[1]]
result_validation_dd <- dd_results[[2]]

# Split EEh-Classifier data
ee_classifier_data <- split_result_data(result_training_eeh)
train_x_eeh <- ee_classifier_data[[1]]
train_y_eeh <- ee_classifier_data[[2]]

ee_classifier_test_data <- split_result_data(result_validation_eeh)
test_x_eeh <- ee_classifier_test_data[[1]]
test_y_eeh <- ee_classifier_test_data[[2]]

# Split EEe-Classifier data
ee_classifier_data <- split_result_data(result_training_eee)
train_x_eee <- ee_classifier_data[[1]]
train_y_eee <- ee_classifier_data[[2]]

ee_classifier_test_data <- split_result_data(result_validation_eee)
test_x_eee <- ee_classifier_test_data[[1]]
test_y_eee <- ee_classifier_test_data[[2]]

# Split DD-Classifier data
dd_classifier_data <- split_result_data(result_training_dd)
train_x_dep <- dd_classifier_data[[1]]
train_y_dep <- dd_classifier_data[[2]]

dd_classifier_test_data <- split_result_data(result_validation_dd)
test_x_dep <- dd_classifier_test_data[[1]]
test_y_dep <- dd_classifier_test_data[[2]]


# Run algorithms using number-fold cross validation
control <- trainControl(method="cv", number=30)
metric <- "Accuracy"

# Models to train
models_ee <- c("lda", "qda", "knn", "svmRadial", "rf")

# Models to train
models_dd <- c("lda", "qda", "knn")

# Results list to store the trained models
model_results_eeh <- list()
model_results_eee <- list()
model_results_dep <- list()

# Loop through each model for EE-Classifier Hipo
for (model in models_ee) {
  modified_model_name <- ifelse(model == 'svmRadial', paste("EEh", 'svm', sep = "_"), 
                      paste("EEh", model, sep = "_"))
  fit <- train_model(model, train_x_eeh, train_y_eeh, metric, control)
  model_results_eeh[[modified_model_name]] <- fit
}

# Loop through each model for EE-Classifier Epi
for (model in models_ee) {
  modified_model_name <- ifelse(model == 'svmRadial', paste("EEe", 'svm', sep = "_"), 
                                paste("EEe", model, sep = "_"))
  fit <- train_model(model, train_x_eee, train_y_eee, metric, control)
  model_results_eee[[modified_model_name]] <- fit
}

# Loop through each model for DD-Classifier
for (model in models_dd) {
  modified_model_name <- paste("DD", model, sep = "_")
  fit <- train_model(model, train_x_dep, train_y_dep, metric, control)
  model_results_dep[[modified_model_name]] <- fit
}

# Combine the results
results_eeh <- resamples(model_results_eeh)
results_eee <- resamples(model_results_eee)
results_dd <- resamples(model_results_dep)

# Combine the summaries
summary_combined <- rbind(summary(results_eeh)$statistics$Accuracy, 
                          summary(results_eee)$statistics$Accuracy, 
                          summary(results_dd)$statistics$Accuracy)
print(summary_combined)

data_plot_name<-"/Users/catalina/Downloads/Results/9_results.eps"
postscript(data_plot_name, 
           width = 6, height = 5, horizontal = FALSE) 

boxplot(t(summary_combined), ylab = "Accuracy")

dev.off() 

