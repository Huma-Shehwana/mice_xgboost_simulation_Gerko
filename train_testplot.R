train_testplot<-function(orig_data, impute_res,prop_values, output_name){
  
  
  
  # orig_data : A list with 10 simulated data sets 
  # impute_res: mice output object. List of list. same dimension as missing_data in main document
  
  #########################################################################################
  #####################     Calculation of testing Error      #############################
  #########################################################################################
 
  test_error_all <- map(impute_res, .f = function(one_p){   # call to testing_error function form mice to calculate test error
    future_map2(orig_data, one_p, .f = function(orig_data_1,imputed_res_1) { 
                mice::testing_error(orig_data_1,imputed_res_1)}, .options = furrr_options(seed = TRUE))
    })
  
  
  #########################################################################################
  ##################             Training and test arrangement       ######################
  #########################################################################################
  
  var_names <- colnames(orig_data[[1]]) # variable names
  
  train_error <- lapply(impute_res, function(inner_list) lapply(inner_list, '[[','trainError')) # extract training error from imputation results
  
  p <- list()
  for(i in var_names){ # for each variable
    
    trainingErrorArranged <- map2(train_error, seq_along(train_error), function(inner_list, outer_index) {
      map2(inner_list, seq_along(inner_list), function(df, inner_index) {
        df[[i]] %>%
          as.data.frame() %>%
          mutate(Percent = prop_values[outer_index],
                 Dataset = inner_index)
      })
    }) # add percentage and dataset index to training error list
    
    trainingError_df <- do.call(rbind, do.call(c, trainingErrorArranged)) # convert to dataframe
    trainingError_df$variable <- i # add variable name
    
    transformed_train <- trainingError_df %>%
      pivot_longer(cols = c("V1","V2","V3","V4","V5"), names_to = "Imputation", values_to = "Error") %>%
      mutate(ErrorType = "TrainingError") #melt data and add trainingError in errortype column
    
    transformed_train$Imputation <- gsub("^V", "", transformed_train$Imputation ) # format imputation column
    
    
    # test error arrangement
    
    testingErrorArranged <- map2(test_error_all, seq_along(test_error_all), function(inner_list, outer_index) {
      map2(inner_list, seq_along(inner_list), function(df, inner_index) {
        df[[i]] %>%
          as.data.frame() %>%
          mutate(Percent = prop_values[outer_index],
                 Dataset = inner_index)
      }) #add percentage and dataset index
    })
    
    testError_df <- do.call(rbind, do.call(c, testingErrorArranged)) # convert to dataframe
    testError_df$variable <- i # add variable name
    
    transformed_test <- testError_df %>%
      pivot_longer(cols = c("V1","V2","V3","V4","V5"), names_to = "Imputation", values_to = "Error") %>%
      mutate(ErrorType = "TestingError") #melt the dataframe
    
    transformed_test$Imputation <- gsub("^V", "", transformed_test$Imputation ) # format imputation column
    
    
    # combine train and test error
    combined_Error<-rbind(transformed_train, transformed_test) #combine train and test error
    
    
    combined_Error$Percent<-paste(combined_Error$Percent, "% missingness", sep = "")
    combined_Error$Dataset <- factor(combined_Error$Dataset, levels = unique(combined_Error$Dataset))
    
    
    color_palette <- c("TrainingError" = "blue", "TestingError" = "red")
    
    # Plot using ggplot2
    p[[i]] <- ggplot(combined_Error, aes(x = Dataset, y = Error, color = ErrorType, group = interaction(ErrorType, Imputation))) +
      geom_line() +
      facet_wrap(~ Percent) +
      scale_color_manual(values = color_palette) +
      labs(x = "Dataset", y = "Error", color = "Error Type") + ggtitle(paste("Variable: ", unique(combined_Error$variable))) +
      theme_minimal()+
      theme(plot.title = element_text(hjust = 0.5)) 
  }
  
  
  ##################################################################
  ###############         Plotting        ##########################
  ##################################################################
  
  
  pdf(output_name)
  print(p)
  dev.off()
  
}

