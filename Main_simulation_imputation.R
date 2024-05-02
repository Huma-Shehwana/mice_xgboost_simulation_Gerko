rm(list = ls())

#source("Simulation_function.R")

library(mice)     # for imputation and amputation
library(purrr)    # for functional programming
library(furrr)    # for functional futures
library(mvtnorm)  # for multivariate normal data
library(magrittr) # for pipes
library(dplyr)    # for data manipulation
library(tibble)   # for tibbles
#library(mixgb) 
library(microbenchmark)
library(ggplot2)


source("train_testplot.R")

seed <- 123
set.seed(seed) 
Num_ds <- 100
prop_values <- c(20, 40, 60, 80)
m <- 5
maxit <- 5


available_cores <- availableCores() - 4
plan(multisession, workers = available_cores)


#dir.create("results_2")

#################################################################################
########################     Step 1: Data generation    #########################
#################################################################################

sigma <- matrix(data = c(1, 0.7, 0.7, 1),     # covariance matrix
                ncol = 2)

#Let's generate 1000 data sets with 1000 entries each. Each data set is saved as a separate element in the list.

simdata <- replicate(n = Num_ds, 
                     expr = mvtnorm::rmvnorm(n = 100, 
                                             mean = c(4, 1), 
                                             sigma = sigma) %>% 
                       as_tibble() %>% # make into a tibble
                       rename(x = V1, z = V2) %>% # rename columns
                       mutate(y = as.numeric(3 * x +  z + 3*I(x*z) + z^2 + rnorm(100, sd = 20))), # add y # mean 0 std of 5 or 10
                     simplify = FALSE) # keep as list of generated sets


true_vals <- c(0,3,1,3,1)

#Regression
simdata %>% 
  map_dfr(~.x %$% # for every simulated set in simdata....
            lm(y ~ x + z + x:z + I(z^2)) %>% # fit model
            coefficients) %>% # extract coefficients
  colMeans() # add all and divide by length (= average)


# mean regression
simdata %>% 
  map(~.x %$% # for every simulated set in simdata....
        summary(lm(y ~ x + z + x:z + I(z^2))) %$% # fit model
        r.squared) %>% # extract coefficients
  unlist() %>%
  mean() # add all and divide by length (= average)



#################################################################################
#############          2. Missing data 
#################################################################################

# make data missing. For each missingness percentage, generate 1000 data sets with given percentage of missing data. 
apply_ampute <- function(simdata, prop_value) {
  simdata %>%
    furrr::future_map(function(x) {
      x %>%
        ampute(prop = prop_value / 100, mech = "MAR", type = "RIGHT") %>%
        .$amp 
    }, .options = furrr_options(seed = TRUE))
}

missing_MAR_list <- map(prop_values, ~ apply_ampute(simdata, .x))
names(missing_MAR_list) <- prop_values

# missing_MAR_list is a list of missingness percentage with a sub-list of 1000 datasets

#NAs_in_data <- map(missing_MAR_list, ~ map_dbl(.x, ~ sum(is.na(.x))))



#################################################################################
####################.           3. Imputation     ###############################
#################################################################################


############################     1 - Default    ##################################

print("starting default imputation")


impute_MAR_default <- missing_MAR_list %>%
  map(function(mat_list) { # for each percentage
    furrr::future_map(mat_list, function(mat) {  # for each dataset
      result <- system.time({
        mice_result <- mice::mice(mat, 
                                  m = m, 
                                  maxit = maxit,
                                  print = FALSE)
      })
      list(mice_result = mice_result, time_taken = result) # Combining both result and time taken
    }, .options = furrr_options(seed = TRUE))
  })


mice_results_default <- map(impute_MAR_default, ~ map(., "mice_result")) #  extract imputation results
time_default <- map(impute_MAR_default, ~ map(., "time_taken"))   # extract time taken for imputation


eval_default <- mice_results_default %>% 
  map(function(mat_list) { # For each percentage
    furrr::future_map(mat_list, function(mat) { # for each dataset
      complete(mat, "all") %>% # create a list of completed data sets
        map(~.x %$% # for every completed data set....
              lm(y ~ x + z + x:z + I(z^2)) # fit  model
        ) %>% 
        pool() %>%  # pool coefficients
        summary(conf.int = TRUE) %>% # summary of coefficients
        mutate(true = c(0, 6, 3, 1, 1), # add true
               cov = conf.low < true & true < conf.high, # coverage
               bias = estimate - true,
               width = conf.high - conf.low) %>% # bias
        column_to_rownames("term")}) %>% # `term` as rownames
      Reduce("+", .) / Num_ds })


save(mice_results_default, file = "results/Default_imputation_nonlinear.RData") # imputation result
save(time_default, file = "results/Default_time_nonlinear.RData") # imputation time
save(eval_default, file = "results/Default_evaluation_nonlinear.RData") #regression results

#rm(mice_results_default,time_default,eval_default,impute_MAR_default)

############################     2 -Imputation using  RF    ##################################

print("starting RF")


impute_MAR_RF <- missing_MAR_list %>% # for each percentage
  map(function(mat_list) {
    furrr::future_map(mat_list, function(mat) { # for each dataset
      result <- system.time({
        mice_result <- mice::mice(mat, 
                                  m = m,
                                  method = "rf",
                                  maxit = maxit,
                                  print = FALSE)
      })
      list(mice_result = mice_result, time_taken = result) # Combining both imputation result and time taken
    }, .options = furrr_options(seed = TRUE))
  })


mice_results_rf <- map(impute_MAR_RF, ~ map(., "mice_result"))  #  extract imputation results   
time_RF <- map(impute_MAR_RF, ~ map(., "time_taken")) #  extract time taken for imputation results


eval_RF <- mice_results_rf %>% 
  map(function(mat_list) {
    furrr::future_map(mat_list, function(mat) {
      complete(mat, "all") %>% # create a list of completed data sets
        map(~.x %$% # for every completed data set....
              lm(y ~ x + z + x:z + I(z^2))
        ) %>% # fit linear model
        pool() %>%  # pool coefficients
        summary(conf.int = TRUE) %>% # summary of coefficients
        mutate(true = c(0, 6, 3, 1, 1), # add true
               cov = conf.low < true & true < conf.high, # coverage
               bias = estimate - true,
               width = conf.high - conf.low) %>% # bias
        column_to_rownames("term")}) %>% # `term` as rownames
      Reduce("+", .) / Num_ds})

save(mice_results_rf, file = "results/RF_imputation_nonlinear.RData") # save imputation results
save(time_RF, file = "results/RF_time_nonlinear.RData") # save time taken for imputation results
save(eval_RF, file = "results/RF_evaluation_nonlinear.RData") # save regression results fitted on imputed data

#rm(mice_results_rf,time_RF,eval_RF,impute_MAR_RF)


###################################     3 - CART   ##############################################

print("starting CART")

impute_MAR_cart <- missing_MAR_list %>% # for each missingness
  map(function(mat_list) {
    furrr::future_map(mat_list, function(mat) {  # for each dataset-
      result <- system.time({
        mice_result <- mice::mice(mat, 
                                  m = m, method = "cart",
                                  maxit = maxit,
                                  print = FALSE)
      })
      list(mice_result = mice_result, time_taken = result) # Combining both result and time taken
    }, .options = furrr_options(seed = TRUE))
  })

mice_results_cart <- map(impute_MAR_cart, ~ map(., "mice_result")) # Extracting imputation results
time_CART <- map(impute_MAR_cart, ~ map(., "time_taken")) # Time taken for imputation

eval_CART <- mice_results_cart %>% 
  map(function(mat_list) {
    furrr::future_map(mat_list, function(mat) {
      complete(mat, "all") %>% # create a list of completed data sets
        map(~.x %$% # for every completed data set....
              lm(y ~ x + z + x:z + I(z^2))) %>% # fit linear model
        pool() %>%  # pool coefficients
        summary(conf.int = TRUE) %>% # summary of coefficients
        mutate(true = c(0, 6, 3, 1, 1), # add true
               cov = conf.low < true & true < conf.high, # coverage
               bias = estimate - true,
               width = conf.high - conf.low) %>% # bias
        column_to_rownames("term")}) %>% # `term` as rownames
      Reduce("+", .) / Num_ds})


save(mice_results_cart, file = "results/CART_imputation_nonlinear.RData") # Save imputation results
save(time_CART, file = "results/CART_time_nonlinear.RData")  # Save time taken for imputation
save(eval_CART, file = "results/CART_evaluation_nonlinear.RData") # save regression results

#rm(mice_results_cart,time_CART,eval_CART,impute_MAR_cart)

###############################    4 - XGBoost - default  ##############################################

impute_MAR_xgb_default <- missing_MAR_list %>% # For each missingness
  map(function(mat_list) {
    furrr::future_map(mat_list, function(mat) {  #For each dataset
      result <- system.time({
        mice_result <- mice::mice(mat, 
                                  m = m, method = "xgb",
                                  maxit = maxit,xgb.params=NULL, # will use default parameters
                                  print = FALSE)
      })
      list(mice_result = mice_result, time_taken = result) # Combining both result and time taken
    }, .options = furrr_options(seed = TRUE))
  })


mice_xgb_default <- map(impute_MAR_xgb_default, ~ map(., "mice_result"))  # imputation results
time_xgb_default <- map(impute_MAR_xgb_default, ~ map(., "time_taken"))           # Time taken


eval_xgb_default <- mice_xgb_default %>% 
  map(function(mat_list) {
    furrr::future_map(mat_list, function(mat) {
      complete(mat, "all") %>% # create a list of completed data sets
        map(~.x %$% # for every completed data set....
              lm(y ~ x + z + x:z + I(z^2))) %>% # fit linear model
        pool() %>%  # pool coefficients
        summary(conf.int = TRUE) %>% # summary of coefficients
        mutate(true = c(0, 6, 3, 1, 1), # add true
               cov = conf.low < true & true < conf.high, # coverage
               bias = estimate - true,
               width = conf.high - conf.low) %>% # bias
        column_to_rownames("term")
    }, .options = furrr_options(seed = TRUE)) %>% # `term` as rownames
      Reduce("+", .) / Num_ds})

save(mice_xgb_default, file = "results/XGB_default_imputation_nonlinear.RData")
save(time_xgb_default, file = "results/XGB_default_time_nonlinear.RData")
save(eval_xgb_default, file = "results/XGB_default_evaluation_nonlinear.RData")

#rm(mice_xgb_default,time_xgb_default,eval_xgb_default,impute_MAR_xgb_default)

#train_testplot(simdata, mice_xgb_default,prop_values, "results/Train_test_xgb_default.pdf")


########################    4 - XGBoost - Random parameter ####################################

random_param_set <- missing_MAR_list %>%
  map(function(mat_list) { 
    furrr::future_map(mat_list, function(mat) {
      result <- system.time({
        params <- xgb_param_calc(mat,response = NULL, select_features=NULL)
      })
      list(params = params, time_taken = result)
    }, .options = furrr_options(seed = TRUE))
  })

random_param_set_parameter <- map(random_param_set, ~ map(., "params"))  # imputation results
random_param_set_time <- map(random_param_set, ~ map(., "time_taken"))           # Time taken


###########     a. Imputation using xgboost prediction 

xgb_random_maxit5 <- future_map2(missing_MAR_list, random_param_set_parameter, 
                                 .f = function(data_inner, params_inner) {
                                   future_map2(data_inner, params_inner, 
                                               .f = function(data_single, params_single) {
                                                 result <- system.time({
                                                   mice_result <- mice(data_single, m = m, method = "xgb", maxit = maxit,xgb.params =  params_single$parameter, print = FALSE)
                                                 })
                                                 list(mice_result = mice_result, time_taken = result)
                                               }, .options = furrr_options(seed = TRUE))
                                 }, .options = furrr_options(seed = TRUE))

mice_xgb_randomParam <- map(xgb_random_maxit5, ~ map(., "mice_result"))  # imputation results
time_xgb_randomParam <- map(xgb_random_maxit5, ~ map(., "time_taken"))           # Time taken


eval_xgb_randomParam <- mice_xgb_randomParam %>% 
  map(function(mat_list) {
    furrr::future_map(mat_list, function(mat) {
      complete(mat, "all") %>% # create a list of completed data sets
        map(~.x %$% # for every completed data set....
              lm(y ~ x + z + x:z + I(z^2))) %>% # fit linear model
        pool() %>%  # pool coefficients
        summary(conf.int = TRUE) %>% # summary of coefficients
        mutate(true = c(0, 6, 3, 1, 1), # add true
               cov = conf.low < true & true < conf.high, # coverage
               bias = estimate - true,
               width = conf.high - conf.low) %>% # bias
        column_to_rownames("term")
    }) %>% # `term` as rownames
      Reduce("+", .) / Num_ds})

save(mice_xgb_randomParam, file = "results/XGB_randomParam_imputation.RData")
save(time_xgb_randomParam, file = "results/XGB_randomParamImputation_time.RData")
save(random_param_set_time, file = "results/XGB_randomParamEstimation_time.RData")
save(eval_xgb_randomParam, file = "results/XGB_randomParam_evaluation.RData")

#train_testplot(simdata,missing_MAR_list, mice_results_xgb_randomParam, "results/Train_test_xgb_randomParam_NL.pdf")


########### b. Imputation using xgboost predicted.observed

xgb_random_maxit5_predictedObserved <- future_map2(missing_MAR_list, random_param_set_parameter, 
                                                   .f = function(data_inner, params_inner) {
                                                     future_map2(data_inner, params_inner, 
                                                                 .f = function(data_single, params_single) {
                                                                   result <- system.time({
                                                                     mice_result <- mice(data_single, m = m, method = "xgb", maxit = maxit,xgb.params =  params_single$parameter, match.type = "predicted.observed", print = FALSE)
                                                                   })
                                                                   list(mice_result = mice_result, time_taken = result)
                                                                 }, .options = furrr_options(seed = TRUE))
                                                   }, .options = furrr_options(seed = TRUE))

mice_results_xgb_randomParam_PO <- map(xgb_random_maxit5_predictedObserved, ~ map(., "mice_result"))  # imputation results
time_xgb_randomParam_PO <- map(xgb_random_maxit5_predictedObserved, ~ map(., "time_taken"))           # Time taken


eval_xgb_randomParam_PO <- mice_results_xgb_randomParam_PO %>% 
  furrr::future_map(function(mat_list) {
    map(mat_list, function(mat) {
      complete(mat, "all") %>% # create a list of completed data sets
        map(~.x %$% # for every completed data set....
              lm(y ~ x + z + x:z + I(z^2))) %>% # fit linear model
        pool() %>%  # pool coefficients
        summary(conf.int = TRUE) %>% # summary of coefficients
        mutate(true = c(0, 6, 3, 1, 1), # add true
               cov = conf.low < true & true < conf.high, # coverage
               bias = estimate - true,
               width = conf.high - conf.low) %>% # bias
        column_to_rownames("term")
    }) %>% # `term` as rownames
      Reduce("+", .) / Num_ds})

save(mice_results_xgb_randomParam_PO, file = "results/XGB_randomParam_PO_imputation_nonlinear.RData")
save(time_xgb_randomParam_PO, file = "results/XGB_randomParam_PO_time_nonlinear.RData")
save(eval_xgb_randomParam_PO, file = "results/XGB_randomParam_PO_evaluation_nonlinear.RData")

#train_testplot(simdata,missing_MAR_list, mice_results_xgb_randomParam_PO, "results/Train_test_xgb_randomParam_NL.pdf")

###############    c. Imputation using xgboost Original Observed

xgb_random_maxit5_OriginalObserved <- future_map2(missing_MAR_list, random_param_set_parameter, 
                                                  .f = function(data_inner, params_inner) {
                                                    future_map2(data_inner, params_inner, 
                                                                .f = function(data_single, params_single) {
                                                                  result <- system.time({
                                                                    mice_result <- mice(data_single, m = m, method = "xgb", maxit = maxit,xgb.params =  params_single$parameter, match.type = "original.observed", print = FALSE)
                                                                  })
                                                                  list(mice_result = mice_result, time_taken = result)
                                                                }, .options = furrr_options(seed = TRUE))
                                                  }, .options = furrr_options(seed = TRUE))

mice_results_xgb_randomParam_OO <- map(xgb_random_maxit5_OriginalObserved, ~ map(., "mice_result"))  # imputation results
time_xgb_randomParam_OO <- map(xgb_random_maxit5_OriginalObserved, ~ map(., "time_taken"))           # Time taken


eval_xgb_randomParam_OO <- mice_results_xgb_randomParam_OO %>% 
  map(function(mat_list) {
    furrr::future_map(mat_list, function(mat) {
      complete(mat, "all") %>% # create a list of completed data sets
        map(~.x %$% # for every completed data set....
              lm(y ~ x + z + x:z + I(z^2))) %>% # fit linear model
        pool() %>%  # pool coefficients
        summary(conf.int = TRUE) %>% # summary of coefficients
        mutate(true = c(0, 6, 3, 1, 1), # add true
               cov = conf.low < true & true < conf.high, # coverage
               bias = estimate - true,
               width = conf.high - conf.low) %>% # bias
        column_to_rownames("term")
    }) %>% # `term` as rownames
      Reduce("+", .) / Num_ds})

save(mice_results_xgb_randomParam_OO, file = "results/XGB_randomParam_OO_imputation_nonlinear.RData")
save(time_xgb_randomParam_OO, file = "results/XGB_randomParam_OO_time_nonlinear.RData")
save(eval_xgb_randomParam_OO, file = "results/XGB_randomParam_OO_evaluation_nonlinear.RData")

#train_testplot(simdata,missing_MAR_list, mice_results_xgb_randomParam_OO, "results/Train_test_xgb_randomParam_NL.pdf")


############################    4 - XGBoost - All parameter #######################################

all_param_set <- missing_MAR_list %>%
  map(function(mat_list) { 
    furrr::future_map(mat_list, function(mat) {
      result <- system.time({
        params <- xgb_param_calc(mat,response = "all", select_features=NULL)
      })
      list(params = params, time_taken = result)
    }, .options = furrr_options(seed = TRUE))
  })


all_param_set_parameter <- map(all_param_set, ~ map(., "params"))  # imputation results
all_param_set_time <- map(all_param_set, ~ map(., "time_taken"))           # Time taken

xgb_param_all_maxit5 <- future_map2(missing_MAR_list, all_param_set_parameter, 
                                    .f = function(data_inner, params_inner) {
                                      future_map2(data_inner, params_inner, 
                                                  .f = function(data_single, params_single) {
                                                    result <- system.time({
                                                      mice_result <- mice(data_single, m = m, method = "xgb", maxit = maxit,xgb.params =  params_single$parameter, print = FALSE)
                                                    })
                                                    list(mice_result = mice_result, time_taken = result)
                                                  }, .options = furrr_options(seed = TRUE))
                                    }, .options = furrr_options(seed = TRUE))


mice_results_xgb_AllParam <- map(xgb_param_all_maxit5, ~ map(., "mice_result"))  # imputation results
time_xgb_AllParam <- map(xgb_param_all_maxit5, ~ map(., "time_taken"))           # Time taken

eval_xgb_param_all_maxit5 <- mice_results_xgb_AllParam %>% 
  map(function(mat_list) {
    furrr::future_map(mat_list, function(mat) {
      complete(mat, "all") %>% # create a list of completed data sets
        map(~.x %$% # for every completed data set....
              lm(y ~ x + z + x:z + I(z^2))) %>% # fit linear model
        pool() %>%  # pool coefficients
        summary(conf.int = TRUE) %>% # summary of coefficients
        mutate(true = c(0, 6, 3, 1, 1), # add true
               cov = conf.low < true & true < conf.high, # coverage
               bias = estimate - true,
               width = conf.high - conf.low) %>% # bias
        column_to_rownames("term")
    }) %>% # `term` as rownames
      Reduce("+", .) / Num_ds})

save(mice_results_xgb_AllParam, file = "results/XGB_AllParam_imputation_nonlinear.RData")
save(time_xgb_AllParam, file = "results/XGB_AllParamImputation_time_nonlinear.RData")
save(all_param_set_time, file = "results/XGB_AllParameterEstimation_time_nonlinear.RData")
save(eval_xgb_param_all_maxit5, file = "results/XGB_AllParam_evaluation_nonlinear.RData")

########################## b. Imputation using predicted.observed

xgb_param_all_maxit5_PO <- future_map2(missing_MAR_list, all_param_set_parameter, 
                                       .f = function(data_inner, params_inner) {
                                         future_map2(data_inner, params_inner, 
                                                     .f = function(data_single, params_single) {
                                                       result <- system.time({
                                                         mice_result <- mice(data_single, m = m, method = "xgb", maxit = maxit,xgb.params =  params_single$parameter,match.type = "predicted.observed", print = FALSE)
                                                       })
                                                       list(mice_result = mice_result, time_taken = result)
                                                     }, .options = furrr_options(seed = TRUE))
                                       }, .options = furrr_options(seed = TRUE))


mice_results_xgb_AllParam_PO <- map(xgb_param_all_maxit5_PO, ~ map(., "mice_result"))  # imputation results
time_xgb_AllParam_PO <- map(xgb_param_all_maxit5_PO, ~ map(., "time_taken"))           # Time taken

eval_xgb_param_all_maxit5_PO <- mice_results_xgb_AllParam_PO %>% 
  map(function(mat_list) {
    furrr::future_map(mat_list, function(mat) {
      complete(mat, "all") %>% # create a list of completed data sets
        map(~.x %$% # for every completed data set....
              lm(y ~ x + z + x:z + I(z^2))) %>% # fit linear model
        pool() %>%  # pool coefficients
        summary(conf.int = TRUE) %>% # summary of coefficients
        mutate(true = c(0, 6, 3, 1, 1), # add true
               cov = conf.low < true & true < conf.high, # coverage
               bias = estimate - true,
               width = conf.high - conf.low) %>% # bias
        column_to_rownames("term")
    }) %>% # `term` as rownames
      Reduce("+", .) / Num_ds})

save(mice_results_xgb_AllParam_PO, file = "results/XGB_AllParam_PO_imputation_nonlinear.RData")
save(time_xgb_AllParam_PO, file = "results/XGB_AllParam_PO_time_nonlinear.RData")
save(eval_xgb_param_all_maxit5_PO, file = "results/XGB_AllParam_PO_evaluation_nonlinear.RData")


####################### c. Imputation using Original.observer

xgb_param_all_maxit5_OO <- future_map2(missing_MAR_list, all_param_set_parameter, 
                                       .f = function(data_inner, params_inner) {
                                         future_map2(data_inner, params_inner, 
                                                     .f = function(data_single, params_single) {
                                                       result <- system.time({
                                                         mice_result <- mice(data_single, m = m, method = "xgb", maxit = maxit,xgb.params =  params_single$parameter, match.type = "original.observed", print = FALSE)
                                                       })
                                                       list(mice_result = mice_result, time_taken = result)
                                                     }, .options = furrr_options(seed = TRUE))
                                       }, .options = furrr_options(seed = TRUE))


mice_results_xgb_AllParam_OO <- map(xgb_param_all_maxit5_OO, ~ map(., "mice_result"))  # imputation results
time_xgb_AllParam_OO <- map(xgb_param_all_maxit5_OO, ~ map(., "time_taken"))           # Time taken

eval_xgb_param_all_maxit5_OO <- mice_results_xgb_AllParam_OO %>% 
  map(function(mat_list) {
    furrr::future_map(mat_list, function(mat) {
      complete(mat, "all") %>% # create a list of completed data sets
        map(~.x %$% # for every completed data set....
              lm(y ~ x + z + x:z + I(z^2))) %>% # fit linear model
        pool() %>%  # pool coefficients
        summary(conf.int = TRUE) %>% # summary of coefficients
        mutate(true = c(0, 6, 3, 1, 1), # add true
               cov = conf.low < true & true < conf.high, # coverage
               bias = estimate - true,
               width = conf.high - conf.low) %>% # bias
        column_to_rownames("term")
    }) %>% # `term` as rownames
      Reduce("+", .) / Num_ds})

save(mice_results_xgb_AllParam_OO, file = "results/XGB_AllParam_OO_imputation_nonlinear.RData")
save(time_xgb_AllParam_OO, file = "results/XGB_AllParam_OO_time_nonlinear.RData")
save(eval_xgb_param_all_maxit5_OO, file = "results/XGB_AllParam_OO_evaluation_nonlinear.RData")


###########################################################################################################
###########################################################################################################
###########################################################################################################




###############################    5 - XGBoost (maxit = 1)  ##############################################

impute_MAR_xgb_default_1 <- missing_MAR_list %>% # For each missingness
  map(function(mat_list) {
    furrr::future_map(mat_list, function(mat) {  #For each dataset
      result <- system.time({
        mice_result <- mice::mice(mat, 
                                  m = m, method = "xgb",
                                  maxit = 1,xgb.params=NULL, # will use default parameters
                                  print = FALSE)
      })
      list(mice_result = mice_result, time_taken = result) # Combining both result and time taken
    }, .options = furrr_options(seed = TRUE))
  })


mice_results_xgb_default_1 <- map(impute_MAR_xgb_default_1, ~ map(., "mice_result"))  # imputation results
time_xgb_default_1 <- map(impute_MAR_xgb_default_1, ~ map(., "time_taken"))           # Time taken


eval_xgb_default_1 <- mice_results_xgb_default_1 %>% 
  map(function(mat_list) {
    furrr::future_map(mat_list, function(mat) {
      complete(mat, "all") %>% # create a list of completed data sets
        map(~.x %$% # for every completed data set....
              lm(y ~ x + z + x:z + I(z^2))) %>% # fit linear model
        pool() %>%  # pool coefficients
        summary(conf.int = TRUE) %>% # summary of coefficients
        mutate(true = c(0, 6, 3, 1, 1), # add true
               cov = conf.low < true & true < conf.high, # coverage
               bias = estimate - true,
               width = conf.high - conf.low) %>% # bias
        column_to_rownames("term")
    }) %>% # `term` as rownames
      Reduce("+", .) / Num_ds})

save(mice_results_xgb_default_1, file = "results/XGB_default_maxit1.RData")
save(time_xgb_default_1, file = "results/XGB_default_maxit1.RData")
save(eval_xgb_default_1, file = "results/XGB_default_maxit1.RData")


#train_testplot(simdata,missing_MAR_list, mice_results_xgb_default_1, "results/Train_test_xgb_default_maxit1.pdf")


########################    5 - XGBoost - Random parameter (maxit = 1) ####################################


xgb_random_maxit1 <- future_map2(missing_MAR_list, random_param_set_parameter, 
                                 .f = function(data_inner, params_inner) {
                                   future_map2(data_inner, params_inner, 
                                               .f = function(data_single, params_single) {
                                                 result <- system.time({
                                                   mice_result <- mice(data_single, m = m, method = "xgb", maxit = 1,xgb.params =  params_single$parameter, print = FALSE)
                                                 })
                                                 list(mice_result = mice_result, time_taken = result)
                                               }, .options = furrr_options(seed = TRUE))
                                 }, .options = furrr_options(seed = TRUE))

mice_results_xgb_randomParam_1 <- map(xgb_random_maxit1, ~ map(., "mice_result"))  # imputation results
time_xgb_randomParam_1 <- map(xgb_random_maxit1, ~ map(., "time_taken"))           # Time taken


eval_xgb_randomParam_1 <- mice_results_xgb_randomParam_1 %>% 
  map(function(mat_list) {
    furrr::future_map(mat_list, function(mat) {
      complete(mat, "all") %>% # create a list of completed data sets
        map(~.x %$% # for every completed data set....
              lm(y ~ x + z + x:z + I(z^2))) %>% # fit linear model
        pool() %>%  # pool coefficients
        summary(conf.int = TRUE) %>% # summary of coefficients
        mutate(true = c(0, 6, 3, 1, 1), # add true
               cov = conf.low < true & true < conf.high, # coverage
               bias = estimate - true,
               width = conf.high - conf.low) %>% # bias
        column_to_rownames("term")
    }) %>% # `term` as rownames
      Reduce("+", .) / Num_ds})

save(mice_results_xgb_randomParam_1, file = "results/XGB_randomParam_maxit1.RData")
save(time_xgb_randomParam, file = "results/XGB_randomParamImputation_time_maxit1.RData")
save(eval_xgb_randomParam, file = "results/XGB_randomParam_evaluation_maxit1.RData")

#train_testplot(simdata,missing_MAR_list, mice_results_xgb_randomParam_1, "results/Train_test_xgb_randomParam_1.pdf")


########### b. Imputation using xgboost predicted.original

xgb_random_maxit1_predictedObserved <- future_map2(missing_MAR_list, random_param_set_parameter, 
                                                   .f = function(data_inner, params_inner) {
                                                     future_map2(data_inner, params_inner, 
                                                                 .f = function(data_single, params_single) {
                                                                   result <- system.time({
                                                                     mice_result <- mice(data_single, m = m, method = "xgb", maxit = 1,xgb.params =  params_single$parameter, match.type = "predicted.observed", print = FALSE)
                                                                   })
                                                                   list(mice_result = mice_result, time_taken = result)
                                                                 }, .options = furrr_options(seed = TRUE))
                                                   }, .options = furrr_options(seed = TRUE))

mice_results_xgb_randomParam_PO_1 <- map(xgb_random_maxit1_predictedObserved, ~ map(., "mice_result"))  # imputation results
time_xgb_randomParam_PO_1 <- map(xgb_random_maxit1_predictedObserved, ~ map(., "time_taken"))           # Time taken


eval_xgb_randomParam_PO_1 <- mice_results_xgb_randomParam_PO_1 %>% 
  furrr::future_map(function(mat_list) {
    map(mat_list, function(mat) {
      complete(mat, "all") %>% # create a list of completed data sets
        map(~.x %$% # for every completed data set....
              lm(y ~ x + z + x:z + I(z^2))) %>% # fit linear model
        pool() %>%  # pool coefficients
        summary(conf.int = TRUE) %>% # summary of coefficients
        mutate(true = c(0, 6, 3, 1, 1), # add true
               cov = conf.low < true & true < conf.high, # coverage
               bias = estimate - true,
               width = conf.high - conf.low) %>% # bias
        column_to_rownames("term")
    }) %>% # `term` as rownames
      Reduce("+", .) / Num_ds})

save(mice_results_xgb_randomParam_PO_1, file = "results/XGB_randomParam_PO_imputation_maxit1.RData")
save(time_xgb_randomParam_PO_1, file = "results/XGB_randomParam_PO_time_maxit1.RData")
save(eval_xgb_randomParam_PO_1, file = "results/XGB_randomParam_PO_eval_maxit1.RData")

#train_testplot(simdata,missing_MAR_list, mice_results_xgb_randomParam_PO_1, "results/Train_test_xgb_randomParam_PO_maxit1.pdf")

###############    c. Imputation using xgboost Original Observed

xgb_random_maxit1_OriginalObserved <- future_map2(missing_MAR_list, random_param_set_parameter, 
                                                  .f = function(data_inner, params_inner) {
                                                    future_map2(data_inner, params_inner, 
                                                                .f = function(data_single, params_single) {
                                                                  result <- system.time({
                                                                    mice_result <- mice(data_single, m = m, method = "xgb", maxit = 1,xgb.params =  params_single$parameter, match.type = "original.observed", print = FALSE)
                                                                  })
                                                                  list(mice_result = mice_result, time_taken = result)
                                                                }, .options = furrr_options(seed = TRUE))
                                                  }, .options = furrr_options(seed = TRUE))

mice_results_xgb_randomParam_OO_1 <- map(xgb_random_maxit1_OriginalObserved, ~ map(., "mice_result"))  # imputation results
time_xgb_randomParam_OO_1 <- map(xgb_random_maxit1_OriginalObserved, ~ map(., "time_taken"))           # Time taken


eval_xgb_randomParam_OO_1 <- mice_results_xgb_randomParam_OO_1 %>% 
  furrr::future_map(function(mat_list) {
    map(mat_list, function(mat) {
      complete(mat, "all") %>% # create a list of completed data sets
        map(~.x %$% # for every completed data set....
              lm(y ~ x + z + x:z + I(z^2))) %>% # fit linear model
        pool() %>%  # pool coefficients
        summary(conf.int = TRUE) %>% # summary of coefficients
        mutate(true = c(0, 6, 3, 1, 1), # add true
               cov = conf.low < true & true < conf.high, # coverage
               bias = estimate - true,
               width = conf.high - conf.low) %>% # bias
        column_to_rownames("term")
    }) %>% # `term` as rownames
      Reduce("+", .) / Num_ds})

save(mice_results_xgb_randomParam_OO_1, file = "results/XGB_randomParam_OO_imputation_maxit1.RData")
save(time_xgb_randomParam_OO_1, file = "results/XGB_randomParam_OO_time_maxit1.RData")
save(eval_xgb_randomParam_OO_1, file = "results/XGB_randomParam_OO_evaluation_maxit1.RData")

#train_testplot(simdata,missing_MAR_list, mice_results_xgb_randomParam_OO, "results/Train_test_xgb_randomParam_NL.pdf")


############################    5 - XGBoost - All parameter (maxit = 1) #######################################


xgb_param_all_maxit1 <- future_map2(missing_MAR_list, all_param_set_parameter, 
                                    .f = function(data_inner, params_inner) {
                                      future_map2(data_inner, params_inner, 
                                                  .f = function(data_single, params_single) {
                                                    result <- system.time({
                                                      mice_result <- mice(data_single, m = m, method = "xgb", maxit = 1,xgb.params =  params_single$parameter, print = FALSE)
                                                    })
                                                    list(mice_result = mice_result, time_taken = result)
                                                  }, .options = furrr_options(seed = TRUE))
                                    }, .options = furrr_options(seed = TRUE))


mice_results_xgb_AllParam_1 <- map(xgb_param_all_maxit1, ~ map(., "mice_result"))  # imputation results
time_xgb_AllParam_1 <- map(xgb_param_all_maxit1, ~ map(., "time_taken"))           # Time taken

eval_xgb_param_all_maxit1 <- mice_results_xgb_AllParam_1 %>% 
  furrr::future_map(function(mat_list) {
    map(mat_list, function(mat) {
      complete(mat, "all") %>% # create a list of completed data sets
        map(~.x %$% # for every completed data set....
              lm(y ~ x + z + x:z + I(z^2))) %>% # fit linear model
        pool() %>%  # pool coefficients
        summary(conf.int = TRUE) %>% # summary of coefficients
        mutate(true = c(0, 6, 3, 1, 1), # add true
               cov = conf.low < true & true < conf.high, # coverage
               bias = estimate - true,
               width = conf.high - conf.low) %>% # bias
        column_to_rownames("term")
    }) %>% # `term` as rownames
      Reduce("+", .) / Num_ds})

save(mice_results_xgb_AllParam_1, file = "results/XGB_AllParam_imputation_maxit1.RData")
save(time_xgb_AllParam_1, file = "results/XGB_AllParamImputation_time_maxit1.RData")
save(eval_xgb_param_all_maxit1, file = "results/XGB_AllParam_evaluation_maxit1.RData")

########################## b. Imputation using predicted.observed

xgb_param_all_maxit1_PO <- future_map2(missing_MAR_list, all_param_set_parameter, 
                                       .f = function(data_inner, params_inner) {
                                         future_map2(data_inner, params_inner, 
                                                     .f = function(data_single, params_single) {
                                                       result <- system.time({
                                                         mice_result <- mice(data_single, m = m, method = "xgb", maxit = 1,xgb.params =  params_single$parameter,match.type = "predicted.observed", print = FALSE)
                                                       })
                                                       list(mice_result = mice_result, time_taken = result)
                                                     }, .options = furrr_options(seed = TRUE))
                                       }, .options = furrr_options(seed = TRUE))


mice_results_xgb_AllParam_PO_1 <- map(xgb_param_all_maxit1_PO, ~ map(., "mice_result"))  # imputation results
time_xgb_AllParam_PO_1 <- map(xgb_param_all_maxit1_PO, ~ map(., "time_taken"))           # Time taken

eval_xgb_param_all_maxit1_PO <- mice_results_xgb_AllParam_PO_1 %>% 
  furrr::future_map(function(mat_list) {
    map(mat_list, function(mat) {
      complete(mat, "all") %>% # create a list of completed data sets
        map(~.x %$% # for every completed data set....
              lm(y ~ x + z + x:z + I(z^2))) %>% # fit linear model
        pool() %>%  # pool coefficients
        summary(conf.int = TRUE) %>% # summary of coefficients
        mutate(true = c(0, 6, 3, 1, 1), # add true
               cov = conf.low < true & true < conf.high, # coverage
               bias = estimate - true,
               width = conf.high - conf.low) %>% # bias
        column_to_rownames("term")
    }) %>% # `term` as rownames
      Reduce("+", .) / Num_ds})

save(mice_results_xgb_AllParam_PO_1, file = "results/XGB_AllParam_PO_imputation_maxit1.RData")
save(time_xgb_AllParam_PO_1, file = "results/XGB_AllParam_PO_time_maxit1.RData")
save(eval_xgb_param_all_maxit1_PO, file = "results/XGB_AllParam_PO_evaluation_maxit1.RData")


####################### c. Imputation using Original.observer

xgb_param_all_maxit1_OO <- future_map2(missing_MAR_list, all_param_set_parameter, 
                                       .f = function(data_inner, params_inner) {
                                         future_map2(data_inner, params_inner, 
                                                     .f = function(data_single, params_single) {
                                                       result <- system.time({
                                                         mice_result <- mice(data_single, m = m, method = "xgb", maxit = 1,xgb.params =  params_single$parameter, match.type = "original.observed", print = FALSE)
                                                       })
                                                       list(mice_result = mice_result, time_taken = result)
                                                     }, .options = furrr_options(seed = TRUE))
                                       }, .options = furrr_options(seed = TRUE))


mice_results_xgb_AllParam_OO_1 <- map(xgb_param_all_maxit1_OO, ~ map(., "mice_result"))  # imputation results
time_xgb_AllParam_OO_1 <- map(xgb_param_all_maxit1_OO, ~ map(., "time_taken"))           # Time taken

eval_xgb_param_all_maxit1_OO <- mice_results_xgb_AllParam_OO_1 %>% 
  furrr::future_map(function(mat_list) {
    map(mat_list, function(mat) {
      complete(mat, "all") %>% # create a list of completed data sets
        map(~.x %$% # for every completed data set....
              lm(y ~ x + z + x:z + I(z^2))) %>% # fit linear model
        pool() %>%  # pool coefficients
        summary(conf.int = TRUE) %>% # summary of coefficients
        mutate(true = c(0, 6, 3, 1, 1), # add true
               cov = conf.low < true & true < conf.high, # coverage
               bias = estimate - true,
               width = conf.high - conf.low) %>% # bias
        column_to_rownames("term")
    }) %>% # `term` as rownames
      Reduce("+", .) / Num_ds})

save(mice_results_xgb_AllParam_OO_1, file = "results/XGB_AllParam_OO_imputation_maxit1.RData")
save(time_xgb_AllParam_OO_1, file = "results/XGB_AllParam_OO_time_maxit1.RData")
save(eval_xgb_param_all_maxit1_OO, file = "results/XGB_AllParam_OO_evaluation_maxit1.RData")


###########################################################################################################
###########################################################################################################
###########################################################################################################



cor_list <- map(simdata, ~ cor(.x)) %>% Reduce("+", .) / length(simdata)

cors_all <- cbind( 
  cor_withY = cor_list[c("y","x","z"),"y"],
  cor_withX =cor_list[c("y","x","z"),"x"],
  cor_withZ = cor_list[c("y","x","z"),"z"]
)


saveRDS(cors_all, "Correlations_XYZ.rds")

