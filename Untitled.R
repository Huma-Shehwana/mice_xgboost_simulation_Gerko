
rm(list = ls())

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
library(tidyr)

source("train_testplot.R")

seed <- 123
Num_ds <- 100
prop_values <- c(20, 40, 60, 80)
m <- 5
maxit <- 1

true_vals <- c(0,3,1,3,1)


set.seed(seed) 

available_cores <- availableCores() - 1
plan(multisession, workers = available_cores)


########################################################################################
##############################     load data     #######################################
########################################################################################

load("results/simdata.RData")
load("results/missing_MAR_list.RData")

load("results/xgb_ParamRandom_maxit5_P_param.RData")
load("results/xgb_ParamAll_maxit5_P_param.RData")


###############################    3.4 - XGBoost - default parameter  ##############################################


########################## match.type = "predicted" - default ############################

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


xgb_ParamDefault_maxit1_P_imp_res <- map(impute_MAR_xgb_default, ~ map(., "mice_result"))  # imputation results
xgb_ParamDefault_maxit1_P_imp_time <- map(impute_MAR_xgb_default, ~ map(., "time_taken"))           # Time taken


xgb_ParamDefault_maxit1_P_eval <- xgb_ParamDefault_maxit1_P_imp_res %>% 
  map(function(mat_list) {
    furrr::future_map(mat_list, function(mat) {
      complete(mat, "all") %>% # create a list of completed data sets
        map(~.x %$% # for every completed data set....
              lm(y ~ x + z + x:z + I(z^2))) %>% # fit linear model
        pool() %>%  # pool coefficients
        summary(conf.int = TRUE) %>% # summary of coefficients
        mutate(true = true_vals, # add true
               cov = conf.low < true & true < conf.high, # coverage
               bias = estimate - true,
               width = conf.high - conf.low) %>% # bias
        column_to_rownames("term")
    }, .options = furrr_options(seed = TRUE)) %>% # `term` as rownames
      Reduce("+", .) / Num_ds})

save(xgb_ParamDefault_maxit1_P_imp_res, file = "results/xgb_ParamDefault_maxit1_P_imp_res.RData")
save(xgb_ParamDefault_maxit1_P_imp_time, file = "results/xgb_ParamDefault_maxit1_P_imp_time.RData")
save(xgb_ParamDefault_maxit1_P_eval, file = "results/xgb_ParamDefault_maxit1_P_eval.RData")

#rm(mice_xgb_default,time_xgb_default,eval_xgb_default,impute_MAR_xgb_default)

train_testplot(simdata, xgb_ParamDefault_maxit1_P_imp_res,prop_values, "results/TT_xgb_ParamDefault_maxit1_P_imp_res.pdf")



########################## match.type = "predicted.observed" - default param ############################
plan(sequential)
set.seed(seed)
plan(multisession, workers = available_cores)

impute_MAR_xgb_default_PO <- missing_MAR_list %>% # For each missingness
  map(function(mat_list) {
    furrr::future_map(mat_list, function(mat) {  #For each dataset
      result <- system.time({
        mice_result <- mice::mice(mat, 
                                  m = m, method = "xgb",
                                  maxit = maxit,xgb.params=NULL, match.type = "predicted.observed", # will use default parameters
                                  print = FALSE)
      })
      list(mice_result = mice_result, time_taken = result) # Combining both result and time taken
    }, .options = furrr_options(seed = TRUE))
  })


xgb_ParamDefault_maxit1_PO_imp_res <- map(impute_MAR_xgb_default_PO, ~ map(., "mice_result"))  # imputation results
xgb_ParamDefault_maxit1_PO_imp_time <- map(impute_MAR_xgb_default_PO, ~ map(., "time_taken"))           # Time taken


xgb_ParamDefault_maxit1_PO_eval <- xgb_ParamDefault_maxit1_PO_imp_res %>% 
  map(function(mat_list) {
    furrr::future_map(mat_list, function(mat) {
      complete(mat, "all") %>% # create a list of completed data sets
        map(~.x %$% # for every completed data set....
              lm(y ~ x + z + x:z + I(z^2))) %>% # fit linear model
        pool() %>%  # pool coefficients
        summary(conf.int = TRUE) %>% # summary of coefficients
        mutate(true = true_vals, # add true
               cov = conf.low < true & true < conf.high, # coverage
               bias = estimate - true,
               width = conf.high - conf.low) %>% # bias
        column_to_rownames("term")
    }, .options = furrr_options(seed = TRUE)) %>% # `term` as rownames
      Reduce("+", .) / Num_ds})

save(xgb_ParamDefault_maxit1_PO_imp_res, file = "results/xgb_ParamDefault_maxit1_PO_imp_res.RData")
save(xgb_ParamDefault_maxit1_PO_imp_time, file = "results/xgb_ParamDefault_maxit1_PO_imp_time.RData")
save(xgb_ParamDefault_maxit1_PO_eval, file = "results/xgb_ParamDefault_maxit1_PO_eval.RData")

#rm(mice_xgb_default,time_xgb_default,eval_xgb_default,impute_MAR_xgb_default)

train_testplot(simdata, xgb_ParamDefault_maxit1_PO_imp_res,prop_values, "results/TT_xgb_ParamDefault_maxit1_PO_imp_res.pdf")


########################## match.type = "original.observed" - default param ############################
plan(sequential)
set.seed(seed)
plan(multisession, workers = available_cores)

impute_MAR_xgb_default_OO <- missing_MAR_list %>% # For each missingness
  map(function(mat_list) {
    furrr::future_map(mat_list, function(mat) {  #For each dataset
      result <- system.time({
        mice_result <- mice::mice(mat, 
                                  m = m, method = "xgb",
                                  maxit = maxit,xgb.params=NULL, match.type = "original.observed", # will use default parameters
                                  print = FALSE)
      })
      list(mice_result = mice_result, time_taken = result) # Combining both result and time taken
    }, .options = furrr_options(seed = TRUE))
  })


xgb_ParamDefault_maxit1_OO_imp_res <- map(impute_MAR_xgb_default_OO, ~ map(., "mice_result"))  # imputation results
xgb_ParamDefault_maxit1_OO_imp_time <- map(impute_MAR_xgb_default_OO, ~ map(., "time_taken"))           # Time taken


xgb_ParamDefault_maxit1_OO_eval <- xgb_ParamDefault_maxit1_OO_imp_res %>% 
  map(function(mat_list) {
    furrr::future_map(mat_list, function(mat) {
      complete(mat, "all") %>% # create a list of completed data sets
        map(~.x %$% # for every completed data set....
              lm(y ~ x + z + x:z + I(z^2))) %>% # fit linear model
        pool() %>%  # pool coefficients
        summary(conf.int = TRUE) %>% # summary of coefficients
        mutate(true = true_vals, # add true
               cov = conf.low < true & true < conf.high, # coverage
               bias = estimate - true,
               width = conf.high - conf.low) %>% # bias
        column_to_rownames("term")
    }, .options = furrr_options(seed = TRUE)) %>% # `term` as rownames
      Reduce("+", .) / Num_ds})

save(xgb_ParamDefault_maxit1_OO_imp_res, file = "results/xgb_ParamDefault_maxit1_OO_imp_res.RData")
save(xgb_ParamDefault_maxit1_OO_imp_time, file = "results/xgb_ParamDefault_maxit1_OO_imp_time.RData")
save(xgb_ParamDefault_maxit1_OO_eval, file = "results/xgb_ParamDefault_maxit1_OO_eval.RData")

#rm(mice_xgb_default,time_xgb_default,eval_xgb_default,impute_MAR_xgb_default)

train_testplot(simdata, xgb_ParamDefault_maxit1_OO_imp_res,prop_values, "results/TT_xgb_ParamDefault_maxit1_OO_imp_res.pdf")



########################    3.5 - XGBoost - Random parameter ####################################



########################## match.type = "predicted" - Random ############################

plan(sequential)
set.seed(seed)
plan(multisession, workers = available_cores)


xgb_ParamRandom_maxit1_P_imp_res_tmp <- future_map2(missing_MAR_list, random_param_set_parameter_maxit5, 
                                                    .f = function(data_inner, params_inner) {
                                                      future_map2(data_inner, params_inner, 
                                                                  .f = function(data_single, params_single) {
                                                                    result <- system.time({
                                                                      mice_result <- mice(data_single, m = m, method = "xgb", maxit = maxit,
                                                                                          xgb.params =  params_single$parameter, print = FALSE)
                                                                    })
                                                                    list(mice_result = mice_result, time_taken = result)
                                                                  }, .options = furrr_options(seed = TRUE))
                                                    }, .options = furrr_options(seed = TRUE))

xgb_ParamRandom_maxit1_P_imp_res <- map(xgb_ParamRandom_maxit1_P_imp_res_tmp, ~ map(., "mice_result"))  # imputation results
xgb_ParamRandom_maxit1_P_imp_time <- map(xgb_ParamRandom_maxit1_P_imp_res_tmp, ~ map(., "time_taken"))           # Time taken


xgb_ParamRandom_maxit1_P_imp_eval <- xgb_ParamRandom_maxit1_P_imp_res %>% 
  map(function(mat_list) {
    furrr::future_map(mat_list, function(mat) {
      complete(mat, "all") %>% # create a list of completed data sets
        map(~.x %$% # for every completed data set....
              lm(y ~ x + z + x:z + I(z^2))) %>% # fit linear model
        pool() %>%  # pool coefficients
        summary(conf.int = TRUE) %>% # summary of coefficients
        mutate(true = true_vals, # add true
               cov = conf.low < true & true < conf.high, # coverage
               bias = estimate - true,
               width = conf.high - conf.low) %>% # bias
        column_to_rownames("term")
    }) %>% # `term` as rownames
      Reduce("+", .) / Num_ds})

save(xgb_ParamRandom_maxit1_P_imp_res, file = "results/xgb_ParamRandom_maxit1_P_imp_res.RData")
save(xgb_ParamRandom_maxit1_P_imp_time, file = "results/xgb_ParamRandom_maxit1_P_imp_time.RData")
save(xgb_ParamRandom_maxit1_P_imp_eval, file = "results/xgb_ParamRandom_maxit1_P_imp_eval.RData")

train_testplot(simdata,xgb_ParamRandom_maxit1_P_imp_res, prop_values, "results/TT_xgb_ParamRandom_maxit1_P_imp_res.pdf")


########################## match.type = "predicted.observed" - Random ############################

plan(sequential)
set.seed(seed)
plan(multisession, workers = available_cores)

xgb_ParamRandom_maxit1_PO_imp_res_tmp <- future_map2(missing_MAR_list, random_param_set_parameter_maxit5, 
                                                     .f = function(data_inner, params_inner) {
                                                       future_map2(data_inner, params_inner, 
                                                                   .f = function(data_single, params_single) {
                                                                     result <- system.time({
                                                                       mice_result <- mice(data_single, m = m, method = "xgb", 
                                                                                           maxit = maxit,xgb.params =  params_single$parameter, 
                                                                                           match.type = "predicted.observed", print = FALSE)
                                                                     })
                                                                     list(mice_result = mice_result, time_taken = result)
                                                                   }, .options = furrr_options(seed = TRUE))
                                                     }, .options = furrr_options(seed = TRUE))

xgb_ParamRandom_maxit1_PO_imp_res <- map(xgb_ParamRandom_maxit1_PO_imp_res_tmp, ~ map(., "mice_result"))  # imputation results
xgb_ParamRandom_maxit1_PO_imp_time <- map(xgb_ParamRandom_maxit1_PO_imp_res_tmp, ~ map(., "time_taken"))  # Time taken


xgb_ParamRandom_maxit1_PO_imp_eval <- xgb_ParamRandom_maxit1_PO_imp_res %>% 
  furrr::future_map(function(mat_list) {
    map(mat_list, function(mat) {
      complete(mat, "all") %>% # create a list of completed data sets
        map(~.x %$% # for every completed data set....
              lm(y ~ x + z + x:z + I(z^2))) %>% # fit linear model
        pool() %>%  # pool coefficients
        summary(conf.int = TRUE) %>% # summary of coefficients
        mutate(true = true_vals, # add true
               cov = conf.low < true & true < conf.high, # coverage
               bias = estimate - true, # bias
               width = conf.high - conf.low) %>% 
        column_to_rownames("term")
    }) %>% # `term` as rownames
      Reduce("+", .) / Num_ds})

save(xgb_ParamRandom_maxit1_PO_imp_res, file = "results/xgb_ParamRandom_maxit1_PO_imp_res.RData")
save(xgb_ParamRandom_maxit1_PO_imp_time, file = "results/xgb_ParamRandom_maxit1_PO_imp_time.RData")
save(xgb_ParamRandom_maxit1_PO_imp_eval, file = "results/xgb_ParamRandom_maxit1_PO_imp_eval.RData")

train_testplot(simdata, xgb_ParamRandom_maxit1_PO_imp_res, prop_values, "results/TT_xgb_ParamRandom_maxit1_PO_imp_res.pdf")

########################## match.type = "original.observed" - Random ############################

plan(sequential)
set.seed(seed)
plan(multisession, workers = available_cores)

xgb_ParamRandom_maxit1_OO_imp_res_tmp <- future_map2(missing_MAR_list, random_param_set_parameter_maxit5, 
                                                     .f = function(data_inner, params_inner) {
                                                       future_map2(data_inner, params_inner, 
                                                                   .f = function(data_single, params_single) {
                                                                     result <- system.time({
                                                                       mice_result <- mice(data_single, m = m, method = "xgb",
                                                                                           maxit = maxit,xgb.params =  params_single$parameter, 
                                                                                           match.type = "original.observed", print = FALSE)
                                                                     })
                                                                     list(mice_result = mice_result, time_taken = result)
                                                                   }, .options = furrr_options(seed = TRUE))
                                                     }, .options = furrr_options(seed = TRUE))

xgb_ParamRandom_maxit1_OO_imp_res <- map(xgb_ParamRandom_maxit1_OO_imp_res_tmp, ~ map(., "mice_result"))  # imputation results
xgb_ParamRandom_maxit1_OO_imp_time <- map(xgb_ParamRandom_maxit1_OO_imp_res_tmp, ~ map(., "time_taken"))           # Time taken


xgb_ParamRandom_maxit1_OO_imp_eval <- xgb_ParamRandom_maxit1_OO_imp_res %>% 
  map(function(mat_list) {
    furrr::future_map(mat_list, function(mat) {
      complete(mat, "all") %>% # create a list of completed data sets
        map(~.x %$% # for every completed data set....
              lm(y ~ x + z + x:z + I(z^2))) %>% # fit linear model
        pool() %>%  # pool coefficients
        summary(conf.int = TRUE) %>% # summary of coefficients
        mutate(true = true_vals, # add true
               cov = conf.low < true & true < conf.high, # coverage
               bias = estimate - true,
               width = conf.high - conf.low) %>% # bias
        column_to_rownames("term")
    }) %>% # `term` as rownames
      Reduce("+", .) / Num_ds})

save(xgb_ParamRandom_maxit1_OO_imp_res, file = "results/xgb_ParamRandom_maxit1_OO_imp_res.RData")
save(xgb_ParamRandom_maxit1_OO_imp_time, file = "results/xgb_ParamRandom_maxit1_OO_imp_time.RData")
save(xgb_ParamRandom_maxit1_OO_imp_eval, file = "results/xgb_ParamRandom_maxit1_OO_imp_eval.RData")

train_testplot(simdata, xgb_ParamRandom_maxit1_OO_imp_res,prop_values, "results/TT_xgb_ParamRandom_maxit1_OO_imp_res.pdf")





########################    3.6 - XGBoost - All parameter ####################################


########################## match.type = "predicted" - All ############################

plan(sequential)
set.seed(seed)
plan(multisession, workers = available_cores)


xgb_ParamAll_maxit1_P_imp_res_tmp <- future_map2(missing_MAR_list, all_param_set_parameter_maxit5, 
                                                 .f = function(data_inner, params_inner) {
                                                   future_map2(data_inner, params_inner, 
                                                               .f = function(data_single, params_single) {
                                                                 result <- system.time({
                                                                   mice_result <- mice(data_single, m = m, method = "xgb", 
                                                                                       maxit = maxit,xgb.params =  params_single$parameter, 
                                                                                       print = FALSE)
                                                                 })
                                                                 list(mice_result = mice_result, time_taken = result)
                                                               }, .options = furrr_options(seed = TRUE))
                                                 }, .options = furrr_options(seed = TRUE))

xgb_ParamAll_maxit1_P_imp_res <- map(xgb_ParamAll_maxit1_P_imp_res_tmp, ~ map(., "mice_result"))  # imputation results
xgb_ParamAll_maxit1_P_imp_time <- map(xgb_ParamAll_maxit1_P_imp_res_tmp, ~ map(., "time_taken"))           # Time taken


xgb_ParamAll_maxit1_P_imp_eval <- xgb_ParamAll_maxit1_P_imp_res %>% 
  map(function(mat_list) {
    furrr::future_map(mat_list, function(mat) {
      complete(mat, "all") %>% # create a list of completed data sets
        map(~.x %$% # for every completed data set....
              lm(y ~ x + z + x:z + I(z^2))) %>% # fit linear model
        pool() %>%  # pool coefficients
        summary(conf.int = TRUE) %>% # summary of coefficients
        mutate(true = true_vals, # add true
               cov = conf.low < true & true < conf.high, # coverage
               bias = estimate - true,
               width = conf.high - conf.low) %>% # bias
        column_to_rownames("term")
    }) %>% # `term` as rownames
      Reduce("+", .) / Num_ds})

save(xgb_ParamAll_maxit1_P_imp_res, file = "results/xgb_ParamAll_maxit1_P_imp_res.RData")
save(xgb_ParamAll_maxit1_P_imp_time, file = "results/xgb_ParamAll_maxit1_P_imp_time.RData")
save(xgb_ParamAll_maxit1_P_imp_eval, file = "results/xgb_ParamAll_maxit1_P_imp_eval.RData")

train_testplot(simdata,xgb_ParamAll_maxit1_P_imp_res, prop_values, "results/TT_xgb_ParamAll_maxit1_P_imp_res.pdf")


########################## match.type = "predicted.observed" - All ############################

plan(sequential)
set.seed(seed)
plan(multisession, workers = available_cores)

xgb_ParamAll_maxit1_PO_imp_res_tmp <- future_map2(missing_MAR_list, all_param_set_parameter_maxit5, 
                                                  .f = function(data_inner, params_inner) {
                                                    future_map2(data_inner, params_inner, 
                                                                .f = function(data_single, params_single) {
                                                                  result <- system.time({
                                                                    mice_result <- mice(data_single, m = m, method = "xgb", 
                                                                                        maxit = maxit,xgb.params =  params_single$parameter, 
                                                                                        match.type = "predicted.observed", print = FALSE)
                                                                  })
                                                                  list(mice_result = mice_result, time_taken = result)
                                                                }, .options = furrr_options(seed = TRUE))
                                                  }, .options = furrr_options(seed = TRUE))

xgb_ParamAll_maxit1_PO_imp_res <- map(xgb_ParamAll_maxit1_PO_imp_res_tmp, ~ map(., "mice_result"))  # imputation results
xgb_ParamAll_maxit1_PO_imp_time <- map(xgb_ParamAll_maxit1_PO_imp_res_tmp, ~ map(., "time_taken"))  # Time taken


xgb_ParamAll_maxit1_PO_imp_eval <- xgb_ParamAll_maxit1_PO_imp_res %>% 
  furrr::future_map(function(mat_list) {
    map(mat_list, function(mat) {
      complete(mat, "all") %>% # create a list of completed data sets
        map(~.x %$% # for every completed data set....
              lm(y ~ x + z + x:z + I(z^2))) %>% # fit linear model
        pool() %>%  # pool coefficients
        summary(conf.int = TRUE) %>% # summary of coefficients
        mutate(true = true_vals, # add true
               cov = conf.low < true & true < conf.high, # coverage
               bias = estimate - true, # bias
               width = conf.high - conf.low) %>% 
        column_to_rownames("term")
    }) %>% # `term` as rownames
      Reduce("+", .) / Num_ds})

save(xgb_ParamAll_maxit1_PO_imp_res, file = "results/xgb_ParamAll_maxit1_PO_imp_res.RData")
save(xgb_ParamAll_maxit1_PO_imp_time, file = "results/xgb_ParamAll_maxit1_PO_imp_time.RData")
save(xgb_ParamAll_maxit1_PO_imp_eval, file = "results/xgb_ParamAll_maxit1_PO_imp_eval.RData")

train_testplot(simdata, xgb_ParamAll_maxit1_PO_imp_res, prop_values, "results/TT_xgb_ParamAll_maxit1_PO_imp_res.pdf")

########################## match.type = "original.observed" - All ############################

plan(sequential)
set.seed(seed)
plan(multisession, workers = available_cores)

xgb_ParamAll_maxit1_OO_imp_res_tmp <- future_map2(missing_MAR_list, all_param_set_parameter_maxit5, 
                                                  .f = function(data_inner, params_inner) {
                                                    future_map2(data_inner, params_inner, 
                                                                .f = function(data_single, params_single) {
                                                                  result <- system.time({
                                                                    mice_result <- mice(data_single, m = m, method = "xgb",
                                                                                        maxit = maxit,xgb.params =  params_single$parameter, 
                                                                                        match.type = "original.observed", print = FALSE)
                                                                  })
                                                                  list(mice_result = mice_result, time_taken = result)
                                                                }, .options = furrr_options(seed = TRUE))
                                                  }, .options = furrr_options(seed = TRUE))

xgb_ParamAll_maxit1_OO_imp_res <- map(xgb_ParamAll_maxit1_OO_imp_res_tmp, ~ map(., "mice_result"))  # imputation results
xgb_ParamAll_maxit1_OO_imp_time <- map(xgb_ParamAll_maxit1_OO_imp_res_tmp, ~ map(., "time_taken"))           # Time taken


xgb_ParamAll_maxit1_OO_imp_eval <- xgb_ParamAll_maxit1_OO_imp_res %>% 
  map(function(mat_list) {
    furrr::future_map(mat_list, function(mat) {
      complete(mat, "all") %>% # create a list of completed data sets
        map(~.x %$% # for every completed data set....
              lm(y ~ x + z + x:z + I(z^2))) %>% # fit linear model
        pool() %>%  # pool coefficients
        summary(conf.int = TRUE) %>% # summary of coefficients
        mutate(true = true_vals, # add true
               cov = conf.low < true & true < conf.high, # coverage
               bias = estimate - true,
               width = conf.high - conf.low) %>% # bias
        column_to_rownames("term")
    }) %>% # `term` as rownames
      Reduce("+", .) / Num_ds})

save(xgb_ParamAll_maxit1_OO_imp_res, file = "results/xgb_ParamAll_maxit1_OO_imp_res.RData")
save(xgb_ParamAll_maxit1_OO_imp_time, file = "results/xgb_ParamAll_maxit1_OO_imp_time.RData")
save(xgb_ParamAll_maxit1_OO_imp_eval, file = "results/xgb_ParamAll_maxit1_OO_imp_eval.RData")

train_testplot(simdata, xgb_ParamAll_maxit1_OO_imp_res,prop_values, "results/TT_xgb_ParamAll_maxit1_OO_imp_res.pdf")



###########################################################################################################
###########################################################################################################
###########################################################################################################



