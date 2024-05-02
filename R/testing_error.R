#' Calculation of testing_error. This function can be used on the output of mice.impute.xgb. 
#'
#'@aliases testing_error
#' @param orig_data_1: original data without missing values
#' @param imputed_res_1: imputation results from mice.impute.xgb
#' @return test error list of all variables in the dataset. Each list contain m values for m different imputations
#' 
#'testing_error
#' @examples
#' testing_error(nhanes, imp_res)
#'
#' @export
#' 

testing_error <- function(orig_data_1 , imputed_res_1){
  
  
  imputed_data_1 <- imputed_res_1$imp
  missing_data_1 <- imputed_res_1$data
  
  if(!(identical(names(imputed_data_1), colnames(orig_data_1)))){
    stop("Original and imputed data has different variables. Testing error can not be calculated")
  }
  
  testError_list = list()
  
  var_names <- colnames(orig_data_1)
  
  
  for(var in var_names){
    missing_index <- which(is.na(missing_data_1[,var]))
    actual_data <- orig_data_1[missing_index, var]
    
    imputed_data_varIndex <- match(var, names(imputed_data_1))
    imputed_data_oneVar <- imputed_data_1[[imputed_data_varIndex]]
    
    
    testing_error<-as.data.frame(apply(imputed_data_oneVar, 2, function(x) {actual_data - x}))
    mean_testing_error <- colMeans(testing_error)
    
    testError_list[[var]] <- matrix(mean_testing_error, nrow = 1)
    
  }
  
  testError_list
}

