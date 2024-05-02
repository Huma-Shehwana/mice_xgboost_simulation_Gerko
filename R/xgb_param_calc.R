#' Creates a \code{xgb_param_calc} argument
#'
#' This helper function creates a valid \code{xgb_param_calc} vector. The
#' \code{xgb_param_calc} vector is an argument to the \code{mice} function that
#' specifies post-processing for a variable after each iteration of imputation.
#' @inheritParams mice
#' @return Character vector of \code{ncol(data)} element
#' @seealso \code{\link{mice}}
#' @examples
#' make.post(nhanes2)
#' @export
#' @import ParamHelpers
#' @import mlrMBO
#' @import xgboost
#' @import ggplot2

xgb_param_calc<-function(data,response = NULL, select_features=NULL, num_cores = 1, iter=30){
  
  
  ##############################################################################
  
  data <- as.data.frame(data)
  num.cc <- sum(complete.cases(data))
  
  if (num.cc == 0) {
    stop("No complete cases in this dataset.")
  }
  
  
  if (num.cc > 0 & num.cc < 20) {
    stop("Less than 20 complete cases. Hyperparameter tuning is not possible")
  }
  
  cc.data <- data[complete.cases(data), ] # only extract complete data. Complete cases will be used for hyperparameter tuning.
  
  Names <- colnames(cc.data) # all feature names
  Types <- as.vector(apply(cc.data,2,class)) #type of each variable
  
  na.col <- which(colSums(is.na(data))!=0) #find features with missing values
  
  if(length(na.col)==0){
    stop("There is no varaible with missing data. Parameter optimization will not be performed")
  }
  
  ################################################################################################
  #################           Target and training features          ##############################
  ################################################################################################
  
  
  # if user has not chosen from the given options, following message will be displayed
  
  if(!is.null(select_features) &!identical(select_features, "all") & !all(select_features %in% seq_along(Names))){
    stop("Please check the format of select_features. It can be either NULL, all or index of columns that you want to include as a feature")
  }
  
  if(!is.null(response) &!identical(response, "all") & !all(response %in% seq_along(Names))){
    stop("Please check the format of response variable. It can be either NULL, all or index of columns that you want to include as a response variable")
  }
  
  # response = NULL, different options for select_features
  
  if (is.null(response)) { 
    
    if(is.null(select_features) | identical(select_features, "all")) {
      r.idx <- sample(na.col, size = 1) #randomly sample one feature from columns with missing data and assign it to response_var in next line
      response_var <- Names[r.idx]
      training_features <- setdiff(Names, response_var) # assign all other features (with/without missing data) to training features set
      
    }else if(is.numeric(select_features)){
      training_features <- Names[select_features] #features are based on user given indexes
      response_var <- sample(setdiff(Names[na.col], training_features), size = 1) # response variable is randomly sampled from missing variables other than select_features
    } 
    
  } else if (is.numeric(response)) {
    if(length(response)>1){ 
      message("You have selected more than one variables for parameter optimization. Please note that only the first index will be used for parameter estimation")
    } 
    
    response_var <- Names[response[1]] # in case of multiple response variables, select first feature from given response variables
    
    if(is.null(select_features) | identical(select_features, "all")) {
      training_features <- setdiff(Names, response_var) # all variables other than response variables
      
    } else if(is.numeric(select_features)){
      if(response %in% select_features){
        message("Response and feature variables should not be same. All variables other than response variable will be used as features")
        training_features <- setdiff(Names, response_var)
      } else {
        training_features <- Names[select_features]
      }
    } 
  } else if(response =="all"){
    response_var <- Names[na.col]
    training_features <- NULL
    message("You have chosen all variables (with missing values) to be iteratively used as a target variable for parameter optimization. Please note that select_features will be automatically chosen for this process.")
  } else {
    stop("Please check the format of response variable It can be either NULL, \"all\" or index of column that you want to use for parameter optimization")
  }
  
  params = list()
  
  if(length(response_var)>1){
    
    response_var <- setNames(response_var, response_var)
    message("Performing bayesian optimization iteratively using each vairable as response variable")
    params <- lapply(response_var,function(x) par_opt(cc.data, setdiff(Names, x), x, iter = iter)) # call to bayesian optimization using assigned features and response variable
    
  } else {
    
    message("Performing bayesian optimization using ", response_var, " as a response variable and ", training_features, " as features")
    params = par_opt(cc.data, training_features, response_var,iter = iter)
  }
  
  
  return(params)
  
}

###########################################################################################################################
########################              Bayesian optimization function            ###########################################
###########################################################################################################################


par_opt <- function(cc.data, training_features, response_var, nfold = 5, iter = iter, early_stopping_rounds = 20, nround = 100){
  
  response_data <- cc.data[,response_var] # data from response variable 
  
  if(class(response_data)=="character"){
    response_data <- as.factor(response_data)
  } # for next steps in code
  
  response_class <- class(response_data)
  
  
  
  ######################################################################################################
  ############             Objective and evaluation function definition             ####################
  ######################################################################################################
  
  var.type <- switch(response_class[1],
                     numeric = "numeric",
                     integer = "numeric",
                     factor = ifelse(nlevels(response_data) == 2, "binary", "multiclass"),
                     ordered = ifelse(nlevels(response_data) == 2, "binary", "multiclass"),
                     stop(paste0("Unsupported data type: ", response_class[1]))
  )
  
  obj.type <- switch(var.type,
                     numeric = "reg:squarederror",
                     binary =  "binary:logistic",
                     multiclass = "multi:softmax",
                     stop(paste0("Unsupported objective function: ", var.type))
  )
  
  eval.metric <- switch(var.type,
                        numeric = "mape",
                        binary = "logloss",
                        multiclass = "mlogloss",
                        stop(paste0("Unsupported data type: ", obj.type))
  )
  
  ####################################################################################################
  #conditions for categorical/multi-level data
  
  continue = TRUE
  
  if(var.type != "numeric"){
    response_data <- as.numeric(unlist(response_data)) - 1 # convert categorical data to integers because xgboost can only deal with numbers and integers as labels
    bin.t <- sort(table(response_data), decreasing = TRUE) # distribution of categories
    
    if (is.na(bin.t[2])) {
      msg <- paste("The variable", var, "in the data only have single class. Optimization algorithm can not run for variable. Thus default aprameters are assigned to ", var)
      message(msg)
      continue = FALSE
    }
  }
  
  ###############################################################################################
  ##########              XGBoost object and cross validation
  ###############################################################################################
  
  p <- length(training_features) + 1
  
  if (p == 2) { # convert training data to sparse matrix
    obs.data <- Matrix::sparse.model.matrix(reformulate(training_features, response_var), data = cc.data)
  } else {
    obs.data <- Matrix::sparse.model.matrix(reformulate(training_features, response_var), data = cc.data)[, -1]
  }
  
  dtrain <- xgboost::xgb.DMatrix(obs.data, label = response_data, missing = NA)
  
  if(var.type=="numeric"){
    
    cv_folds = rBayesianOptimization::KFold(response_data, # cross validation 
                                            nfolds= nfold,
                                            seed= TRUE)
    
  } else {
    
    cv_folds = rBayesianOptimization::KFold(response_data, 
                                            nfolds= nfold, stratified = TRUE, #so that distribution of different classes is similar in all folds
                                            seed= TRUE)
    
  }
  
  if(any(lapply(cv_folds, length)<5)){
    msg <- paste("Hyperparameter tuning is not possible for variable", response_var, 
                 " because it does not have sufficient data-points for cross-validation. Default parameters are being assigned to ", response_var)
    message(msg)
    continue = FALSE
  }
  
  ###############################################################################################
  ####################          Single objective function       ################################
  ###############################################################################################
  
  obj.fun  <- smoof::makeSingleObjectiveFunction(
    name = "xgb_cv_bayes",
    fn =   function(x){
      if(var.type!="multiclass"){
        cv <- xgb.cv(params = list(
          booster = "gbtree",
          eta = x["eta"],
          max_depth = x["max_depth"],
          min_child_weight = x["min_child_weight"],
          gamma = x["gamma"],
          alpha = x["alpha"],
          lambda = x["lambda"],
          objective = obj.type, 
          eval_metric = eval.metric),
          data = dtrain,
          nround = nround,
          early_stopping_rounds = early_stopping_rounds,
          folds = cv_folds,
          prediction = FALSE,
          maximize = FALSE,
          showsd = FALSE,
          verbose = 0)
      } else {
        N <- length(levels(as.factor(cc.data[,response_var]))) # number of classes
        cv <- xgb.cv(params = list(
          booster = "gbtree",
          eta = x["eta"],
          max_depth = x["max_depth"],
          num_class = N,     
          min_child_weight = x["min_child_weight"],
          gamma = x["gamma"],
          alpha = x["alpha"],
          lambda = x["lambda"],
          objective = obj.type, 
          eval_metric = eval.metric),
          data = dtrain,
          nround = nround,
          early_stopping_rounds = early_stopping_rounds,
          folds = cv_folds,
          prediction = FALSE,
          maximize = FALSE,
          showsd = FALSE,
          verbose = 0)
      }
      eval_column <- names(cv$evaluation_log)[grep("test.*mean", names(cv$evaluation_log))] #extract column containing test metrics
      min(cv$evaluation_log[, eval_column]) #minimum test metrics
    },
    par.set = makeParamSet(  #setting bounds of parameter which are being optimized
      makeNumericParam("eta",lower = 0.0001, upper = 1),
      makeNumericParam("alpha",lower = 0, upper = 100),
      makeNumericParam("lambda",lower = 0, upper = 100),
      makeNumericParam("gamma",lower = 0, upper = 10),
      makeIntegerParam("max_depth",lower= 1,upper = 20),
      makeIntegerParam("min_child_weight", lower= 0,upper = 50)
    ),
    minimize = TRUE # minimize the evaluation metric
  ) #end of single objective function
  
  
  #################################################################################################################
  ###################          Call to objective function using bayesian optimization      ########################
  #################################################################################################################
  
  control = makeMBOControl()
  
  des = generateDesign(n=10L,
                       par.set = getParamSet(obj.fun), 
                       fun = lhs::randomLHS) 
  
  control = setMBOControlTermination(control, iters = iter)
  
  
  if(continue==TRUE){
    run = mbo(fun = obj.fun, 
              control = control, 
              design = des)
    
    # plot of iterations on x-axis and evaluation metric on x-axis
    plot_fig<-run$opt.path$env$path  %>% 
      mutate(Round = row_number()) %>% ggplot(aes(x= Round, y= y)) + 
      geom_point() +
      labs(title = sprintf("Response Variable: %s , Features: %s", response_var, paste(training_features,collapse = ",")))+
      ylab(sprintf("Evaluation score %s", eval.metric)) + theme(plot.title = element_text(hjust = 0.5))
    
    
    best_solution <- run$opt.path$env$path[which.min(run$opt.path$env$path$y),]
    best_eval <- best_solution$y
    best_parameters <- run$x
  } else {
    
    best_parameters <- list("eta" = 0.3, "lambda" = 1, "alpha" = 0, "gamma" = 0, "max_depth" = 3, "min_child_weight" = 1)
    plot_fig <- NA
    
  }
  
  list("parameter" = best_parameters, "fig" = plot_fig)
} # end of bayesian optimization function






