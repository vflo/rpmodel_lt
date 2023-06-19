#' General method for trainnig NN models. 
#' 
#' @param data A data frame containing observational data for all the 
#' predictors and training variables with all NAs removed.
#' 
#' @param predictors A character string defining which variables 
#' (column name in \code{df}) are used as predictors variables.
#' 
#' @param target A character string defining which variable 
#' (column name in \code{df}) is to be used as target variable.
#' 
sfn_predict_nn <- function( data, predictors, nam_target, weights = NULL, nn = NULL, 
                        threshold = 0.03, do_predict = TRUE, do_modobs = FALSE, 
                        trainfrac = 1.0, package="nnet", 
                        seed = 1, hidden = NULL ){
  
  if (package=="nnet"){
    
    require( nnet )
    require( caret )
    
    if (is.null(nn)){
      
      forml  <- as.formula(  paste( nam_target, "~", paste( predictors, collapse=" + " ) ) )
      
      preprocessParams <- caret::preProcess( data, method=c("range") )
      traincotrlParams <- caret::trainControl( method="repeatedcv", 
                                               number=4, 
                                               repeats=4, 
                                               verboseIter=FALSE, 
                                               p=0.75 ) ## take best of 10 repetitions of 
      ## training with 75% used for training 
      ## (25% for testing)
      
      if (is.null(hidden)){
        tune_grid <- expand.grid( .decay = c(0.1), .size = seq(4,20,2))
      } else {
        tune_grid <- expand.grid( .decay = c(0.1), .size = c(hidden) )
      }
      
      set.seed(seed)
      nn <- caret::train(
        form        = forml,
        data        = data, #training,
        weights     = weights,
        method      = "nnet",
        linout      = TRUE,
        # tuneGrid    = tune_grid,
        preProcess  = 'range',
        trControl   = traincotrlParams,
        trace       = TRUE
      )
      
    }
    
    if (do_predict){
      vals <- as.vector( predict( nn, data ) )  # try( predict( nn, newdata=testing ) )
    } else {
      vals <- rep( NA, nrow(data) )
    }
    
  } else {
    
    rlang::abort("predict_nn(): No other training package implemented than nnet.")
    
  }
  
  return( list( nn=nn, vals=vals, hidden_best=nn$bestTune$size ) )
  
}
