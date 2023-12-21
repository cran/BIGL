#' Functions for scaling, and rescaling residuals. May lead to unstable behaviour in practice
#' @details Residuals are calculated with respect to the average observation on
#' the off-axis point, so replicates are required!
#' @param sampling_errors A vector of raw  residuals
#' @param ... passed on to predictVar
scaleResids = function(sampling_errors, ...){
    predVar = predictVar(...)
    sampling_errors/sqrt(predVar)
}
#' Backscale residuals
#' @param scaledResids scaled residuals
#' @inheritParams scaleResids
backscaleResids = function(scaledResids, ...){
    predVar = predictVar(...)
    scaledResids*sqrt(predVar)
}
#'Predict variance
#' @param means a vector of means
#' @inheritParams generateData
predictVar = function(means, model, invTransFun){
    predVar = invTransFun(model[1] + model[2]*means)
    if(model["min"] == 0){
        predVar[predVar<=0] = 0.000001 #Correct for negative variances
    } else {
        predVar[predVar<=0] = model["min"] #Correct for negative variances
    }
    predVar[predVar > model["max"]] = model["max"] #Upper bound
    predVar
}


#' Add residuals by adding to mean effects
#' @inheritParams scaleResids
#' @inheritParams predictVar
addResids = function(means, ...){
    means + sampleResids(means, ...)
}
#' Sample residuals according to a new model
#' @inheritParams fitSurface
#' @inheritParams predictVar
#' @inheritParams scaleResids
#' @return sampled residuals
sampleResids = function(means, sampling_errors, method, rescaleResids,...){
    if(method %in% c("equal", "unequal")){
        return(sample(sampling_errors, size = length(means), replace = TRUE))
    } else if(method == "model"){
        resids = if(rescaleResids){
            scaledResids = scaleResids(sampling_errors, means, ...)
            sampledResids = sample(scaledResids, replace = TRUE)
            backscaleResids(sampledResids, means, ...)
        } else{
            rnorm(length(means), sd = sqrt(predictVar(means, ...)))
        }
        return(resids)
    } else{}
}

#' Sample residuals according to a new model
#' @inheritParams fitSurface
#' @inheritParams predictVar
#' @inheritParams scaleResids
#' @inheritParams bootConfInt
#' @importFrom stats rgamma
#' @return sampled residuals
wildbootAddResids <- function(means, sampling_errors, method, rescaleResids, model, invTransFun, wild_bootstrap, wild_bootType,...){
  if(wild_bootstrap){
      errors = switch(wild_bootType,
                      
                      # Rademacher
                      "rademacher" = {sampling_errors*(2*rbinom(length(means), size = 1, prob = 0.5)-1)},          # Rademacher distribution
                      "gamma"  =  {sampling_errors*(rgamma(length(means),shape = 4, scale = 0.5)-2)},              # Gamma distribution
                      
                      # Normal
                     "normal" = {noff <- length(means)
                      mu1 <- 0.5*(sqrt(17/6)+sqrt(1/6))
                      mu2 <- 0.5*(sqrt(17/6)-sqrt(1/6))
                      W1 = rnorm(noff,mu1,sqrt(0.5))
                      Z1 = rnorm(noff,mu2,sqrt(0.5))
                      sampling_errors*(W1*Z1-mu1*mu2)},                                                             # Normal distribution
                      
                      # Two-point distribution by Mammen (1993)
                     
                      "two-point" = {vals <- c(-(sqrt(5)-1)/2, (sqrt(5)+1)/2)
                      probs <- rev(abs(vals)/sqrt(5))
                      sampling_errors*sample(vals, size = length(means), replace = TRUE, prob = probs)})
  }else{
      errors <-  sampleResids(means, sampling_errors, method, rescaleResids, model, invTransFun,...)
     }
  means+errors
}

