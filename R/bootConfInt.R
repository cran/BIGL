#' Obtain confidence intervals for the raw effect sizes on every off-axis point and overall
#'
#' @param Total data frame with all effects and mean effects
#' @inheritParams fitSurface
#' @inheritParams meanR
#' @inheritParams generateData
#' @param posEffect a boolean, are effects restricted to be positive
#' @param respS the observed response surface
#' @return A list with components
#' \item{offAxis}{The off-axis bootstrapped confidence intervals}
#' \item{single}{A mean effect and percentile and studentized boostrap intervals}
bootConfInt = function(Total, idUnique, bootStraps,
    transforms, respS, B.B, method,
    CP, reps, n1, cutoff, R, fitResult,
    bootRS, data_off, posEffect = all(Total$effect >= 0),
    transFun, invTransFun, model, rescaleResids, wild_bootstrap, wild_bootType, control, digits, ...) {
  Total <- Total[Total$d1 & Total$d2, ]
  sampling_errors <- Total$effect - Total$meaneffect								# sampling errors for off Axis points
  A <- getA(data_off, fitResult, method, CP, reps, n1, transFun = transFun,
      invTransFun = invTransFun)
  bootEffectSizesList <- lapply(bootStraps, function(bb) {
        #Do use bootstrapped response surface for complete mimicry of variability
        dat_off_resam <- within(Total, {
              effect <- wildbootAddResids(Total$meaneffect, sampling_errors, method,
                          rescaleResids, model, invTransFun, wild_bootstrap, wild_bootType) 
              #Sample with replacement
              if (posEffect) {
                effect <- abs(effect)
              }
            })
        #Transforms have occurred in Total already
        bootR <- getR(data = dat_off_resam, idUnique = dat_off_resam$d1d2,
            transforms = NULL, respS = if(bootRS) bb$respS else respS)
        bootA <- getA(dat_off_resam, bb$simFit, method, CP, reps, n1, transFun = transFun,
            invTransFun = invTransFun)
        list("R" = bootR, "A" = bootA)
      })
  bootEffectSizes <- vapply(bootEffectSizesList, FUN.VALUE = c(R), function(x) x$R)
  bootAs <- vapply(bootEffectSizesList, FUN.VALUE = diag(A), function(x) sqrt(diag(x$A)))
  
  # specify here two ways of constructing confidence intervals
  # 1) simultaneous CI; i.e. controls FWER
  # 2) False coverage rate CI, i.e. controls false coverage proportion
  
  
  
  #Off axis confidence interval
  if (control == "FCR") {                                                                            # Control False coverage rate
    tgrid <- seq(0, 10, by = 0.01) 
    
    bootEffectSizesStand <- abs(bootEffectSizes-c(R))/bootAs                                     # Tb - standardized statistic
    tcount <- c()
    for (i in 1:ncol(bootEffectSizesStand)){
      count <- c()
      for (j in  1:length(tgrid)){
        tt <- sum(bootEffectSizesStand[,i]>tgrid[j])
        count <- c(count,tt)
      }
      tcount <- rbind(tcount,count)
    }
    prob <- apply(tcount,2,sum)/(length(R)*ncol(bootEffectSizesStand))
    
    id <- min(which(prob <= 1-cutoff))
    effectSizeQuant <- tgrid[id]                                                                    # t-alpha  we need for FCR
    
    confInt <- c(R) + outer(effectSizeQuant*sqrt(diag(A)),
                            c("lower" = -1, "upper" = 1))
    R <- round(R, digits = digits)                                                                     # round the result t two decimal places
    confInt <- round(confInt, digits = digits)
    rownames(confInt) <- rownames(bootEffectSizesStand)
    
  } else if (control == "dFCR"){                                                              # controlling directional false coverage rate
    tgrid <- seq(0, 5, by = 0.01) 
    bootEffectSizesStand <- abs(bootEffectSizes-c(R))/bootAs                             # Tb - standardized statistic
    EE.med<-median(abs(R))                                                               # estimated effect median value
    zz <- (abs(R)<EE.med)
    EE <- R
    EE[zz] <- 0                                                                            # Estimated effect size set to zero if smaller than median effect size
    Nb <- matrix(nrow = ncol(bootEffectSizes),ncol= length(tgrid))
    
    for (j in 1:ncol(bootEffectSizes)){
      thresh <- lapply(tgrid, function(i){
        low <- bootEffectSizes[,j]-i*bootAs[,j]
        upp <- bootEffectSizes[,j]+ i*bootAs[,j]
        res <- ((low>0)&(upp>0)&(EE<=0))|((low<0)&(upp<0)&(EE>=0))|((low<0)&(upp>0)&((EE<low)|(EE>upp)))          # directional false coverage 
        list("res" = res)
      })
      threshtot <-  sapply(thresh, function(y) y[["res"]])
      Nb[j,]<- apply(threshtot, 2, mean)
    }
    
    E.Nb <- apply(Nb,2,mean)
    t.alpha <- tgrid[which.min(abs(E.Nb-(1-cutoff)))]
    confInt <- c(R) + outer(t.alpha*sqrt(diag(A)),
                            c("lower" = -1, "upper" = 1))
    R <- round(R, digits = digits)                                                                   # round the result to two decimal places
    confInt <- round(confInt, digits = digits)
    rownames(confInt) <- rownames(bootEffectSizesStand) 
    
  } else {                                                                                         # default controls FWER
    
    # Off axis confidence interval, control FWER
    bootEffectSizesStand <- abs(bootEffectSizes-c(R))/bootAs
    maxEffectSizes <- apply(bootEffectSizesStand, 2, max)
    effectSizeQuant <- quantile(maxEffectSizes, cutoff, na.rm = TRUE)
    confInt <- c(R) + outer(effectSizeQuant*sqrt(diag(A)),
                            c("lower" = -1, "upper" = 1))
    R <- round(R, digits = digits)                                                                     # round the result t two decimal places
    confInt <- round(confInt, digits = digits)
    
    rownames(confInt) <- rownames(bootEffectSizesStand)
  }
  
  coefFit <- fitResult$coef
  eq  <- coefFit["m1"] == coefFit["b"] && coefFit["m2"] == coefFit["b"]
  inc <- coefFit["m1"] >= coefFit["b"] && coefFit["m2"] >= coefFit["b"]
  dec <- coefFit["m1"] <= coefFit["b"] && coefFit["m2"] <= coefFit["b"]
  
  call <- rep("None", length(R))
  call[confInt[, "lower"] > 0] <- if (eq) {
        "Undefined"
      } else if (inc) {
        "Syn"
      } else if (dec) {
        "Ant"
      } else 
        "Undefined"
  call[confInt[, "upper"] < 0] <- if (eq) {
        "Undefined"
      } else if (inc) {
        "Ant"
      } else if (dec) {
        "Syn"
      } else 
        "Undefined"
  confInt <- data.frame("estimate" = R, confInt, "call" = call)
  
  #Single measure of effect size
  singleMeasure <- mean(R)
  bootR <- colMeans(bootEffectSizes)
  bootRstand <- (bootR-singleMeasure)/vapply(bootEffectSizesList, 
      FUN.VALUE = double(1), function(x) mean(x$A))
  sdA <- mean(A)
  studentizedCI <- singleMeasure + sdA*quantile(bootRstand, 
      c((1-cutoff)/2, (1+cutoff)/2), na.rm = TRUE)
  names(studentizedCI) <- c("lower", "upper")
  overallCall <- if (eq || any(is.na(studentizedCI))) {
        "Undefined"
      } else {
        if (studentizedCI["lower"] > 0) {
          if (inc) {
            "Syn"
          } else if (dec) {
            "Ant"
          } else
            "Undefined"
        } else if (studentizedCI["upper"] < 0) {
          if (inc) {
            "Ant"
          } else if (dec) {
            "Syn"
          } else 
            "Undefined"
        } else {
          "None"
        }
      }
  
  ans <- list("offAxis" = confInt,
      "single" = list("meanEffect" = singleMeasure,
          "confIntMeanEffect" = studentizedCI,
          "Call" = overallCall),
      "cutoff" = cutoff)
  class(ans) <- append("BIGLconfInt", class(ans))
  return(ans)
}