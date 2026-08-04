### This function is for tuning the number of variables and/or components for a mint.block.spls model

############
# function for finding optimal number of variables per block per component
############

# The user has to enter all blocks in the list called X and use indY to indicate which is treated as Y
# This keeps things simpler
# will need to spell this out as the mint.block.(s)pls functions allow Y to be entered separately

tune.mint.block.spls <- function(
    X,
    indY,
    study,
    ncomp = 2,
    mode = c('regression')[1],
    test.keepX = NULL,
    already.tested.X = NULL,
    test.measure = c("covRatio", "covPred")[1],
    use.ttest = FALSE,
    threshold.ttest = 0.05,
    design = "full",
    scale = TRUE,
    tol = 1e-06,
    max.iter = 100,
    BPPARAM = SerialParam(),
    seed = NULL
){
  
  BPPARAM$RNGseed <- seed
  

  #########
  # If there are no variable selection options, run the function for assessing the number of components
  #########

  
  if (is.null(test.keepX)) {
    print("test.keepX is set to NULL, tuning only for number of components...")
    loocv.ncomp <- .tune.mint.block.spls.ncomp(
      X = X, indY = indY, study = study, ncomp = ncomp, mode = mode,
      design = design, scale = scale, tol = tol, max.iter = max.iter, 
      BPPARAM = BPPARAM, seed = seed
    )
    loocv.result <- list(
      loocv.ncomp = loocv.ncomp
    )
    return(loocv.result)
  }
  

  #########  
  # If there are variable selection options then try them
  #########
  
  # fit initial model to use as source for tuning functions
  fit_initial <- mint.block.spls(
    X = X, indY = indY, study = study, ncomp = ncomp, mode = mode,
    design = design, scale = scale, tol = tol, max.iter = max.iter
  )
  
  # extract pre-processed data from the initial model
  X <- fit_initial$X
  indY <- fit_initial$indY
  mode <- fit_initial$mode
  
  study <- fit_initial$study
  study.names <- levels(study)
  n <- length(study)
  design.mat <- fit_initial$design
  
  # check mode
  if (!is.element(mode,c('regression'))) {
    stop("tuning is only available for modes: 'regression'", call. = FALSE)
  }
  
  # check the specified measure
  test.measure <- test.measure[1]
  if(!is.element(test.measure,c("covRatio", "covPred"))){
    stop("Select either 'covRatio' or 'covPred' for the test.measure")
  }
  
  # checks on variable selection options list
  if (paste(names(test.keepX),collapse=",") != paste(names(X),collapse=",")){
    stop(
      paste0("test.keepX should be a list with element names: ",paste(names(X), collapse = ", "))
    )
  }

  if (any(sapply(test.keepX, function(x) {!is.vector(x)})) == TRUE){
    stop("Each entry of 'test.keepX' must be a vector of keepX values")
  }
  
  if (any(sapply(test.keepX, function(x) {is.list(x)})) == TRUE){
    stop("Each entry of 'test.keepX' must be a vector of keepX values")
  }
  
  if (any(sapply(test.keepX, function(x) {length(x) == 0})) == TRUE){
    stop("Each entry of 'test.keepX' must be a vector of keepX values")
  }
  
  # check for elements of length = 0?
  
  # create a list for storing the optimal number of variables
  if (is.null(already.tested.X)){
    choice.keepX <- sapply(names(test.keepX),function(x) {NULL}, simplify=FALSE, USE.NAMES=TRUE)
  } else {
    if (paste(names(already.tested.X),collapse=",") != paste(names(X),collapse=",")){
      stop(
        paste0("already.tested.X should be a list with element names: ",paste(names(X), collapse = ", "))
      )
    }
    # we require the same number of already tuned components on each block
    if (length(unique(sapply(already.tested.X, length))) > 1){
      stop("Each entry of 'already.tested.X' must have the same number of keepX values")
    }
    # check each element is a non-empty vector and not a list
    if (any(sapply(already.tested.X, function(x) {!is.vector(x)})) == TRUE){
      stop("Each entry of 'already.tested.X' must be a vector of keepX values")
    }
    if (any(sapply(already.tested.X, function(x) {is.list(x)})) == TRUE){
      stop("Each entry of 'already.tested.X' must be a vector of keepX values")
    }
    if (any(sapply(already.tested.X, function(x) {length(x) == 0})) == TRUE){
      stop("Each entry of 'already.tested.X' must be a vector of keepX values")
    }
    # there needs to be at least 1 component to tune variable selection counts for
    if (any(sapply(already.tested.X, function(x) {length(x) >= ncomp})) == TRUE){
      stop("Each entry of 'already.tested.X' must have less keepX values than 'ncomp'")
    }
    choice.keepX <- already.tested.X
  }
  
  # create list of variable selection combinations to try
  keep.vals <- expand.grid(test.keepX)
  keep.vals <- keep.vals[order(rowSums(keep.vals)),, drop = FALSE] # ordered by total number of variables
  keep.vals <- lapply(seq(nrow(keep.vals)),function(x){
    sapply(colnames(keep.vals),function(y){
      keep.vals[x,y]
    },simplify=FALSE, USE.NAMES=TRUE)
  })
  
  
  ############### Loop over components to perform variable selection for each component
  
  # empty list for storage
  check.comp <- vector("list")
  compareVals.comp <- vector("list")
  compareVals.sum.comp <- vector("list")
  
  # only components for which we are tuning the number of variables
  comp.real <- (length(choice.keepX[[1]])+1):ncomp
  for (comp in comp.real){
    
    # get fitted variates for each variable selection option
    fit.values <- bplapply(
      keep.vals, FUN = function(keep.val){ 
        # combine chosen and candidate variable counts
        keepX <- sapply(names(choice.keepX), function(block){
          c(choice.keepX[[block]], keep.val[[block]])
        },simplify=FALSE, USE.NAMES=TRUE)
        # perform logocv
        bplapply(
          study.names, FUN = function(study.name){
            # test set data
            InTest <- study==study.name
            X.test <- lapply(X, function(block){block[InTest,, drop = FALSE]})
            study.test <- factor(as.character(study[InTest]))
            # get the values
            modelVals <- .tune.mint.block.spls.model.values(
              X.train = X, study.train = study, 
              X.test = X.test, study.test = study.test,
              indY = indY, ncomp = comp, mode = mode, keepX = keepX,
              design = design, scale = scale, tol = tol, max.iter = max.iter)$variates
            return(modelVals)
          }, BPPARAM = BPPARAM )
      }, BPPARAM = BPPARAM)
    fit.values <- lapply(fit.values, function(x){
      names(x) <- study.names
      return(x)
    })
    
    # get predicted variates for each variable selection option
    pred.values <- bplapply(
      keep.vals, FUN = function(keep.val){ 
        # combine chosen and candidate variable counts
        keepX <- sapply(names(choice.keepX), function(block){
          c(choice.keepX[[block]], keep.val[[block]])
        },simplify=FALSE, USE.NAMES=TRUE)
        # perform logocv
        bplapply(
          study.names, FUN = function(study.name){
            # test set data
            InTest <- study==study.name
            X.test <- lapply(X, function(block){block[InTest,, drop = FALSE]})
            study.test <- factor(as.character(study[InTest]))
            # training dataset
            InTrain <- !InTest
            X.train <- lapply(X, function(block){block[InTrain,, drop = FALSE]})
            study.train <- factor(as.character(study[InTrain]))
            # get the values
            modelVals <- .tune.mint.block.spls.model.values(
              X.train = X.train, study.train = study.train, 
              X.test = X.test, study.test = study.test,
              indY = indY, ncomp = comp, mode = mode, keepX = keepX,
              design = design, scale = scale, tol = tol, max.iter = max.iter)$variates
            return(modelVals)
          }, BPPARAM = BPPARAM )
      }, BPPARAM = BPPARAM)
    pred.values <- lapply(pred.values, function(x){
      names(x) <- study.names
      return(x)
    })
    
    # calculate covariance by study and variable selection option
    compareVals.study <- sapply(1:length(keep.vals), function(keep.val){
      sapply(study.names, function(study.name){
        # covariance of fitted variates values
        fit.variates.comp <- lapply(fit.values[[keep.val]][[study.name]], function(block){
          block[,comp]
        })
        fit.variates.comp <- do.call(cbind, fit.variates.comp)
        covFit <- cov(fit.variates.comp)
        covFit <- sum(upper.tri(covFit, diag = FALSE) * design.mat * covFit)
        # covariance of predicted variates vlaues
        pred.variates.comp <- lapply(pred.values[[keep.val]][[study.name]], function(block){
          block[,comp]
        })
        pred.variates.comp <- do.call(cbind, pred.variates.comp)
        covPred <- cov(pred.variates.comp)
        covPred <- sum(upper.tri(covPred, diag = FALSE) * design.mat * covPred)
        # values for evaluation
        Vals <- data.frame(
          covPred = covPred,
          covRatio = covPred/covFit
        )
        return(Vals)
      }, simplify = FALSE, USE.NAMES=TRUE)
    }, simplify = FALSE)
    
    # concatenate values for each study to produce a data frame for each variable selection option
    compareVals <- lapply(1:length(compareVals.study), function(x){
      do.call(rbind, compareVals.study[[x]])
    })
    
    # try variable selection options to find optimal values
    # first loop to find variable selection that gives the best mean value overall
    opt.index <- 1
    opt.compareVal <- compareVals[[opt.index]]
    for (i in (2:length(compareVals))){
      compareVal <- compareVals[[i]]
      if (mean(compareVal[,test.measure]) > mean(opt.compareVal[,test.measure])){
        opt.compareVal <- compareVal
        opt.index <- i
      } 
    }
    
    # if t-tests option chosen then second loop over variable selection options, which are in ascending order
    # select lowest number of variables whose study values are not statistically different from those of the global optimal
    if (use.ttest & (length(opt.compareVal[,test.measure])>=3) & (opt.index>1)){
      for (i in 1:(opt.index-1)){
        compareVal <- compareVals[[i]]
        if (mean(unlist(keep.vals[[i]]) <= unlist(keep.vals[[opt.index]]))==1){
          test.result <- t.test(x = opt.compareVal[,test.measure], y = compareVal[,test.measure], alternative = "greater", paired = TRUE)
          if (test.result$p.value > threshold.ttest){ # non significant result so can go with the lower number of variables
            opt.compareVal <- compareVal
            opt.index <- i
            break
          }
        }
      }
    }
    
    # Update choice.keepX with the the optimal values
    choice.keepX.comp <- keep.vals[[opt.index]]
    choice.keepX <- sapply(names(choice.keepX), function(block){
      c(choice.keepX[[block]], choice.keepX.comp[[block]])
    },simplify=FALSE, USE.NAMES=TRUE)
    
    # aggregate the values for visualisation and add in associated keep.vals options
    compareVals.sum <- lapply(compareVals.study, function(compareVal.study){
      sapply(study.names, function(study.name){
        return(compareVal.study[[study.name]] * (1/length(study.names)))
      }, simplify = FALSE)
    })
    
    compareVals.sum <- sapply(1:length(compareVals.sum), function(x) {
      compareVals.df <- Reduce("+",compareVals.sum[[x]])
      keep.vals.df <- as.data.frame(keep.vals[[x]])
      colnames(keep.vals.df) <- paste0("keep_",colnames(keep.vals.df))
      compareVals.sum.df <- merge(compareVals.df,keep.vals.df)
      return(compareVals.sum.df)
    }, simplify = FALSE)
    compareVals.sum = do.call(rbind, compareVals.sum)
    
    
    # save results for the component
    
    compareVals.comp[[paste0("comp",comp)]] <- compareVals
    compareVals.sum.comp[[paste0("comp",comp)]] <- compareVals.sum


  }
  
  ############### End of loop over components
  
  # list of variable selection tuning results
  loocv.variable.selection <- list(
    keep.vals = keep.vals, 
    study.names = study.names,
    study = study,
    mode = mode,
    design.mat = design.mat,
    compareVals.comp = compareVals.comp,
    compareVals.sum.comp = compareVals.sum.comp,
    test.measure = test.measure,
    choice.keepX = choice.keepX
  )
  
  # get Q2 for the selected number of components
  loocv.ncomp <- .tune.mint.block.spls.ncomp(
    X = X, indY = indY, study = study, ncomp = ncomp, 
    # keepX = choice.keepX,
    mode = mode, design = design, scale = scale, tol = tol, max.iter = max.iter, 
    BPPARAM = BPPARAM, seed = seed
  )
  
  # return combined results
  loocv.result <- list(
    loocv.variable.selection = loocv.variable.selection,
    loocv.ncomp = loocv.ncomp
  )
  return(loocv.result)
  
}

################
# function for assessing the number of components for given number of variables per component
# can be used for spls as well as pls mint.block models
################

.tune.mint.block.spls.ncomp <- function(
    X,
    indY,
    study,
    ncomp = 2,
    mode = c('regression')[1],
    keepX,
    design = "full",
    scale = TRUE,
    tol = 1e-06,
    max.iter = 100,
    BPPARAM = SerialParam(),
    seed = NULL
){
  
  BPPARAM$RNGseed <- seed
  
  # hardcode Q2 cutoff
  limQ2 <- 0.0975
  
  # fit initial model to use as source for tuning functions
  fit_initial <- mint.block.spls(
    X = X, indY = indY, study = study, ncomp = ncomp, mode = mode, keepX = keepX,
    design = design, scale = scale, tol = tol, max.iter = max.iter
  )
  
  # extract pre-processed data from the initial model
  X <- fit_initial$X
  indY <- fit_initial$indY
  mode <- fit_initial$mode
  if (!is.element(mode,c('regression'))) {
    stop("tuning is only available for modes: 'regression'", call. = FALSE)
  }
  study <- fit_initial$study
  study.names <- levels(study)
  n <- length(study)
  
  # create list of keepX values - needed for the PRESS function to work
  keepX <- do.call(rbind,fit_initial$keepX)
  keepX[,names(X)[indY]] <- as.vector(unlist(fit_initial$keepY))
  keepX <- keepX[,names(X)]
  keepX <- as.list(keepX)
  
  # fit model with all data to calculate total and study-wise RSS
  fit.values <- .tune.mint.block.spls.model.values(
    X.train = X, study.train = study, 
    X.test = X, study.test = study,
    indY = indY, ncomp = ncomp, mode = mode, keepX = keepX,
    design = design, scale = scale, tol = tol, max.iter = max.iter)
  # total RSS for Q2 calculation
  RSS <- fit.values$SS
  block.names <- names(RSS)
  RSS.Q2 <- lapply(RSS, function(x){rbind(rep(n - 1, ncol(x)), x[-nrow(x),])})
  # RSS and RSS for Q2 calculation by study
  fit.squares <- fit.values$squares
  RSS.study <- sapply(study.names, function(study.name){
    sapply(block.names, function(block.name){
      do.call(rbind, lapply(1:ncomp, function(comp){
        colSums(fit.squares[[block.name]][[comp]][study==study.name,], na.rm=TRUE)
        }))
    }, simplify=FALSE, USE.NAMES=TRUE)
  }, simplify=FALSE, USE.NAMES=TRUE)
  RSS.Q2.study <- sapply(study.names, function(study.name){
    lapply(RSS.study[[study.name]], function(block){rbind(rep(sum(study==study.name) - 1, ncol(block)), block[-nrow(block),])})
  }, simplify=FALSE, USE.NAMES=TRUE)
  
  
  # Now perform LOGOCV to calculate total and study-wise PRESS
  PRESS.study <- bplapply(study.names, FUN = function(study.name){
    # test set
    InTest <- study==study.name
    X.test <- lapply(X, function(block){block[InTest,, drop = FALSE]})
    study.test <- factor(as.character(study[InTest]))
    # training dataset
    InTrain <- !InTest
    X.train <- lapply(X, function(block){block[InTrain,, drop = FALSE]})
    study.train <- factor(as.character(study[InTrain]))
    # cross validation
    .tune.mint.block.spls.model.values(
      X.train = X.train, study.train = study.train,
      X.test = X.test, study.test = study.test,
      indY = indY, ncomp = ncomp, mode = mode, keepX = keepX,
      design = design, scale = scale, tol = tol, max.iter = max.iter)$SS
  }, BPPARAM = BPPARAM)
  names(PRESS.study) <- study.names
  # sum PRESS over studies by block
  PRESS <- sapply(block.names, function(block.name){
    Reduce("+", lapply(PRESS.study, function(fold)(fold[[block.name]])))
    }, simplify=FALSE, USE.NAMES=TRUE)
  
  
  # compute Q2 criterion using total PRESS and RSS values  
  Q2 <- sapply(block.names, function(block.name){
    1 - rowSums(PRESS[[block.name]])/rowSums(RSS.Q2[[block.name]])
    }, simplify=FALSE, USE.NAMES=TRUE)
  
  # Q2.variable <- sapply(block.names, function(block.name){
  #   1 - PRESS[[block.name]]/RSS.Q2[[block.name]]
  #   }, simplify=FALSE, USE.NAMES=TRUE)
  
  # select the optimal number of components based on Q2
  choice.ncomp <- sapply(block.names, function(block.name){
    excess.comp <- which(Q2[[block.name]] < limQ2)
    if(length(excess.comp)==0){
      return(ncomp)
    } else {
      return(max(min(excess.comp)-1,1))
    }
  }, USE.NAMES=TRUE)
  
  # calculate Q2 by study
  Q2.study <- sapply(block.names, function(block.name){
    sapply(study.names, function(study.name){
      1 - rowSums(PRESS.study[[study.name]][[block.name]])/rowSums(RSS.Q2.study[[study.name]][[block.name]])
    }, simplify=FALSE, USE.NAMES=TRUE)
  }, simplify=FALSE, USE.NAMES=TRUE)
  
  # choose optimal components by study
  choice.ncomp.study <- lapply(Q2.study,function(block){
    lapply(block, function(fold){
      excess.comp <- which(fold < limQ2)
      if(length(excess.comp)==0){
        return(ncomp)
      } else {
        return(max(min(excess.comp)-1,1))
      }
    })
  })
  
  # put Q2.study into long format for plotting
  Q2.study <- sapply(block.names, function(block.name){
    Q2.block <- sapply(study.names, function(study.name){
      data.frame(
        study = study.name,
        comp = 1:length(Q2.study[[block.name]][[study.name]]),
        Q2 = Q2.study[[block.name]][[study.name]]
      )
    }, simplify=FALSE, USE.NAMES=TRUE)
    do.call(rbind, Q2.block)
  }, simplify=FALSE, USE.NAMES=TRUE)
  
  # put PRESS.study into long format for plotting and calculate means
  MPSE.study <- sapply(block.names, function(block.name){
    MPSE.block <- sapply(study.names, function(study.name){
      data.frame(
        study = study.name,
        comp = 1:nrow(PRESS.study[[study.name]][[block.name]]),
        MPSE = rowSums(PRESS.study[[study.name]][[block.name]]) / sum(study==study.name)
      )
    }, simplify=FALSE, USE.NAMES=TRUE)
    do.call(rbind, MPSE.block)
  }, simplify=FALSE, USE.NAMES=TRUE)
  
  # return the result
  result <- list(
    limQ2 = limQ2,
    mode = mode,
    RSS.study = RSS.study,
    RSS.Q2.study = RSS.Q2.study,
    RSS = RSS,
    RSS.Q2 = RSS.Q2,
    PRESS.study = PRESS.study,
    MPSE.study = MPSE.study,
    PRESS = PRESS,
    Q2.study = Q2.study,
    choice.ncomp.study = choice.ncomp.study,
    Q2 = Q2,
    choice.ncomp = choice.ncomp,
    indY = indY
  )
  
}

####################
## single function that retrieves either fitted or predicted variates and squared errors
####################

.tune.mint.block.spls.model.values <- function(
    X.train, 
    X.test,
    study.train, 
    study.test,
    indY,
    ncomp,
    mode,
    keepX,
    design,
    scale,
    tol,
    max.iter
){
  
  # fit and predict
  fit <- mint.block.spls(
    X = X.train, indY = indY, study = study.train, ncomp = ncomp, mode = mode, keepX = keepX, 
    design = design, scale = scale, tol = tol, max.iter = max.iter
  )
  X.test.X <- X.test[-indY]
  X.test.Y <- X.test[[indY]]
  pred <- predict(object = fit, newdata = X.test.X, study.test = study.test)
  
  # calculate variates and sums of squares contingent on the mode
  if (mode == "regression"){
    # predicted Y variates: calculated using loadings from model fit and Y test actual and predicted
    Ymat <- X.test.Y
    Rmat <- replace(Ymat, is.na(Ymat), 0)
    Umat <- matrix(0, nrow=nrow(Ymat), ncol = ncomp)
    for (comp in 1:ncomp){
      Umat[,comp] <- Rmat %*% fit$loadings[[indY]][,comp]
      Rmat <- Ymat - pred$WeightedPredict[,,comp] 
      Rmat <- replace(Rmat, is.na(Rmat), 0)
    }
    # squared residuals or errors
    squares <- vector("list")
    squares[[names(X.train)[indY]]] <- lapply(1:ncomp, function(x){
      (X.test.Y - pred$WeightedPredict[,,x])^2
    })
    # sum of squares
    SS <- vector("list")
    SS[[names(X.train)[indY]]] <- do.call(rbind, lapply(1:ncomp, function(x){
      colSums((X.test.Y - pred$WeightedPredict[,,x])^2, na.rm=TRUE)
      }))
  } else {
    stop("tuning is only available for modes: 'regression'", call. = FALSE)
  }
  
  # make sure all variates have column and row names
  rownames(Umat) <- rownames(Ymat)
  colnames(Umat) <- colnames(pred$variates[[1]])
  # list of variates with names and order consistent with input data
  variates <- pred$variates
  variates[[names(X.train)[indY]]] <- Umat
  variates <- variates[names(X.train)]
  
  # return values
  modelVals <- list(
    squares = squares,
    SS = SS,
    variates = variates
  )
  return(modelVals)
  
}
