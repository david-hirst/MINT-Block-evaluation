### This function is for tuning the number of variables and/or components for a mint.block.splsda model

############
# function for finding optimal number of variables per block per component
############

# The user can currently either enter Y as a factor vector or as an indicator matrix in the list X, with indY used to specify the position

tune.mint.block.splsda <- function(
    X,
    Y,
    indY,
    study,
    ncomp = 2,
    test.keepX = NULL,
    already.tested.X = NULL,
    prediction.distance = c("max.dist", "centroids.dist", "mahalanobis.dist"),
    prediction.type = c("AveragedPredict.class", "WeightedPredict.class", "MajorityVote", "WeightedVote"),
    error.type = c("ER", "BER"),
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
  
  ##########
  # If there are no variable selection options, run the function for assessing the number of components
  #########
  
  if (is.null(test.keepX)) {
    print("test.keepX is set to NULL, tuning only for number of components...")
    loocv.ncomp <- .tune.mint.block.splsda.ncomp(
      X = X, Y = Y, indY = indY, study = study, ncomp = ncomp, 
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
  
  # ensure only single prediction type and distance, and that they are compatible
  prediction.distance <- prediction.distance[1]
  prediction.type <- prediction.type[1]
  if (is.element(prediction.distance,c("centroids.dist", "mahalanobis.dist")) &
      is.element(prediction.type,c("AveragedPredict.class", "WeightedPredict.class"))){
    stop(paste0(prediction.distance,
                " can only be used with either prediction type MajorityVote or WeightedVote. ", 
                "Only max.dist can be used with prediction type ",
                prediction.type
                )
         )
  }
  
  # check the error types
  error.type <- error.type[1]
  if(!is.element(error.type,c("ER","BER"))){
    stop("Select either 'ER' or 'BER' for the error.type")
  }
  
  # fit initial model to use as source for tuning functions
  fit_initial <- mint.block.splsda(
    X = X, Y = Y, indY = indY, study = study, ncomp = ncomp,  
    design = design, scale = scale, tol = tol, max.iter = max.iter
  )
  
  # extract pre-processed data from the initial model
  X <- fit_initial$X
  Y <- fit_initial$Y
  names(Y) <- rownames(X[[1]])
  study <- fit_initial$study
  study.names <- levels(study)
  
  # checks on variable selection options list
  # check alignment with blocks in X
  if (paste(names(test.keepX),collapse=",") != paste(names(X),collapse=",")){
    stop(
      paste0("test.keepX should be a list with element names: ",paste(names(X), collapse = ", "))
    )
  }
  # check each element is a vector and not a list
  if (any(sapply(test.keepX, function(x) {!is.vector(x)})) == TRUE){
    stop("Each entry of 'test.keepX' must be a vector of keepX values")
  }
  if (any(sapply(test.keepX, function(x) {is.list(x)})) == TRUE){
    stop("Each entry of 'test.keepX' must be a vector of keepX values")
  }
  if (any(sapply(test.keepX, function(x) {length(x) == 0})) == TRUE){
    stop("Each entry of 'test.keepX' must be a vector of keepX values")
  }
  
  # create a list, called choice.keepX, for storing optimal number of variables
  # if already.tested.X is NULL then create an empty list 
  # otherwise perform some checks and initialise choice.keepX with already.tested.X
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
    # check each element is a vector and not a list
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
    if (length(already.tested.X[[1]]) >= ncomp){
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
  
  # empty list to store aggregate results for each component
  loocv.global.comp <- vector("list")
  loocv.study.comp <- vector("list")
  
  ############### Loop over components to perform variable selection for each component
  
  # only components for which we are tuning the number of variables
  comp.real <- (length(choice.keepX[[1]])+1):ncomp
  for (comp in comp.real){
    
    # perform loocv for each combination of variable selection counts
    loocv <- bplapply(keep.vals, FUN = function(keep.val){ 
        # combine choosen and candidate variable counts
        keepX <- sapply(names(choice.keepX), function(block){
          c(choice.keepX[[block]], keep.val[[block]])
        },simplify=FALSE, USE.NAMES=TRUE)
        # Now perform LOOCV
        loocv <- .tune.mint.block.splsda.loocv(
          study.names = study.names,
          X = X, Y = Y, study = study, ncomp = comp,  keepX = keepX,
          prediction.distance = prediction.distance, prediction.type = prediction.type,
          design = design, scale = scale, tol = tol, max.iter = max.iter,
          BPPARAM = BPPARAM, seed = seed
        )
        # extract the error rates for each study
        loocv.study <- lapply(loocv$loocv.study, function(study){
          unname(study$pred.error[[prediction.type]][[prediction.distance]][[error.type]][comp])
        })
        loocv.study <- unlist(loocv.study)
        # extract the global error rates
        loocv.global <- loocv$loocv.global[[prediction.type]][[prediction.distance]][[error.type]][["mean"]][comp]
        # return the result
        loocv <- list(
          keepX = keepX,
          loocv.study = loocv.study,
          loocv.global = loocv.global
        )
        return(loocv)
      }, BPPARAM = BPPARAM)
    
    # find variable selection that gives the lowest global error rate
    loocv.global <- unlist(lapply(loocv,function(x){x$loocv.global}))
    optimal.index <- min(which(loocv.global==min(loocv.global)))
    optimal.loocv <- loocv[[optimal.index]]
    
    # if required then run t-tests to find option with lowest number of variables
    # that has study-wise error rates not statistically significantly higher 
    # than the option with the lowest global error rate
    if (use.ttest & (length(optimal.loocv$loocv.study)>=3) & (optimal.index>1)){
      for (i in 1:(optimal.index-1)){
        candidate.loocv <- loocv[[i]]
        if (mean(unlist(candidate.loocv$keepX) <= unlist(optimal.loocv$keepX))==1){
          test.result <- t.test(x = candidate.loocv$loocv.study, y = optimal.loocv$loocv.study, 
                               alternative = "greater", paired = TRUE)
          if (test.result$p.value > threshold.ttest){ 
            optimal.loocv <- candidate.loocv # non significant result so can go with the lower number of variables
            optimal.index <- i
            break
          }
        }
      }
    }
    
    # update the cumulative optimal variable selection counts
    choice.keepX <- optimal.loocv$keepX
    
    # save all global errors for the component
    loocv.keepX <- sapply(names(choice.keepX),function(block){
      keepX.block <- lapply(loocv, function(candidate){
        candidate$keepX[[block]][comp]
      })
      keepX.block <- unlist(keepX.block)
    }, simplify = FALSE, USE.NAMES = TRUE)
    loocv.keepX <- as.data.frame(do.call(cbind,loocv.keepX))
    loocv.global.comp[[paste0("comp",comp)]] <- data.frame(
      loocv.keepX, 
      loocv.error = loocv.global
      )
    
    # save all study-wise errors for the component
    loocv.study <- lapply(loocv, function(candidate){
      loocv.errors <- data.frame(
        study = names(candidate$loocv.study),
        loocv.error = unname(candidate$loocv.study)
      )
      loocv.keepX <- lapply(candidate$keepX, function(block){
        block[comp]
      })
      loocv.keepX <- as.data.frame(do.call(cbind,loocv.keepX))
      loocv.study <- data.frame(loocv.keepX, loocv.errors)
    })
    loocv.study <- do.call(rbind,loocv.study)
    loocv.study.comp[[paste0("comp",comp)]] <- loocv.study
    
    # end of loop
    
  }
  
  # list of variable selection tuning results
  loocv.variable.selection <- list(
    prediction.type = prediction.type,
    prediction.distance = prediction.distance,
    error.type = error.type,
    choice.keepX = choice.keepX,
    loocv.global.comp = loocv.global.comp,
    loocv.study.comp = loocv.study.comp
  )
  
  # check the number of components for the sparse model
  loocv.ncomp <- .tune.mint.block.splsda.ncomp(
    X = X, Y = Y, study = study, ncomp = ncomp, keepX = choice.keepX,
    design = design, scale = scale, tol = tol, max.iter = max.iter, 
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
# can be used for splsda as well as plsda models
################

.tune.mint.block.splsda.ncomp <- function(
    X,
    Y,
    indY,
    study,
    ncomp = 2,
    keepX = NULL,
    design = "full",
    scale = TRUE,
    tol = 1e-06,
    max.iter = 100,
    BPPARAM = SerialParam(),
    seed = NULL
) {
  
  
  BPPARAM$RNGseed <- seed
  
  # fit initial model to use as source for tuning functions
  fit_initial <- mint.block.splsda(
    X = X, Y = Y, indY = indY, study = study, ncomp = ncomp,  keepX = keepX,  
    design = design, scale = scale, tol = tol, max.iter = max.iter
    )
  
  # extract pre-processed data from the initial model
  X <- fit_initial$X
  Y <- fit_initial$Y
  names(Y) <- rownames(X[[1]])
  study <- fit_initial$study
  study.names <- levels(study)
  
  # Now perform LOOCV
  loocv <- .tune.mint.block.splsda.loocv(
    study.names = study.names,
    X = X, Y = Y, study = study, ncomp = ncomp,  keepX = keepX,  
    design = design, scale = scale, tol = tol, max.iter = max.iter,
    BPPARAM = BPPARAM, seed = seed
  )
    
  # Select optimal number of components based on global average
  n.comp.global <- lapply(loocv$loocv.global, function(type){
    n.comp.dist <- lapply(type, function(dist){
      n.comp.rate <- lapply(dist, function(rate){
        min(which(rate[["mean"]]==min(rate[["mean"]])))
      })
      unlist(n.comp.rate)
    })
    do.call(rbind, n.comp.dist)
  })
  
  # return global and study-wise error rates in long format for plotting
  error.types <- names(loocv$loocv.global[[1]][[1]])
  summary.types <- names(loocv$loocv.global[[1]][[1]][[1]])
  prediction.types <- names(loocv$loocv.global)
  # obtain long format data frame of errors per error type
  errors.global <- sapply(error.types, function(error.type){
    global.error.summary <- sapply(summary.types, function(summary.type){
      global.error.type <- sapply(prediction.types, function(prediction.type){
        global.error.dist <- lapply(loocv$loocv.global[[prediction.type]], function(dist){
          dist[[error.type]][[summary.type]]
        })
        global.error.type <- do.call(rbind, global.error.dist)
        rownames(global.error.type) <- paste0(rownames(global.error.type),"_",prediction.type)
        colnames(global.error.type) <- paste0("comp",1:ncol(global.error.type))
        return(global.error.type)
      }, simplify = FALSE, USE.NAMES = TRUE)
      # single matrix of errors
      global.error <- do.call(rbind, global.error.type)
      # put into long format for plotting
      global.error.comps <- colnames(global.error)
      global.error <- as.data.frame(global.error)
      global.error$prediction.method <- rownames(global.error)
      global.error <- reshape(global.error,
                              varying = global.error.comps,
                              v.names = summary.type,
                              timevar = "ncomp",
                              idvar =  "prediction.method",
                              direction = "long")
      
    }, simplify = FALSE, USE.NAMES = TRUE)
    # merge into a single dataframe for the given error type
    errors.global <- Reduce(
      f = function(x, y){
        merge(x, y, by = c("prediction.method","ncomp"), all = TRUE)
      }, 
      x = global.error.summary
    )
  }, simplify = FALSE, USE.NAMES = TRUE)
  
  # obtain long format data frame of study-wise errors per error type
  errors.study <- sapply(study.names, function(study.name){
    errors <- loocv$loocv.study[[study.name]]$pred.error
    errors <- sapply(prediction.types, function(prediction.type){
      errors <- errors[[prediction.type]]
      distances <- names(errors)
      errors <- sapply(distances, function(dist){
        errors <- errors[[dist]]
        errors <- sapply(error.types, function(error.type){
          errors <- data.frame(
            study.name = study.name,
            prediction.method = paste0(dist,"_",prediction.type),
            error.type = error.type,
            error = errors[[error.type]],
            ncomp = 1:length(errors[[error.type]])
          )
        }, simplify = FALSE)
        errors <- do.call(rbind,errors)
      }, simplify = FALSE)
      errors <- do.call(rbind, errors)
    }, simplify = FALSE)
    errors <- do.call(rbind, errors)
  }, simplify = FALSE)
  errors.study <- do.call(rbind, errors.study)
  # split by error type
  errors.study <- sapply(error.types, function(error.type){
    errors.study[errors.study$error.type==error.type,]
    }, simplify = FALSE, USE.NAMES = TRUE)
  
  # final list to return
  result <- list(
    # loocv.error.rates = loocv,
    n.comp.global = n.comp.global,
    errors.global = errors.global,
    errors.study = errors.study
  )
  return(result)

}

#############
#### Function to perform LOOCV to assess predictions
#############

.tune.mint.block.splsda.loocv <- function(
    study.names,
    X,
    Y,
    study,
    ncomp = 2,
    keepX = NULL,
    prediction.distance = c("max.dist", "centroids.dist", "mahalanobis.dist"),
    prediction.type = c("AveragedPredict.class", "WeightedPredict.class", "MajorityVote", "WeightedVote"),
    design = "full",
    scale = TRUE,
    tol = 1e-06,
    max.iter = 100,
    BPPARAM = SerialParam(),
    seed = NULL
){
  # This function performs loocv with each study taking a turn at being the test set
  
  BPPARAM$RNGseed <- seed
  
  # first get errors and error rates by study
  loocv.study <- bplapply(study.names, function(study.name){
    
    # split data into test and train based on study name
    InTest <- study==study.name
    X.test <- lapply(X, function(x){x[InTest,, drop = FALSE]})
    Y.Test <- Y[InTest]
    study.test <- factor(as.character(study[InTest]))
    InTrain <- !InTest
    X.train <- lapply(X, function(x){x[InTrain,, drop = FALSE]})
    Y.train <- Y[InTrain]
    study.train <- factor(as.character(study[InTrain]))
    
    # fit a model with the training set
    fit.train <- mint.block.splsda(
      X = X.train, Y = Y.train, study = study.train, ncomp = ncomp,  keepX = keepX,
      design = design, scale = scale, tol = tol, max.iter = max.iter
    )
    
    # make predictions for the test set
    pred.test <- predict(
      object = fit.train, newdata = X.test, study.test = study.test
    )
    # subset the prediction list object to only include the class predictions
    # each element is now a set of class predictions for a given prediction type
    # each sub element is a matrix of class predictions for a given distance
    # the rows of the matrix are samples and the columns ar the dimensionality
    pred.test <- pred.test[prediction.type]
    
    # calculate prediction errors
    pred.error <- lapply(pred.test, function(type){
      distances <- names(type)[is.element(names(type), prediction.distance)]
      sapply(distances, function(dist){
        # errors matrix: TRUE if class prediction is wrong, otherwise FALSE 
        # rows are samples and columns are components (cumulative impact)
        errors <- type[[dist]] != Y.Test
        # If the prediction is NA then this is treated as an error
        errors[is.na(errors)] <- TRUE
        # get overall error rate
        ER <- apply(errors, 2, mean, na.rm = TRUE)
        # get errors by response factor level 
        LER <- sapply(levels(Y.Test), function(level){
          level.errors <- errors[Y.Test==level,, drop = FALSE]
          LER <- apply(level.errors, 2, mean, na.rm = TRUE)
        }, simplify = FALSE, USE.NAMES = TRUE)
        LER <- do.call(rbind, LER)
        # get balanced error rate
        BER <- apply(LER, 2, mean, na.rm = TRUE)
        pred.error <- list(
          errors = errors,
          ER = ER,
          BER = BER
        )
      }, simplify = FALSE, USE.NAMES = TRUE)
    })
    
    # return as a list
    result <- list(
      Y.Test = Y.Test,
      pred.error = pred.error
    )
    
  }, BPPARAM = BPPARAM)
  names(loocv.study) = study.names
  
  # calculate averages and standard deviations of study-wise error rates
  loocv.global <- sapply(prediction.type, function(type){
    distances <- names(loocv.study[[1]]$pred.error[[type]])
    sapply(distances, function(dist){
      # get the overall Error Rate averages and standard deviations
      ER <- lapply(loocv.study, function(study){
        study$pred.error[[type]][[dist]]$ER
      })
      ER <- do.call(rbind, ER)
      # ER <- apply(ER, 2, mean, na.rm = TRUE)
      ER <- list(
        mean = apply(ER, 2, mean, na.rm = TRUE),
        sd = apply(ER, 2, sd, na.rm = TRUE),
        se = apply(ER, 2, FUN = function(x){sd(x,na.rm=TRUE)/sqrt(length(x[!is.na(x)]))})
      )
      # get the Balanced Error Rate averages and standard deviations
      BER <- lapply(loocv.study, function(study){
        study$pred.error[[type]][[dist]]$BER
      })
      BER <- do.call(rbind, BER)
      # BER <- apply(BER, 2, mean, na.rm = TRUE)
      BER <- list(
        mean = apply(BER, 2, mean, na.rm = TRUE),
        sd = apply(BER, 2, sd, na.rm = TRUE),
        se = apply(BER, 2, FUN = function(x){sd(x,na.rm=TRUE)/sqrt(length(x[!is.na(x)]))})
      )
      # return results
      pred.error <- list(
        ER = ER,
        BER = BER
      )
    }, simplify = FALSE, USE.NAMES = TRUE)
  }, simplify = FALSE, USE.NAMES = TRUE) 
  
  # return global and study results
  loocv <- list(
    loocv.study = loocv.study,
    loocv.global = loocv.global
  )
  return(loocv)
  
}


####################
#### archive
####################

# alternatively calculate global errors from raw studywise errors?
# loocv.global = sapply(prediction.type, function(type){
#   distances = names(loocv.study[[1]]$pred.error[[type]])
#   sapply(distances, function(dist){
#     # get the errors
#     errors = lapply(loocv.study, function(study){
#       study$pred.error[[type]][[dist]]$errors
#     })
#     errors = do.call(rbind, errors)
#     # get the true values
#     Y.Test = lapply(loocv.study, function(study){
#       study$Y.Test
#     })
#     Y.Test = unlist(Y.Test)
#     # calculate error rate
#     ER = apply(errors, 2, mean, na.rm = TRUE)
#     # error rate per factor level
#     LER = sapply(levels(Y.Test), function(level){
#       level.errors = errors[Y.Test==level,, drop = FALSE] # keep as matrix if only 1 component
#       LER = apply(level.errors, 2, mean, na.rm = TRUE)
#     }, simplify = FALSE, USE.NAMES = TRUE)
#     LER = do.call(rbind, LER)
#     # balanced error rate
#     BER = apply(LER, 2, mean, na.rm = TRUE)
#     # return results
#     pred.error = list(
#       ER = ER,
#       BER = BER
#     )
#   }, simplify = FALSE, USE.NAMES = TRUE)
# }, simplify = FALSE, USE.NAMES = TRUE)

