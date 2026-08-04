# Benchmarking analysis functions

###################
# function to calculate r2 values by study
# based on lm fits with a variate as the response and the study as an explanatory factor
# returns data frame to use for plots
###################

extract.fit.r2 <- function(fits, blocks = c("taxa", "GO", "continuous_env"), comps = seq(ncomp), studies = sample.study){
  
  r2.all.fits <- lapply(fits, function(fit){
    sapply(blocks, function(block){
      sapply(comps, function(comp){
        r2 <- summary(lm(fit$variates[[block]][,comp] ~ as.factor(studies)))$r.squared
      }, simplify = FALSE, USE.NAMES = TRUE)
    }, simplify = FALSE, USE.NAMES = TRUE)
  })
  
  r2.df <- do.call(rbind, lapply(names(r2.all.fits), function(fit){
    do.call(rbind, lapply(names(r2.all.fits[[fit]]), function(block){
      data.frame(
        fit = fit,
        block = block,
        comp = paste0("comp ", comps),
        r2 = unlist(r2.all.fits[[fit]][[block]])
      )
    }))
  }))
  
  return(r2.df)
  
}

###################
# function to extract silhouette coefficients
# calculated per block per fit based on variates values
###################

extract.fit.sc <- function(fits, 
                          blocks = c("taxa", "GO", "continuous_env"),
                          studies = sample.study
){
  
  study.fctr <- as.factor(studies)
  study.num <- as.numeric(study.fctr)
  
  sc.all.fits <- lapply(fits, function(fit){
    sapply(blocks, function(block){
      d <- stats::dist(x = fit$variates[[block]], method = "euclidean")
      sc <- silhouette(x = study.num, dist = d)
      sc <- data.frame(
        study = sc[,1],
        width = sc[,3]
      ) %>%
        group_by(study) %>%
        summarise(ave_width = mean(width)) %>%
        as.data.frame()
    }, simplify = FALSE, USE.NAMES = TRUE)
  })
  
  sc.df <- do.call(rbind, lapply(names(sc.all.fits), function(fit){
    do.call(rbind, lapply(names(sc.all.fits[[fit]]), function(block){
      data.frame(
        fit = fit,
        block = block,
        study = sc.all.fits[[fit]][[block]]$study,
        ave_width = sc.all.fits[[fit]][[block]]$ave_width
      )
    }))
  }))
  
  sc.df$study.name <- levels(study.fctr)[sc.df$study]
  
  return(sc.df)
}

###################
# function to extract ARI values based on kmeans clustering
# calculated per block per fit based on variates values
###################

extract.fit.ari <- function(fits, blocks = c("taxa", "GO", "continuous_env"), 
                           studies = sample.study, iter.max = 100, nstart = 20){
  
  ari.all.fits <- lapply(fits, function(fit){
    sapply(blocks, function(block){
      km.clusters <- kmeans(x = fit$variates[[block]], centers = length(unique(studies)), iter.max = iter.max, nstart = nstart)$cluster
      ari <- mclust::adjustedRandIndex(x = km.clusters, y = studies)
    }, simplify = FALSE, USE.NAMES = TRUE)
  })
  
  # convert to a dataframe
  ari.all.fits <- do.call(rbind, lapply(names(ari.all.fits), function(fit.name){
    do.call(rbind, lapply(names(ari.all.fits[[fit.name]]), function(block.name){
      data.frame(
        fit = fit.name,
        block = block.name,
        ari = ari.all.fits[[fit.name]][[block.name]]
      )
    }))
  }))
  
  return(ari.all.fits)
  
}

###################
# function to extract features by block, fit and component
###################

extract.fit.comp.features <- function(blocks = c("taxa", "GO"), fits, comps = seq(ncomp)){
  fit.features <- sapply(blocks, function(block){
    ls <- lapply(fits, function(fit){
      sapply(comps, function(comp){
        rownames(fit$loadings[[block]])[fit$loadings[[block]][,comp] != 0]
      }, simplify = FALSE)
    })
    ls <- unlist(ls,recursive=FALSE)
  }, simplify = FALSE, USE.NAMES = TRUE)
}

###################
# function to extract feature selection across subsets
###################

extract.fit.features <- function(
    blocks, 
    fits, 
    selection_frequency_th = 1,
    begin_comp = 1,
    end_comp = 1
){
  sapply(blocks, function(block){
    lapply(fits, function(fit){
      ftrs <- sapply(names(fit), function(study){
        rownames(fit[[study]]$loadings[[block]])[
          rowSums(fit[[study]]$loadings[[block]][,begin_comp:end_comp,drop = FALSE] != 0)>0
          ]
      }, simplify = FALSE)
      ftrs <- unlist(ftrs)
      ftr_freqs <- table(ftrs)
      ftr_freqs <- data.frame(
        feature = names(ftr_freqs),
        selection_frequency = c(ftr_freqs)
        )
      ftr_freqs$selection_frequency_group <- "less"
      ftr_freqs$selection_frequency_group[ftr_freqs$selection_frequency > selection_frequency_th] <- "more"
      return(ftr_freqs)
    })
  }, simplify = FALSE, USE.NAMES = TRUE)
}

###################
# function to extract feature selection frequency proportions across subsets
###################

extract.fit.feature.props <- function(
    blocks = c("taxa", "GO"), 
    fits, 
    selection_frequency_th = 1,
    begin_comp = 1,
    end_comp = 1
    ){
  
  freq_props <- sapply(blocks, function(block){
    lapply(fits, function(fit){
      ftrs <- sapply(names(fit), function(study){
        rownames(fit[[study]]$loadings[[block]])[rowSums(fit[[study]]$loadings[[block]][,begin_comp:end_comp,drop = FALSE] != 0)>0]
      }, simplify = FALSE)
      ftrs <- unlist(ftrs)
      ftr_freqs <- table(ftrs)
      freq_freqs <- table(ftr_freqs)
      freq_props <- prop.table(freq_freqs)
      data.frame(
        selection_frequency = as.numeric(names(freq_props)),
        proportion_of_features = c(freq_props)
      )
    })
  }, simplify = FALSE, USE.NAMES = TRUE)
  
  # convert freq_props to a data frame for plotting
  freq_props <- do.call(rbind, lapply(names(freq_props), function(block){
    do.call(rbind, lapply(names(freq_props[[block]]), function(fit){
      data.frame(
        block = block,
        fit = fit,
        selection_frequency = freq_props[[block]][[fit]]$selection_frequency,
        proportion_of_features = freq_props[[block]][[fit]]$proportion_of_features,
        wghtd_freq = freq_props[[block]][[fit]]$selection_frequency * freq_props[[block]][[fit]]$proportion_of_features
      )
    }))
  }))
  
  freq_props$selection_frequency_group <- "less"
  freq_props$selection_frequency_group[freq_props$selection_frequency > selection_frequency_th] <- "more"
  
  return(freq_props)
  
}


######################
##### function to extract predictions and actuals for holdout samples for (s)plsda models
#####################

extract.classifications <- function(predictions, pred.method, pred.distance, dim = ncomp){
  # get classification for all models and studies
  classifications <- lapply(names(predictions), function(prediction){
    lapply(names(predictions[[prediction]]), function(study){
      data.frame(
        study.predictions = predictions[[prediction]][[study]]$pred[[pred.method]][[pred.distance]][,dim],
        study.actuals = predictions[[prediction]][[study]]$Y.test,
        model.name = prediction,
        study.name = study
      )
    })
  })
  # convert to a dataframe
  classifications <- do.call("rbind",do.call("rbind",classifications))
  classifications$IfCorrect <- classifications$study.predictions==classifications$study.actuals
  return(classifications)
}


