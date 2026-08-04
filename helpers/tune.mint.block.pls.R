### This function is for tuning the number of components for a mint.block.pls model

############
# function for finding optimal number of components
############

# The user has to enter all blocks in the list called X and use indY to indicate which is treated as Y
# This keeps things simpler
# will need to spell this out as the mint.block.(s)pls functions allow Y to be entered separately

tune.mint.block.pls <- function(
    X,
    indY,
    study,
    ncomp = 2,
    mode = c('regression')[1],
    design = "full",
    scale = TRUE,
    tol = 1e-06,
    max.iter = 100,
    BPPARAM = SerialParam(),
    seed = NULL
){
  
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
