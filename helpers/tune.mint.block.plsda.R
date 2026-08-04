### This function is for tuning the number of components for a mint.block.plsda model

tune.mint.block.plsda <- function(
    X,
    Y,
    indY,
    study,
    ncomp = 2,
    design = "full",
    scale = TRUE,
    tol = 1e-06,
    max.iter = 100,
    BPPARAM = SerialParam(),
    seed = NULL
){
  
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
