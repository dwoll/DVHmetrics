## function for cubic local polynomial smoothing
## dose, dose rel, volume, volume rel, N DVH nodes
getKSmooth <- function(d, dR, v, vR, nodes=NULL, rangeD=NULL) {
    nodes <- max(nodes, length(d))
    
    ## smooth
    bwDV  <- try(KernSmooth::dpill(d, v))
    bwDVR <- try(KernSmooth::dpill(d, vR))
    smDV  <- try(KernSmooth::locpoly(d, v,  bandwidth=bwDV,
                                     gridsize=nodes, degree=3))
    smDVR <- try(KernSmooth::locpoly(d, vR, bandwidth=bwDVR,
                                     gridsize=nodes, degree=3))
    
    ## dose
    dose <- if(!inherits(smDV, "try-error")) {
        smDV$x
    } else if(!inherits(smDVR, "try-error")) {
        smDVR$x
    } else {
        NA_real_
    }

    ## dose rel -> just use equally spaced grid points as is done in 
    doseRel <- if(!any(is.na(dR))) {
        rangeD <- range(dR)
        seq(rangeD[1], rangeD[2], length.out=nodes)
    } else {
        NA_real_
    }

    volume <- if(!inherits(smDV, "try-error")) {
        smDV$y
    } else { NA_real_ }

    volumeRel <- if(!inherits(smDVR, "try-error")) {
        smDVR$y
    } else { NA_real_ }

    return(cbind(dose=dose, doseRel=doseRel, volume=volume, volumeRel=volumeRel))
}

## function for cubic local polynomial smoothing
## dose, dose rel, volume, volume rel, N DVH nodes
getSmoothSpl <- function(d, dR, v, vR, nodes=NULL, rangeD=NULL) {
    nodes <- max(nodes, length(d))
    if(is.null(rangeD)) { rangeD  <- range(d) }

    ## smooth
    smDV  <- try(smooth.spline(d, v))
    smDVR <- try(smooth.spline(d, vR))

    ## dose, dose rel -> just use equally spaced grid points
    dose    <- seq(rangeD[1],  rangeD[2], length.out=nodes)
    doseRel <- if(!any(is.na(dR))) {
        rangeDR <- range(dR)
        seq(rangeDR[1], rangeDR[2], length.out=nodes)
    } else {
        NA_real_
    }

    volume <- if(!inherits(smDV, "try-error")) {
        predict(smDV,  dose)$y
    } else { NA_real_ }

    volumeRel <- if(!inherits(smDVR, "try-error")) {
        predict(smDVR, dose)$y
    }

    return(cbind(dose=dose, doseRel=doseRel, volume=volume, volumeRel=volumeRel))
}

## function for cubic spline interpolation
## dose, dose rel, volume, volume rel, N DVH nodes
getInterpSpl <- function(d, dR, v, vR, nodes=NULL, rangeD=NULL) {
    nodes <- max(nodes, length(d))
    if(is.null(rangeD)) { rangeD  <- range(d) }
    
    ## interpolation
    smDV  <- try(splinefun(d, v,  method="fmm"))
    smDVR <- try(splinefun(d, vR, method="fmm"))
    
    ## dose, dose rel -> just use equally spaced grid points
    dose    <- seq(rangeD[1],  rangeD[2], length.out=nodes)
    doseRel <- if(!any(is.na(dR))) {
        rangeDR <- range(dR)
        seq(rangeDR[1], rangeDR[2], length.out=nodes)
    } else {
        NA_real_
    }

    volume <- if(!inherits(smDV, "try-error")) {
        smDV(dose)
    } else {
        NA_real_
    }

    volumeRel <- if(!inherits(smDVR, "try-error")) {
        smDVR(dose)
    } else {
        NA_real_
    }

    return(cbind(dose=dose, doseRel=doseRel, volume=volume, volumeRel=volumeRel))
}

## function for linear interpolation
## dose, dose rel, volume, volume rel, N DVH nodes
getInterpLin <- function(d, dR, v, vR, nodes=NULL, rangeD=NULL) {
    nodes <- max(nodes, length(d))
    if(is.null(rangeD)) { rangeD  <- range(d) }
    
    ## interpolation
    smDV  <- try(approxfun(d, v,  method="linear", rule=2, ties=max))
    smDVR <- try(approxfun(d, vR, method="linear", rule=2, ties=max))
    
    ## dose, dose rel -> just use equally spaced grid points
    dose    <- seq(rangeD[1],  rangeD[2], length.out=nodes)
    doseRel <- if(!any(is.na(dR))) {
        rangeDR <- range(dR)
        seq(rangeDR[1], rangeDR[2], length.out=nodes)
    } else {
        NA_real_
    }

    volume <- if(!inherits(smDV, "try-error")) {
        smDV(dose)
    } else {
        NA_real_
    }

    volumeRel <- if(!inherits(smDVR, "try-error")) {
        smDVR(dose)
    } else {
        NA_real_
    }

    return(cbind(dose=dose, doseRel=doseRel, volume=volume, volumeRel=volumeRel))
}
