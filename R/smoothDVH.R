## function for cubic local polynomial smoothing
## dose, dose rel, volume, volume rel, N DVH nodes
getKSmooth <- function(d, dR, v, vR, nodes=NULL) {
    nodes <- max(nodes, length(d))
    bwDV  <- KernSmooth::dpill(d, v)
    bwDVR <- KernSmooth::dpill(d, vR)
    smDV  <- KernSmooth::locpoly(d, v,  bandwidth=bwDV,
                                 gridsize=nodes, degree=3)
    smDVR <- KernSmooth::locpoly(d, vR, bandwidth=bwDVR,
                                 gridsize=nodes, degree=3)
    
    ## dose rel -> just use equally spaced grid points
    doseRel <- if(!any(is.na(dR))) {
        rangeDR <- range(dR)
        seq(rangeDR[1], rangeDR[2], length.out=nodes)
    } else {
        NA_real_
    }

    return(cbind(dose=smDV$x, doseRel=doseRel, volume=smDV$y, volumeRel=smDVR$y))
}

## function for cubic local polynomial smoothing
## dose, dose rel, volume, volume rel, N DVH nodes
getSmoothSpl <- function(d, dR, v, vR, nodes=NULL) {
    nodes <- max(nodes, length(d))

    ## smooth
    smDV  <- smooth.spline(d, v)
    smDVR <- smooth.spline(d, vR)

    ## dose, dose rel -> just use equally spaced grid points
    rangeD  <- range(d)
    dose    <- seq(rangeD[1],  rangeD[2],  length.out=nodes)
    doseRel <- if(!any(is.na(dR))) {
        rangeDR <- range(dR)
        seq(rangeDR[1], rangeDR[2], length.out=nodes)
    } else {
        NA_real_
    }
    volume    <- predict(smDV,  dose)$y
    volumeRel <- predict(smDVR, dose)$y

    return(cbind(dose=dose, doseRel=doseRel, volume=volume, volumeRel=volumeRel))
}

## function for cubic spline interpolation
## dose, dose rel, volume, volume rel, N DVH nodes
getInterpSpl <- function(d, dR, v, vR, nodes=NULL) {
    nodes <- max(nodes, length(d))
    
    ## interpolation
    smDV  <- splinefun(d, v,  method="fmm")
    smDVR <- splinefun(d, vR, method="fmm")
    
    ## dose, dose rel -> just use equally spaced grid points
    rangeD  <- range(d)
    dose    <- seq(rangeD[1],  rangeD[2],  length.out=nodes)
    doseRel <- if(!any(is.na(dR))) {
        rangeDR <- range(dR)
        seq(rangeDR[1], rangeDR[2], length.out=nodes)
    } else {
        NA_real_
    }
    volume    <- smDV(dose)
    volumeRel <- smDVR(dose)

    return(cbind(dose=dose, doseRel=doseRel, volume=volume, volumeRel=volumeRel))
}
