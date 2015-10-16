## S3 generic method
getMeanDVH <-
function(x, cumul=TRUE, purge=3, byPat=TRUE, patID=NULL, structure=NULL,
         interp=c("linear", "spline", "smoothSpl"), fixed=TRUE) {
    UseMethod("getMeanDVH")
}

## for completeness sake - "mean" of just 1 DVH
getMeanDVH.DVHs <-
function(x, cumul=TRUE, purge=3, byPat=TRUE, patID=NULL, structure=NULL,
         interp=c("linear", "spline", "smoothSpl"), fixed=TRUE) {
    interp <- match.arg(interp)

    x <- if(byPat) {
        setNames(list(x), x$structure)
    } else {
        setNames(list(x), x$patID)
    }

    class(x) <- "DVHLst"
    attr(x, which="byPat") <- byPat

    getMeanDVH.DVHLst(x, cumul=cumul, purge=purge, byPat=byPat, patID=patID,
                      structure=structure, interp=interp, fixed=fixed)
}

getMeanDVH.DVHLst <-
function(x, cumul=TRUE, purge=3, byPat=TRUE, patID=NULL, structure=NULL,
         interp=c("linear", "spline", "smoothSpl"), fixed=TRUE) {
    interp <- match.arg(interp)

    ## make sure DVH list is organized as required for byPat
    if(is.null(attributes(x)$byPat) || attributes(x)$byPat != byPat) {
        stop(c("DVH list organization by-patient / by-structure ",
               "either could not be determined or is different from byPat"))
    }

    ## if patIDs are selected, filter them here -> strips DVHLst class
    if(!is.null(patID)) {
        p <- trimWS(patID, side="both")
        x <- Filter(function(y) { 
            if(fixed) {
                any(y$patID %in% p)
            } else {
                any(grepl(paste(p, collapse="|"), y$patID))
            }
        }, x)

        if(length(x) < 1L) { stop("No selected patient found") }
    }

    ## if structures are selected, filter them here
    if(!is.null(structure)) {
        s <- trimWS(structure, side="both")
        x <- Filter(function(y) {
            if(fixed) {
                any(y$structure %in% s)
            } else {
                any(grepl(paste(s, collapse="|"), y$structure))
            }
        }, x)

        if(length(x) < 1L) { stop("No selected structure found") }
    }

    ## get dose range and number of dose nodes
    rangeD <- c(0, max(vapply(x, function(z) {    max(z$dvh[ , "dose"]) }, numeric(1))))
    nodes  <-      max(vapply(x, function(z) { length(z$dvh[ , "dose"]) }, numeric(1)))
        
    ## coarser dose grid for M+SD but with at least 100 nodes
    nodes <- max(100, ceiling(nodes/purge))

    ## average cumulative or differential DVH?
    x <- if(cumul) {
        ## cumulative with linear interpolation
        convertDVH.DVHLst(x, toType="asis", interp=interp,
                          nodes=nodes, rangeD=rangeD, perDose=TRUE)
    } else {
        ## differential with linear interpolation
        convertDVH.DVHLst(x, toType="differential", interp=interp,
                          nodes=nodes, rangeD=rangeD, perDose=TRUE)
    }

    ## extract actual DVH from each DVHs object
    dvhDFL <- if(cumul) {
        ## cumulative DVH
        lapply(x, function(y) {
            data.frame(y$dvh, patID=y$patID, structure=y$structure,
                       stringsAsFactors=FALSE)
        })
    } else {
        ## differential DVH
        lapply(x, function(y) {
            data.frame(y$dvhDiff, patID=y$patID, structure=y$structure,
                       stringsAsFactors=FALSE)
        })
    }
    
    ## combine list to data frame
    dvhDF <- do.call("rbind", dvhDFL)

    ## generate point-wise mean/sd for dose -> aggregate over dose
    if(byPat) {
        dfM  <- aggregate(cbind(volume, volumeRel) ~ patID + dose,
                          data=dvhDF, FUN=mean, na.action=na.pass)
        dfSD <- aggregate(cbind(volume, volumeRel) ~ patID + dose,
                          data=dvhDF, FUN=sd,   na.action=na.pass)
    } else {
        dfM  <- aggregate(cbind(volume, volumeRel) ~ structure + dose,
                          data=dvhDF, FUN=mean, na.action=na.pass)
        dfSD <- aggregate(cbind(volume, volumeRel) ~ structure + dose,
                          data=dvhDF, FUN=sd,   na.action=na.pass)
    }

    ## rename columns in aggregated data frames
    ## dfM stays as is because "volPlot" is used in top layer of the plot
    namesDFM  <- names(dfM)
    namesDFSD <- names(dfSD)
    namesDFM[ namesDFM  == "volume"]    <- "volumeM"
    namesDFSD[namesDFSD == "volume"]    <- "volumeSD"
    namesDFM[ namesDFM  == "volumeRel"] <- "volumeRelM"
    namesDFSD[namesDFSD == "volumeRel"] <- "volumeRelSD"
    names(dfM)  <- namesDFM
    names(dfSD) <- namesDFSD
    
    ## combine mean and sd
    dfMSD <- merge(dfM, dfSD)
    
    ## add M +/- 1SD, 2SD
    dfMSD$volLo1SD <- dfMSD$volumeM -   dfMSD$volumeSD
    dfMSD$volHi1SD <- dfMSD$volumeM +   dfMSD$volumeSD
    dfMSD$volLo2SD <- dfMSD$volumeM - 2*dfMSD$volumeSD
    dfMSD$volHi2SD <- dfMSD$volumeM + 2*dfMSD$volumeSD

    dfMSD$volRelLo1SD <- dfMSD$volumeRelM -   dfMSD$volumeRelSD
    dfMSD$volRelHi1SD <- dfMSD$volumeRelM +   dfMSD$volumeRelSD
    dfMSD$volRelLo2SD <- dfMSD$volumeRelM - 2*dfMSD$volumeRelSD
    dfMSD$volRelHi2SD <- dfMSD$volumeRelM + 2*dfMSD$volumeRelSD

    rownames(dfMSD) <- NULL
    dfMSD
}

## x is a DVH list (1 per id or 1 per structure) of lists
## plots many DVH files
## either for many patients   -> multiple structures per DVH
## or     for many structures -> multiple patients   per DVH
getMeanDVH.DVHLstLst <-
function(x, cumul=TRUE, purge=3, byPat=TRUE, patID=NULL, structure=NULL,
         interp=c("linear", "spline", "smoothSpl"), fixed=TRUE) {
    interp <- match.arg(interp)

    ## re-organize x into by-patient or by-structure form if necessary
    isByPat <- attributes(x)$byPat
    if(is.null(isByPat) || (isByPat != byPat)) {
        x <- reorgByPat(x, byPat=byPat)
    }

    ## if byPat=TRUE and patIDs are selected, filter them here
    if(byPat && !is.null(patID)) {
        p <- trimWS(patID, side="both")
        x <- Filter(function(y) {
            if(fixed) {
                any(y[[1]]$patID %in% p)
            } else {
                any(grepl(paste(p, collapse="|"), y[[1]]$patID))
            }
        }, x)

        if(length(x) < 1L) { stop("No selected patient found") }
        patID <- NULL
    }

    ## if byPat=FALSE and structures are selected, filter them here
    if(!byPat && !is.null(structure)) {
        s <- trimWS(structure, side="both")
        x <- Filter(function(y) {
            if(fixed) {
                any(y[[1]]$structure %in% s)
            } else {
                any(grepl(paste(s, collapse="|"), y[[1]]$structure))
            }
        }, x)

        if(length(x) < 1L) { stop("No selected structure found") }
        structure <- NULL
    }

    resDFL <- Map(getMeanDVH, x, cumul=cumul, purge=purge, byPat=byPat,
                  patID=list(patID), structure=list(structure),
                  interp=interp, fixed=fixed)

    resDF <- do.call("rbind", resDFL)
    rownames(resDF) <- NULL
    resDF
}
