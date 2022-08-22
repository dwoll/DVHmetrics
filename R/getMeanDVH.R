## S3 generic method
getMeanDVH <-
function(x, fun=list(mean=mean, median=median, sd=sd),
         cumul=TRUE, thin=1, byPat=TRUE, patID=NULL, structure=NULL,
         fixed=TRUE, returnDVHObj=FALSE) {
    UseMethod("getMeanDVH")
}

## for completeness' sake - "mean" of just 1 DVH
getMeanDVH.DVHs <-
function(x, fun=list(mean=mean, median=median, sd=sd),
         cumul=TRUE, thin=1, byPat=TRUE, patID=NULL, structure=NULL,
         fixed=TRUE, returnDVHObj=FALSE) {

    x <- if(byPat) {
        setNames(list(x), x$structure)
    } else {
        setNames(list(x), x$patID)
    }

    class(x) <- "DVHLst"
    attr(x, which="byPat") <- byPat

    getMeanDVH.DVHLst(x, fun=fun, cumul=cumul, thin, byPat=byPat,
                      patID=patID, structure=structure, fixed=fixed,
                      returnDVHObj=returnDVHObj)
}

getMeanDVH.DVHLst <-
function(x, fun=list(mean=mean, median=median, sd=sd),
         cumul=TRUE, thin=1, byPat=TRUE, patID=NULL, structure=NULL,
         fixed=TRUE, returnDVHObj=FALSE) {

    extract_info <- function(comp) {
        vals <- sapply(x, function(z) { z[[comp]] })
        paste(unique(unname(vals)), collapse="_")
    }

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
    nodes <- max(100, ceiling(nodes/thin))

    ## average cumulative or differential DVH?
    x <- if(cumul) {
        ## cumulative with linear interpolation
        convertDVH.DVHLst(x, toType="asis", interp="linear",
                          nodes=nodes, rangeD=rangeD, perDose=TRUE)
    } else {
        ## differential with linear interpolation
        convertDVH.DVHLst(x, toType="differential", interp="linear",
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
    getAggr <- function(fun, symbol) {
        dat <- if(byPat) {
            aggregate(cbind(volume, volumeRel) ~ patID + dose,
                      data=dvhDF, FUN=fun, na.action=na.pass)
        } else {
            aggregate(cbind(volume, volumeRel) ~ structure + dose,
                      data=dvhDF, FUN=fun, na.action=na.pass)
        }

        ## rename columns in aggregated data frames
        namesDat <- names(dat)
        namesDat[namesDat == "volume"]    <- paste0("volume",    symbol)
        namesDat[namesDat == "volumeRel"] <- paste0("volumeRel", symbol)
        names(dat) <- namesDat
        cbind(dat, doseRel=(dat$dose / max(dat$dose)))
    }

    ## get all point-wise estimates
    if(!returnDVHObj) {
        dfL <- Map(getAggr, fun, toupper(names(fun)))
        ## combine point-wise estimates

        dfMSD <- Reduce(merge, dfL)
        rownames(dfMSD) <- NULL

        ## add information about original patIDs / structures
        if(byPat) {
            dfMSD$structure <- abbreviate(paste(sort(unique(dvhDF$structure)), collapse="_"),
                                          minlength=20)
        } else {
            dfMSD$patID     <- abbreviate(paste(sort(unique(dvhDF$patID)),     collapse="_"),
                                          minlength=20)
        }

        dfMSD
    } else {
        if(length(fun) > 1L) {
            warning("Only first element of 'fun' will be used")
        }

        dvh  <- getAggr(fun[[1L]], "")
        DVHs <- if(byPat) {
            list(dvh       =data.matrix(dvh[ , c("dose", "doseRel", "volume", "volumeRel")]),
                 patName   =extract_info("patName"),
                 patID     =extract_info("patID"),
                 date      =extract_info("date"),
                 DVHtype   =extract_info("DVHtype"),
                 plan      =extract_info("plan"),
                 structure =abbreviate(paste(sort(unique(dvhDF$structure)), collapse="_"),
                                       minlength=20),
                 structVol =as.numeric(extract_info("structVol")),
                 doseUnit  =extract_info("doseUnit"),
                 volumeUnit=extract_info("volumeUnit"),
                 doseMin   =NA_real_,
                 doseMax   =NA_real_,
                 doseRx    =NA_real_,
                 isoDoseRx =NA_real_,
                 doseAvg   =NA_real_,
                 doseMed   =NA_real_,
                 doseMode  =NA_real_,
                 doseSD    =NA_real_)
        } else {
            list(dvh       =data.matrix(dvh[ , c("dose", "doseRel", "volume", "volumeRel")]),
                 patName   =extract_info("patName"),
                 patID     =abbreviate(paste(sort(unique(dvhDF$patID)),     collapse="_"),
                                       minlength=20),
                 date      =extract_info("date"),
                 DVHtype   =extract_info("DVHtype"),
                 plan      =extract_info("plan"),
                 structure =extract_info("structure"),
                 structVol =as.numeric(extract_info("structVol")),
                 doseUnit  =extract_info("doseUnit"),
                 volumeUnit=extract_info("volumeUnit"),
                 doseMin   =NA_real_,
                 doseMax   =NA_real_,
                 doseRx    =NA_real_,
                 isoDoseRx =NA_real_,
                 doseAvg   =NA_real_,
                 doseMed   =NA_real_,
                 doseMode  =NA_real_,
                 doseSD    =NA_real_)
        }

        rownames(DVHs$dvh) <- NULL
        ## set class
        class(DVHs) <- "DVHs"
        DVHs
    }

}

## x is a DVH list (1 per id or 1 per structure) of lists
## plots many DVH files
## either for many patients   -> multiple structures per DVH
## or     for many structures -> multiple patients   per DVH
getMeanDVH.DVHLstLst <-
function(x, fun=list(mean=mean, median=median, sd=sd),
         cumul=TRUE, thin=1, byPat=TRUE, patID=NULL, structure=NULL,
         fixed=TRUE, returnDVHObj=FALSE) {

    extract_info <- function(comp) {
        vals <- sapply(x, function(z) { z[[comp]] })
        paste(unique(unname(vals)), collapse="_")
    }

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
                any(y[[1L]]$patID %in% p)
            } else {
                any(grepl(paste(p, collapse="|"), y[[1L]]$patID))
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
                any(y[[1L]]$structure %in% s)
            } else {
                any(grepl(paste(s, collapse="|"), y[[1L]]$structure))
            }
        }, x)

        if(length(x) < 1L) { stop("No selected structure found") }
        structure <- NULL
    }

    resDFL <- Map(getMeanDVH, x, fun=list(fun), cumul=cumul, thin,
                  byPat=byPat, patID=list(patID), structure=list(structure),
                  fixed=fixed, returnDVHObj=returnDVHObj)

    if(!returnDVHObj) {
        resDF <- do.call("rbind", resDFL)
        rownames(resDF) <- NULL
        resDF
    } else {
        class(resDFL) <- "DVHLst"
        attr(resDFL, which="byPat") <- !byPat
        resDFL
    }
}

