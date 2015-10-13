## S3 generic method
showDVH <-
function(x, cumul=TRUE, byPat=TRUE, patID=NULL, structure=NULL,
         rel=TRUE, guessX=TRUE, thresh=1, addMSD=FALSE, show=TRUE, ...) {
    UseMethod("showDVH")
}

## plots 1 DVH file for 1 id and 1 structure
showDVH.DVHs <-
function(x, cumul=TRUE, byPat=TRUE, patID=NULL, structure=NULL,
         rel=TRUE, guessX=TRUE, thresh=1, addMSD=FALSE, show=TRUE, ...) {
    x <- if(byPat) {
        setNames(list(x), x$structure)
    } else {
        setNames(list(x), x$patID)
    }

    class(x) <- "DVHLst"
    attr(x, which="byPat") <- byPat

    showDVH.DVHLst(x, cumul=cumul, byPat=byPat, patID=patID, structure=structure,
                   rel=rel, guessX=guessX, thresh=thresh, addMSD=addMSD, show=show)
}

## plots 1 list of DVH objects
## for byPat=TRUE:  1 patient   -> multiple structures
## for byPat=FALSE: 1 structure -> multiple patients
showDVH.DVHLst <-
function(x, cumul=TRUE, byPat=TRUE, patID=NULL, structure=NULL,
         rel=TRUE, guessX=TRUE, thresh=1, addMSD=FALSE, show=TRUE, ...) {
    ## make sure DVH list is organized as required for byPat
    if(is.null(attributes(x)$byPat) || attributes(x)$byPat != byPat) {
        stop(c("DVH list organization by-patient / by-structure ",
               "either could not be determined or is different from byPat"))
    }

    guessX <- as.numeric(guessX)

    ## if patIDs are selected, filter them here
    if(!is.null(patID)) {
        p <- trimWS(patID, side="both")
        x <- Filter(function(y) any(grepl(paste(p, collapse="|"), y$patID, ...)), x)
        if(length(x) < 1L) { stop("No selected patient found") }
    }

    ## if structures are selected, filter them here
    if(!is.null(structure)) {
        s <- trimWS(structure, side="both")
        x <- Filter(function(y) any(grepl(paste(s, collapse="|"), y$structure, ...)), x)
        if(length(x) < 1L) { stop("No selected structure found") }
    }

    ## combine all data frames
    dvhDFL <- if(cumul) {
        ## cumulative DVH
        lapply(x, function(y) {
            data.frame(y$dvh, patID=y$patID, structure=y$structure,
                       stringsAsFactors=FALSE)
        })
    } else {
        lapply(x, function(y) {
            ## differential DVH - create first
            y$dvhDiff <- convertDVH(y$dvh, toType="differential",
                                    toDoseUnit="asis", perDose=TRUE)

            data.frame(y$dvhDiff, patID=y$patID, structure=y$structure,
                       stringsAsFactors=FALSE)
        })
    }
    
    dvhDF <- do.call("rbind", dvhDFL)

    ## choose upper x-axis limit (dose) - add 10%
    ## TODO: up to first dose for which all structures reach threshold
    xMax <- if(guessX == 1L) {
        volGEQ <- lapply(x, function(y) {
            y$dvh[ , "dose"][y$dvh[ , "volumeRel"] >= thresh] })
        1.1*max(unlist(volGEQ), na.rm=TRUE)
    } else {
        1.1*max(c(guessX, dvhDF$dose))
    }
    
    ## choose upper y-axis limit (volume) - add 3%
    yMaxAbs <- 1.03*max(dvhDF$volume)
    yMaxRel <- 1.03*max(dvhDF$volumeRel)

    ## title string
    strTitle <- if(byPat) {
        paste0("patient ",   x[[1]]$patID)
    } else {
        paste0("structure ", x[[1]]$structure)
    }
    
    ## check if relative volume is available if requested
    yMax <- if(rel) {
        yMaxRel
    } else {
        yMaxAbs
    }

    if(rel && all(is.na(dvhDF$volumeRel))) {
        warning("All relative volumes are missing, will try to show absolute volume")
        yMax <- yMaxAbs
        rel  <- FALSE
    } else if(!rel && all(is.na(dvhDF$volume))) {
        warning("All absolute volumes are missing, will try to show relative volume")
        yMax <- yMaxRel
        rel  <- TRUE
    }

    ## check if absolute dose is available
    isDoseRel <- if(all(is.na(dvhDF$dose))) {
        warning("All absolute doses are missing, will try to show relative dose")
        dvhDF$dose <- dvhDF$doseRel
        TRUE
    } else {
        FALSE
    }

    ## check if point-wise volume mean and sd should be plotted
    if(addMSD) {
        ## generate point-wise mean/sd and merge back to data
        if(byPat) {
            dfM  <- aggregate(cbind(volume, volumeRel) ~ patID + dose, data=dvhDF, FUN=mean)
            dfSD <- aggregate(cbind(volume, volumeRel) ~ patID + dose, data=dvhDF, FUN=sd)
            names(dfM)  <- c("patID", "dose", "volM",  "volRelM")
            names(dfSD) <- c("patID", "dose", "volSD", "volRelSD")
        } else {
            dfM  <- aggregate(cbind(volume, volumeRel) ~ structure + dose, data=dvhDF, FUN=mean)
            dfSD <- aggregate(cbind(volume, volumeRel) ~ structure + dose, data=dvhDF, FUN=sd)
            names(dfM)  <- c("structure", "dose", "volM",  "volRelM")
            names(dfSD) <- c("structure", "dose", "volSD", "volRelSD")
        }

        dfMSD <- merge(dfM, dfSD)
        dfMSD$volLo1SD    <- dfMSD$volM    -   dfMSD$volSD
        dfMSD$volHi1SD    <- dfMSD$volM    +   dfMSD$volSD
        dfMSD$volLo2SD    <- dfMSD$volM    - 2*dfMSD$volSD
        dfMSD$volHi2SD    <- dfMSD$volM    + 2*dfMSD$volSD
        dfMSD$volRelLo1SD <- dfMSD$volRelM -   dfMSD$volRelSD
        dfMSD$volRelHi1SD <- dfMSD$volRelM +   dfMSD$volRelSD
        dfMSD$volRelLo2SD <- dfMSD$volRelM - 2*dfMSD$volRelSD
        dfMSD$volRelHi2SD <- dfMSD$volRelM + 2*dfMSD$volRelSD
        #dfMSD$doseRel     <- dfMSD$dose
        #dfMSD$volume      <- dfMSD$volM
        #dfMSD$volumeRel   <- dfMSD$volRelM
        #dvhDF <- merge(dvhDF, dfMSD, all.x=TRUE)
    }

    volUnit <- if(rel) {
        "%"
    } else {
        x[[1]]$volumeUnit
    }

    doseUnit <- if(isDoseRel) {
        "%"
    } else {
        x[[1]]$doseUnit
    }

    ## plot - empty base layer
    diag <- ggplot()

    ## point-wise 1-SD, 2-SD shaded areas?
    diag <- if(addMSD) {
        if(rel) {
            diag +
            geom_ribbon(data=dfMSD,
                        aes_string(x="dose", ymin="volRelLo1SD", ymax="volRelHi1SD"),
                        alpha=0.25, linetype="blank") +
            geom_ribbon(data=dfMSD,
                        aes_string(x="dose", ymin="volRelLo2SD", ymax="volRelHi2SD"),
                        alpha=0.2, linetype="blank")
        } else {
            diag +
            geom_ribbon(data=dfMSD,
                        aes_string(x="dose", ymin="volLo1SD", ymax="volHi1SD"),
                        alpha=0.25, linetype="blank") +
            geom_ribbon(data=dfMSD,
                        aes_string(x="dose", ymin="volLo2SD", ymax="volHi2SD"),
                        alpha=0.3, linetype="blank")
        }
    } else {
        diag
    }
    
    ## actual DVHs
    diag <- if(rel) {                    # relative volume
        if(byPat) {
            diag + geom_line(data=dvhDF,
                             aes_string(x="dose", y="volumeRel", color="structure"),
                             size=1.2)
        } else {
            diag + geom_line(data=dvhDF,
                             aes_string(x="dose", y="volumeRel", color="patID"),
                             size=1.2)
        }
    } else {
        if(byPat) {
            diag + geom_line(data=dvhDF,
                             aes_string(x="dose", y="volume", color="structure"),
                             size=1.2)
        } else {
            diag + geom_line(data=dvhDF,
                             aes_string(x="dose", y="volume", color="patID"),
                             size=1.2)
        }
    }

    ## rescale x-axis
    diag <- if(is.finite(xMax)) {
        diag + coord_cartesian(xlim=c(0, xMax)) + expand_limits(y=0)
    } else {
        diag
    }
    
    ## rescale y-axsis
    diag <- if(is.finite(yMax)) {
        diag + coord_cartesian(ylim=c(0, yMax))
    } else {
        diag
    }

    ## point-wise mean DVH?
    diag <- if(addMSD) {
        if(rel) {
            diag + geom_line(data=dfMSD,
                             aes_string(x="dose", y="volRelM"),
                             color="black", size=1.2)
        } else {
            diag + geom_line(data=dfMSD,
                             aes_string(x="dose", y="volM"),
                             color="black", size=1.2)
        }
    } else {
        diag
    }

    diag <- diag + ggtitle(strTitle) +
        theme_bw() +
        scale_y_continuous(expand=c(0, 0.6)) +
        xlab(paste0("Dose [",   doseUnit, "]")) +
        ylab(paste0("Volume [", volUnit,  "]"))

    if(show) {
        print(diag)
    }

    return(invisible(diag))
}

## x is a DVH list (1 per id or 1 per structure) of lists
## plots many DVH files
## either for many patients   -> multiple structures per DVH
## or     for many structures -> multiple patients   per DVH
showDVH.DVHLstLst <-
function(x, cumul=TRUE, byPat=TRUE, patID=NULL, structure=NULL,
         rel=TRUE, guessX=TRUE, thresh=1, addMSD=FALSE, show=TRUE, ...) {
    ## re-organize x into by-patient or by-structure form if necessary
    isByPat <- attributes(x)$byPat
    if(is.null(isByPat) || (isByPat != byPat)) {
        x <- reorgByPat(x, byPat=byPat)
    }

    ## if byPat=TRUE and patIDs are selected, filter them here
    if(byPat && !is.null(patID)) {
        p <- trimWS(patID, side="both")
        x <- Filter(function(y) any(grepl(paste(p, collapse="|"), y[[1]]$patID, ...)), x)
        if(length(x) < 1L) { stop("No selected patient found") }
        patID <- NULL
    }

    ## if byPat=FALSE and structures are selected, filter them here
    if(!byPat && !is.null(structure)) {
        s <- trimWS(structure, side="both")
        x <- Filter(function(y) any(grepl(paste(s, collapse="|"), y[[1]]$structure, ...)), x)
        if(length(x) < 1L) { stop("No selected structure found") }
        structure <- NULL
    }

    diagL <- Map(showDVH, x, cumul=cumul, byPat=byPat, rel=rel,
                 patID=list(patID), structure=list(structure),
                 guessX=guessX, thresh=thresh, addMSD=addMSD, show=show, ...)

    return(invisible(diagL))
}
