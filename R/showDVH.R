## S3 generic method
showDVH <-
function(x, cumul=TRUE, byPat=TRUE, patID=NULL, structure=NULL,
         rel=TRUE, guessX=TRUE, thresh=1, show=TRUE, ...) {
    UseMethod("showDVH")
}

## plots 1 DVH file for 1 id and 1 structure
showDVH.DVHs <-
function(x, cumul=TRUE, byPat=TRUE, patID=NULL, structure=NULL,
         rel=TRUE, guessX=TRUE, thresh=1, show=TRUE, ...) {
    x <- if(byPat) {
        setNames(list(x), x$structure)
    } else {
        setNames(list(x), x$patID)
    }

    class(x) <- "DVHLst"
    attr(x, which="byPat") <- byPat

    #NextMethod("showDVH")
    showDVH.DVHLst(x, cumul=cumul, byPat=byPat, patID=patID, structure=structure,
                   rel=rel, guessX=guessX, thresh=thresh, show=show)
}

## plots 1 list of DVH objects
## for byPat=TRUE:  1 patient   -> multiple structures
## for byPat=FALSE: 1 structure -> multiple patients
showDVH.DVHLst <-
function(x, cumul=TRUE, byPat=TRUE, patID=NULL, structure=NULL,
         rel=TRUE, guessX=TRUE, thresh=1, show=TRUE, ...) {
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

    ## choose upper x-axis limit (dose) - add 10%
    ## TODO: up to first dose for which all structures reach threshold
    xMax <- if(guessX == 1L) {
        volGEQ <- lapply(x, function(y) {
            y$dvh[ , "dose"][y$dvh[ , "volumeRel"] >= thresh] })
        1.1*max(unlist(volGEQ), na.rm=TRUE)
    } else {
        1.1*c(guessX,
              max(vapply(x, function(y) {
                  max(y$dvh[ , "dose"], na.rm=TRUE) }, numeric(1))))
    }
    
    ## title string
    strTitle <- if(byPat) {
        paste0("patient ",   x[[1]]$patID)
    } else {
        paste0("structure ", x[[1]]$structure)
    }

    ## combine all data frames
    dvhDFL <- if(cumul) {
        ## cumulative DVH
        lapply(x, function(y) {
            if(is.null(y$dvh)) {
                y$dvh <- convertDVH(y$dvhDiff, toType="cumulative",
                                    toDoseUnit="asis")
            }

            data.frame(y$dvh, patID=y$patID, structure=y$structure,
                       stringsAsFactors=FALSE)
        })
    } else {
        lapply(x, function(y) {
            ## differential DVH - create if not yet present
            if(is.null(y$dvhDiff)) {
                y$dvhDiff <- convertDVH(y$dvh, toType="differential",
                                        toDoseUnit="asis")
            }

            data.frame(y$dvhDiff, patID=y$patID, structure=y$structure,
                       stringsAsFactors=FALSE)
        })
    }
    
    dvhDF <- do.call("rbind", dvhDFL)
    
    ## check if relative volume is available if requested
    if(rel && all(is.na(dvhDF$volumeRel))) {
        warning("All relative volumes are missing, will try to show absolute volume")
        rel <- FALSE
    } else if(!rel && all(is.na(dvhDF$volume))) {
        warning("All absolute volumes are missing, will try to show relative volume")
        rel <- FALSE
    }

    ## check if absolute dose is available
    isDoseRel <- if(all(is.na(dvhDF$dose))) {
        warning("All absolute doses are missing, will try to show relative dose")
        dvhDF$dose <- dvhDF$doseRel
        TRUE
    } else {
        FALSE
    }

    diag0 <- if(rel) {                    # relative volume
        if(byPat) {
            ggplot(dvhDF, aes_string(x="dose", y="volumeRel",
                                     colour="structure"))
        } else {
            ggplot(dvhDF, aes_string(x="dose", y="volumeRel",
                                     colour="patID"))
        }
    } else {                             # absolute volume
        if(byPat) {
            ggplot(dvhDF, aes_string(x="dose", y="volume",
                                     colour="structure"))
        } else {
            ggplot(dvhDF, aes_string(x="dose", y="volume",
                                     colour="patID"))
        }
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

    diag1 <- diag0 + geom_line(size=1.5)
    # ggplot2::theme(legend.justification=c(1, 1), legend.position=c(1, 1))
    diag2 <- if(is.finite(xMax)) {
        diag1 + coord_cartesian(xlim=c(0, xMax)) + expand_limits(y=0)
    } else {
        diag1
    }
    
    diag3 <- diag2 + ggtitle(strTitle) +
        theme_bw() +
        scale_y_continuous(expand=c(0, 0.6)) +
        xlab(paste0("Dose [",   doseUnit, "]")) +
        ylab(paste0("Volume [", volUnit,  "]"))

    if(show) {
        print(diag3)
    }

    return(invisible(diag3))
}

## x is a DVH list (1 per id or 1 per structure) of lists
## plots many DVH files
## either for many patients   -> multiple structures per DVH
## or     for many structures -> multiple patients   per DVH
showDVH.DVHLstLst <-
function(x, cumul=TRUE, byPat=TRUE, patID=NULL, structure=NULL,
         rel=TRUE, guessX=TRUE, thresh=1, show=TRUE, ...) {
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
                 guessX=guessX, thresh=thresh, show=show, ...)

    return(invisible(diagL))
}
