## plots 1 list of average DVHs +/1 SD
showMeanDVH <-
function(x, byPat=TRUE, patID=NULL, structure=NULL, rel=TRUE,
         guessX=TRUE, thresh=1, show=TRUE, fixed=TRUE, showSD=TRUE, color=TRUE) {
    guessX <- as.numeric(guessX)

    ## combine all mean DVHs in x - might a list or just 1 data frame
    dvhDF <- if(is.data.frame(x)) {
        x
    } else if(is.list(x)) {
        if(all(vapply(x, is.data.frame, logical(1)))) {
            do.call("rbind", x)
        } else {
            stop("x must be a data frame or a list of data frames")
        }
    } else {
        stop("x must be a data frame or a list of data frames")
    }

    ## if patIDs are selected, filter them here -> strips DVHLst class
    if(!is.null(patID)) {
        p <- trimWS(patID, side="both")
        dvhDF <- if(fixed) {
            dvhDF[dvhDF$patID %in% p, ]
        } else {
            dvhDF[grepl(paste(p, collapse="|"), dvhDF$patID), ]
        }

        if(nrow(dvhDF) < 1L) { stop("No selected patient found") }
    }

    ## if structures are selected, filter them here
    if(!is.null(structure)) {
        s <- trimWS(structure, side="both")
        dvhDF <- if(fixed) {
           dvhDF[dvhDF$structure %in% s, ]
        } else {
           dvhDF[grepl(paste(s, collapse="|"), dvhDF$structure), ]
        }

        if(nrow(dvhDF) < 1L) { stop("No selected structure found") }
    }

    ## choose upper x-axis limit (dose) - add 10%
    ## TODO: up to first dose for which all structures reach threshold
    xMax <- if((guessX == 1L) && !all(is.na(dvhDF[ , "volumeRelMEAN"]))) {
        volGEQ <- dvhDF[ , "dose"][dvhDF[ , "volumeRelMEAN"] >= thresh]
        1.1*max(volGEQ, na.rm=TRUE)
    } else {
        1.1*max(c(guessX, dvhDF$dose))
    }
    
    ## set plot volume to absolute or relative
    ## and choose upper y-axis limit (volume) - add 3%
    if(rel) {
        ## relative volume - check if available
        if(all(is.na(dvhDF$volumeRelMEAN))) {
            warning("All relative volumes are missing, will try to show absolute volume")
            yMax <- 1.03*max(dvhDF$volumeMEAN, na.rm=TRUE)
            rel  <- FALSE
            dvhDF$volPlot <- dvhDF$volumeMEAN
            dvhDF$volSD   <- dvhDF$volumeSD
        } else {
            yMax <- min(c(103, 1.03*max(dvhDF$volumeRelMEAN, na.rm=TRUE)))
            dvhDF$volPlot <- dvhDF$volumeRelMEAN
            dvhDF$volSD   <- dvhDF$volumeRelSD
        }
    } else {
        ## absolute volume - check if available
        if(all(is.na(dvhDF$volumeMEAN))) {
            warning("All absolute volumes are missing, will try to show relative volume")
            yMax <- min(c(103, 1.03*max(dvhDF$volumeRelMEAN, na.rm=TRUE)))
            rel  <- TRUE
            dvhDF$volPlot <- dvhDF$volumeRelMEAN
            dvhDF$volSD   <- dvhDF$volumeRelSD
        } else {
            yMax <- 1.03*max(dvhDF$volumeMEAN, na.rm=TRUE)
            dvhDF$volPlot <- dvhDF$volumeMEAN
            dvhDF$volSD   <- dvhDF$volumeSD
        }
    }

    ## add +/- 1 SD, 2SD areas
    dvhDF$volLo1SD <- dvhDF$volPlot -   dvhDF$volSD
    dvhDF$volHi1SD <- dvhDF$volPlot +   dvhDF$volSD
    dvhDF$volLo2SD <- dvhDF$volPlot - 2*dvhDF$volSD
    dvhDF$volHi2SD <- dvhDF$volPlot + 2*dvhDF$volSD

    ## check if absolute dose is available
    isDoseRel <- if(all(is.na(dvhDF$dose))) {
        warning("All absolute doses are missing, will try to show relative dose")
        dvhDF$dose <- dvhDF$doseRel
        TRUE
    } else {
        FALSE
    }

    ## title string
    strTitle <- if(byPat) {
        paste0("structures ", paste(sort(unique(dvhDF$structure)), collapse=", "))
    } else {
        paste0("patients ",   paste(sort(unique(dvhDF$patID)),     collapse=", "))
    }

    ## do we have many structure/ID categories? -> more legend columns
    nCateg <- if(byPat) {
        length(unique(dvhDF$patID))
    } else {
        length(unique(dvhDF$structure))
    }

    ## 15 entries per legend column
    nLegendCols <- (nCateg %/% 15) + 1            # number of categories

    ## do the actual plotting
    diag <- if(byPat) {
        ggplot(dvhDF, aes_string(x="dose", y="volPlot",
                                 group="patID", color="patID", fill="patID")) +
            facet_grid(as.formula("structure ~ ."))
    } else {
        ggplot(dvhDF, aes_string(x="dose", y="volPlot",
                                 group="structure", color="structure", fill="structure")) +
            facet_grid(as.formula("patID ~ ."))
    }

    ## rescale x-axis, y-axis
    diag <- if(is.finite(xMax)) {
        if(is.finite(yMax)) {
            diag + coord_cartesian(xlim=c(0, xMax), ylim=c(0, yMax))
        } else {
            diag + coord_cartesian(xlim=c(0, xMax))
        }
    } else {
        if(is.finite(yMax)) {
            diag + coord_cartesian(ylim=c(0, yMax))
        } else {
            diag
        }
    }
    
    ## point-wise 1-SD, 2-SD shaded areas
    diag <- if(showSD) {
        diag +
        geom_ribbon(aes_string(x="dose", ymin="volLo1SD", ymax="volHi1SD"),
                    alpha=0.25, linetype="blank") +
        geom_ribbon(aes_string(x="dose", ymin="volLo2SD", ymax="volHi2SD"),
                    alpha=0.2, linetype="blank")
    } else {
        diag
    }
    
    ## im grey -> linetype
    diag <- if(color) {
        diag +
            geom_line(aes_string(x="dose", y="volPlot"), size=1.2)
    } else {
        if(byPat) {
            diag +
                geom_line(aes_string(x="dose", y="volPlot", linetype="patID"), size=1.2) +
                scale_color_grey(start=0, end=0.5) +
                scale_fill_grey()
        } else {
            diag +
                geom_line(aes_string(x="dose", y="volPlot", linetype="structure"), size=1.2) +
                scale_color_grey(start=0, end=0.5) +
                scale_fill_grey()
        }
    }

    ## final diagram
    diag <- diag + theme_bw() +
        expand_limits(y=0) +                      # make sure 0 is included
        scale_y_continuous(expand=c(0, 0.6)) +    # no space below y=0
        xlab("Dose") +
        ylab("Volume") +
        guides(color=guide_legend(ncol=nLegendCols))

    if(show) {
        print(diag)
    }

    return(invisible(diag))
}
