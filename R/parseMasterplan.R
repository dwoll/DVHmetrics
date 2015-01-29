#####---------------------------------------------------------------------------
## parse character vector from Eclipse DVH file
parseMasterplan <- function(x, planInfo=FALSE) {
    planInfo <- as.character(planInfo)

    ## function to extract one information element from a number of lines
    ## make sure only first : is matched -> not greedy
    getElem <- function(pattern, ll, trim=TRUE, iCase=FALSE, collWS=TRUE) {
        line <- ll[grep(pattern, ll)]
        elem <- sub("^.+?:[[:blank:]]+([[:alnum:][:punct:][:blank:]]+$)", "\\1",
                    line, ignore.case=iCase, perl=TRUE)
        elem <- if(trim) {
            trimWS(elem, side="both")
        } else {
            elem
        }

        if(collWS) {
            collWS(elem)
        } else {
            elem
        }
    }

    getDose <- function(pattern, ll, doseRx) {
        line <- ll[grep(pattern, ll)]
        elem <- sub("^.+?:[[:blank:]]+([[:alnum:][:punct:]]+[[:blank:]]*$)", "\\1", line, perl=TRUE)
        num  <- trimWS(elem)
        if(grepl("\\[%\\]", line)) {
            ## relative dose
            doseRx * as.numeric(num)/100
        } else {
            as.numeric(num)
        }
    }

    getDoseUnit <- function(ll) {
        line <- ll[grep("^(Absolute|Relative) dose - ", ll, ignore.case=TRUE)]
        elem <- sub("^.+\\[(GY|CGY)\\]", "\\1", line, perl=TRUE, ignore.case=TRUE)
        toupper(trimWS(elem))
    }

    getVolUnit <- function(ll) {
        line <- ll[grep("^(Absolute|Relative) volume - ", ll, ignore.case=TRUE)]
        elem <- sub("^.+\\[(ccm|%)\\]", "\\1", line, perl=TRUE, ignore.case=TRUE)
        toupper(trimWS(elem))
    }

    getDVHtype <- function(ll) {
        line <- ll[grep("mode", ll)]
        elem <- sub("^(Cumulative|Differential) mode", "\\1", line, perl=TRUE, ignore.case=TRUE)
        tolower(trimWS(elem))
    }

    sStart <- grep("^Dose: [[:alnum:]]", x)   # start of sections
    sLen   <- diff(c(sStart, length(x)+1))         # length of sections
    if((length(sLen) < 1L) || all(sLen < 1L)) {
        stop("No structures found")
    }

    structList <- Map(function(start, len) x[start:(start+len-1)], sStart, sLen)

    ## extract file header and header info
    header    <- x[seq_len(sStart[1]-1)]      # header
    patName   <- NA_character_
    patID     <- getElem("^Case:",  header)   # patient id
    plan      <- getElem("^Plan:",  header)   # treatment plan
    quadrant  <- NA_character_
    DVHdate   <- NA_character_
    DVHtype   <- getDVHtype(header)
    isoDoseRx <- 100
    doseRx    <- NA_real_
    doseUnit  <- getDoseUnit(header)
    if(!grepl("^(GY|CGY)$", doseUnit)) {
        warning("Could not determine dose measurement unit")
        doseUnit <- NA_character_
    }

    volumeUnit <- getVolUnit(header)
    volumeUnit <- if(grepl("^CCM", volumeUnit)) {
        "CC"
    } else if(grepl("^%", volumeUnit)) {
        "PERCENT"
    } else {
        warning("Could not determine volume measurement unit")
        NA_character_
    }

    ## extract DVH from one structure section and store in a list
    ## with DVH itself as a data frame
    getDVH <- function(strct, info) {
        ## extract information from info list
        doseRx    <- info$doseRx
        isoDoseRx <- info$isoDoseRx

        ## extract structure, volume, dose min, max, mean, median and sd
        structure <- getElem("^ROI*:", strct)
        structVol <- NA_real_

        ## find DVH
        ## DVH column headers
        colHead  <- grep("DOSE.+VOLUME", strct, ignore.case=TRUE, perl=TRUE)
        dvhStart <- colHead+1                 # first numeric line of DVH
        dvhLen   <- length(strct) - dvhStart
        if((length(dvhLen) < 1L) || dvhLen < 1L) {
            stop("No DVH data found")
        }

        ## column headers
        vars1 <- unlist(strsplit(strct[colHead], split="[[:blank:]]+", fixed=FALSE, perl=TRUE))

        ## remove leading and trailing white space
        vars2 <- tolower(trimWS(vars1))

        ## make sure we recognize all columns in the DVH
        patDose <- "^dose"
        patVol  <- "^volume"
        patBin  <- "^bin"
        hits <- sum(c(grepl(patDose, vars2),
                      grepl(patVol,  vars2),
                      grepl(patBin,  vars2)))
        if(length(vars2) != hits) {
            stop(c("Could not identify all DVH columns"),
        		 paste(vars2, collapse=", "))
        }

        ## final column headers
        ## rename volume to relative volume if necessary
        vars3 <- vars2
        if(info$volumeUnit == "PERCENT") {
            vars3[grepl(patVol, vars3)] <- "volumeRel"
        }

        ## extract DVH as a matrix and store preceding information
        ## read line length(strct) for cases where file does not end with a
        ## blank line -> this will then be last DVH line, otherwise blank
        ## check if dvh is all blank -> no data
        if(all(nchar(strct[dvhStart:length(strct)]) == 0L)) {
            return(NULL)
        }
        
        con <- textConnection(strct[dvhStart:length(strct)])
        dvh <- data.matrix(read.table(con,
                                      header=FALSE, stringsAsFactors=FALSE,
                                      colClasses=rep("numeric", length(vars3)),
                                      comment.char="", nrows=dvhLen))
        close(con)

        ## remove bin numbers
        colnames(dvh) <- vars3
        dvh <- dvh[ , -grep(patBin, vars3, ignore.case=TRUE)]

        ## structVol is NA in Masterplan -> 100% or maximum of DVH volumes
        structVol <- if(info$volumeUnit == "PERCENT") {
            warning(c("No information on absolute structure volume available"))
            NA_real_
        } else {
            if(info$DVHtype == "cumulative") {
                max(dvh[ , "volume"])
            } else {
                ## reconstruct volumes -> volume is per gray -> mult with bin width
                volBin <- dvh[ , "volume"]*diff(c(-dvh[1, "dose"], dvh[ , "dose"]))
                sum(volBin)
            }
        }

        ## add information we don't have yet: relative/absolute volume
        if((       "volumeRel" %in% vars3) && !("volume"    %in% vars3)) {
            dvh <- cbind(dvh, volume=structVol*(dvh[ , "volumeRel"]/100))
        } else if(("volume"    %in% vars3) && !("volumeRel" %in% vars3)) {
            dvh <- cbind(dvh, volumeRel=100*(dvh[ , "volume"]/structVol))
        }

        ## add information we don't have yet: relative/absolute dose
        ## considering isoDoseRx
        if((    "doseRel" %in% vars3) && !("dose"    %in% vars3)) {
            dvh <- cbind(dvh, dose=dvh[ , "doseRel"]*doseRx / isoDoseRx)
            # (doseRx/(isoDoseRx/100))*(dvh$doseRel/100)
        } else if(("dose" %in% vars3) && !("doseRel" %in% vars3)) {
            dvh <- cbind(dvh, doseRel=dvh[ , "dose"]*isoDoseRx / doseRx)
            # 100*(dvh$dose/(doseRx/(isoDoseRx/100)))
        }

        DVH <- list(dvh=dvh,
                    patName=info$patName,
                    patID=info$patID,
                    date=info$date,
                    DVHtype=info$DVHtype,
                    plan=info$plan,
                    quadrant=info$quadrant,
                    structure=structure,
                    structVol=structVol,
                    doseUnit=info$doseUnit,
                    volumeUnit=info$volumeUnit,
                    doseRx=doseRx,
                    isoDoseRx=isoDoseRx)

        ## convert differential DVH to cumulative
        ## and add differential DVH separately
        if(info$DVHtype == "differential") {
            DVH$dvh     <- dvhConvert(dvh, toType="cumulative", toDoseUnit="asis")
            DVH$dvhDiff <- dvh
            dvhDiff     <- dvh
        } else {
            ## need differential DVH for mean dose etc.
            dvhDiff <- dvhConvert(dvh, toType="differential", toDoseUnit="asis")
        }

        ## Masterplan does not export mean/min/max ...
        ## dose category mid-points
        N       <- nrow(dvh)
        doseMid <- dvh[-N, "dose"] + diff(dvh[ , "dose"])/2
        volume  <- if(!all(is.na(dvh[ , "volume"]))) {
            dvh[ , "volume"]
        } else {
            dvh[ , "volumeRel"]
        }

        DVH$doseMin <- min(dvh[-N, "dose"][-diff(volume) > 0])
        DVH$doseMax <- max(doseMid[             -diff(volume) > 0])
        DVH$doseMed <- approx(dvh[ , "volumeRel"], dvh[ , "dose"], 50, method="linear", rule=1)$y
        DVH$doseAvg <- sum(doseMid*(-1)*diff(dvh[, "volumeRel"]/100))
        DVH$doseSD  <- sqrt(sum(doseMid^2*(-1)*diff(dvh[, "volumeRel"]/100)) - DVH$doseAvg^2)
        
        ## for mode, abs or rel volume does not matter
        DVH$doseMode <- if(!all(is.na(dvhDiff[ , "volume"]))) {
            dvhDiff[which.max(dvhDiff[ , "volume"]),    "dose"]
        } else {
            dvhDiff[which.max(dvhDiff[ , "volumeRel"]), "dose"]
        }

        ## set class
        class(DVH) <- "DVHs"
        return(DVH)
    }

    ## list of DVH data frames with component name = structure
    info <- list(patID=patID, patName=patName, date=DVHdate,
                 DVHtype=DVHtype, plan=plan, quadrant=quadrant,
                 doseRx=doseRx, isoDoseRx=isoDoseRx, doseUnit=doseUnit,
                 volumeUnit=volumeUnit)
    dvhL <- lapply(structList, getDVH, info=info)
    dvhL <- Filter(function(y) !is.null(y), dvhL)
    names(dvhL) <- sapply(dvhL, function(y) y$structure)
    if(length(unique(names(dvhL))) < length(dvhL)) {
        warning("Some structures have the same name - this can lead to problems")
    }

    class(dvhL) <- "DVHLst"
    attr(dvhL, which="byPat") <- TRUE

    return(dvhL)
}
