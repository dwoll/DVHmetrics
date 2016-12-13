#####---------------------------------------------------------------------------
## parse character vector from Tomo HiArt DVH file
parsePRIMO <- function(x, planInfo=FALSE, courseAsID=FALSE) {
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

    getDVHtype <- function(ll) {
        line <- ll[grep("^# Mode:", ll)]
        elem <- sub("^# Mode:[[:blank:]]+(Cumulative|Differential)", "\\1",
                    line, perl=TRUE, ignore.case=TRUE)
        tolower(trimWS(elem))
    }

    sStart <- grep("^# _________________", x) # start of DVH matrix
    sLen   <- diff(c(sStart, length(x)+1))  # length of sections

    ## extract file header and header info
    header     <- x[seq_len(sStart[1]-1)]                        # header
    patName    <- gsub("[^a-z0-9]", "\\1", tempfile(pattern="", tmpdir=""))
    patID      <- getElem("^# Project:",  header)   # patient id
    plan       <- NA_character_
    doseRx     <- NA_real_
    isoDoseRx  <- NA_real_
    DVHdate    <- NA_character_
    DVHtype    <- getDVHtype(header)
    doseRx     <- NA_real_

    ## get dose and volume units
    varDose  <- sub("^.+Dose \\((.+?)\\)[[:blank:]]+.+$",
                    "\\1", header[grep("Dose \\(.+?\\)",   header)])
    varVol   <- sub("^.+Volume \\((.+?)\\)$",
                    "\\1", header[grep("Volume \\(.+?\\)", header)])
    if(grepl("%", varDose, ignore.case=TRUE)) {
        isDoseRel <- TRUE
        doseUnit  <- "PERCENT"
    } else {
        ## TODO: need example file for this
        warning("PRIMO files with absolute dose are not implemented")
        isDoseRel <- FALSE
        doseUnit  <- NA_character_
    }

    if(grepl("%", varVol, ignore.case=TRUE)) {
        isVolRel   <- TRUE
        volumeUnit <- "PERCENT"
    } else {
        warning("PRIMO files with absolute volume are not implemented")
        ## TODO: need example file for this
        isVolRel   <- FALSE
        volumeUnit <- NA_character_
    }

    ## find columns for structure, dose, volume
    varTxt <- x[sStart+1]
    varTxt <- gsub("^#[[:blank:]]+(.+)$", "\\1", varTxt)
    vars   <- as.matrix(read.table(text=varTxt, header=FALSE, sep="\t",
                                   stringsAsFactors=FALSE, comment.char="")[1, ])

    nStructs  <- length(vars) - 4
    structIdx <- 2:(nStructs+1)
    doseIdx   <- rep(1, nStructs)
    volumeIdx <- 3:(nStructs+2)

    ## read all data
    datAll <- data.matrix(read.table(text=x[(sStart[1]+3):(sStart[1]+sLen-1)],
                                     header=FALSE, dec=",", sep="\t",
                                     stringsAsFactors=FALSE, comment.char=""))

    ## extract DVH from one structure section and store in a list
    ## with DVH itself as a matrix
    getDVH <- function(strIdx, dIdx, vIdx, info) {
        structure <- vars[[strIdx]][1]

        ## extract DVH as a matrix and set variable names
        dvh <- datAll[ , c(dIdx, vIdx)]
        haveVars <- if(isVolRel) {
            if(isDoseRel) {
                c("doseRel", "volumeRel")
            } else {
                c("dose",    "volumeRel")
            }
        } else {
            if(isDoseRel) {
                c("doseRel", "volume")
            } else {
                c("dose",    "volume")
            }
        }
        
        colnames(dvh) <- haveVars

        ## add information we don't have yet
        ## relative/absolute volume/dose
        if(!("volume" %in% haveVars)) {
            dvh <- cbind(dvh, volume=NA_real_)
        }

        if(!("volumeRel" %in% haveVars)) {
            dvh <- cbind(dvh, volumeRel=NA_real_)
        }

        if(!("dose" %in% haveVars)) {
            dvh <- cbind(dvh, dose=NA_real_)
        }

        if(!("doseRel" %in% haveVars)) {
            dvh <- cbind(dvh, doseRel=NA_real_)
        }

        ## check if dose is increasing
        stopifnot(isIncreasing(dvh))

        ## differential or cumulative DVH
        DVHtype <- dvhType(dvh)

        DVH <- list(dvh=dvh,
                    patName=info$patName,
                    patID=info$patID,
                    date=info$date,
                    DVHtype=DVHtype,
                    plan=info$plan,
                    structure=structure,
                    structVol=NA_real_,
                    doseUnit=info$doseUnit,
                    volumeUnit=info$volumeUnit,
                    doseRx=doseRx,
                    isoDoseRx=isoDoseRx,
                    doseMin=NA_real_,
                    doseMax=NA_real_,
                    doseAvg=NA_real_,
                    doseMed=NA_real_,
                    doseMode=NA_real_,
                    doseSD=NA_real_)

        ## convert differential DVH to cumulative
        ## and add differential DVH separately
        if(DVHtype == "differential") {
            warning("I assume differential DVH is per unit dose\nbut I have no information on this")
            DVH$dvh <- convertDVH(dvh, toType="cumulative",
                                  toDoseUnit="asis", perDose=TRUE)
            DVH$dvhDiff <- dvh
        }

        ## set class
        class(DVH) <- "DVHs"
        return(DVH)
    }

    ## list of DVH data frames with component name = structure
    info <- list(patID=patID, patName=patName, date=DVHdate,
                 plan=plan, doseRx=doseRx, isoDoseRx=isoDoseRx,
                 doseUnit=doseUnit, volumeUnit=volumeUnit)
    dvhL <- Map(getDVH, structIdx, doseIdx, volumeIdx, info=list(info))
    dvhL <- Filter(Negate(is.null), dvhL)
    names(dvhL) <- sapply(dvhL, function(y) y$structure)
    if(length(unique(names(dvhL))) < length(dvhL)) {
        warning("Some structures have the same name - this can lead to problems")
    }

    class(dvhL) <- "DVHLst"
    attr(dvhL, which="byPat") <- TRUE

    return(dvhL)
}
