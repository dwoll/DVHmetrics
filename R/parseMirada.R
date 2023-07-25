#####---------------------------------------------------------------------------
## parse character vector from Mirada file
## courseAsID ignored
parseMirada <- function(x, planInfo=FALSE, courseAsID=FALSE, ...) {
    planInfo <- as.character(planInfo)
    dots     <- list(...)
    
    ## extract file header and header info
    ## Mirada file does not have a PatID
    header1  <- strsplit(x[1], "\t")[[1]]
    patName  <- header1[2]
    DVHdate  <- tryCatch(as.Date(strptime(header1[3], "%m/%d/%Y")),
                         error=function(e) { NA_character_ })
    
    patID      <- gsub("[^a-z0-9]", "\\1", tempfile(pattern="", tmpdir="")) # as.character(trunc(runif(1, 0, 10000001)))
    doseUnit   <- NA_character_
    volumeUnit <- NA_character_
    doseRx     <- NA_real_

    ## extract actual DVH columns
    idx_sep   <- which(!nzchar(trimws(x)))[1]
    DVHspan   <- x[(idx_sep+1):(length(x)-1)]
    ## construct proper DVH names
    DVHnames0 <- strsplit(DVHspan[1], "\t")[[1]]
    dose_cols <- seq(1L, length(DVHnames0)-1L, by=2L)
    vol_cols  <- dose_cols+1L
    isVolRel  <- all(grepl("^% of volume$", DVHnames0[vol_cols]))
    struct_names0 <- DVHnames0[dose_cols]
    struct_names  <- trimws(gsub("^([[:alnum:]_-]+)[[:blank:]]+\\(.+\\)$", "\\1", struct_names0))
    struct_vols   <- as.numeric(gsub("^[[:alnum:]_-]+[[:blank:]]+\\(([[:digit:].]+).+\\)$", "\\1", struct_names0))
    volumeUnit0   <- gsub("^[[:alnum:]_-]+[[:blank:]]+\\([[:digit:].]+(.+)\\)$", "\\1", struct_names0)
    volumeUnit1   <- unique(substr(volumeUnit0, 1L, 2L))
    
    stopifnot(length(volumeUnit1) == 1L)
    volumeUnit <- if(volumeUnit1 == "cm") {
        "CC"
    } else {
        stop(paste0("Volume unit cm3 expected, but found ", volumeUnit1))
    }

    DVHall <- read.table(text=DVHspan[-1], header=FALSE, sep="\t",
                         stringsAsFactors=FALSE, comment.char="")
    
    ## DVH structures are given in wide format -> split
    getDVH <- function(idx_dose, idx_vol, structName, structVol, info) {
        dvh <- data.matrix(DVHall[ , idx_dose:idx_vol])
        ## set names
        colnames(dvh) <- if(info$isVolRel) {
            c("dose", "volumeRel")
        } else {
            c("dose", "volume")
        }
        
        ## check if dose is increasing
        stopifnot(isIncreasing(dvh))
        
        ## add absolute/relative volume
        dvh <- if(info$isVolRel) {
            cbind(dvh,
                  volume=(1/100)*dvh[ , "volumeRel"]*structVol,
                  doseRel=NA_real_)
        } else {
            ## check if structure volume should be assumed
            ## to be equal to max given volume in DVH
            structVol <- if(hasName(dots, "volume_from_dvh")) {
                if(dots[["volume_from_dvh"]]) {
                    max(dvh[ , "volume"])
                }
            } else {
                NA_real_
            }
            
            cbind(dvh,
                  volumeRel=100*(dvh[ , "volume"] / structVol),
                  doseRel=NA_real_)
        }
        
        ## differential or cumulative DVH
        DVHtype <- dvhType(dvh)
        
        DVH <- list(dvh=dvh,
                    patName=info$patName,
                    patID=info$patID,
                    date=info$date,
                    DVHtype=DVHtype,
                    structure=structName,
                    structVol=structVol,
                    doseUnit=info$doseUnit,
                    volumeUnit=info$volumeUnit,
                    doseRx=info$doseRx,
                    doseMin=NA_real_,
                    doseMax=NA_real_,
                    doseAvg=NA_real_,
                    doseMed=NA_real_,
                    doseMode=NA_real_,
                    doseSD=NA_real_)
        
        ## convert differential DVH (not per unit dose!) to cumulative
        ## and add differential DVH separately
        if(DVHtype == "differential") {
            DVH$dvh <- convertDVH(dvh, toType="cumulative",
                                  toDoseUnit="asis", perDose=FALSE)
            DVH$dvhDiff <- dvh
        }
        
        ## set class
        class(DVH) <- "DVHs"
        return(DVH)
    }

    ## list of DVH data frames with component name = structure
    info <- list(patID=patID,
                 patName=patName,
                 date=DVHdate,
                 doseRx=doseRx,
                 doseUnit=doseUnit,
                 volumeUnit=volumeUnit,
                 isVolRel=isVolRel)
    
    dvhL <- Map(getDVH,
                dose_cols, vol_cols,
                struct_names,
                struct_vols,
                list(info))

    dvhL <- Filter(Negate(is.null), dvhL)
    names(dvhL) <- vapply(dvhL, function(y) y$structure, character(1))
    if(length(unique(names(dvhL))) < length(dvhL)) {
        warning("Some structures have the same name - this can lead to problems")
    }

    class(dvhL) <- "DVHLst"
    attr(dvhL, which="byPat") <- TRUE

    return(dvhL)
}
