## we need ... because getMetric() will also pass parameters
## intended for other functions through ...
getTCP <-
function(x,
         TCPtcd50=NULL, TCPm=NULL, TCPn=NULL, TCPgamma50=NULL, NTCPs=NULL,
         EUDa=NULL, EUDfn=NULL, EUDab=NULL,
         TCPtype=c("probit", "logit", "poisson", "relative_seriality"), ...) {
    out <- getNTCP(x=x,
                   NTCPtd50=TCPtcd50,
                   NTCPm=TCPm,
                   NTCPn=TCPn,
                   NTCPgamma50=TCPgamma50,
                   NTCPs=NTCPs,
                   EUDa=EUDa,
                   EUDfn=EUDfn,
                   EUDab=EUDab,
                   NTCPtype=TCPtype)
    outNames <- names(out)
    outNames[outNames == "NTCP"] <- "TCP"
    setNames(out, outNames)
}

getNTCP <-
function(x,
         NTCPtd50=NULL, NTCPm=NULL, NTCPn=NULL, NTCPgamma50=NULL, NTCPs=NULL,
         EUDa=NULL, EUDfn=NULL, EUDab=NULL,
         NTCPtype=c("probit", "logit", "poisson", "relative_seriality"), ...) {
    UseMethod("getNTCP")
}

getNTCP.DVHs <-
function(x,
         NTCPtd50=NULL, NTCPm=NULL, NTCPn=NULL, NTCPgamma50=NULL, NTCPs=NULL,
         EUDa=NULL, EUDfn=NULL, EUDab=NULL,
         NTCPtype=c("probit", "logit", "poisson", "relative_seriality"), ...) {
    NTCPtype <- tolower(NTCPtype)
    NTCPtype <- match.arg(NTCPtype)
    stopifnot(!is.null(NTCPtd50))

    if(length(NTCPtd50) > 1L) {
    	warning(paste0("Will only use NTCPtd50=", NTCPtd50[1]))
    	NTCPtd50 <- NTCPtd50[1]
    }

    if(NTCPtd50 <= 0) {
        warning("'NTCPtd50' must be > 0")
    	return(NA_real_)
    }

    if(is.null(NTCPm)) {
        stopifnot(!is.null(NTCPgamma50))
        
        if(length(NTCPgamma50) > 1L) {
    		warning(paste0("Will only use NTCPgamma50=", NTCPgamma50[1]))
    		NTCPgamma50 <- NTCPgamma50[1]
        }

        if(NTCPgamma50 <= 0) {
            warning("'NTCPgamma50' must be > 0")
    		return(NA_real_)
        }

        NTCPm <- 1 / (NTCPgamma50*sqrt(2*pi))
    } else {
        if(length(NTCPm) > 1L) {
    		warning(paste0("Will only use NTCPm=", NTCPm[1]))
    		NTCPm <- NTCPm[1]
        }

        if(NTCPm <= 0) {
            warning("'NTCPm' must be > 0")
    		return(NA_real_)
        }

        NTCPgamma50 <- 1 / (NTCPm*sqrt(2*pi))
    }

    if(NTCPtype %in% c("probit", "logit", "poisson")) {
        if(!is.null(EUDa)) {
            if(length(EUDa) > 1L) {
                warning(paste0("Will only use EUDa=", EUDa[1]))
                EUDa <- EUDa[1]
            }
            
            if(isTRUE(all.equal(EUDa, 0))) {
                warning("'EUDa' must not be zero")
                return(NA_real_)
            }
            
            NTCPn <- 1/EUDa
        } else {
            stopifnot(!is.null(NTCPn))
            
            if(length(NTCPn) > 1L) {
                warning(paste0("Will only use NTCPn=", NTCPn[1]))
                NTCPn <- NTCPn[1]
            }
            
            if(is.infinite(NTCPn)) {
                warning("'NTCPn' must not be infinite")
                return(NA_real_)
            }
            
            EUDa <- 1/NTCPn
        }
        
        EUD <- getEUD(x, EUDa=EUDa, EUDfn=EUDfn, EUDab=EUDab)$EUD
    }
    
    NTCP <- if(NTCPtype == "probit") {
        ## Lyman probit model based on EUD
        ## quantile at which to evaluate standard normal cdf
        ## equivalent to
        ## 0.5*(1+pracma::erf((EUD-NTCPtd50)/(sqrt(2)*NTCPm*NTCPtd50)))
        qq <- (EUD-NTCPtd50) / (NTCPm*NTCPtd50)
        pnorm(qq, mean=0, sd=1, lower.tail=TRUE)
    } else if(NTCPtype == "logit") {
        ## Niemierko logit model based on EUD
        1 / (1 + ((NTCPtd50/EUD)^(4*NTCPgamma50)))
    } else if(NTCPtype == "poisson") {
        ## Kaellman Poisson model based on EUD
        2^(-exp(exp(1)*NTCPgamma50*(1-(EUD/NTCPtd50))))
    } else if(NTCPtype == "relative_seriality") {
        ## relative seriality model
        stopifnot(!is.null(NTCPs))
        if(length(NTCPs) > 1L) {
            warning(paste0("Will only use NTCPs=", NTCPs[1]))
            NTCPs <- NTCPs[1]
        }
        
        if(NTCPs <= 0) {
            warning("'NTCPs' must be > 0")
            return(NA_real_)
        }

        ## get differential DVH
        x_diff <- convertDVH(x, toType="differential", perDose=FALSE)
        Di     <- x_diff$dvhDiff[ , "dose"]
        Vi_rel <- x_diff$dvhDiff[ , "volumeRel"] / 100 # volumeRel is in %
        ## per-voxel NTCP
        # ntcp_voxels <- exp(-exp(exp(1) * NTCPgamma50 - ((Di / NTCPtd50) * (exp(1) * NTCPgamma50 - log(log(2.0))))))
        # ntcps       <- (1.0 - ntcp_voxels^NTCPs)^Vi_rel
        # (1.0 - prod(ntcps))^(1.0 / NTCPs)
        P_Di <- 2^(-exp(exp(1)*NTCPgamma50*(1 - (Di/NTCPtd50))))
        (1 - prod((1 - P_Di^NTCPs)^Vi_rel))^(1/NTCPs)
    }

    data.frame(NTCP=NTCP,
               patID=x$patID,
               structure=x$structure,
               stringsAsFactors=FALSE)
}

getNTCP.DVHLst <-
function(x,
         NTCPtd50=NULL, NTCPm=NULL, NTCPn=NULL, NTCPgamma50=NULL, NTCPs=NULL,
         EUDa=NULL, EUDfn=NULL, EUDab=NULL,
         NTCPtype=c("probit", "logit", "poisson", "relative_seriality"), ...) {
    NTCPl <- Map(getNTCP, x,
                 NTCPtd50=list(NTCPtd50),
                 NTCPm=list(NTCPm),
                 NTCPn=list(NTCPn),
                 NTCPgamma50=list(NTCPgamma50),
                 NTCPs=list(NTCPs),
                 EUDa=list(EUDa),
                 EUDfn=list(EUDfn),
                 EUDab=list(EUDab),
                 NTCPtype=list(NTCPtype))
    NTCPdf <- do.call("rbind", NTCPl)
    rownames(NTCPdf) <- NULL
    NTCPdf
}

getNTCP.DVHLstLst <-
function(x,
         NTCPtd50=NULL, NTCPm=NULL, NTCPn=NULL, NTCPgamma50=NULL, NTCPs=NULL,
         EUDa=NULL, EUDfn=NULL, EUDab=NULL,
         NTCPtype=c("probit", "logit", "poisson", "relative_seriality"), ...) {
    NTCPl <- Map(getNTCP, x,
                 NTCPtd50=list(NTCPtd50),
                 NTCPm=list(NTCPm),
                 NTCPn=list(NTCPn),
                 NTCPgamma50=list(NTCPgamma50),
                 NTCPs=list(NTCPs),
                 EUDa=list(EUDa),
                 EUDfn=list(EUDfn),
                 EUDab=list(EUDab),
                 NTCPtype=list(NTCPtype))
    NTCPdf <- do.call("rbind", NTCPl)
    rownames(NTCPdf) <- NULL
    NTCPdf
}
