get_name_pair <- function(x, y, sep=" <-> ") {
    paste0(x, sep, y) 
}

get_name_elem <- function(x, which=1L, sep=" <-> ") {
    stopifnot(which %in% c(1L, 2L))
    stopifnot(grepl(sep, x))
    if(which == 1) {
        gsub(paste0("^(.+)", sep, ".+$"), "\\1", x)
    } else {
        gsub(paste0("^.+", sep, "(.+$)"), "\\1", x)
    }
}

read_mesh_one <- function(x, name,
                          reconstruct=c("No", "AFS", "Poisson"),
                          spacing=1) {
    reconstruct <- match.arg(tolower(reconstruct),
                             choices=c("no", "afs", "poisson"))
    
    mesh_name <- if(missing(name)) {
        basename(tools::file_path_sans_ext(x))
    } else {
        basename(tools::file_path_sans_ext(name))
    }
    
    mesh_0 <- readMeshFile(x)
    mesh_r <- if(reconstruct == "no") {
        toRGL(Mesh(vertices=mesh_0$vertices, faces=mesh_0$faces))
    } else if(reconstruct == "afs") {
        AFSreconstruction(mesh_0$vertices)
    } else if(reconstruct == "poisson") {
        PoissonReconstruction(mesh_0$vertices, spacing=spacing)
    }
    
    vol_0 <- try(meshVolume(mesh_r))
    ctr_0 <- try(meshCentroid(mesh_r))
    vol <- if(!inherits(vol_0, "try-error")) {
        vol_0
    } else {
        NA_real_
    }
    
    ctr <- if(!inherits(ctr_0, "try-error")) {
        ctr_0
    } else {
        rep(NA_real_, 3)
    }
    
    list(name    =mesh_name,
         mesh    =mesh_r,
         volume  =vol,
         centroid=ctr)
}

read_mesh_obs <- function(x, name,
                          reconstruct=c("No", "AFS", "Poisson"),
                          spacing=1) {
    reconstruct <- match.arg(tolower(reconstruct),
                             choices=c("no", "afs", "poisson"))
    
    mesh_names <- if(missing(name)) {
        basename(tools::file_path_sans_ext(x))
    } else {
        basename(tools::file_path_sans_ext(name))
    }
    
    meshL <- lapply(seq_along(x), function(i) {
        read_mesh_one(x[i],
                      name=mesh_names[i],
                      reconstruct=reconstruct,
                      spacing=spacing)
    })
    
    setNames(meshL, mesh_names)
}

read_mesh <- function(x, name,
                      reconstruct=c("No", "AFS", "Poisson"),
                      spacing=1) {
    reconstruct <- match.arg(tolower(reconstruct),
                             choices=c("no", "afs", "poisson"))
    
    obs_names <- if(missing(name)) {
        names(x)
    } else {
        names(name)
    }
    
    mesh_names <- if(missing(name)) {
        lapply(x, function(fv) { basename(tools::file_path_sans_ext(fv)) })
    } else {
        name
    }
    
    Map(read_mesh_obs,
        setNames(x, obs_names),
        mesh_names,
        reconstruct=reconstruct,
        spacing=spacing)
}

print_mesh_one <- function(x) {
    vol_fmt_str <- if(is.na(x$volume))        { "%s" }           else { "%.2f" }
    ctr_fmt_str <- if(any(is.na(x$centroid))) { "[%s, %s, %s]" } else { "[%.2f, %.2f, %.2f]" }
    cat(sprintf("Mesh: %s", x$name),
        "\n",
        capture.output(print(x$mesh)),
        "\n",
        sprintf(paste0("Volume: ", vol_fmt_str), x$volume),
        "\n",
        sprintf(paste0("Centroid: ", ctr_fmt_str),
                x$centroid[1],
                x$centroid[2],
                x$centroid[3]),
        "\n",
        "\n", sep="")
}

print_mesh <- function(x) {
    invisible(Map(print_mesh_one, unlist(x, recursive=FALSE)))
}

## starting from list of observers, each with a list of meshes
## generate all observer-pairs for each corresponding mesh-list entry
get_mesh_pairs <- function(x, sep=" <-> ", names_only=FALSE) {
    if(length(x) <= 1L) { stop("Need more than 1 mesh for comparisons") }
    
    ## number of meshes per observer
    n_obs_meshes <- lengths(x) # may be different per observer
    n_meshes     <- max(n_obs_meshes)
    
    ## for given pair, put corresponding meshes in a list
    get_pair_mesh <- function(pair, mesh) {
        idx1  <- pairs_idx[pair, 1] # index observer 1
        idx2  <- pairs_idx[pair, 2] # index observer 2
        obs1  <- x[[idx1]]          # observer 1
        obs2  <- x[[idx2]]          # observer 2
        ## do both observers have the mesh?
        if((length(obs1) >= mesh) && (length(obs2) >= mesh)) {
            mesh1 <- x[[idx1]][[mesh]]
            mesh2 <- x[[idx2]][[mesh]]
            
            ## add group information -> same structure
            if(names_only) {
                list(name=get_name_pair(mesh1$name, mesh2$name, sep=sep),
                     group=sprintf("strct_%.3d", mesh))
            } else {
                list(name=get_name_pair(mesh1$name, mesh2$name, sep=sep),
                     mesh1=mesh1,
                     mesh2=mesh2,
                     group=sprintf("strct_%.3d", mesh))
            }
        } else {
            NULL
        }
    }
    
    pairs_idx <- if(length(x) >= 2L) {
        t(combn(seq_along(x), 2))
    } else {
        matrix(c(1, 1), ncol=2)
    }
    
    ll_outer <- lapply(seq_len(n_meshes), function(idx_mesh) {
        lapply(seq_len(nrow(pairs_idx)), function(idx_pair) { get_pair_mesh(idx_pair, idx_mesh) })
    })
    
    ## weed out NULL components
    ll <- Filter(Negate(is.null), unlist(ll_outer, recursive=FALSE))
    pair_names <- lapply(ll, function(x) { x$name })
    setNames(ll, pair_names)
}

## union and intersection for list of two meshes x
get_mesh_ui_pair <- function(x) {
    union_0     <- try(MeshesUnion(       list(x$mesh1$mesh, x$mesh2$mesh), clean=TRUE))
    intersect_0 <- try(MeshesIntersection(list(x$mesh1$mesh, x$mesh2$mesh), clean=TRUE))
    
    if(!inherits(union_0, "try-error") && !inherits(intersect_0, "try-error")) {
        union     <- toRGL(union_0)
        intersect <- toRGL(intersect_0)
    } else {
        union     <- NULL
        intersect <- NULL
    }
    
    list(name=get_name_pair(x$mesh1$name, x$mesh2$name),
         group=x$group,
         union=union,
         intersection=intersect)
}

get_mesh_ui <- function(x) {
    pairL <- get_mesh_pairs(x)
    Map(get_mesh_ui_pair, pairL)
}

get_mesh_metro_pair <- function(x, chop=TRUE, ...) {
    metro <- vcgMetro(x$mesh1$mesh,
                      x$mesh2$mesh,
                      ...)
    
    if(chop) {
        metro[["distances1"]]    <- NULL
        metro[["distances2"]]    <- NULL
        metro[["forward_hist"]]  <- NULL
        metro[["backward_hist"]] <- NULL
    }
    
    metro[["name"]]  <- get_name_pair(x$mesh1$name, x$mesh2$name)
    metro[["group"]] <- x$group
    metro
}

get_mesh_metro <- function(x, chop=TRUE, ...) {
    pairL <- get_mesh_pairs(x)
    Map(get_mesh_metro_pair, pairL, chop=chop, ...)
}

## distance measures, union, intersection for each mesh pair
get_mesh_agree_pair <- function(x, metro, ui, chop=TRUE, ...) {
    ## distance-based measures
    if(missing(metro)) {
        metro <- get_mesh_metro_pair(x, chop=chop, ...)
    }
    
    DCOM <- sqrt(sum((x$mesh2$centroid - x$mesh1$centroid)^2))
    HD_forward  <- metro$ForwardSampling$maxdist
    HD_backward <- metro$BackwardSampling$maxdist
    
    if(is.finite(HD_forward) && is.finite(HD_backward)) {
        HD_max <- max(c(HD_forward, HD_backward))
        HD_avg <- (HD_forward + HD_backward) / 2
    } else {
        HD_max <- NA_real_
        HD_avg <- NA_real_
    }
    
    ## average surface distance based on weighted average of sampled distances
    ## not on actual vertex distances as stored in distances1, distances2
    n1 <- metro$ForwardSampling$nsamples
    n2 <- metro$BackwardSampling$nsamples
    if((n1 > 0L) && (n2 > 0L)) {
        w1   <- n1 / (n1+n2)
        w2   <- n2 / (n1+n2)
        ASD  <-      w1* metro$ForwardSampling$meandist   + w2* metro$BackwardSampling$meandist
        RMSD <- sqrt(w1*(metro$ForwardSampling$RMSdist^2) + w2*(metro$BackwardSampling$RMSdist^2))
    } else {
        ASD  <- NA_real_
        RMSD <- NA_real_
    }
    
    ## volume-overlap-based measures
    ## check if union/intersection are supplied
    if(missing(ui)) {
        ui <- get_mesh_ui_pair(x)
    }
    
    if(!is.null(ui$union) && !is.null(ui$intersection)) {
        vol_mesh1     <- x$mesh1$volume
        vol_mesh2     <- x$mesh2$volume
        vol_union     <- meshVolume(ui$union)
        vol_intersect <- meshVolume(ui$intersection)
        
        JSC <- vol_intersect   / vol_union
        DSC <- 2*vol_intersect / (vol_mesh1 + vol_mesh2)
    } else {
        JSC <- NA_real_
        DSC <- NA_real_
    }
    
    data.frame(mesh1=x$mesh1$name,
               mesh2=x$mesh2$name,
               group=x$group,
               DCOM=DCOM,
               HD_max=HD_max,
               HD_avg=HD_avg,
               ASD=ASD,
               RMSD=RMSD,
               JSC=JSC,
               DSC=DSC)
}

get_mesh_agree <- function(x, chop=TRUE, ...) {
    pairL  <- get_mesh_pairs(x)
    uiL    <- Map(get_mesh_ui_pair,    pairL)
    metroL <- Map(get_mesh_metro_pair, pairL, chop=chop, ...)
    agreeL <- Map(get_mesh_agree_pair,
                  pairL,
                  metro=metroL,
                  ui=uiL,
                  chop=chop)
    
    d <- do.call("rbind", agreeL)
    rownames(d) <- NULL
    d
}

get_mesh_agree_long <- function(x) {
    vars_varying <- c("DCOM",
                      "HD_max", "HD_avg", "ASD", "RMSD",
                      "JSC", "DSC")
    
    vars_id <- names(x)[!(names(x) %in% vars_varying)]
    
    dL <- reshape(x,
                  direction="long",
                  idvar=vars_id,
                  varying=vars_varying,
                  v.names="observed",
                  timevar="metric")
    
    rownames(dL) <- NULL
    
    dL[["metric"]] <- factor(dL[["metric"]],
                             levels=seq_along(vars_varying),
                             labels=vars_varying)
    
    dL
}

get_mesh_agree_aggr <- function(x, na.rm=FALSE) {
    d_agreeL <- get_mesh_agree_long(x)
    d_agreeL[["observed_ln"]] <- log(d_agreeL[["observed"]])
    
    d_mean   <- aggregate(observed    ~ group + metric, FUN=mean,   data=d_agreeL, na.rm=na.rm)
    d_median <- aggregate(observed    ~ group + metric, FUN=median, data=d_agreeL, na.rm=na.rm)
    d_sd     <- aggregate(observed    ~ group + metric, FUN=sd,     data=d_agreeL, na.rm=na.rm)
    d_var    <- aggregate(observed    ~ group + metric, FUN=var,    data=d_agreeL, na.rm=na.rm)
    d_varlog <- aggregate(observed_ln ~ group + metric, FUN=var,    data=d_agreeL, na.rm=na.rm)
    
    d_aggr <- Reduce(function(x, y) { suppressWarnings(merge(x, y, by=c("group", "metric"))) },
                     list(d_mean, d_median, d_sd, d_var, d_varlog))
    
    names(d_aggr)       <- c("group", "metric", "Mean", "Median", "SD", "VAR", "VAR_log")
    d_aggr[["CV"]]      <- d_aggr[["SD"]] / d_aggr[["Mean"]]
    d_aggr[["CV_ln"]]   <- sqrt(exp(d_aggr[["VAR_log"]]) - 1)
    d_aggr[["VAR"]]     <- NULL
    d_aggr[["VAR_log"]] <- NULL
    
    d_aggr
}

get_mesh_agree_aggr_long <- function(x) {
    vars_varying <- c("Mean", "Median", "SD", "CV", "CV_ln")
    vars_id      <- names(x)[!(names(x) %in% vars_varying)]
    
    dL <- reshape(x,
                  direction="long",
                  idvar=vars_id,
                  varying=vars_varying,
                  v.names="observed",
                  timevar="statistic")
    
    rownames(dL) <- NULL
    
    dL[["statistic"]] <- factor(dL[["statistic"]],
                                levels=seq_along(vars_varying),
                                labels=vars_varying)
    
    dL
}

mesh_list_to_observer_list <- function(x) {
    ll <- Map(function(i, name) {
        setNames(list(i), name)
    }, x, names(x))
    
    setNames(ll, sprintf("Observer_%.2d", seq_along(ll)))
}
