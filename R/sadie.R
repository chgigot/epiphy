#------------------------------------------------------------------------------#
#' Spatial Analysis by Distance IndicEs (SADIE)
#'
#' The SADIE approach is implemented in this package. Refer to the documentation
#' of the function \code{\link{sadie}} for more information and references about
#' this approach.
#' The \code{sadie} function estimates the parameters of a SADIE approach based
#' on the number of permutations specified (\code{nperm}). It also computes the
#' Ia indice and Pa probability.
#'
#' @param data A data frame or a matrix with the counts for the different sampling units. The data
#'     set must decribe a complete (no missing data) squared or rectangular
#'     data set.
#'     or a matrix (x-y).
#' @param cost A \code{length(count)} by \code{length(count)} matrix specifing
#'     the cost of transporting single unit between the corresponding source and
#'     destination sampling unit. If not specified, the distance between each
#'     couple of adjacent sampling units (rook's sense) is assumed to be equal
#'     to 1. If the cost matrix is not given, the two first columns of data must
#'     correspond to "x" and "y" coordinates, and an euclidean distance between
#'     all the recorded points will be computed.
#' @param index The name of the index to use: "Perry", "Li-Madden-Xu" or "all".
#'     By default, only Perry's index is computed.
#' @param nperm Number of permutations for each sampling unit.
#' @param rseed Randomisation seed. Unseful for checking pursposes.
#'
#' For historical reasons, the inflow is negative and the outflow is positive.... not true parce ce aue on parle de cost
#' Two transportation algorithms are available: The Shortlist Method (Gottschlich C. and Schuhmacher D.)
#' and the revised simplex algorithm (Luenberger and Ye (2008, Section 6.4)). The first one is substantially faster than the second one.
#' Both of these algorithms are available in the R package \code{transportation}. There were integrated and
#' fully C++-coded in this package in order to minimize R-code for fast reasons.
#'
#' @references
#'
#' Perry JN. 1995. Spatial analysis by distance indices. Journal of Animal
#' Ecology 64, 303–314. \href{http://dx.doi.org/10.2307/5892}{doi:10.2307/5892}
#'
#' Perry JN, Winder L, Holland JM, Alston RD. 1999. Red–blue plots for detecting
#' clusters in count data. Ecology Letters 2, 106–113.
#' \href{http://dx.doi.org/10.1046/j.1461-0248.1999.22057.x}{doi:10.1046/j.1461-0248.1999.22057.x}
#'
#' Li B, Madden LV, Xu X. 2012. Spatial analysis by distance indices: an
#' alternative local clustering index for studying spatial patterns. Methods in
#' Ecology and Evolution 3, 368–377.
#' \href{http://dx.doi.org/10.1111/j.2041-210X.2011.00165.x}{doi:10.1111/j.2041-210X.2011.00165.x}
#'
#' @examples
#' # Possible datasets:
#' #aphids
#' #arthropods
#' #codling_moths
#' res <- sadie(aphids)
#' summary(res)
#' plot(res)
#'
#' @name sadie
#' @export
#------------------------------------------------------------------------------#
sadie <- function(data, ...) {
    UseMethod("sadie")
}

#------------------------------------------------------------------------------#
#' @rdname sadie
#' @method sadie data.frame
#' @export
#------------------------------------------------------------------------------#
sadie.data.frame <- function(data, index = c("Perry", "Li-Madden-Xu", "all"),
                             nperm = 100, rseed = TRUE, seed = 12345, cost,
                             threads = 1, ..., method = "shortsimplex") {

    index <- match.arg(index)
    n_col <- ncol(data)
    stopifnot(n_col == 3)

    if (missing(cost)) {
        # Assumption: the 1st and 2nd columns of the input data frame correspond
        # to "x" and "y" coordinates, respectively.
        cost <- as.matrix(dist(data[1:2])) ### ATTENTION pour plus compliqué données eg arthropodes
    }
    stopifnot(is.matrix(cost))
    # More checking and data wrangling...

    # More interesting
    if (!rseed) set.seed(seed)

    N <- nrow(data)
    data[[n_col]] <- as.numeric(data[[n_col]]) # Just in case it's integer, to have the same nature for data and fdata

    # Reshape long data into a wide data frame (with standrad R)
    ##cn   <- colnames(data)
    ##n_cn <- length(cn)
    ##frm  <- as.formula(paste0(cn[n_cn], "~",
    ##                          paste0(cn[1:(n_cn - 1)], collapse = "+")))
    ##data_array <- unclass(xtabs(frm, data)) # unclass, to get a "naked" matrix/array
    ##attr(data_array, "call") <- NULL # Just to clean the object

    # Create a homogeneous (flatten) version of data_array

    ##fdata_array <- array(data = rep(mean(data_array), N),
    ##                     dim = dim(data_array),
    ##                     dimnames = dimnames(data_array))

    fcount <- rep(mean(data[[n_col]]), N)

    # Use optimal transportation algorithm
    opt_transport <- wrap_transport(data[[n_col]], fcount, cost, method)
    opt_transport <- as.matrix(opt_transport, N)

    # Computation the costs
    cost_flows <- vapply(1:N, function(x) {
        costToti(x, opt_transport, cost, type = "both")
    }, numeric(1L))
    # Due to the convention "in => neg value", only negative values are added
    # up. Da = total observed distance (or total cost).
    Da <- sum(cost_flows[cost_flows < 0])

    # Deal with indices
    if (index == "all") index <- c("Perry", "Li-Madden-Xu")
    info_P   <- list(clust = NA_real_, Ea = NA_real_)
    info_LMX <- list(clust = NA_real_, Dis_all = matrix(NA_real_))
    new_prob <- NA_real_
    if (any(index == "Perry")) {
        cat("Computation of Perry's indices:\n")
        info_P <- clust_P(opt_transport, cost, N, nperm,
                          start = data[[n_col]], end = fcount,
                          method, cost_flows, threads)
    }
    if (any(index == "Li-Madden-Xu")) {
        cat("Computation of Li-Madden-Xu's indices:\n")
        info_LMX <- clust_LMX(opt_transport, cost, N, nperm,
                              start = data[[n_col]], end = fcount,
                              method, cost_flows, threads)
        # Test:
        new_prob <- abs(info_LMX[["Dis_all"]]) > abs(cost_flows)
        new_prob <- rowSums(new_prob)
        new_prob <- (new_prob + 1) / (nperm + 1) # res2 + 1 => rank of Di(s)
    }

    # Compute other outputs (if possible)
    # Ci-dessous, ça veut dire que l'on calcule toujours le perry's
    if (any(is.na(info_P[["Ea"]]))) {
        warning("Ia and Pa cannot be computed if Perry's indices are not.")
        Ia <- NA_real_
        Pa <- NA_real_
    } else {
        Ia <- Da / mean(info_P[["Ea"]])
        # On met des moins pour mémoire  à cause de cette $%@*&^ de convention inversée !!!
        Pa <- sum((-info_P[["Ea"]]) > (-Da)) / length(info_P[["Ea"]])
    }
    info_clust <- data.frame(data,
                             cost_flows = cost_flows,
                             idx_P      = info_P[["clust"]],
                             idx_LMX    = info_LMX[["clust"]],
                             prob       = new_prob)

    # Summary indices
    summary_idx <- data.frame(overall = c(mean(abs(info_clust$idx_P)),
                                          mean(abs(info_clust$idx_LMX))),
                              inflow  = c(mean(info_clust$idx_P[info_clust$idx_P < 0]),
                                          mean(info_clust$idx_LMX[info_clust$idx_LMX < 0])),
                              outflow = c(mean(info_clust$idx_P[info_clust$idx_P > 0]),
                                          mean(info_clust$idx_LMX[info_clust$idx_LMX > 0])))
    rownames(summary_idx) <- c("Perry's index", "Li-Madden-Xu's index")

    #    return(Sadie())
    res <- list(info_clust = info_clust,
                Da = Da,
                Ia = Ia,
                Pa = Pa,
                Ea_perry = info_P[["Ea"]],
                Ea_li = info_LMX[["Dis_all"]],
                summary_idx = summary_idx,
                nperm = nperm,
                rseed = rseed,
                seed = NA_real_)
    attr(res, "class") <- "sadie"
    attr(res, "call")  <- match.call()
    res
}

#------------------------------------------------------------------------------#
#' @rdname sadie
#' @export
#------------------------------------------------------------------------------#
sadie.matrix <- function(data, index = c("Perry", "Li-Madden-Xu", "all"),
                         nperm = 100, rseed = TRUE, seed = 12345, cost,
                         threads = 1) {
    sadie.data.frame(as.data.frame(data), index, nperm, rseed, seed, cost, threads)
}

#------------------------------------------------------------------------------#
#' @rdname sadie
#' @export
#------------------------------------------------------------------------------#
sadie.count <- function(data, index = c("Perry", "Li-Madden-Xu", "all"),
                        nperm = 100, rseed = TRUE, seed = 12345, cost,
                        threads = 1) {
    mapped_data <- map_data(data)
    sadie.data.frame(mapped_data, index, nperm, rseed, seed, cost, threads)
}


#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
print.sadie <- function(x, ...) {
    cat("\nCall:\n")
    print(attr(x, "call"))
    cat("\nIa: ", format(x$Ia, digits = 1, nsmall = 4),
        " (Pa = ", format.pval(x$Pa), ")\n\n", sep = "")
}

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
summary.sadie <- function(object, ...) {
    res <- object
    class(res) <- "summary.sadie"
    res
}

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
print.summary.sadie <- function(x, ...) {
    cat("\nCall:\n")
    print(attr(x, "call"))
    n <- 6L
    cat("\nFirst ", n, " rows of clustering indices:\n", sep = "") # Cf Image pkg for a nice display
    print(head(x$info_clust, n = n))
    #printCoefmat(x$info_clust)
    #cat("number of permutations: ", object$nperm, "\n",
    #    "random seed: ", object$rseed, "\n", sep = "")
    cat("\nSummary indices:\n")
    print(x$summary_idx)
    cat("\nMain outputs:")
    cat("\nIa: ", format(x$Ia, digits = 1, nsmall = 4),
        " (Pa = ", format.pval(x$Pa), ")\n", sep = "")
    cat("\n'Total cost': ", x$Da, "\n",
        "Number of permutations: ", x$nperm, "\n",
        "Random Seed?: ", x$rseed, "(", x$seed, ")\n\n", sep = "")
}

#------------------------------------------------------------------------------#
# Utilities
#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
wrap_transport <- function(start, end, cost, method = "shortsimplex",
                           control = list(), ...) {
    # As we always have named rows and columns for end, we use this one
    # for the generator. But in theory, start should be used too.

    ##gen <- lapply(dimnames(end), as.numeric)
    ##res <- transport::transport(transport::pgrid(start, generator = gen),
    ##                            transport::pgrid(end, generator = gen),
    ##                            p = p, method = method, control = control, ...)
    ##class(res) <- c("transport", class(res))
    ##res

    res <- transport::transport(start, end, costm = cost, method = method,
                                control = control, ...)
    res <- res[res$from != res$to, ] # We do not take into account what happens on the diagonal
    class(res) <- c("transport", class(res))
    res
}

#------------------------------------------------------------------------------#
#' @method as.matrix transport
#' @export
#------------------------------------------------------------------------------#
as.matrix.transport <- function(x, dim_mat, dimnames = NULL, ...) { # N: Number of sampling units
    res <- matrix(rep(0, dim_mat^2), nrow = dim_mat, ncol = dim_mat,
                  dimnames = dimnames, ...)
    lapply(1:nrow(x), function(i1) {
        res[x[i1, 1], x[i1, 2]] <<- x[i1, 3]
    })
    res
}

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
costToti <- function(i, flow, cost, type = "both", average = FALSE) {
    stopifnot(type == "in" | type == "out" | type == "both")
    #    sorties             entrées
    if ((sum(flow[i, ]) <= sum(flow[, i])) &
        ((type == "in") | (type == "both"))) {
        res <- -sum(flow[, i] * cost[, i])
        if (average) res <- res / sum(flow[, i]) # sûr que i est à droite?
        return(res) # / sum(solution[-i, i])
        #       sorties            entrées
    } else if ((sum(flow[i, ]) > sum(flow[, i])) &
               ((type == "out") | (type == "both"))) {
        res <- sum(flow[i, ] * cost[i, ])
        if (average) res <- res / sum(flow[i, ])
        return(res)# /
    } else {
        return(0)
    }
}

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
costTot <- function(flow, cost) {
    sum(vapply(1:nrow(flow),
               function(x) costToti(x, flow, cost, type = "in"),
               numeric(1)))
}

#------------------------------------------------------------------------------#
# standardised and dimensionless clustering index (nui) for a donor
# Perry
#' @export
#------------------------------------------------------------------------------#
clust_P <- function(flow, cost, dim_mat, nperm, start, end,
                    method = "shortsimplex", cost_flows, threads) {
    idx <- seq_len(nrow(flow)) # or: seq_len(length(start))
    randomisations <- pbapply::pblapply(seq_len(nperm), function(i1) {
        # Yi & Yc
        rand_idx      <- sample(idx)
        new_start     <- start[rand_idx] # new_start n'a pas de nom de col et row contrairement à start (????)
        opt_transport <- wrap_transport(new_start, end, cost, method)
        opt_transport <- as.matrix(opt_transport, dim_mat) # Replace n with N

        # ~ Assez long ci-dessous
        idx_same_count <- lapply(idx, function(i2) which(new_start == start[i2]))
        # res_Yci <- rep(NA, times = length(idx_same_count))
        # lapply(seq_len(length(res_Yci)), function(j) {
        #     res_Yci[j] <<- mean(vapply(seq_len(length(idx_same_count[[j]])),
        #                                function(x) {
        #                                    if (cost_flows[[j]] < 0) type <- "in"
        #                                    else                     type <- "out"
        #                                    #costToti(idx_same_count[[j]][x], opt_transport, cost, type = type)
        #                                    costTotiCPP(idx_same_count[[j]][x], opt_transport, cost, type)
        #                                },
        #                                numeric(1L)))
        # })
        res_Yci <- vapply(seq_len(length(idx_same_count)), function(i2) {
            mean(vapply(seq_len(length(idx_same_count[[i2]])),
                        function(x) {
                            if (cost_flows[[i2]] < 0) type <- "in"
                            else                      type <- "out"
                            #costToti(idx_same_count[[i2]][x], opt_transport, cost, type = type)
                            costTotiCPP(idx_same_count[[i2]][x], opt_transport, cost, type)
                        },
                        numeric(1L)))
        }, numeric(1L))
        res_Yii <- vapply(idx,
                          function(x) {
                              if (cost_flows[[x]] < 0) type <- "in"
                              else                     type <- "out"
                              #costToti(x, opt_transport, cost, type = type)
                              costTotiCPP(x, opt_transport, cost, type)
                          },
                          numeric(1L))

        list(flow = opt_transport,
             idx  = rand_idx,
             Yci = res_Yci,
             Yii = res_Yii,
             #costTot = costTot(opt_transport, cost))) # Pour calculer Ia
             costTot = costTotCPP(opt_transport, cost)) # To compute Ia
    }, cl = threads)

    Ycs <- simplify2array(lapply(randomisations, function(x) x[["Yci"]]))
    Ycs <- abs(rowMeans(Ycs))
    Yis <- simplify2array(lapply(randomisations, function(x) x[["Yii"]]))
    Yis <- abs(rowMeans(Yis))
    Y0  <- mean(Ycs) / 2
    Yi  <- vapply(idx,
                  function(x) {
                      if (cost_flows[[x]] < 0) type <- "in"
                      else                     type <- "out"
                      costTotiCPP(x, flow, cost, type)
                  },
                  numeric(1L))

    list(clust = ((Yi * Y0) / (Yis * Ycs)),
         Ea = unlist(lapply(randomisations, function(x) x[["costTot"]])))
}

#------------------------------------------------------------------------------#
# standardised and dimensionless clustering index (nui) for a donor
#' @export
#------------------------------------------------------------------------------#
clust_LMX <- function(flow, cost, dim_mat, nperm, start, end,
                             method = "shortsimplex", cost_flows,
                             threads) {
    idx            <- seq_len(nrow(flow)) # or: seq_len(length(start))
    #base           <- list(rep(NA, nrow(flow)))
    #randomisations <- rep(base, nperm) # NON : + 1 pour ajout du "Di" (le vrai !) à la fin

    # New idx : Li-Madden-Xu
    #for (i in seq_len(nperm)) { ## NEW
    #lapply(seq_len(nperm), function(i1) {
    randomisations <- pbapply::pblapply(seq_len(nperm), function(i1) {
        #for (j in idx) { ## NEW
        vapply(idx, function(i2) {
            sub_idx   <- idx[idx != i2] ## NEW
            rand_idx  <- sample(sub_idx)
            ##TO DOUBLE CHECK## rand_idx <- insert(rand_idx, i2, i2) ##NEW : insert is in pkgg R.utils à récupérer uniquement fn intéressante en cpp après
            rand_idx  <- append(rand_idx, i2, i2)
            new_start <- start[rand_idx]
            #dimnames(new_start) <- list(rownames(start), colnames(start))
            opt_transport <- wrap_transport(new_start, end, cost, method)
            opt_transport <- as.matrix(opt_transport, dim_mat)

            # imp!
            if (cost_flows[[i2]] < 0) type <- "in"  # NEW NEW
            else                      type <- "out" # NEW NEW
            costTotiCPP(i2, opt_transport, cost, type)
        }, numeric(1L))

        # lapply(idx, function(j) {
        #     sub_idx   <- idx[idx != j] ## NEW
        #     rand_idx  <- sample(sub_idx)
        #     ##TO DOUBLE CHECK## rand_idx <- insert(rand_idx, j, j) ##NEW : insert is in pkgg R.utils à récupérer uniquement fn intéressante en cpp après
        #     rand_idx  <- append(rand_idx, j, j)
        #     new_start <- start[rand_idx]
        #     #dimnames(new_start) <- list(rownames(start), colnames(start))
        #     opt_transport <- wrap_transport(new_start, end, cost, method)
        #     opt_transport <- as.matrix(opt_transport, dim_mat)
        #
        #     # imp!
        #     if (cost_flows[[j]] < 0) type <- "in"  # NEW NEW
        #     else                     type <- "out" # NEW NEW
        #     randomisations[[i]][[j]] <<- costTotiCPP(j, opt_transport, cost, type = type)
        # })
    }, cl = threads)

    # Di, le vrai! ... NOT NEED EN THEORIE CAR DEJA CALCULER avec la var cost_flows
    #for (j in idx) {
    #lapply(idx, function(j) { ### Déjà présent dans le "main.R" // redondant
    # pbapply::pblapply(idx, function(j) { ### Déjà présent dans le "main.R" // redondant
    #     #idx        <- idx
    #     #new_start <- matrix(start[idx], nrow = nrow(start))
    #     new_start <- start ### A simplifier
    #     opt_transport <- wrap_transport(new_start, end, cost, method)
    #     opt_transport <- as.matrix(opt_transport, dim_mat)
    #
    #     # imp!
    #     if (cost_flows[[j]] < 0) type <- "in"  # NEW NEW
    #     else                     type <- "out" # NEW NEW
    #     randomisations[[length(randomisations)]][[j]] <<- costTotiCPP(j, opt_transport, cost, type = type)
    # }, cl = threads)

    randomisations <- simplify2array(randomisations)
    Dis_bar        <- rowMeans(cbind(randomisations, cost_flows)) # Négatif tous pour l'instant. cost_flows = Dis
    Dis_bar        <- abs(Dis_bar)
    return(list(clust = (cost_flows / Dis_bar), # positif pour l'instant par construction
                Dis_all = randomisations))
}

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
Ia <- function(Da, Ea) {return(Da/Ea)} # Da = observed, Ea = mean of randomized distances
# Perry_1995, explication Pa et Ia



#========



#------------------------------------------------------------------------------#
#' @include sadie.R
#' @export
#------------------------------------------------------------------------------#
plot.sadie <- function(x, y, ..., isoclines = FALSE, onlySignificant = FALSE,
                       resolution = rep(100, 2)) { # + specify significancy we want (0.95,...)

    data1 <- x$info_clust

    dataLoess <- stats::loess(idx_P ~ x * y, data = data1,
                              degree = 2, span = 0.2)
    x <- seq(min(data1$x), max(data1$x), length = resolution[1]) # xResolution
    y <- seq(min(data1$y), max(data1$y), length = resolution[2]) # yResolution
    interpolated <- predict(dataLoess, expand.grid(x = x, y = y))

    data2 <- data.frame(expand.grid(x = x, y = y), z = as.vector(interpolated))

    data3 <- data.frame(lapply(split(data2, data2$y), function(XX) XX$z))

    thres <- c(-1.5, 1.5)
    cpt <- function(x) {
        if (x < 0) {
            if (x < thres[1]) return("black") # low
            else              return("white")  # no
        } else {
            if (x > thres[2]) return("red") # high
            else              return("white")   # no
        }
    }
    data1$col <- vapply(data1$idx_P, cpt, character(1))
    #data1$col <- "no-se"

    #g <- ggplot(data2, aes(x = x, y = y, z = z)) +
    #    geom_raster(aes(fill = z)) + # tile ou rectangle plut^ot, car pas toujours des carrés !!!
    #    geom_contour(color = "white", size = 0.25) +
    #    geom_point(data = data1, inherit.aes = FALSE,
    #               aes(x = x, y = y,
    #                   size = abs(idx_P)),
    #               #                                 fill = as.factor(col)),
    #               colour = "black", pch = 21) +
    #    scale_fill_gradient(low = "white", high = "black")
    ##scale_fill_distiller(palette = "Spectral")
    ##theme_bw()
    ##geom_tile() +
    ##    stat_contour(color="white", size=0.25) +
    ##    viridis::scale_fill_viridis() +
    ##    ggthemes::ggthemestheme_tufte(base_family="Helvetica")
    #print(g)



    if (isoclines) {
        filled.contour(x = unique(data2$x), y = unique(data2$y), z = as.matrix(data3),
                       color.palette = terrain.colors,
                       plot.axes = {
                           with(data1, points(x, y, pch = 21,
                                              cex = abs(idx_P),
                                              col = "black", bg = col));
                           contour(x = unique(data2$x), y = unique(data2$y),
                                   z = as.matrix(data3),
                                   col = "black", lty = "solid", add = TRUE,
                                   vfont = c("sans serif", "plain"))
                       })
    } else {
        with(data1, plot(x, y, pch = 21, cex = abs(idx_P),
                         col = "black", bg = col))
    }
}
















