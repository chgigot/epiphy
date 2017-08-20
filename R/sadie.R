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
#'     or a matrix (x-y)
#' @param cost A \code{length(count)} by \code{length(count)} matrix specifing
#'     the cost of transporting single unit between the corresponding source and
#'     destination sampling unit. If not specified, the distance between each
#'     couple of adjacent sampling units (rook's sense) is assumed to be equal
#'     to 1.
#' @param index The name of the index to use: "perry", "li-madden-xu (LMX)" or "all".
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
#' #aphid_counts
#' #arthropods_counts
#' #codling_moth_counts
#' res <- sadie(aphid_counts)
#' summary(res)
#' plot(res)
#'
#' @name sadie
#' @export
#------------------------------------------------------------------------------#
sadie <- function(data, index, nperm,
                  rseed, seed, cost) {
    UseMethod("sadie")
}

#------------------------------------------------------------------------------#
#' @rdname sadie
#' @method sadie data.frame
#' @export
#------------------------------------------------------------------------------#
sadie.data.frame <- function(data, index = c("perry", "LMX", "all"), # needed to specify possible indices here cause of match.arg
                             nperm = 100, rseed = TRUE, seed = 12345, cost) { # Apparement en S3, il faut spécifier les arg par défaut ici, à verfier

    index <- match.arg(index)
    if (missing(cost)) {
        cost <- as.matrix(dist(data[, c("x", "y")]))#, method = "manhattan"))
    }

    # IMPORTANT: count is a matrix, now

    # Checking and data wrangling

    # More interesting
    if (!rseed) set.seed(seed)

    # A faire plus propre / generic ci-dessous
    count <- tidyr::spread(data, y, d) # Ne plus utiliser tidyr
    row.names(count) <- count$x
    count <- as.matrix(dplyr::select(count, -x)) # Ne plus utiliser dplyr

    n <- prod(dim(count))
    fcount <- matrix(data = rep(mean(count), n),
                     nrow = nrow(count))
    dimnames(fcount) <- list(rownames(count), colnames(count)) # Useful????

    #optTransDT  <- wrapTransport(count, fcount)
    #optTransMat <- optTransDT2Mat(optTransDT)
    optTransport <- optTransDT2Mat(wrapTransport(count, fcount), n) # Bizarre que l'on prend pas cost en compte ici, sans doute basé sur les noms des row/col, à vérfier (et amélirer en reprennant le code C++ par exemple)
    ##optTransport <- optTransDT2Mat(wrapTransportEMD(count, fcount), n)


    Da <- costTot(optTransport, cost) # C'est quoi Da, le cost total, OK, trouver un nom plus explicite, ou bien le préciser en commentaire

    cost_of_flow <- sapply(1:n, function(x) costToti(x, optTransport, cost, type = "both"))

    if (index == "all") index <- c("perry", "LMX")
    perry_ <- list(clustering = NA_real_, Ea = NA_real_)
    li_madden_xu <- list(clustering = NA_real_, Dis_all = matrix(NA_real_))
    new_prob <- NA_real_
    if (any(index == "perry")) {
        perry_ <- clusteringIdx(optTransport, cost, n, nperm,
                                dataBeg = count, dataEnd = fcount,
                                method = "shortsimplex", cost_of_flow)
    }
    if (any(index == "LMX")) {
        li_madden_xu <- clusteringIdxNew(optTransport, cost, n, nperm,
                                         dataBeg = count, dataEnd = fcount,
                                         method = "shortsimplex", cost_of_flow)

        tmp <- li_madden_xu[["Dis_all"]][, 1:nperm] # On ne prend pas le dernier = le "vrai"
        Dis <- li_madden_xu[["Dis_all"]][, (nperm + 1 )]## un peu idiot, déjà calculer avec Da <- costTot(optTransMat, cost)
        res2 <- abs(tmp) > abs(Dis)
        res2 <- rowSums(res2)
        res2 <- (res2 + 1) / (nperm + 1) # res2 + 1 => rank of Di(s)
        new_prob <- res2

    }
    # Ci-dessous, ça veut dire que l'on calcule toujours le perry's
    Ia <- Da / mean(perry_[["Ea"]])
    Pa <- sum((-perry_[["Ea"]]) > (-Da)) / # On met des moins pour mémoire  à cause de cette $%@*&^ de convention inversée !!!
        length(perry_[["Ea"]])

    clusteringIndices <- data.frame(data,
                                    cost.of.outflow   = cost_of_flow,
                                    original.index = perry_[["clustering"]],
                                    new.index      = li_madden_xu[["clustering"]],
                                    prob           = new_prob)

    # Summary indices
    summary_indices <- data.frame(overall = c(mean(abs(clusteringIndices$original.index)), NA),
                                  inflow  = c(mean(clusteringIndices$original.index[clusteringIndices$original.index < 0]), NA),
                                  outflow = c(mean(clusteringIndices$original.index[clusteringIndices$original.index > 0]), NA))
    rownames(summary_indices) <- c("Perry's index", "Li-Madden-Xu's index")

    #    return(Sadie())
    res <- list(clusteringIndices = clusteringIndices,
                Da = Da,
                Ia = Ia,
                Pa = Pa,
                Ea_perry = perry_[["Ea"]],
                Ea_li = li_madden_xu[["Dis_all"]],
                summary_indices = summary_indices,
                nperm = nperm,
                rseed = rseed,
                seed = NA_real_)
    attr(res, "class") <- "sadie"
    attr(res, "call")  <- match.call()
    res
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
    print(head(x$clusteringIndices, n = n))
    #printCoefmat(x$clusteringIndices)
    #cat("number of permutations: ", object$nperm, "\n",
    #    "random seed: ", object$rseed, "\n", sep = "")
    cat("\nSummary indices:\n")
    print(x$summary_indices)
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
wrapTransport <- function(dataBeg, dataEnd, method = "shortsimplex") {
    stdout <- vector('character')
    con    <- textConnection('stdout', 'wr', local = TRUE)
    sink(con)
    gen <- list(as.numeric(rownames(dataEnd)), as.numeric(colnames(dataEnd))) # Parceque on a toujours les noms des colonnes et lignes dans dataEnd, penser à les ajouter dans dataBegbis pour calcul idx old et new
    res <- transport::transport(transport::pgrid(dataBeg, generator = gen),
                                transport::pgrid(dataEnd, generator = gen),
                                method = method)
    sink()
    close(con)
    return(res)
}

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
wrapTransportEMD <- function(dataBeg, dataEnd) {
    emdist::emd2d(dataBeg, dataEnd) # remplacer par la version C de vecteurs (et non de matrices)... d'aillerus meme chose pour transport (eviter ggrid)
}

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
optTransDT2Mat <- function(x, suLen) {
    res <- matrix(rep(0, suLen^2), nrow = suLen, ncol = suLen)
    lapply(1:nrow(x), function(j) {
        res[x[j, 1], x[j, 2]] <<- x[j, 3]
    })
    return(res)
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
clusteringIdx <- function(flow, cost, n, nperm, dataBeg, dataEnd, method = "shortsimplex", cost_of_flow) {
    seqNrowFlow <- seq_len(nrow(flow))
    progress       <- txtProgressBar(min = 0, max = nperm, style = 3)
    #randomisations <- parallel::mclapply(seq_len(nperm), function(i) {
    randomisations <- lapply(seq_len(nperm), function(i) {
        # Yi & Yc
        idx       <- seqNrowFlow
        randIdx   <- sample(idx)
        dataBegbis <- matrix(dataBeg[randIdx], nrow = nrow(dataBeg))
        #dimnames(dataBegbis) <- list(rownames(dataBeg), colnames(dataBeg))
        optTransDT <- wrapTransport(dataBegbis, dataEnd, method)
        ##optTransDT <- wrapTransportEMD(dataBegbis, dataEnd)

        # ~ Assez long ci-dessous
        optTransMat <- optTransDT2Mat(optTransDT, n)
        idxSameCounts <- lapply(idx, function(j) which(dataBegbis == dataBeg[j]))
        tmp <- vector("list", length(idxSameCounts))
        lapply(seq_len(length(tmp)), function(j) {
            tmp[[j]] <<- mean(vapply(seq_len(length(idxSameCounts[[j]])),
                                     function(x) {
                                         if (cost_of_flow[[j]] < 0) type <- "in"
                                         else                       type <- "out"
                                         costToti(idxSameCounts[[j]][x], optTransMat, cost, type = type)
                                     },
                                     numeric(1)))
        })
        resYci <- unlist(tmp) ## Vraiment pertinent que tmp soit liste ???
        resYii <- vapply(idx,
                         function(x) {
                             if (cost_of_flow[[x]] < 0) type <- "in"
                             else                       type <- "out"
                             costToti(x, optTransMat, cost, type = type)
                         },
                         numeric(1))
        setTxtProgressBar(progress, value = i)

        return(list(flow = optTransMat,
                    idx  = randIdx,
                    Yci = resYci,
                    Yii = resYii,
                    costTot = costTot(optTransMat, cost))) # Pour calculer Ia
    })
    #}, mc.cores = 8)

    Ycs <- rowMeans(simplify2array(lapply(randomisations, function(x) x[["Yci"]])))
    Ycs <- abs(Ycs)
    Yis <- rowMeans(simplify2array(lapply(randomisations, function(x) x[["Yii"]])))
    Yis <- abs(Yis)
    # Pas sûr du tout que les Yii soient bien calculés


    Y0  <- mean(Ycs) / 2
    #Y0  <- mean(Ycs)
    Yi  <- vapply(seqNrowFlow,
                  function(x) {
                      if (cost_of_flow[[x]] < 0) type <- "in"
                      else                       type <- "out"
                      #type <- "in"
                      costToti(x, flow, cost, type = type)
                  },
                  numeric(1))

    return(list(clustering = ((Yi * Y0) / (Yis * Ycs)),
                Ea = unlist(lapply(randomisations, function(x) x[["costTot"]]))))
}

#------------------------------------------------------------------------------#
# standardised and dimensionless clustering index (nui) for a donor
#' @export
#------------------------------------------------------------------------------#
clusteringIdxNew <- function(flow, cost, n, nperm, dataBeg, dataEnd, method = "shortsimplex", cost_of_flow) {
    seqNrowFlow <- seq_len(nrow(flow))
    base           <- list(rep(NA, nrow(flow)))
    randomisations <- rep(base, (nperm + 1)) # + 1 pour ajout du "Di" (le vrai !) à la fin
    progress       <- txtProgressBar(min = 0, max = nperm, style = 3)

    # New idx : LMX
    #for (i in seq_len(nperm)) { ## NEW
    lapply(seq_len(nperm), function(i) {
        #for (j in seqNrowFlow) { ## NEW
        lapply(seqNrowFlow, function(j) {
            idx     <- seqNrowFlow
            idx     <- idx[idx != j] ## NEW
            randIdx <- sample(idx)
            ##TO DOUBLE CHECK## randIdx <- insert(randIdx, j, j) ##NEW : insert is in pkgg R.utils à récupérer uniquement fn intéressante en cpp après
            randIdx <- append(randIdx, j, j)
            dataBegbis <- matrix(dataBeg[randIdx], nrow = nrow(dataBeg))
            #dimnames(dataBegbis) <- list(rownames(dataBeg), colnames(dataBeg))
            optTransDT <- wrapTransport(dataBegbis, dataEnd, method)

            # ~ Assez long ci-dessous
            optTransMat <- optTransDT2Mat(optTransDT, n)

            # imp!
            if (cost_of_flow[[j]] < 0) type <- "in"  # NEW NEW
            else                       type <- "out" # NEW NEW
            #type <- "in"
            randomisations[[i]][[j]] <<- costToti(j, optTransMat, cost, type = type)
        })
        setTxtProgressBar(progress, value = i) ## Ajout du +1 après si integration de la boucle ci dessous
    })

    # Di, le vrai!
    #for (j in seqNrowFlow) {
    lapply(seqNrowFlow, function(j) { ### Déjà présent dans le "main.R" // redondant
        idx        <- seqNrowFlow
        #dataBegbis <- matrix(dataBeg[idx], nrow = nrow(dataBeg))
        dataBegbis <- dataBeg ### A simplifier
        optTransDT <- wrapTransport(dataBegbis, dataEnd, method)

        # ~ Assez long ci-dessous
        optTransMat <- optTransDT2Mat(optTransDT, n)

        # imp!
        if (cost_of_flow[[j]] < 0) type <- "in"  # NEW NEW
        else                       type <- "out" # NEW NEW
        #type <- "in"
        randomisations[[length(randomisations)]][[j]] <<- costToti(j, optTransMat, cost, type = type)
    })

    randomisations <- simplify2array(randomisations)
    Dis_bar        <- rowMeans(randomisations) # Négatif tous pour l'instant
    Dis_bar        <- abs(Dis_bar)
    Dis            <- randomisations[, (nperm + 1)] # Négatif tous pour l'instant
    return(list(clustering = (Dis / Dis_bar), # positif pour l'instant par construction
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
plot.sadie <- function(x, y, ..., onlySignificant = FALSE, resolution = rep(100, 2)) { # + specify significancy we want (0.95,...)

    data1 <- x$clusteringIndices

    dataLoess <- stats::loess(original.index ~ x * y, data = data1,
                              degree = 2, span = 0.2)
    x <- seq(min(data1$x), max(data1$x), length = resolution[1]) # xResolution
    y <- seq(min(data1$y), max(data1$y), length = resolution[2]) # yResolution
    interpolated <- predict(dataLoess, expand.grid(x = x, y = y))

    data2 <- data.frame(expand.grid(x = x, y = y), z = as.vector(interpolated))

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
    data1$col <- vapply(data1$original.index, cpt, character(1))
    #data1$col <- "no-se"

    #g <- ggplot(data2, aes(x = x, y = y, z = z)) +
    #    geom_raster(aes(fill = z)) + # tile ou rectangle plut^ot, car pas toujours des carrés !!!
    #    geom_contour(color = "white", size = 0.25) +
    #    geom_point(data = data1, inherit.aes = FALSE,
    #               aes(x = x, y = y,
    #                   size = abs(original.index)),
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



    with(data1, plot(x, y, pch = 21, cex = abs(original.index),
                     col = "black", bg = col))

}
















