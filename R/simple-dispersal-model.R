#------------------------------------------------------------------------------#
#' Simple dispersal model
#'
#' Even if this model does not fall into the scope stricto sensu of the package
#' \code{epiphy}, it is provided to perform some assessments and test
#' tmplemented statistical methods.
#'
#' The model is not directly made available to the user. One may import all the
#' utilities related to this simple dispersal model using:
#' \code{invisible(list2env(epiphy:::simple_model, environment()))}.
#'
#' @param nfoci The number of initial infections.
#' @param xrate Multiplication rate, the number of new infections produced by
#'     each pre-existing infection at each generation.
#' @param ngen The number of generations.
#' @param lambda The rate parameter for the exponential distribution.
#'
#' @references
#'
#' Esker PD, Sparks AH, Antony G, Bates M, Dall' Acqua W, Frank EE, Huebel L,
#' Segovia V, Garrett KA, 2007. Ecology and Epidemiology in R: Modeling
#' dispersal gradients. The Plant Health Instructor.
#' \href{http://dx.doi.org/10.1094/PHI-A-2007-1226-03}{doi:10.1094/PHI-A-2007-1226-03}
#'
#' @examples
#' # To use the "hiden" functions of this simple model:
#' invisible(list2env(epiphy:::simple_model, environment()))
#'
#' set.seed(12345)
#' foci <- disperse(nfoci = 1, xrate = 3, ngen = 10, lambda = 20, ngen_active = 1)
#' plot(foci, cex = 0.1)
#' quad <- quadrat(surf_dim = c(x = 1, y = 1), nint = c(x = 90, y = 90))
#' collection <- collect(foci, quad)
#'
#' # To perform analyses:
#' collection$n <- 1
#' my_data <- incidence(collection, mapping(r = incidence))
#' my_data <- clump(my_data, unit_size = c(x = 3, y = 3)) # /!\ avec clump(...), on retire les colonnes non mappées
#' res <- fit_distr(my_data)
#' res
#' summary(res)
#' plot(res)
#'
#' lambdas <- c(5, 10, 20, 40, 80)
#'
#' # Multi-data set analyses
#' set.seed(12345)
#' quad <- quadrat(surf_dim = c(x = 1, y = 1), nint = c(x = 90, y = 90))
#' my_data <- list()
#' for (i in 1:100) {
#'     nfoci <- sample(1:100, 1)
#'     xrate <- sample(5:15, 1)
#'     print(paste0("set ", i, " to be done (nfoci = ", nfoci, "; xrate = ", xrate,")"))
#'     foci <- disperse(nfoci = nfoci, xrate = xrate, ngen = 2, lambda = 20)
#'     my_data[[i]] <- collect(foci, quad)
#' }
#' # saveRDS(my_data, "simple_model_data_final_5-15.rds")
#' # res <- readRDS("simple_model_data_final_5-15.rds")
#'
#' ## Power law
#' my_data2 <- lapply(my_data, function(x) {
#'     x$n <- 1
#'     x <- incidence(x, mapping(r = incidence))
#'     clump(x, unit_size = c(x = 3, y = 3))
#' })
#'
#' res <- power_law(my_data2)
#' res
#' summary(res)
#' plot(res)
#'
#' ## Spatial hierarchy
#' low <- my_data2
#' high <- lapply(my_data2, level_up)
#'
#' res <- spatial_hier(low, high)
#' res
#' summary(res)
#' plot(res)
#'
#' @keywords internal
#------------------------------------------------------------------------------#
simple_model <- list()

simple_model$disperse <- function(nfoci, xrate, lambda, ngen, ngen_active = 0) {
    res <- as.data.frame(dispersalCPP(nfoci, xrate, lambda, ngen, ngen_active))
    colnames(res) <- c("x", "y", "gen")
    structure(res, class = c("disperse", "data.frame"))
}

simple_model$plot.disperse <- function(x, ...) {
    plot(x$x, x$y, pch = 19, col = "red", xlim = c(0, 1), ylim = c(0, 1), ...)
}

simple_model$quadrat <- function(surf_dim = c(x = 1, y = 1),
                                 nint = c(x = 90, y = 90)) {
    xquad <- surf_dim[["x"]] / nint[["x"]]
    yquad <- surf_dim[["y"]] / nint[["y"]]

    prex <- head(seq(0, surf_dim[["x"]], by = xquad), -1)
    prey <- head(seq(0, surf_dim[["y"]], by = yquad), -1)

    prex <- cbind.data.frame(x = seq_len(length(prex)), x1 = prex)
    prey <- cbind.data.frame(y = seq_len(length(prey)), y1 = prey)

    quad <- expand.grid(x1 = prex[["x1"]], y1 = prey[["y1"]])
    quad <- merge(merge(quad, prex, by = "x1"), prey, by = "y1")
    quad <- quad[c("x", "y", "x1", "y1")]
    quad <- quad[with(quad, order(x, y)), ]
    quad <- cbind(quad, x2 = quad$x1 + xquad, y2 = quad$y1 + yquad)
    structure(quad, class = c("quadrat", "data.frame"))
}

simple_model$collect_ <- function(disperse, quadrat) {
    res <- as.data.frame(collectCPP(as.matrix(disperse), as.matrix(quadrat)))
    res <- setNames(res, c("x", "y", "x1", "y1", "x2", "y2",
                           "count", "incidence"))
    structure(res, class = c("collect", "quadrat", "data.frame"))
}

simple_model$collect <- function(disperse, quadrat) {
    quadrat$count     <- NA
    quadrat$incidence <- NA
    lapply(seq_len(nrow(quadrat)), function(i1) {
        inx <- (quadrat[i1, "x1"] < disperse[, "x"]) & (disperse[, "x"] <= quadrat[i1, "x2"])
        iny <- (quadrat[i1, "y1"] < disperse[, "y"]) & (disperse[, "y"] <= quadrat[i1, "y2"])
        quadrat[i1, "count"]     <<- sum(inx * iny) # Nice trick!
        quadrat[i1, "incidence"] <<- ifelse(quadrat[i1, "count"] == 0, 0, 1)
    })
    structure(quadrat, class = c("collect", "quadrat", "data.frame"))
}

simple_model$print2.collect <- function(x, ..., type = c("incidence", "count")) {
    type <- match.arg(type)
    #if (!add) {
    #plot(0, type = "n", xlim = xlim, ylim = ylim, xlab = "x", ylab = "y",
    #     asp = asp, ...)
    plot(0, type = "n", xlim = c(0, 1), ylim = xlim, ...)
    #}
    x$coord <- I(lapply(seq_len(nrow(x)), function(i1) {
        coord <- rbind(unlist(unname(x[i1, c("x1", "y1")])),
              unlist(unname(x[i1, c("x1", "y2")])),
              unlist(unname(x[i1, c("x2", "y2")])),
              unlist(unname(x[i1, c("x2", "y1")])))
        dimnames(coord) <- list(NULL, c("x", "y"))
        polygon(coord, border = "grey", col = x[i1, "incidence"])
        # We do not need to close the box, it's supposed by "polygon"
    }))
    #x$col <- pal(length(0:round(maxi)))[round(x[[type]])]
    lapply(x$coord, function(xxx) polygon(xxx, border = "grey"))
    polygon(x$coord, border = "grey", col = x$incidence)# x$col)
}







# Figures
#simple_model$plot.disperse <- function(x, ...) {
#    plot(x, pch = 19, col = "red", xlim = c(0, 1), ylim = c(0, 1), ...)
#}

#simple_model$plot.quadrat <- ...

# simple_model$plot.collect <- function(x, ..., xlim = c(0, 1), ylim = xlim,
#                                       asp = 1) {


# }




