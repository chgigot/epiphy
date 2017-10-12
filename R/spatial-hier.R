#------------------------------------------------------------------------------#
#' @include utils.R
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
spatialHier <- function(low, high) {

    # Initial checks
    if (missing(low) || missing(high))
        stop("low and high must be specified.")
    # V'erifier type de x et y, afin d'associer la bonne m'ethode
    # automatiquement
    if (length(low) != length(high))
        stop("'low' and 'high' lengths differ.")
    # penser à verif partout pareil, et pas seulement [[1]] !!!
    nLow  <- unique(low[[1]]$obs$n) # pas très propre
    if (length(nLow) != 1) stop("Function only work for same n. so far.")
    n <- nLow

    nHigh <- unique(high[[1]]$obs$n) # pas très propre
    if (length(nHigh) != 1) stop("Function only work for same n. so far.")
    if (nHigh != 1)         stop("n of high must be equal to 1.")

    pLow  <- vapply(low,  function(x) {sum(x$obs$d) / sum(x$obs$n)}, numeric(1))
    pHigh <- vapply(high, function(x) {sum(x$obs$d) / sum(x$obs$n)}, numeric(1))

    ###
    # Modif. here?
    ###


    ##x <- log10(1 - pLow)
    ##y <- log10(1 - pHigh)

    xy <- data.frame(x = cloglog(pLow),
                     y = cloglog(pHigh))

    # Keep rows with only finite values
    xy <- xy[is.finite(xy$x) & is.finite(xy$y), ]

    ##x.tmp <- x
    ##x <- x[!is.infinite(x) & !is.infinite(y)] # Peut-^etre travailler avec une df , plus prorpe?
    ##y <- y[!is.infinite(x.tmp) & !is.infinite(y)]

    model    <- lm(y ~ offset(x), data = xy)
    # ... is the same as:
    #model    <- lm((y - x) ~ 1, data = xy) # Just an intercept.
    ##logNu    <- unname(coef(model))
    nu       <- unname(exp(coef(model)))
    coordObs <- data.frame(x = pLow, y = pHigh)
    ##yThe     <- logNu + xy$x
    ##coordThe <- data.frame(x = cloglog(xy$x, rev = TRUE), y = cloglog(yThe, rev = TRUE))
    ##rownames(coordThe) <- rownames(xy) # gérer les noms row plus proprement / xy et pas coordObs a cause des infinite supprimé
    # Other poss:
    coordThe <- data.frame(x = pLow, y = 1 - (1 - pLow)^nu)

    ##model    <- lm(y ~ x + 0)
    ##yThe     <- predict(model, data.frame(x), type="response")
    ##coordObs <- data.frame(x = (1 - 10^x), y = (1 - 10^y))
    ##coordObs <- data.frame(x = pLow, y = pHigh)
    ##coordThe <- data.frame(x = (1 - 10^x), y = (1 - 10^yThe))

    #model2 <- glm(pHigh ~ offset(I(cloglog(pLow))), family = quasibinomial(link = "cloglog"))
    ###----------------------------------------------------------------------###
    ### Corrected case:

    ##nHigh   <- vapply(high, function(x) sum(x@obs$n), numeric(1)) # Attention, il y a un nHigh au dessus !!! Garder celui-là et trouver autre chose pour celui du dessus

    ##xy2 <- data.frame(x = cloglog(pLow),
    ##                  y = pHigh,
    ##                  w = nHigh)

    # Keep rows with only finite values
    ##xy2 <- xy2[is.finite(xy2$x), ]

    # GLM model
    ##model2 <- glm(y ~ offset(x), data = xy2, family = binomial(link = "cloglog"), weights = w)
    #res <- predict(model2, type = "response", newdata = data.frame(pLow = 0), se.fit = TRUE)
    #nu <- unname(exp(res$fit))

    ##nu <- unname(exp(coef(model2)))

    ##coordThe <- data.frame(x = pLow,
    ##                       y = 1 - (1 - pLow)^nu)

    structure(list(call = match.call(),
                   model = model,
                   nu = nu,
                   n = n,
                   coordObs = coordObs,
                   coordThe = coordThe),
              class = "spatialHier")
}

#------------------------------------------------------------------------------#
# setMethod("initialize",
#           signature(.Object = "Relationship"),
#           definition = function(.Object, x, y) {
#               if (missing(x) || missing(y))
#                   stop("x and y must be specified.")
#               # V'erifier type de x et y, afin d'associer la bonne m'ethode
#               # automatiquement
#               if (length(x) != length(y))
#                   stop("'x' and 'y' lengths differ.")
#
#               ## Only for incidence spatial hierarchies now
#               res <- spatialHier2(x, y)
#
#               .Object@n        <- res$n
#               .Object@model    <- res$model
#               .Object@coordObs <- res$coordObs
#               .Object@coordThe <- res$coordThe
#
#               return(.Object)
#           }
# )


#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
summary.spatialHier <- function(object, ...) {
    # Retrieve result matrice (to which we will add extra estimates)
    param <- coef(summary(object$model))

    baseLog <- exp(1)
    nu <- estimateCoef(object$model, bquote(.(baseLog)^x1))

    param <- rbind(param, unlist(nu))
    rownames(param) <- c("log_base(nu)", "nu")

    structure(list(coefficients = param),
              class = "summary.spatialHier")
}

#------------------------------------------------------------------------------#
#' @method print summary.spatialHier
#' @export
#------------------------------------------------------------------------------#
print.summary.spatialHier <- function(x, ...) {
    cat("Coefficients:\n")
    printCoefmat(x$coefficients)
}


#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
plot.spatialHier <- function(x, y, type = c("regular", "cloglog"),
                             size = 2, col = "black", ...) {

    type <- match.arg(type)

    if (type == "regular") {
        # Observed
        g <- ggplot(data = x$coordObs, aes(x = x, y = y)) +
            scale_x_continuous(limits = c(0,1)) +
            scale_y_continuous(limits = c(0,1)) +
            geom_point(color=col, size=size)

         # Theo
         g <- g + geom_line(data = x$coordThe, aes(x = x, y = y), color = col)

        # Binomial theorical
        xBin <- seq(0, 1, by = 0.01) # attention si on utilise x ici (au lieu de xBin), ća écrase le parametre d'entree
        yBin <- 1 - (1 - xBin)^x$n
        theo <- data.frame(x = xBin, y = yBin)
        g <- g + geom_line(data = theo, aes(x = x, y = y), linetype = 2) +
            labs(x = expression(p[low]), y = expression(p[high]))

        ## Offrir cette possibilité dans les options
        g <- g + theme_bw()
        print(g)

    } else if (type == "cloglog") {
        # Observed
        g <- ggplot(data = cloglog(x$coordObs), aes(x = x, y = y)) +
            # scale_x_continuous(limits = c(0,1)) +
            # scale_y_continuous(limits = c(0,1)) +
            geom_point(color=col, size=size)

        # Theo
        g <- g + geom_line(data = cloglog(x$coordThe), aes(x = x, y = y), color = col)

        # Binomial theorical
        xBin <- seq(0, 1, by = 0.001) # attention si on utilise x ici (au lieu de xBin), ća écrase le parametre d'entree
        yBin <- 1 - (1 - xBin)^x$n
        theo <- data.frame(x = xBin, y = yBin)
        g <- g + geom_line(data = cloglog(theo), aes(x = x, y = y), linetype = 2) +
            labs(x = expression(p[low]), y = expression(p[high]))
        print(g)

    } else stop("type must be 'regular' or 'log-log'.")

}







