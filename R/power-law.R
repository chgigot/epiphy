#------------------------------------------------------------------------------#
#' @include utils.R
#------------------------------------------------------------------------------#
NULL

#==============================================================================#
#' Taylor's and Binomial Power Laws
#'
#' These functions allow to...
#'
#' \deqn{log_{10}(s_y^{2}) = log_{10}(a) + b log_{10}(Y)}
#' with \eqn{Y} is equivalent to \eqn{\bar{y}(1 - \bar{y})} in case of incidence
#' data, and \eqn{\bar{y}} in case of count data. \eqn{y} correspond to the
#' proportion of diseased plants or plant parts in each sampling unit, whereas
#' \eqn{y} corresponds to the absolute value in case of data count.
#'
#' @export
#------------------------------------------------------------------------------#
powerLaw <- function(x, baseLog = 10,...) {

    if (missing(x)) stop("x must be specified.")

    n       <- unique(x[[1]]$n) ### TRES TRES SAL !!! En particulier si on travaille avec du Poisson !!!

    # if(!is.list(x)) => error
    # len <- length(x)
    # if(len == 0) => error
    # if(len == 1) => warning
    # type <- is(x[[1]])
    # switch(type,
    #        Incidence = {}, type Incidence à checher pour tous éléments de la liste
    #        Count = {}, type Count à checher pour tous éléments de la liste
    #        error("Aie"))
    #### Attention, aucune protection ici !!!!!
    #if (!is.IncidenceGroup(x) && !is.CountGroup(x)) {
    #stop("x must be a valid IncidenceGroup or CountGroup object.")
    #} else {

    classObjs <- class(x[[1]])

    powerLawFn <- function(list, type, baseLog, ...) {
        if (length(list) == 1) {
            stop("Only 1 point is not enough to perform a linear regression.")
        }
        switch(type,
               "incidence" = {
                   datas <- lapply(list, function(x) data.frame(freq = (x@obs$d / x@obs$n), n = x@obs$n))
                   y     <- vapply(datas, function(x) var(x$freq), numeric(1)) # y before x to avoid to crach x too early... find another name rather than x for the argument (object)
                   x     <- vapply(datas, function(x) (mean(x$freq) * (1 - mean(x$freq))) / x$n[1], numeric(1)) # Work only if n the same everywheerer
               },
               "count" = {
                   datas <- lapply(list, function(x) x@obs$d) # Peut ^etre choisir autre chose que x, dans x@obs$d ?
                   y     <- vapply(datas, function(x) var(x), numeric(1))
                   x     <- vapply(datas, function(x) mean(x), numeric(1))
               }
        )

        modelFormula <- as.formula(substitute(
            log(y, base = baseLog) ~ log(x, base = baseLog),
            list(baseLog = baseLog)))
        model    <- lm(modelFormula, ...)
        coordObs <- data.frame(x = x, y = y)
        yThe     <- predict(model, data.frame(log10(x)), type="response")
        coordThe <- data.frame(x = x, y = baseLog^(yThe))
        return(list(model = model,
                    coordObs = coordObs,
                    coordThe = coordThe))
    }

    x <- powerLawFn(x, classObjs, baseLog, ...) # x[[1]] car une liste d'objet de meme type (a verifier avant)

    # Warning message:
    # 'newdata' had 2 rows but variables found have 6 rows
    #return(.Object)

    baseLogAd <- b <- list()
    baseLogAd$est  <- coefficients(x[[1]])[[1]]
    b$est          <- coefficients(x[[1]])[[2]]

    # Retrieve result matrice (to which we will add extra estimates)
    param <- coef(summary(x$model))

    if (classObjs == "incidence") {

        Ad <- estimateCoef(x$model, bquote(.(baseLog)^x1))
        ad <- estimateCoef(x$model, bquote(.(baseLog)^x1 * .(n)^(-x2)))
        AD <- estimateCoef(x$model, bquote(.(baseLog)^x1 * .(n)^(2 * (1 - x2))))
        aD <- estimateCoef(x$model, bquote(.(baseLog)^x1 * .(n)^(2 - x2)))

        param <- rbind(param, unlist(Ad), unlist(ad), unlist(AD), unlist(aD))
        rownames(param) <- c("log_base(Ap)", "b", "Ap", "ap", "An", "an")

    }

    structure(list(call = match.call(),
                   model = x[[1]],
                   par = param,
                   n = n,
                   baseLog = baseLog,
                   coordObs = x[[2]],
                   coordThe = x[[3]]),
              class = "powerLaw")
}

### Restructuration à prévoir ici: Plus de IncidenceGroup ou CountGroup... berk!
### Un truc du genre :
#verif <- function(list) {
#    if(!is.list(list)) stop("Err:/")
#    switch(is(list[[1]]),
#           Incidence={type <- "Incidence"},
#           Count={type <- "Count"},
#           stop("Err:/")
#    )
#    if (!any(sapply(a, is.Incidence))) stop("Err:/")
#}
########## ETC

#==============================================================================#
#' @export
#------------------------------------------------------------------------------#
plot.powerLaw <- function(x, y, col = "black", size = 2, observed = TRUE,
                          model = TRUE, bisector = TRUE, print = TRUE, type, ...) { # Pas très propre il me semble ?

    if (missing(type)) type <- "log"

    if (type == "log") {

        baseLog <- x$baseLog
        if (baseLog == exp(1)) nameBaseLog <- "e" else nameBaseLog <- baseLog
        g <- ggplot(data = log(x$coordObs, base = baseLog), aes(x = x, y = y))
        minxy <- log(min(min(x$coordObs$x), min(x$coordObs$y)), base = baseLog)
        maxxy <- log(max(max(x$coordObs$x), max(x$coordObs$y)), base = baseLog)
        g <- g + scale_x_continuous(limits = c(minxy, maxxy))
        g <- g + scale_y_continuous(limits = c(minxy, maxxy))
        g <- g + labs(x = bquote(log[.(nameBaseLog)] ~ "(binomial variance)"),
                      y = bquote(log[.(nameBaseLog)] ~ "(observed variance)"))
        if (observed) g <- g + geom_point(color=col, size=size)
        if (model)    g <- g + geom_line(data = log(x$coordThe, base = baseLog), aes(x = x, y = y), color=col) ## Not necessary aes
        if (bisector) g <- g + geom_abline(intercept = 0, slope = 1, linetype=2, color=col)

        ## Offrir cette possibilité dans les options
        g <- g + theme_bw()
        return(g) # can only print apparently .... I did not succeed in return a ggplot object !!
        ## SI SI maintenant ća marche, il faut retourner g, en non print(g) !!!!

    } else if (type == "regular") {

        g <- ggplot(data = x$coordObs, aes(x=x, y=y))
        #minxy <- min(min(x@coordObs$x), min(x@coordObs$y))
        #maxxy <- max(max(x@coordObs$x), max(x@coordObs$y))
        #g <- g + scale_x_continuous(limits = c(minxy, maxxy))
        #g <- g + scale_y_continuous(limits = c(minxy, maxxy))
        g <- g + labs(x = expression(s[bin]^2), y = expression(s[obs]^2))
        if (observed) g <- g + geom_point(color=col, size=size)
        if (model)    g <- g + geom_line(data = x$coordThe, aes(x = x, y = y), color=col) ## Not necessary aes
        #if (bisector) g <- g + geom_abline(intercept = 0, slope = 1, linetype = 2, color = col)
        return(g) # can only print apparently .... I did not succeed in return a ggplot object !!

    } else stop("type must be 'regular' or 'log'.")
}

#==============================================================================#
#' @export
#------------------------------------------------------------------------------#
print.powerLaw <- function(x, ...) {
    cat("\nPower Law Analysis:\n")
    printCoefmat(x$par)
}

# For count



