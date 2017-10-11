#------------------------------------------------------------------------------#
#' @include intensity-classes.R
#------------------------------------------------------------------------------#
NULL

#------------------------------------------------------------------------------#
#' Taylor\'s and binary power laws
#'
#' Assesses the overall degree of heterogeneity in a collection of data sets at
#' the sampling-unit scale.
#'
#' The power law describes the relationship between the observed variance of
#' individuals within a data set (\code{s^2}) and the corresponding variance
#' under the assumption of no aggregation (\code{s\'^2}). It can be expressed
#' under its logarithmic form as: \code{log(s^2) = log(a) + b log(Y)}, with:
#' \itemize{
#'     \item \code{Y = p} in the case of count data (Taylor\'s power law).
#'     \item \code{Y = p(1 - p)} in the case of incidence data (binary power law).
#' }
#' \code{p} corresponds to the mean proportion of recorded individuals in case
#' of incidence data, and the absolute value in case of count data.
#'
#' @param x A list of \code{intensity} objects (\code{count} or
#'     \code{incidence} objects).
#' @param log_base Logarithm base to be used.
#' @param ... Not yet implemented.
#'
#' @return A \code{power_law} object.
#'
#' @examples
#' require(magrittr)
#' my_data <- do.call(c, lapply(citrus_ctv, function(citrus_field) {
#'    incidence(citrus_field) %>%
#'        clump(unit_size = c(x = 3, y = 3)) %>%
#'        split(by = "t")
#' }))
#' # my_data is a list of incidence object, each one corresponding to a given
#' # given time at a given location.
#' my_power_law <- power_law(my_data)
#' summary(my_power_law)
#' plot(my_power_law)
#'
#' @references
#'
#' Taylor LR. 1961. Aggregation, variance and the mean. Nature 189: 732–35.
#'
#' Hughes G, Madden LV. 1992. Aggregation and incidence of disease. Plant
#' Pathology 41 (6): 657–660.
#' \href{http://dx.doi.org/10.1111/j.1365-3059.1992.tb02549.x}{doi:10.1111/j.1365-3059.1992.tb02549.x}
#'
#' Madden LV, Hughes G, van den Bosch F. 2007. Spatial aspects of epidemics -
#' III: Patterns of plant disease. In: The study of plant disease epidemics,
#' 235–78. American Phytopathological Society, St Paul, MN.
#'
#' @export
#------------------------------------------------------------------------------#
power_law <- function(x, log_base = exp(1), ...) {

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

    powerLawFn <- function(list, type, log_base, ...) {
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
            log(y, base = log_base) ~ log(x, base = log_base),
            list(log_base = log_base)))
        model    <- lm(modelFormula, ...)
        coordObs <- data.frame(x = x, y = y)
        yThe     <- predict(model, data.frame(log10(x)), type="response")
        coordThe <- data.frame(x = x, y = log_base^(yThe))
        return(list(model = model,
                    coordObs = coordObs,
                    coordThe = coordThe))
    }

    x <- powerLawFn(x, classObjs, log_base, ...) # x[[1]] car une liste d'objet de meme type (a verifier avant)

    # Warning message:
    # 'newdata' had 2 rows but variables found have 6 rows
    #return(.Object)

    baseLogAd <- b <- list()
    baseLogAd$est  <- coefficients(x[[1]])[[1]]
    b$est          <- coefficients(x[[1]])[[2]]

    # Retrieve result matrice (to which we will add extra estimates)
    param <- coef(summary(x$model))

    if (classObjs == "incidence") {

        Ad <- estimateCoef(x$model, bquote(.(log_base)^x1))
        ad <- estimateCoef(x$model, bquote(.(log_base)^x1 * .(n)^(-x2)))
        AD <- estimateCoef(x$model, bquote(.(log_base)^x1 * .(n)^(2 * (1 - x2))))
        aD <- estimateCoef(x$model, bquote(.(log_base)^x1 * .(n)^(2 - x2)))

        param <- rbind(param, unlist(Ad), unlist(ad), unlist(AD), unlist(aD))
        rownames(param) <- c("log_base(Ap)", "b", "Ap", "ap", "An", "an")

    }

    structure(list(call = match.call(),
                   model = x[[1]],
                   par = param,
                   n = n,
                   log_base = log_base,
                   coordObs = x[[2]],
                   coordThe = x[[3]]),
              class = "power_law")
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

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
plot.power_law <- function(x, y, col = "black", size = 2, observed = TRUE,
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

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
print.power_law <- function(x, ...) {
    cat("\nPower Law Analysis:\n")
    printCoefmat(x$par)
}

# For count


#------------------------------------------------------------------------------#
#' aaaaa: A wAy to pAinlessly switch between different power LAw formulAtions
#'
#' \code{aaaaa} was designed to avoid headaches that are likely to occur when
#' working with different formulations of the binomial power law analysis.
#'
#' The binomial power law can be expressed as: \eqn{s_y^2 = (intercept)(s_{bin}^2)^b}.
#' But different forms of (intercept) are possible depending on the formulation of the
#' binomial power law.
#' \tabular{ccccc}{
#'       \tab Ad         \tab ad      \tab AD         \tab aD      \cr
#'    Ad \tab 1          \tab n^b     \tab n^{2(b-1)} \tab n^{b-2} \cr
#'    ad \tab n^{-b}     \tab 1       \tab n^{b-2}    \tab n^{-2}  \cr
#'    AD \tab n^{2(1-b)} \tab n^{2-b} \tab 1          \tab n^{-b}  \cr
#'    aD \tab n^{2-b}    \tab n^2     \tab n^b        \tab 1       \cr
#' }
#'
#' @param intercept Intercept parameter to be converted.
#' @param from Kind of the input intercept parameter.
#' @param to Desired kind for the ouput intercept parameter.
#' @param slope Slope parameter.
#' @param n Number of individuals per sampling unit.
#'
#' @examples
#'
#' aaaaa(from = , to = , n = , b = )
#' aaaaa(to = , data = <powerLaw object>)
#'
#' @export
#------------------------------------------------------------------------------#
aaaaa <- function(intercept, from = c("Ad", "ad", "AD", "aD"),
                  to = c("Ad", "ad", "AD", "aD"), slope, n) {
    from <- match.arg(from)
    to   <- match.arg(to)
    b    <- slope
    dico <- expand.grid(from = c("Ad", "ad", "AD", "aD"),
                        to = c("Ad", "ad", "AD", "aD"),
                        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    #             | col Ad     |col ad  | col AD     | col aD
    dico$coef <- c(1,           n^b,     n^(2*(b-1)), n^(b-2), # row Ad
                   n^(-b),      1,       n^(b-2),     n^(-2),  # row ad
                   n^(2*(1-b)), n^(2-b), 1,           n^(-b),  # row AD
                   n^(2-b),     n^2,     n^b,         1)       # row aD
    item <- dico[which(dico$from == from & dico$to == to), ]
    res <- intercept * item[["coef"]]
    attr(res, "params") <- c(intercept = intercept, coef = item[["coef"]],
                             slope = b, n = n)
    res
}

#------------------------------------------------------------------------------#
# @export
#------------------------------------------------------------------------------#
#aaaaa.powerLaw <- function(obj, to)




