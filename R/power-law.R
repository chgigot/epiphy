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
#' my_power_law
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
power_law <- function(list, log_base = exp(1), na.rm = FALSE, ...) { # S'OCCUPER DES NA

    # Checks:
    stopifnot(is.list(list))
    if (length(list) < 2) {
        stop("Less than 2 points is not enough to perform linear regressions.")
    }
    object_class <- unique(vapply(list, function(x) class(x)[[1]], character(1)))
    stopifnot(length(object_class) == 1)
    stopifnot(object_class %in% c("count", "incidence"))

    # Perform power law analysis:
    switch(object_class,
           "count" = {
               data <- lapply(list, function(obj) map_data(obj)$r)
               x    <- vapply(data, function(obj) mean(obj, na.rm = na.rm), numeric(1L))
               y    <- vapply(data, function(obj) var(obj, na.rm = na.rm), numeric(1L))
           },
           "incidence" = {
               data <- lapply(list, function(obj) {
                   mapped_data <- map_data(obj)
                   data.frame(p = (mapped_data$r / mapped_data$n),
                              n = mapped_data$n)
               })
               # For incidence data as proportions:
               # v_t = p(1 - p)/n (Madden & Hughes, 1995)
               x    <- vapply(data, function(obj) {
                   with(obj, (mean(p, na.rm = na.rm) * (1 - mean(p, na.rm = na.rm))) / mean(n, na.rm = na.rm))
               }, numeric(1L))
               y    <- vapply(data, function(obj) var(obj$p, na.rm = na.rm), numeric(1L))
           }
    )
    coord_obs <- data.frame(x = x, y = y)
    model_formula <- as.formula(bquote(
        log(y, base = .(log_base)) ~ log(x, base = .(log_base))
    ))
    model     <- lm(model_formula, ...)
    y_the     <- predict(model, type = "response")
    coord_the <- data.frame(x = x, y = log_base^(y_the))

    # Retrieve summary matrice of coefficients, and eventually add some extra
    # estimates:
    par <- coef(summary(model))
    switch (object_class,
        "count" = {
            # Nothing to do.
        },
        "incidence" = {
            n  <- mean(vapply(data, function(obj) mean(obj$n, na.rm = na.rm), numeric(1L)), na.rm = na.rm) ### PAS TOP
            Ar <- estimateCoef(model, bquote(.(log_base)^x1))
            ar <- estimateCoef(model, bquote(.(log_base)^x1 * .(n)^(-x2)))
            AR <- estimateCoef(model, bquote(.(log_base)^x1 * .(n)^(2 * (1 - x2))))
            aR <- estimateCoef(model, bquote(.(log_base)^x1 * .(n)^(2 - x2)))
            par <- rbind(par, unlist(Ar), unlist(ar), unlist(AR), unlist(aR))
            rownames(par) <- c("log_base(Ar)", "b", "Ar", "ar", "AR", "aR")
        }
    )

    # Return the following object:
    structure(list(call      = match.call(),
                   data      = data,
                   model     = model,
                   par       = par,
                   log_base  = log_base,
                   coord_obs = coord_obs,
                   coord_the = coord_the),
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
    cat("# Power law analysis:\n")
    printCoefmat(x$par)
}

# For count


#------------------------------------------------------------------------------#
#' a2a: A wAy to pAinlessly switch between different power LAw formulAtions
#'
#' \code{a2a} was designed to avoid headaches that are likely to occur when
#' working with different formulations of the binomial power law analysis.
#'
#' The binomial power law can be expressed as: \eqn{s_y^2 = (intercept)(s_{bin}^2)^b}.
#' But different forms of (intercept) are possible depending on the formulation of the
#' binomial power law.
#' \tabular{ccccc}{
#'       \tab Ar         \tab ar      \tab AR         \tab aR      \cr
#'    Ar \tab 1          \tab n^b     \tab n^{2(b-1)} \tab n^{b-2} \cr
#'    ar \tab n^{-b}     \tab 1       \tab n^{b-2}    \tab n^{-2}  \cr
#'    AR \tab n^{2(1-b)} \tab n^{2-b} \tab 1          \tab n^{-b}  \cr
#'    aR \tab n^{2-b}    \tab n^2     \tab n^b        \tab 1       \cr
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
#' a2a(from = , to = , n = , b = )
#' a2a(to = , data = <powerLaw object>)
#'
#' @export
#------------------------------------------------------------------------------#
a2a <- function(intercept, from = c("Ar", "ar", "AR", "aR"),
                  to = c("Ar", "ar", "AR", "aR"), slope, n) {
    from <- match.arg(from)
    to   <- match.arg(to)
    b    <- slope
    dico <- expand.grid(from = c("Ar", "ar", "AR", "aR"),
                        to = c("Ar", "ar", "AR", "aR"),
                        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    #-----------------------------------------------------------------#
    #             | col Ar     |col ar  | col AR     | col aR
    dico$coef <- c(1,           n^b,     n^(2*(b-1)), n^(b-2), # row Ar
                   n^(-b),      1,       n^(b-2),     n^(-2),  # row ar
                   n^(2*(1-b)), n^(2-b), 1,           n^(-b),  # row AR
                   n^(2-b),     n^2,     n^b,         1)       # row aR
    #-----------------------------------------------------------------#
    item <- dico[which(dico$from == from & dico$to == to), ]
    res  <- intercept * item[["coef"]]
    attr(res, "par") <- c(intercept = intercept,
                          coef      = item[["coef"]],
                          slope     = b,
                          n         = n)
    attr(res, "class") <- c("a2a", "numeric")
    res
}

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
print.a2a <- function(x, ...) cat(x, "\n", sep = "")

#------------------------------------------------------------------------------#
# @export
#------------------------------------------------------------------------------#
#aaaaa.powerLaw <- function(obj, to)




