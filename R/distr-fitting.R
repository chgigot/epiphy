#------------------------------------------------------------------------------#
#' @include betabinom.R
#' @include intensity-classes.R
#------------------------------------------------------------------------------#

optimBetaBinom <- function(inits, x, n) {
    # log-likelihood function to minimize
    lbetabin <- function(inits, x, n) {
        sum(-dbetabinom(x = x, size = n, shape1 = inits[1], shape2 = inits[2], log = TRUE))
    }

    epsilon <- 1e-7
    model <- NULL
    try(model <- optim(par = inits,
                       fn = lbetabin,  ############# ATTENTION, meme si try(..) encore des warnings pour chisq.test... par la suite.
                       x = x, n = n,
                       method = "L-BFGS-B",
                       hessian = TRUE,
                       lower = c(epsilon, epsilon)),
        silent = TRUE)
    ##### +++ control=list(maxit=300,trace=TRUE,REPORT=1)) ############################
    structure(c(model, list(nobs = length(x))),
              class = c("optim", "list")) # return "model"
}

# avec bbmle::mle2 plut besoin de ce logLik là
#' @export
logLik.optim <- function(object, ...) {
    res <- (- object$value) # Do not forget the minus
    attributes(res) <- list(nobs = object$nobs, df = NA)
    structure(res, class = "logLik")
}

epsilon <- 1e-7

### Compute expected beta-binomial probabilities using approach in book from
# Madden et al, 2007
# pbbinom <- function( q, size, prob, theta ) { # q, size belongs to Int
#     for (x in q) {
#         alpha <- prob / theta
#         beta  <- (1-prob)/theta
#         y     <- 1
#         for (i in 0:(size-1))
#             y <- y * ((1-prob+i*theta) / (1+i*theta))
#         if (x!=0) {
#             k <- y
#             for (j in 1:x) {
#                 k <- k * ((size+1-j)/j) * ((alpha-1+j)/(size+beta-j))
#                 y <- y + k
#             }
#         }
#         if ( !exists("val", inherits = FALSE) ) {val <- y} else {val <- c(val, y)} # SI NON EXIST VAL DANS L'ENVIRONNEMENT DE LA FONCTION !!! TRES IMPORTANT
#     }
#     return(val)
# }
#
# dbetabinom2 <- function(x, size, prob, theta) {
#     values <- pbbinom(x, size, prob, theta)
#     (values - c(0, values[1:max(x)])) # Pas propre // temporaire
# }

#==============================================================================#
# @export
#bbmle::summary

# @export
#bbmle::coef

# @export
#bbmle::logLik

# @export
#bbmle::vcov


#------------------------------------------------------------------------------#
#' Factory for Maximum Likelihood Estimation (MLE) functions
#'
#' Creates functions with the signature \code{function(data, start, ...)} in
#' order to performe maximum likelihood estimation. It's based on the
#' \code{\link[bbmle]{mle2}} (package \code{bbmle}) function to perform the
#' optimization. The default optimizer handle by \code{mle2} is the built-in
#' \code{optim}, but other optimizer are possible (see \code{mle2} documentation
#' for more details).
#'
#' Two functions are already implemented for convenience reasons: \code{mleBetabinom} and
#' \code{mleNbinom}. But, using the tool \code{fitDiffDistr}, one should not handle such
#' functions. Nevertheless, it should be usefull if you want to use an alternative
#' function for the aggregation case inside the \code{fitDiffDistr} function.
#'
#' @param .negLogLik Negative log-likelihood function.
#' @param .data Data set to work with. Must be a named list (A data frame is a
#'   named list.)
#' @param .start Starting values for the parameters of interest. Must be a named
#'   list.
#' @param .preEstimate Function to get initial values to work with. Could be a
#'   function computing the parameters specified in \code{.start} based on
#'   method of moment estimation (MME).
#' @param ... Extra parameters which will becomes default parameters for the
#'   created function.
#'
#' @seealso fitDiffDistr
#'
#' @examples
#' rand <- glm(d/n ~ 1, family = binomial, data = dataMadden1987, weights = n)
#' p <- estimateCoef(rand, quote(1/(1+exp(-x1))))
#' # p <- estimateCoef(rand, quote(logistic(x1)))
#' params <- rbind(coef(summary(rand)), unlist(p))
#' rownames(params) <- c("logit(p)", "p")
#' printCoefmat(params)
#'
#'
#' agg <- mleBetabinom(data = dataMadden1987, lower = list(alpha = -Inf, beta = -Inf))
#' # coef(summary(agg))
#'
#' # to compute p, one can use the function estimateCoef
#' # see ?estimateCoef for more information
#' p     <- estimateCoef(agg, quote(x1 / (x1 + x2)), student = FALSE)
#' theta <- estimateCoef(agg, quote(1  / (x1 + x2)), student = FALSE)
#' rho   <- estimateCoef(agg, quote(1  / (x1 + x2 + 1)), student = FALSE)
#' params <- rbind(coef(summary(agg)), unlist(p), unlist(theta), unlist(rho))
#' rownames(params) <- c("alpha", "beta", "p", "theta", "rho")
#' printCoefmat(params)
#'
#' # Computing agg, is the same as doing the following hard way :
#' agg2 <- mleFactory(
#'     .negLogLik = function(d, n, alpha, beta) {
#'         sum(-dbetabinom(x = d, size = n, shape1 = alpha, shape2 = beta,
#'         log = TRUE)) },
#'     .data = c("d", "n"),
#'     .start = c("alpha", "beta"),
#'     .preEstimate = epiphy:::mmeBetabinom,
#'     method = "L-BFGS-B",
#'     lower = c(alpha = 1e-7, beta = 1e-7))
#' agg2 <- agg2(data = dataMadden1987)
#' coef(summary(agg2))
#'
#' @export
#------------------------------------------------------------------------------#
mleFactory <- function(.negLogLik, .data, .start, .preEstimate = NULL, ...) {
    # Checks
    argNames <- names(formals(.negLogLik))
    stopifnot(!missing(.data) || !missing(.start))
    stopifnot(is.character(.data) && is.character(.start))
    stopifnot(all(argNames %in% c(.data, .start)))
    # Empty names ("") not allowed... to check

    # Extra paramters to become defaults for the created function
    defaultParams <- list(...)

    function(data, start, ...) {

        # Initial checks
        stopifnot(!missing(data))
        stopifnot(is.list(data))
        stopifnot(all(names(data) %in% .data))

        if (missing(start)) start <- list()
        stopifnot(is.list(start))
        lapply(seq_len(length(.start)), function(i) {
            if (is.null(start[[.start[i]]])) { # Access to a non-existing element returns NULL
                stopifnot(!is.null(.preEstimate))
                argPreEstimate <- names(formals(.preEstimate))
                stopifnot(all(argPreEstimate %in% c("data", "name")))
                start[[.start[i]]] <<- .preEstimate(data, name = .start[i])
            }
        })

        extraParams <- list(...)
        baseParams <- list(minuslogl = .negLogLik,
                           data = data,
                           start = start)

        params <- modifyList(defaultParams, extraParams)
        params <- modifyList(params, baseParams)

        # MLE computation
        model <- NULL
        try(model <- do.call(bbmle::mle2, params), silent = TRUE)
        #if (!is.null(model)) model <- as(model, "est")
        # Astuce pour eviter de tout coller la sortie construite par do.call
        if (!is.null(model)) model@call <- match.call()
        if (!is.null(model)) model@call.orig <- match.call()
        model
    }
}

# negative log-likelihood function to minimize
nllBetabinom <- function(d, n, alpha, beta) {
    sum(-dbetabinom(x = d, size = n, shape1 = alpha, shape2 = beta, log = TRUE))
}

#mmeBetabinom <- function(x, n = max(x)) {
### DOUBLE-check MME !!!!
mmeBetabinom <- function(data, name) {
    x <- data$d
    n <- data$n[[1]] # if pas de n, alors n = max(x) /// et c'est pas tres propre
    m1 <- mean(x)
    m2 <- mean(x^2)
    s2 <- var(x) # x != proportion, but number of diseased individual per sampling unit.
    switch (name,
            "prob"  = m1 / n,
            "theta" = (s2 - n * m1 * (1 - m1)) / (n^2 * m1 * (1 - m1) - s2),
            "alpha" = (n * m1 - m2) / (n * (m2 / m1 - m1 - 1) + m1),
            "beta"  = (n - m1) * (n - m2 / m1) / (n * (m2 / m1 - m1 - 1) + m1)
    )
}


#------------------------------------------------------------------------------#
#' @rdname mleFactory
#' @export
#------------------------------------------------------------------------------#
mleBetabinom <- mleFactory(.negLogLik = nllBetabinom,
                           .data = c("d", "n"),
                           .start = c("alpha", "beta"),
                           .preEstimate = mmeBetabinom,
                           method = "L-BFGS-B",
                           lower = c(alpha = epsilon, beta = epsilon))

#------------------------------------------------------------------------------#
# Estimation of additional parameters
#
# Based on the delta method (msm::deltamethod())
#
#------------------------------------------------------------------------------#

#extraCoefFactory()

#extraCoefBinom()
#extraCoefBetabinom()

#extraCoefPoisson()
#extraCoefNbinom()

# mleBetabinom <- mleFactory(.negLogLik = nllBetabinom,
#                            .data = c("d", "n"),
#                            .start = c("alpha", "beta"),
#                            .preEstimate = mmeBetabinom,
#                            .extraCoef = list(
#                                p     = quote(x1 / (x1 + x2)),
#                                theta = quote(1  / (x1 + x2)),
#                                rho   = quote(1  / (x1 + x2 + 1)),
#                            )
#                            method = "L-BFGS-B",
#                            lower = c(alpha = epsilon, beta = epsilon))


## if (colnames(coef(summary(agg)))[3] == "z value")

#myfunc <- function(v1) {
#    deparse(substitute(v1))
#}
#
#myfunc(foo)
#[1] "foo"

# log-likelihood function to minimize
nllNbinom <- function(d, mu, k) {
    sum(-dnbinom(x = d, mu = mu, size = k, log = TRUE))
}

# Method of Moments Estimation
mmeNbinom <- function(data, name) {
    x <- data$d
    m1 <- mean(x)
    m2 <- mean(x^2)
    s2 <- var(x)
    switch (name,
            "mu" = m1,
            "k" = m1^2 / (s2 - m1) # Madden et al, 2007, p. 249
    )
}

#------------------------------------------------------------------------------#
#' @rdname mleFactory
#' @export
#------------------------------------------------------------------------------#
mleNbinom <- mleFactory(.negLogLik = nllNbinom,
                        .data = c("d"),
                        .start = c("mu", "k"),
                        .preEstimate = mmeNbinom,
                        method = "L-BFGS-B",
                        lower = c(mu = epsilon, k = epsilon))


### !!!!!!!!!!!!!!!
# lower non pris en compte via cette manière
# D'où des resultatsd legerement differentes entre :
# mleNbinom(x=dataMadden1987$d)
# ... et ...
# mleNbinom2(list(x=dataMadden1987$d))
# Doit faire ceci pour corriger :
# mleNbinom2(list(x=dataMadden1987$d), lower = c(mu=epsilon, k=epsilon))


#mleNbinom(data)

# glmFactory(formula, family, ...) {
#     function(data) {
#         glm(formula = formula, family = family,..., data = data)
#     }
# }
#
# glmBinom(data)
# glmPoisson(data)

#==============================================================================#
#==============================================================================#
# mu1 <- function(x) {
#     stopifnot(is.vector(x))
#     mean(x)
# }
#
# mu2 <- function(x) {
#     stopifnot(is.vector(x))
#     mean(x^2)
# }

#------------------------------------------------------------------------------#
# Do not export
#------------------------------------------------------------------------------#
#mme <- list()

#------------------------------------------------------------------------------#
# Maximum likelihood
#
# @param x A vector of data.
# @param n A vector of data (size > 1 supported)
# @param start Named list with inital parameter values.
# @param epsilon Limit.
#
# @name mle
#------------------------------------------------------------------------------#
#NULL

#------------------------------------------------------------------------------#
# Beta-binomial distribution
#------------------------------------------------------------------------------#
# Method of Moments Estimation

#### TO DOUBLE CHECK below !!!!!!!!!!!!!!! ######
# mme$betabinom <- function(x, n = max(x)) {
#     m1 <- mu1(x)
#     m2 <- mu2(x)
#     s2 <- var(x) # x != proportion, but number of diseased individual per sampling unit.
#     list(prob  = mu1(x) / n,
#          theta = (s2 - n * m1 * (1 - m1)) / (n^2 * m1 * (1 - m1) - s2),
#          alpha = (n * m1 - m2) / (n * (m2 / m1 - m1 - 1) + m1),
#          beta  = (n - m1) * (n - m2 / m1) / (n * (m2 / m1 - m1 - 1) + m1))
# }

#------------------------------------------------------------------------------#
# @rdname mle
# @export
#------------------------------------------------------------------------------#
# mleBetabinom <- function(x, n = max(x), start, epsilon = 1e-7) {
#     # log-likelihood function to minimize
#     llbetabinom <- function(x, n, alpha, beta) {
#         sum(-dbetabinom(x = x, size = n, shape1 = alpha, shape2 = beta, log = TRUE))
#     }
#
#     if (missing(start)) start <- list()
#     stopifnot(is.list(start))
#     if (is.null(start$alpha)) start$alpha <- mme$betabinom(x = x, n = n)$alpha
#     if (is.null(start$beta))  start$beta  <- mme$betabinom(x = x, n = n)$beta
#
#     model <- NULL
#     try(model <- bbmle::mle2(llbetabinom,
#                              data = list(x = x, n = n),
#                              start = list(alpha = start$alpha, beta = start$beta),
#                              lower = c(alpha = epsilon, beta = epsilon),
#                              method = "L-BFGS-B"), # alpha > 0 and beta > 0
#         silent = TRUE)
#     model
# }

#------------------------------------------------------------------------------#
# Negative binomial distribution
#
# The function \code{\link[MASS]{glm.nb}} can be used as well. But as the returned
# object is not "standard", in the sense that for example theta is not reachable
# through the coefficient table (with all the nice extra information: standard
# error, p-values,...). So to deal with that, a consistent interface is used in
# this package, based on bbmle::mle2() function which offers a nice and
# consistent interface.
#------------------------------------------------------------------------------#
# Method of Moments Estimation
# mme$nbinom <- function(x) {
#     list(mu = mu1(x),
#          k = mu1(x)^2 / (var(x) - mu1(x)))
# }

#------------------------------------------------------------------------------#
# @rdname mle
# @export
#------------------------------------------------------------------------------#
# mleNbinom <- function(x, start, epsilon = 1e-7) {
#     # log-likelihood function to minimize
#     llnbinom <- function(x, mu, k) {
#         # Alternative parametrization, with mu and k (see ?dnbinom)
#         sum(-dnbinom(x = x, mu = mu, size = k, log = TRUE))
#     }
#
#     if (missing(start)) start <- list()
#     stopifnot(is.list(start))
#     if (is.null(start$mu)) start$mu <- mme$nbinom(x = x)$mu
#     if (is.null(start$k))  start$k  <- mme$nbinom(x = x)$k
#
#     model <- NULL
#     try(model <- bbmle::mle2(llnbinom,
#                              data = list(x = x),
#                              start = list(mu = start$mu, k = start$k),
#                              lower = c(mu = epsilon, k = epsilon),
#                              method = "L-BFGS-B"), # mu > 0 and k > 0
#         silent = TRUE)
#     model
# }


#------------------------------------------------------------------------------#
# Notes
#------------------------------------------------------------------------------#
# exp((AICmin−AICi)/2) peut être compris comme la probabilité pour que le ième candidat modèle minimise l'estimation de la perte d'information
# ref:Burnham et Anderson 2002, §6.4.5
#exp((118.9-152.1)/2)
#------------------------------------------------------------------------------#
# Christensen (1990 *) said: " The likelihood ratio test statistic is G2 = −2[log L(m0) − log L( m)], where ˆm0 is the MLE of m under the assumption that H0 is true and m is the MLE under the “unrestricted” model ".
# * Christensen, R. (1990)Log-Linear Models and Logistic Regression. Springer-Verlag New York, Inc., pp. 332-336.
#------------------------------------------------------------------------------#
# http://stackoverflow.com/questions/1826519/function-returning-more-than-one-value
# https://stat.ethz.ch/pipermail/r-help/2004-June/053343.html
# From: G. Grothendieck

# Some refs:
# - About negbinom (tutorial):
# http://rstudio-pubs-static.s3.amazonaws.com/2516_d204dc109a2c44d58899afb418ae3885.html
# http://data.princeton.edu/wws509/r/overdispersion.html
# http://stats.stackexchange.com/questions/10419/what-is-theta-in-a-negative-binomial-regression-fitted-with-r
# - About Hessian and Co.
# http://stats.stackexchange.com/questions/27033/in-r-given-an-output-from-optim-with-a-hessian-matrix-how-to-calculate-paramet
# http://stats.stackexchange.com/questions/68080/basic-question-about-fisher-information-matrix-and-relationship-to-hessian-and-s






#----------------------------------------------




#------------------------------------------------------------------------------#
# Distribution fitting
# These function allow to...
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#' Try to fit distributions to data sets
#'
#' Depending of the nature of the data submitted, different distributions are
#' tryied to be fitted to the data.
#'
#' Two distributions are used. One under the assumption of no aggregation and
#' the other one under the assumption of aggregation.
#'
#' @param x A data of the class \code{Intensity}.
#' @param type The kind of distribution to fit: "random", "aggregated" or "both".
#'
#' @return An object of class \code{DistFitting}, with the following components:
#'
#' \tabular{ll}{
#'     \code{from} \tab Kind of the \code{Intensity} object. \cr
#'     \code{random} \tab Theoretical data under random assumption. \cr
#'     \code{aggregated} \tab Theoretical data under aggregated assumption. \cr
#'     \code{freqTable} \tab Entry 4 \cr
#' }
#'
#' @references
#'
#' @examples
#' \dontrun{
#' data <- Incidence(Cochran1936)
#' fitDistr(data)
#' }
#'
#' @name fitDistr
#' @export
#------------------------------------------------------------------------------#
fitDistr <- function(object, ...) {
    UseMethod("fitDistr")
}

#------------------------------------------------------------------------------#
#' @describeIn fitDistr
#' For incidence data, the two distributions fitted are binomial and beta-binomial
#' distributions.
#'
#' @export
#------------------------------------------------------------------------------#
fitDistr.incidence <- function(object, ..., progress = TRUE, simulatePValue) { # Extra arguments useful?
    # Add a parameter? gof = c("Chisq", "LLR", "no") ???
    #
    # ...
    #
    # if (!requireNamespace("XNomial", quietly = TRUE) && gof == "LLR") {
    #     warning("Package 'XNomial' needed to use the option 'LLR'. Auto-switch to 'Chisq'.")
    #     gof <- "Chisq"
    # }
    #
    # ...
    #

    # Initial checks and data preparation
    d <- object$obs$d # All the raw data
    N <- length(d)
    n <- unique(object$obs$n)
    if (length(n) != 1) stop(paste0("Current implementation only deals ",
                                    "with equal size sampling units."))

    # Summary table (with observed and theoretical frequencies)
    freq <- as.data.frame(table(d))
    names(freq) <- c("category","observed")
    freq <- merge(x = data.frame(category = 0:n),
                  y = freq,
                  by = "category", all = TRUE)
    freq[is.na(freq)] <- 0
    freq[] <- lapply(freq, as.numeric) # Force each column to be numeric (and only numeric!)

    # Binomial distribution
    rand <- list()
    rand$model <- glm(d/n ~ 1, family = binomial, weights = rep(n, N))

    lp <- predict(rand$model, type = "link", se.fit = TRUE)
    ## Ci-dessous, on utilise le fait que lp a retourné une liste
    ## pourrait faire : lp[[residual.scale]] <- NULL
    ## Et aussi modifyList(...)
    lp$fit    <- lp$fit[[1]]
    lp$se     <- lp$se.fit[[1]]
    lp$zvalue <- lp$fit / lp$se
    lp$pvalue <- 2 * pnorm(-abs(lp$zvalue))

    p <- predict(rand$model, type = "response", se.fit = TRUE)
    p$fit    <- p$fit[[1]] #equivalent to former: pBinom <- fitted(modelBinom)[[1]]
    p$se     <- p$se.fit[[1]] ## Vérifier la pertinence de cela
    p$zvalue <- p$fit / p$se
    p$pvalue <- 2 * pnorm(-abs(p$zvalue))

    mat <- cbind(c(lp$fit, p$fit),
                 c(lp$se, p$se),
                 c(lp$zvalue, p$zvalue),
                 c(lp$pvalue, p$pvalue))
    rownames(mat) <- c("logit(p)", "p")
    colnames(mat) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    rand$par <- mat

    # z value (Pr(>|z|)
    # 2*pnorm(Estimate / Std. Error) = Pr(>|z|)
    # Pourquoi pas t-value??? (contrairement à lm??)
    #pvalue <- 2 * pt(-abs(tvalue), df.r)
    #coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
    #dimnames(coef.table) <- list(names(coef.p), c(dn, "t value", "Pr(>|t|)"))

    freq$binomial <- N * dbinom(x = 0:n, size = n, prob = p$fit)

    # Beta-binomial distribution
    inits <- c(p$fit, 1 - p$fit)
    modelBBD <- optimBetaBinom(inits, d, n) # !!!!! si retour NULL, pas d'assignement list possible
    # use stats4::mle() instead
    agg <- list()

    if (!is.null(modelBBD)) { # If the optimization procedure has been successful
        agg$model <- modelBBD
        alpha <- beta <- p <- theta <- rho <- list()
        alpha$fit <- agg$model$par[1]
        beta$fit  <- agg$model$par[2]
        #varcov <- solve(model$hessian) # var-cov matrix
        #se     <- sqrt(diag(varcov))   # standard errors
        p$fit     <- alpha$fit / (alpha$fit + beta$fit)
        theta$fit <- 1 / (alpha$fit + beta$fit)
        rho$fit <- 1 / (alpha$fit + beta$fit + 1)
        # theta: index of aggregation.
        # Aggregation increases with increasing theta

        # Work with the hessian
        estVcov  <- solve(agg$model$hessian) # estimator of the asymptotic covariance matrix
        estSe    <- sqrt(diag(estVcov))       # solve(mat): compute the inverse of the matrix mat
        alpha$se <- estSe[[1]]
        beta$se  <- estSe[[2]]
        p$se     <- msm::deltamethod(~ x1 / (x1 + x2), agg$model$par, estVcov)
        theta$se <- msm::deltamethod(~ 1 / (x1 + x2), agg$model$par, estVcov)
        rho$se   <- msm::deltamethod(~ 1 / (x1 + x2 + 1), agg$model$par, estVcov)

        p$zvalue <- p$fit / p$se
        p$pvalue <- 2 * pnorm(-abs(p$zvalue))

        theta$zvalue <- theta$fit / theta$se
        theta$pvalue <- 2 * pnorm(-abs(theta$zvalue))

        alpha$zvalue <- alpha$fit / alpha$se
        alpha$pvalue <- 2 * pnorm(-abs(alpha$zvalue))

        beta$zvalue <- beta$fit / beta$se
        beta$pvalue <- 2 * pnorm(-abs(beta$zvalue))

        rho$zvalue <- rho$fit / rho$se
        rho$pvalue <- 2 * pnorm(-abs(rho$zvalue))

        # # Below should work with bbmle::mle2
        # # Retrieve result matrice (to which we will add extra estimates)
        # param <- coef(summary(modelBBD)) # alpha and beta
        # #mylist <- list(x1 = coef(modelBBD)[[1]], # alpha
        # #               x2 = coef(modelBBD)[[2]]) # beta
        #
        # p     <- estimateCoef(modelBBD, quote(x1 / (x1 + x2)), student = FALSE)
        # theta <- estimateCoef(modelBBD, quote(1  / (x1 + x2)), student = FALSE)
        # rho   <- estimateCoef(modelBBD, quote(1  / (x1 + x2 + 1)), student = FALSE)
        #
        # param <- rbind(param, unlist(p), unlist(theta), unlist(rho))
        # rownames(param) <- c("alpha", "beta", "p", "theta", "rho")


        freq$betaBinomial <- N * dbetabinom(x = 0:n, size = n, # remplacer par dbetabinom après
                                            prob = p$fit, theta = theta$fit)

        ## cf. emdbook for 1 specific function
        ## emdbook::dbetabinom(x, prob, size,  theta, shape1, shape2, log = FALSE)
        mat <- cbind(c(NA, p$fit, theta$fit, rho$fit, alpha$fit, beta$fit),
                     c(NA, p$se, theta$se, rho$se, alpha$se, beta$se),
                     c(NA, p$zvalue, theta$zvalue, rho$zvalue, alpha$zvalue, beta$zvalue),
                     c(NA, p$pvalue, theta$pvalue, rho$pvalue, alpha$pvalue, beta$pvalue))
        rownames(mat) <- c("logit(p)", "p", "theta", "rho", "alpha", "beta")
        colnames(mat) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
        agg$par <- mat

    } else {
        agg$model <- list(NULL)
        freq$betaBinomial <- NA
        agg$par <- list(NULL)
        agg$test <- list(NULL)
    }

    # Et si optim n'a pas convergé ? Gérer ce cas ci-dessous.
    # Et puis, plein de warnings!!!
    # gof
    #testB <- with(freq, chisq.test(x = observed, p = binomial, rescale.p = TRUE))
    #testBBD <- with(freq, chisq.test(x = observed, p = betaBinomial, rescale.p = TRUE))
    # Add simulate.p.value = TRUE ???

    ####### Apply the Cochran rule regarding using Xchi2
    # Cochran : 80 % des classes doivent satisfaire la règle des cinq éléments tandis que les autres doivent être non vides.

    # Returns:

    structure(list(call = match.call(),
                   random = list(model = rand$model,
                                 par = rand$par,
                                 test = NULL),
                   aggregated = list(model = agg$model,
                                     par = agg$par,
                                     test = NULL),
                   frequency = freq,
                   test = llrtest(rand$model, agg$model)),
              class = "fitDistr") # "cmpDistr" should be better

}

#------------------------------------------------------------------------------
#' @export
summary.fitDistr <- function(object, ...) {
    cat("\nCall:\n")
    print(object$call)
    cat("\n----------\nRandom:\n")
    printCoefmat(res$random$par)
    cat("\n----------\nAggregated:\n")
    printCoefmat(res$aggregated$par)
    cat("\n----------\n")
    print(object$test)
}

#' @export
print.fitDistr <- function(x, ...) {
    cat("\nCall:\n")
    print(x$call)
    #cat("\nSource: <", x, ">\n", sep = "")
    cat("\nParameters:\n")
    cat("Random (BD distr.):\n")
    cat("  p     = ", formatC(x$random$par["p", 1]),
        " (\u00b1 ", formatC(x$random$par["p", 2]), ")\n", sep = "" )
    cat("Aggregated (BBD distr.):\n")
    cat("  p     = ", formatC(x$aggregated$par["p", 1]),
        " (\u00b1 ", formatC(x$aggregated$par["p", 2]), ")\n", sep = "")
    cat("  theta = ", formatC(x$aggregated$par["theta", 1]),
        " (\u00b1 ", formatC(x$aggregated$par["theta", 2]), ")\n", sep = "")
    cat("  rho   = ", formatC(x$aggregated$par["rho", 1]),
        " (\u00b1 ", formatC(x$aggregated$par["rho", 2]), ")\n", sep = "")
    cat("  alpha = ", formatC(x$aggregated$par["alpha", 1]),
        " (\u00b1 ", formatC(x$aggregated$par["alpha", 2]), ")\n", sep = "")
    cat("  beta  = ", formatC(x$aggregated$par["beta", 1]),
        " (\u00b1 ", formatC(x$aggregated$par["beta", 2]), ")\n", sep = "")
}

#' @export
plot.fitDistr <- function(x, y, ..., plot = TRUE) {
    g <- list(
        #geom_bar(data = (x$frequency %>% tidyr::gather("key", "value", -category)),
        #         aes(x = category, y = value, fill = key),
        #         stat = "identity", position = "dodge")
        geom_bar(data = (x$frequency %>% tidyr::gather("key", "value", -category)),
                 aes(x = category, y = value, fill = key),
                 stat = "identity", position = "dodge", colour = "black", size = 0.5),
        scale_fill_manual(breaks = c("observed", "betaBinomial", "binomial"),
                          values = c(observed = "black", betaBinomial = "grey", binomial = "white"))#,
        #                          guide = FALSE)
        #scale_fill_discrete(guide = FALSE)

        #    ggplot(data = (x$freqTable %>% tidyr::gather("key", "value", -category)),
        #           aes(x = category, y = value, fill = key)) +
        #        geom_bar(stat = "identity", position = "dodge", colour = "black", size = 0.5) +
        #        scale_fill_manual(breaks = c("observed", "betaBinomial", "binomial"),
        #                          values = c(observed = "black", betaBinomial = "grey", binomial = "white"),
        #                          guide = FALSE)
        #scale_fill_discrete(guide = FALSE)
    )
    if (!plot) g else ggplot() + g
}


#plot.fitDistr <- function(x, y, ...,
#                          type = c("all", "distribution", "parameters"),
#                          parameters = c("p-theta", "p-rho", "alpha-beta"),
#                          trend = TRUE) {}







