#------------------------------------------------------------------------------#
# Utilities
#' @keywords internal
#------------------------------------------------------------------------------#
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
    return(abs(x - round(x)) < tol)
}

#------------------------------------------------------------------------------#
#' Some link functions
#'
#' Logit, probit and cloglog functions are available.
#' The logit and the logistic (with rev = TRUE), i.e. the inverse-logit functions.
#' Probit is a wrapper around \code{qnorm} (for \eqn{probit}) and \code{pnorm} (for \eqn{probit^{-1}})
#' Complementary log-log transformation.
#'
#' @name link
#' @export
#------------------------------------------------------------------------------#
logit <- function(x, rev = FALSE) {
    if (!rev) log(x / (1 - x))
    else      1 / (1 + exp(-x))
}

#------------------------------------------------------------------------------#
#' @rdname link
#' @export
#------------------------------------------------------------------------------#
probit <- function(x, rev = FALSE) {
    if (!rev) qnorm(x)
    else      pnorm(x)
}

#------------------------------------------------------------------------------#
#' @rdname link
#' @export
#------------------------------------------------------------------------------#
cloglog <- function(x, rev = FALSE) {
    if (!rev) log(-log(1 - x))
    else      1 - exp(-exp(x))
}

#------------------------------------------------------------------------------#

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

#------------------------------------------------------------------------------#
# Log-likelihood ratio test / likelihood ratio test statistic (LRS)
#
# idea from  lmtest::lrtest()
#
# Performs a Log-likelihood ratio test
#
# retravailler cette fonction pour qu'elle exploite directemenent lmtest::lrtest()
#
# not exported currently
#------------------------------------------------------------------------------#
llrtest <- function(modRand, modAgg) {
    llModRand <- as.numeric(logLik(modRand))
    llModAgg  <- as.numeric(logLik(modAgg))
    value <- 2 * (llModAgg - llModRand)
    df <- 1 #2 * (N - 1)
    pvalue <- 1 - pchisq(value, df) # use that : "lower.tail = FALSE)" instead
    chisq <- qchisq(pvalue, df)

    tab <- matrix(rep(NA, 8), nrow = 2)
    dimnames(tab) <- list(c("random :", "aggregated :"),
                          c("LogLik", "Df", "Chisq", "Pr(>Chisq)"))
    tab[, 1]  <- c(llModRand, llModAgg)
    tab[2, 2] <- df
    tab[2, 3] <- value
    tab[2, 4] <- pvalue

    title <- "Likelihood ratio test\n"
    #topnote <- paste0("Model (random)\nModel (agg)")
    #heading = c(title, topnote)
    structure(as.data.frame(tab), heading = title,
              class = c("anova", "data.frame"))
}


#------------------------------------------------------------------------------#
# Estimate extra coef, based on deltamethod
#' @export
#------------------------------------------------------------------------------#
estimateCoef <- function(model, expr, envir = NULL, type = c("t", "norm")) { # List uniquement pour estimate mean / trouver autre nom
    # NEW !!!
    type <- match.arg(type)
    resFormula <- paste0(deparse(expr), collapse = "") # Utilise si très longue formule !!!!!, car expr sur plusieurs lignes si longue formule... pas top
    resFormula <- as.formula(paste0("~", resFormula))
    # END NEW !!!
    # resFormula <- as.formula(paste0("~", deparse(expr)))
    res        <- list()
    if (is.null(envir)) {
        coefs <- coef(model)
        envir <- setNames(as.list(coefs), paste0("x", 1:length(coefs)))
    }
    res$est    <- eval(expr, envir)
    res$se     <- msm::deltamethod(resFormula, coef(model), vcov(model))
    res$tvalue <- res$est / res$se
    res$pvalue <- ifelse((type == "t"), ## faire une condition plus intelligente (t and norm)
                         2 * pt(abs(res$tvalue), df.residual(model), lower.tail = FALSE),
                         2 * pnorm(abs(res$tvalue), lower.tail = FALSE))

    # renvoyer les bons noms : par exempel:
    # Estimate Std. Error   z value        Pr(z)
    # ou
    # Estimate Std. Error t value Pr(>|t|)
    res
}


#' @export
extraCoef <- function(model) {
    #myfunc <- function(v1) {
    #    deparse(substitute(v1))
    #}
    #
    #myfunc(foo)
    #[1] "foo"

    # Retrieve result matrice (to which we will add extra estimates)
    coefs   <- coef(stats4::summary(model)) # NON pas normal de devoir mettre bbmle ici, pas besoin de le faire dans estimateCoef
    # not need to do that forcoef(model) ... why?

    # pour betabinom uniquement
    extraCoefs = list(
        p     = quote(x1 / (x1 + x2)),
        theta = quote(1  / (x1 + x2)),
        rho   = quote(1  / (x1 + x2 + 1))
    )

    newCoefs <- do.call(rbind, lapply(extraCoefs, function(i) {
        unlist(estimateCoef(model, i, type = "norm")) # faire plus intelligent ici avec le "norm"
    }))
    rbind(coefs, newCoefs)
}

