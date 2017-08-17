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
#' @param index The name of the index to use: "perry", "li-madden-xu" or "all".
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
#' @name sadie
#' @export
#------------------------------------------------------------------------------#
sadie <- function(data, index = c("perry", "LMX", "all"), nperm = 100,
                  rseed = TRUE, seed = 12345, cost) {
    UseMethod("sadie")
}

#------------------------------------------------------------------------------#
#' @rdname sadie
#' @method sadie data.frame
#' @export
#------------------------------------------------------------------------------#
sadie.data.frame <- function(data, index, nperm, rseed, seed, cost) {
    return(NULL)
}



