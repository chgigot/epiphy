#------------------------------------------------------------------------------#
#' \code{epiphy}: An R package to analyze plant disease epidemics.
#'
#' This package provides a common framework for spatialized plant disease
#' intensity data collected at one or more time points. Many statistical methods
#' developed over the last decades to describe and quantify plant disease
#' epidemics are implemented (e.g., aggregation indices, Taylor and binary power
#' laws, distribution fitting, cross-product approaches). A bundle of historical
#' data sets that were mainly published in plant disease epidemiology literature
#' is also made available.
#'
#' @author
#'
#' \strong{Maintainer:} Christophe Gigot <ch.gigot@gmail.com>
#'
#' @seealso
#'
#' Useful references:
#'
#' Gosme M. 2008. Comment analyser la structure spatiale et modéliser le
#' développement spatio-temporel des épiphyties? Canadian Journal of Plant
#' Pathology, 30:4-23.
#'
#' Madden LV, Hughes G, van den Bosch F. 2007. Spatial aspects of epidemics -
#' III: Patterns of plant disease. In: The study of plant disease epidemics,
#' 235–278. American Phytopathological Society, St Paul, MN.
#'
#' @keywords internal
#'
#' @docType package
#' @name epiphy
#' @useDynLib epiphy
#' @importFrom Rcpp sourceCpp
#' @import ggplot2
#------------------------------------------------------------------------------#
NULL

# TODO: Needed for checking process...
utils::globalVariables(c("Number per sampling unit", "Frequency", "key",
                         "x", "y", "z", "i", "phs"))
