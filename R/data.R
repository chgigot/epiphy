#------------------------------------------------------------------------------#
#' Incidence of tomato spotted wilt virus (TSWV) disease in field trials.
#'
#' Intensively mapped TSWV incidence data reported by Cochran (1936) and Bald
#' (1937). The disease assessments were performed in field trials at the Waite
#' Institute (Australia) in 1928 and 1929. TSWV is a virus disease spread by
#' thrips.
#'
#' The data set \code{field_1928}, reported by Bald (1937), was a set of four
#' plots. Each plot consisted of 14 rows containing 33 plants each, so that
#' there were 462 plants in each plot. The tomato variety Early Dwarf Red was
#' used in two plots, and the variety Burwood Prize in the other two. The
#' tomatoes were planted out on 15th October 1928. The two plots dedicated to a
#' given variety experienced different irrigation practices, using either
#' overhead sprays or trenches. Otherwise, all were treated alike. Weekly
#' records of TSWV incidence were performed from 6th November to 12th December.
#'
#' The data set \code{field_1929}, reported by Cochran (1936), was a field of 24
#' rows containing 60 plants each, so that there were 1440 plants. The tomatoes
#' were planted out in 26th November 1929. TSWV incidence records made on 18th
#' December 1929, 31st December 1929 and 22nd January 1930 are reported in this
#' data set.
#'
#' @format
#'
#' There are two data frames:
#'
#' \code{field_1928}: A data frame with 11088 rows and 8 variables:
#' \tabular{rll}{
#'     [, 1]   \tab plot       \tab Plot id. \cr
#'     [, 2]   \tab variety    \tab Variety name. \cr
#'     [, 3]   \tab irrigation \tab Irrigation system. \cr
#'     [, 4:5] \tab x,y        \tab Grid spatial coordinates. \cr
#'     [, 6]   \tab t          \tab Date of disease assessments. 1: 6 Nov, 2: 14
#'                                  Nov, 3: 21 Nov, 4: 28-29 Nov, 5: 5 Dec, 6:
#'                                  12 Dec 1928. \cr
#'     [, 7]   \tab i          \tab Disease incidence. 0: Healthy, 1: Diseased. \cr
#'     [, 8]   \tab n          \tab Sampling unit size. n = 1 means that the
#'                                  sampling unit size is the plant. \cr
#' }
#'
#' \code{field_1929}: A data frame with 4320 rows and 5 variables:
#' \tabular{rll}{
#'     [, 1:2] \tab x,y \tab Grid spatial coordinates. \cr
#'     [, 3]   \tab t   \tab Date of disease assessments. 1: 18 Dec, 2: 31 Dec
#'                           1929, 3: 22 Jan 1930. \cr
#'     [, 4]   \tab i   \tab Disease incidence. 0: Healthy, 1: Diseased. \cr
#'     [, 5]   \tab n   \tab Sampling unit size. n = 1 means that the sampling
#'                           unit size is the plant. \cr
#' }
#'
#' @source
#'
#' Cochran WG. 1936. The statistical analysis of field counts of
#'     diseased plants. Supplement to the Journal of the Royal Statistical
#'     Society 3, 49–67. \href{http://dx.doi.org/10.2307/2983677}{doi:10.2307/2983677}
#'
#' Bald JG. 1937. Investigations on "spotted wilt" of tomatoes. III.
#'     Infection in field plots. Bulletin 106. Melbourne, Australia: Council for
#'     Scientific and Industrial Research.
#------------------------------------------------------------------------------#
"tomato_tswv"

#------------------------------------------------------------------------------#
#' Incidence of citrus tristeza virus (CTV) disease in three fields.
#'
#' CTV incidence data for three orchards in eastern Spain reported for
#' consecutive years.
#'
#' Both \code{IVI3and4} and \code{IVI6and7} orchards consisted of 216 trees each
#' of Washington navel orange on Troyer citrange planted in 1978 on a 2 x 6-m
#' spacing. \code{El_Realengo} orchard consisted of 400 Marsh seedless
#' grapefruit on Troyer citrange planted in 1973 on a 5.5 x 5.5-m spacing.
#'
#' @format
#'
#' There are three data frames:
#' \itemize{
#'     \item \code{IVI3and4}: A data frame with 864 rows and 5 variables.
#'     \item \code{IVI6and7}: A data frame with 648 rows and 5 variables.
#'     \item \code{El_Realengo}: A data frame with 2000 rows and 5 variables.
#' }
#'
#' The structure is the same for all the data frames:
#'
#' \tabular{rll}{
#'     [, 1:2] \tab x,y \tab Grid spatial coordinates. \cr
#     [, 3:4] \tab X,Y \tab Metric spatial coordinates. \cr
#'     [, 3]   \tab t   \tab Year of disease assessments. \cr
#'     [, 4]   \tab i   \tab Disease incidence. 0: Healthy, 1: Diseased. \cr
#'     [, 5]   \tab n   \tab Sampling unit size. n = 1 means that the sampling
#'                           unit size is the plant. \cr
#' }
#'
#' @source
#'
#' Gottwald TR, Cambra M, Moreno P, Camarasa E, Piquer J. 1996. Spatial
#'     and temporal analyses of citrus tristeza virus in eastern Spain.
#'     Phytopathology 86, 45–55.
#'
#' Gibson GJ. 1997. Investigating mechanisms of spatiotemporal epidemic
#'     spread using stochastic models. Phytopathology 87, 139-46.
#'     \href{http://dx.doi.org/10.1094/PHYTO.1997.87.2.139}{doi:10.1094/PHYTO.1997.87.2.139}
#------------------------------------------------------------------------------#
"citrus_ctv"

#------------------------------------------------------------------------------#
#' Incidence of three viruses in an Australian hop garden.
#'
#' Three viruses, i.e. Hop latent virus (HpLV), Hop mosaic virus (HpMV), and
#' Apple mosaic virus (ApMV), were monitored in an Australian hop garden for two
#' consecutive years (1996 and 1997). The hop garden was established in 1989
#' with the variety Victoria in a commercial hop farm at Bushy Park, Tasmania,
#' Australia. It consisted of 25 rows containing 51 plants each, so that there
#' were 1275 hop plants in total. There were 2.1 m between rows, and 1.8 m
#' between plants within rows.
#'
#' @format
#'
#' There are three data frames, one for each virus (\code{HpLV}, \code{HpMV} and
#' \code{ApMV}). Each data frame consists of 2550 rows and 7 variables:
#' \tabular{lll}{
#'     [, 1:2] \tab x,y   \tab Grid spatial coordonates. \cr
#'     [, 3:4] \tab xm,ym \tab Metric spatial coordinates. \cr
#'     [, 5]   \tab t     \tab Year of disease assessments. \cr
#'     [, 6]   \tab i     \tab Incidence. 0: Healthy, 1: Diseased. \cr
#'     [, 7]   \tab n     \tab Sampling unit size. n = 1 means that the sampling
#'                             unit size is the plant. \cr
#' }
#'
#' @source
#'
#' Pethybridge SJ, Madden LV. 2003. Analysis of spatiotemporal dynamics
#'     of virus spread in an Australian hop garden by stochastic modeling. Plant
#'     Disease 87:56-62.
#'     \href{http://dx.doi.org/10.1094/PDIS.2003.87.1.56}{doi:10.1094/PDIS.2003.87.1.56}
#------------------------------------------------------------------------------#
"hop_viruses"

#------------------------------------------------------------------------------#
#' Incidence of bacterial blight of onion.
#'
#' Assessments of bacterial blight of onion at two dates. The experimental plot
#' was sown with naturally X. axonopodis pv. allii-contaminated onion (A. cepa
#' L. cv. Chateau-vieux) seed lot, with a contamination rate of about 0.04\%.
#'
#' @format
#'
#' A data frame with 1134 rows and 5 variables:
#' \tabular{rll}{
#'     [, 1:2] \tab x,y \tab Grid spatial coordinates. \cr
#'     [, 3]   \tab t   \tab Date of disease assessments. \cr
#'     [, 4]   \tab i   \tab Disease incidence. 0: Healthy, 1: Diseased. \cr
#'     [, 5]   \tab n   \tab Sampling unit size. n = 1 means that the sampling
#'                           unit size is the plant. \cr
#' }
#'
#' @source
#'
#' Roumagnac P, Pruvost O, Chiroleu F, Hughes G. 2004. Spatial and
#'     temporal analyses of bacterial blight of onion caused by Xanthomonas
#'     axonopodis pv. allii. Phytopathology 94, 138–146.
#'     \href{http://dx.doi.org/10.1094/PHYTO.2004.94.2.138}{doi:10.1094/PHYTO.2004.94.2.138}
#------------------------------------------------------------------------------#
"onion_bacterial_blight"

#------------------------------------------------------------------------------#
#' Examples of simulated epidemic data.
#'
#' Epidemics were generated using the stochastic simulator from Xu and Madden
#' (2004). The data consist of the numbers of diseased plants per sampling
#' unit (out of a total of n = 100 plants in each sampling unit). N = 144
#' sampling units, and different values for the parameters \code{pattern} and
#' \code{mu} were used for the simulations.
#'
#' @format A data frame with 864 rows and 6 variables:
#' \tabular{rll}{
#'     [, 1]   \tab pattern \tab Either clumped (i.e. aggregated), random or
#'                               regular. \cr
#'     [, 2]   \tab mu      \tab Median spore dispersal parameter. \cr
#'     [, 3:4] \tab x,y     \tab Grid spatial coordinates. \cr
#'     [, 5]   \tab i       \tab Number of diseased plants (from 0 to 100). \cr
#'     [, 6]   \tab n       \tab Sampling unit size. Here, n = 100 plants per
#'                               sampling unit. \cr
#' }
#'
#' @source
#'
#' Xu XM, Madden LV. 2004. Use of SADIE statistics to study spatial
#'   dynamics of plant disease epidemics. Plant Pathology 53, 38–49.
#'   \href{http://dx.doi.org/10.1111/j.1365-3059.2004.00949.x}{doi:10.1111/j.1365-3059.2004.00949.x}
#------------------------------------------------------------------------------#
"simulated_epidemics"

#------------------------------------------------------------------------------#
#' Incidence of tobacco plants infected with viruses.
#'
#' Experimental plot consisted of 75 sampling units with 40 tobacco plants in
#' each one.
#'
#' @format A data frame with 75 rows and 2 variables:
#' \tabular{rll}{
#'     [, 1] \tab i \tab Number of diseased plants (from 0 to 40). \cr
#'     [, 2] \tab n \tab Sampling unit size. Here, n = 40 plants per sampling
#'                       unit. \cr
#' }
#'
#' @source
#'
#' Madden LV, Pirone TP, Raccah B. 1987. Analysis of spatial patterns of
#'     virus-diseased tobacco plants. Phytopathology 77, 1409–1417.
#------------------------------------------------------------------------------#
"tobacco_viruses"

#------------------------------------------------------------------------------#
#' Incidence of dogwood anthracnose.
#'
#' Incidence data from the Dogwood Anthracnose Impact Assessment Program for
#' 1990 and 1991, in the Southeast of the USA, reported by Zarnoch et al (1995).
#' Only plots with exactly n = 10 dogwood trees are present in the data set (168
#' and 161 plots in 1990 and 1991, respectively).
#'
#' @format A data frame with 329 rows and 3 variables:
#' \tabular{rll}{
#'     [, 1] \tab t \tab Year of disease assessments (1990 or 1991).. \cr
#'     [, 2] \tab i \tab Number of diseased plants (from 0 to 10). \cr
#'     [, 3] \tab n \tab Sampling unit size. Here, n = 10 plants per sampling
#'                       unit (or plot). \cr
#' }
#'
#' @source
#'
#' Zarnoch SJ, Anderson RL, Sheffield RM. 1995. Using the \eqn{\beta}-binomial
#'     distribution to characterize forest health. Canadian journal of forest
#'     research 25, 462–469.
#------------------------------------------------------------------------------#
"dogwood_anthracnose"

#------------------------------------------------------------------------------#
#' Incidence of ray blight disease of pyrethrum.
#'
#' An assessment of the incidence of ray blight disease of pyrethrum in 62
#' sampling units, containing 6 plants each.
#'
#' @format A data frame with 62 rows and 2 variables:
#' \tabular{rll}{
#'     [, 1] \tab i \tab Number of diseased plants (from 0 to 6). \cr
#'     [, 2] \tab n \tab Sampling unit size. Here, n = 6 plants per sampling
#'                       unit. \cr
#' }
#'
#' @source
#'
#' Pethybridge SJ, Esker P, Hay F, Wilson C, Nutter FW. 2005.
#'     Spatiotemporal description of epidemics caused by Phoma ligulicola in
#'     Tasmanian pyrethrum fields. \emph{Phytopathology} 95, 648–658.
#'     \href{http://dx.doi.org/10.1094/PHYTO-95-0648}{doi:10.1094/PHYTO-95-0648}
#------------------------------------------------------------------------------#
"pyrethrum_ray_blight"

#------------------------------------------------------------------------------#
#' Counts of aphids.
#'
#' Counts of 554 aphids of the species Sitobion avenae, sampled on 28 June 1996
#' in a 250 x 180-m field of winter wheat near Wimborne, Dorset, UK. The 63
#' sampling units, made of the inspection of five tillers each, were located on
#' a 9 x 7 rectangular grid at intervals of 30 m.
#'
#' @format A data frame with 63 rows and 3 variables:
#' \tabular{lll}{
#'     [, 1:2] \tab x,y   \tab Grid spatial coordinates. \cr
#'     [, 3:4] \tab xm,ym \tab Metric spatial coordonates. \cr
#'     [, 5]   \tab i     \tab Counts of aphids. \cr
#' }
#'
#' @source
#'
#' Perry JN, Winder L, Holland JM, Alston RD. 1999. Red-blue plots for
#'     detecting clusters in count data. Ecology Letters 2, 106-13.
#'     \href{http://dx.doi.org/10.1046/j.1461-0248.1999.22057.x}{doi:10.1046/j.1461-0248.1999.22057.x}
#------------------------------------------------------------------------------#
"aphids"

#------------------------------------------------------------------------------#
#' Counts of arthropods.
#'
#' A sampling unit was made of a pitfall to collect arthropods in a field of
#' organic winter wheat, near Wimborne, Dorset, UK in 1996. The sampling units
#' were located on a 9 x 7 rectangular grid at intervals of 30 m. There were
#' six sampling dates.
#'
#' @format A data frame with 378 rows and 4 variables:
#' \tabular{lll}{
#'     [, 1:2] \tab x,y   \tab Grid spatial coordinates. \cr
#'     [, 3:4] \tab xm,ym \tab Metric spatial coordonates. \cr
#'     [, 5]   \tab t     \tab Sampling date. 1: 7 Jun, 2: 14 Jun, 3: 21 Jun, 4:
#'                             28 Jun, 5: 5 Jul, 6: 12 Jul 1996. \cr
#'     [, 6]   \tab i     \tab Counts of arthropods. \cr
#' }
#'
#' @source
#'
#' Holland JM, Winder L, Perry JN. 1999. Arthropod prey of farmland
#'     birds: Their spatial distribution within a sprayed field with and without
#'     buffer zones. Aspects of Applied Biology 54: 53-60.
#------------------------------------------------------------------------------#
"arthropods"

#------------------------------------------------------------------------------#
#' Count of codling moth larvae.
#'
#' Codling moth diapausing larvae were collected in an apple orchard in
#' south-eastern France. Larvae were caught on strip traps wrapped around tree
#' trunks in July 2008 and collected the following October. 30 traps were used.
#'
#' @format A data frame with 30 rows and 3 variables:
#' \tabular{lll}{
#'     [, 1:2] \tab x,y \tab Metric spatial coordonates. \cr
#'     [, 3]   \tab i   \tab Counts of codling moth larvae. \cr
#' }
#'
#' @source
#'
#' Lavigne C, Ricci B, Franck P, Senoussi R. 2010. Spatial analyses of
#'     ecological count data: A density map comparison approach. Basic and
#'     Applied Ecology 11: 734-42.
#'     \href{http://dx.doi.org/10.1016/j.baae.2010.08.011}{doi:10.1016/j.baae.2010.08.011}
#------------------------------------------------------------------------------#
"codling_moths"

#------------------------------------------------------------------------------#
#' Offspring survival of rats experiencing different diets.
#'
#' Results of an experiment where two groups of 16 female rats were fed
#' different diets during pregnancy and lactation periods. One group's diet
#' contained a chemical under review, and the other one was a control. For each
#' litter, the number of pups alive at 4 days, and the number of pups weaned
#' (i.e. that survived the 21-day lactation period) were recorded.
#'
#' @format A data frame with 32 rows and 3 variables:
#' \tabular{rll}{
#'     [, 1]   \tab group \tab Either control or treated group. \cr
#'     [, 2]   \tab i     \tab Pups weaned. \cr
#'     [, 3]   \tab n     \tab Pups alive at 4 days. \cr
#' }
#'
#' @source
#'
#' Weil CS. 1970. Selection of the valid number of sampling units and a
#'     consideration of their combination in toxicological studies involving
#'     reproduction, teratogenesis or carcinogenesis. Food and Cosmetics
#'     Toxicology 8: 177-182.
#------------------------------------------------------------------------------#
"offspring_survival"




