#------------------------------------------------------------------------------#
#' Incidence of TSWV in tomato fields (2 varieties + 2 irrigation systems) and Incidence of TSWV recorded in a tomato field over time.
#'
#' Records1: Tomato plants with symptoms of Tomato spotted wilt virus disease. The assessments of
#' disease incidence were made in 1928, and reported by Bald (1937). Four rectangular plots,
#' each containing ??? plants (?? rows and ?? columns). Combinations of two cultivars (Burwood Prize and Early Dwarf
#' Red) and two irrigation systems (overhead sprays and trenches) were studied.
#'
#' Two varieties were used: Early Dwarf Red and Burwood Prize. There were set in the field on 15th October, 1928.
#' There were two blocks of each variety. After the seedlings were rooted, there were watered by overhead sprays, or by trenches? Otherwise
#' all were treated alike. The plants were pruned and staked, and weekly infection records were taken from 6th November to 12th December.
#'
#' Note that this plot of tomatoes was subjectd to the severest early-season epidemic of spotted wilt recorded at the Waite Intityte, and the severity
#' of the disease was correlated witht the invasion of thrips in exceptional numbers.
#'
#' Field plan showing weekly records of infection in four plots of tomatoes.
#'
#'
#'
#' Regarding the time, here are the meanings:
#' 0 = Healthy plants
#' Below: diseased plants
#' 1 = 6 November 1928
#' 2 = 14 November 1928
#' 3 = 21 November 1928
#' 4 = 28-29 November 1928
#' 5 = 5 December 1928
#' 6 = 12 December 1928
#'
#'
#' @section Records2:
#'
#' Intensively mapped disease incidence data reported by Cochran (1936). This virus disease is carried by a species of thrips. Data were obtained from field trislsmade
#' made at the Waite Institute, Australia. It's an examination by Bald of the spread of this disease.
#' There were 1440 plants assessed in the field. The tomatoes were planted out in November 26, 1929.
#' Counts of diseased plants at the first count (December 18), the second count (December 31) and the
#' last count (January 22) are recorded.
#'
#' Epidemiological data reported by Cochran (1936). A field map of Tomato spotted wilt virus disease incidence.
#' Sowing date was 26 November 1929. Intensively mapped disease data.
#'
#'
#' @format records1: A data frame with 11088 rows and 7 variables:
#' \tabular{rll}{
#'     [, 1]   \tab variety    \tab \cr
#'     [, 2]   \tab irrigation \tab \cr
#'     [, 3:4] \tab x,y \tab Spatial coordinates. \cr
#'     [, 5]   \tab t   \tab Time of the disease assessments. \cr
#'     [, 6]   \tab d   \tab Disease incidence. 1: Diseased, 0: Healthy. \cr
#'     [, 7]   \tab n   \tab Sampling unit size. n = 1 means that the sampling unit size is the plant. \cr
#' }
#'
#' **records2: A data frame with 4320 rows and 5 variables:
#' \tabular{rll}{
#'     [, 1:2] \tab x,y \tab Spatial coordinates. \cr
#'     [, 3]   \tab t   \tab Time of the disease assessments. 1: 18 December 1929, 2: 31 December 1929, 3: 22 January 1930. \cr
#'     [, 4]   \tab d   \tab Disease incidence. 1: Diseased, 0: Healthy. \cr
#'     [, 5]   \tab n   \tab Sampling unit size. n = 1 means that the sampling unit size is the plant. \cr
#' }
#'
#' @source Cochran WG. 1936. The statistical analysis of field counts of
#'     diseased plants. Supplement to the Journal of the Royal Statistical
#'     Society 3, 49–67. \href{http://dx.doi.org/10.2307/2983677}{doi:10.2307/2983677}
#'
#' @source Bald JG. 1937. Investigations on "spotted wilt" of tomatoes. III.
#'   Infection in field plots. Bulletin 106. Melbourne, Australia: Council for
#'   Scientific and Industrial Research.
#------------------------------------------------------------------------------#
"tomato_tswv"

#------------------------------------------------------------------------------#
#' CTV disease incidence at 3 assessment dates.
#'
#' Map of Citrus tristeza virus disease incidence at three assessment years.
#' 21 x 21 grid???
#'
#' Epidemic of citrus tristeza virus (CTV) observed in eastern Spain (epidemic IVIA3&4 [13]). Black squares denote positions of infections observed at
#' previous observation times, and white squares represent the new infections recorded in each disease map. From left to right, mappings represent infections dis-
#' covered up to May 1991 (45 infections), May 1992 (77 infections), and May 1993 (127 infections), respectively.
#'
#' A, Epidemic of citrus tristeza virus (CTV) observed in eastern Spain (epidemic IVIA3&4 [13]). Black squares denote positions of infections observed at
#' previous observation times, and white squares represent the new infections recorded in each disease map. From left to right, mappings represent infections dis-
#' covered up to May 1991 (45 infections), May 1992 (77 infections), and May 1993 (127 infections), respectively. B, Corresponding information for epidemic
#' IVIA6&7 (13). From left to right, mappings represent infections discovered up to May 1990 (60 infections), May 1991 (77 infections), and May 1992 (127
#' infections), respectively.
#'
#' Plots IVIA3&4 and IVIA6&7, consisted of 216 trees each of Washington navel orange on Troyer citrange planted in 1978
#' on a 2 x 6-m spacing. Both plots were located on the grounds of the Instituto Valenciano de Investigaciones
#' Agrarias, near Moncada, Valencia, Spain and were assayed yearly from 1981 to 1994.
#'
#' @format There 3 data sets here:
#'
#' Field IVIA3&4: A data frame with ... rows and 7 variables:
#' \tabular{rll}{
#'     [, 1:2] \tab x,y \tab Spatial coordinates. \cr
#'     [, 3:4] \tab X,Y \tab Metric spatial coordinates. \cr
#'     [, 5]   \tab t   \tab Time/Year of the disease assessments (1991, 1992, 1993 or 1994). \cr
#'     [, 6]   \tab r   \tab Disease incidence records. 1: Diseased, 0: Healthy. \cr
#'     [, 7]   \tab n   \tab Sampling unit size. n = 1 means that the sampling unit size is the plant. \cr
#' }
#'
#' Field IVIA6&7: A data frame with ... rows and 7 variables:
#' \tabular{rll}{
#'     [, 1:2] \tab x,y \tab Spatial coordinates. \cr
#'     [, 3:4] \tab X,Y \tab Metric spatial coordinates. \cr
#'     [, 5]   \tab t   \tab Time/Year of the disease assessments (1990, 1991, or 1992). \cr
#'     [, 6]   \tab r   \tab Disease incidence records. 1: Diseased, 0: Healthy. \cr
#'     [, 7]   \tab n   \tab Sampling unit size. n = 1 means that the sampling unit size is the plant. \cr
#' }
#'
#' Field El-Realengo: A data frame with ... rows and 7 variables:
#' \tabular{rll}{
#'     [, 1:2] \tab x,y \tab Spatial coordinates. \cr
#'     [, 3:4] \tab X,Y \tab Metric spatial coordinates. \cr
#'     [, 5]   \tab t   \tab Time/Year of the disease assessments (1981, 1982, 1984, 1985, 1990). \cr
#'     [, 6]   \tab r   \tab Disease incidence records. 1: Diseased, 0: Healthy. \cr
#'     [, 7]   \tab n   \tab Sampling unit size. n = 1 means that the sampling unit size is the plant. \cr
#' }
#'
#' @source Gottwald TR, Cambra M, Moreno P, Camarasa E, Piquer J. 1996. Spatial
#'     and temporal analyses of citrus tristeza virus in eastern Spain.
#'     Phytopathology 86, 45–55.
#'
#' @source Gibson GJ. 1997. Investigating mechanisms of spatiotemporal epidemic
#'     spread using stochastic models. Phytopathology 87, 139-46.
#'     \href{http://dx.doi.org/10.1094/PHYTO.1997.87.2.139}{doi:10.1094/PHYTO.1997.87.2.139}
#------------------------------------------------------------------------------#
"citrus_ctv"

#------------------------------------------------------------------------------#
#' Bacterial blight of onion.
#'
#' Map of bacterial blight of onion incidence at two assessment dates.
#' 21 x 21 grid???
#'
#' Experimental plots were sown with the same naturally X. axonopodis pv. allii-contaminated onion
#' (A. cepa L. cv. Chateau-vieux) seed lot. The contamination rate of seed was about 0.04%.
#'
#' @format A data frame with 1134 rows and 5 variables (x, y, t, d, n) d- record
#' range(x) = 1:27
#' range(y) = 1:21
#'
#' @source Roumagnac P, Pruvost O, Chiroleu F, Hughes G. 2004. Spatial and
#'     temporal analyses of bacterial blight of onion caused by Xanthomonas
#'     axonopodis pv. allii. Phytopathology 94, 138–146.
#'     \href{http://dx.doi.org/10.1094/PHYTO.2004.94.2.138}{doi:10.1094/PHYTO.2004.94.2.138}
#------------------------------------------------------------------------------#
"onion_bacterial_blight"

#------------------------------------------------------------------------------#
#' Examples of simulated epidemic data.
#'
#' Epidemics generated using the stochastic simulator of Xu and Madden (2004).
#' The data are number of diseased plants per quadrat (out of a total of n = 100
#' plant in each). N = 144 quadrats. Differents parameters (pattern and mu).
#'
#' @format A data frame with 864 rows and 6 variables:
#' \tabular{rll}{
#'     [, 1]   \tab pattern \tab Either clumped (i.e. aggregated), random or regular. \cr
#'     [, 2]   \tab mu \tab Median spore dispersal parameter. \cr
#'     [, 3:4] \tab x,y \tab Spatial coordinates. \cr
#'     [, 5]   \tab d   \tab Diseased individuals (from 0 to 100). \cr
#'     [, 6]   \tab n   \tab Sampling unit size. Here, n = 100 individuals per sampling unit. \cr
#' }
#'
#' @source Xu XM, Madden LV. 2004. Use of SADIE statistics to study spatial
#'   dynamics of plant disease epidemics. Plant Pathology 53, 38–49.
#'   \href{http://dx.doi.org/10.1111/j.1365-3059.2004.00949.x}{doi:10.1111/j.1365-3059.2004.00949.x}
#------------------------------------------------------------------------------#
"simulated_epidemics"

#------------------------------------------------------------------------------#
#' Incidence of tobacco plants infected with viruses
#'
#' Field with N = 75 sampling units with n = 40 plants in each one.
#'
#' @format A data frame with 75 rows and 2 variables:
#' \tabular{rll}{
#'     [, 1] \tab d \tab Disease incidence. 1: Diseased, 0: Healthy. \cr
#'     [, 2] \tab n \tab Sampling unit size. n = 1 means that the sampling unit size is the plant. \cr
#' }
#'
#' @source Madden LV, Pirone TP, Raccah B. 1987. Analysis of spatial patterns of
#'     virus-diseased tobacco plants. Phytopathology 77, 1409–1417.
#------------------------------------------------------------------------------#
"tobacco_viruses"

#------------------------------------------------------------------------------#
#' Dogwood anthracnose data.
#'
#' Dogwood trees with symptoms of dogwood anthracnose. The assessments of disease
#' incidence were based on plots containing 10 dogwood trees and reported by Zarnoch et al (1995).
#'
#' @source Zarnoch SJ, Anderson RL, Sheffield RM. 1995. Using the β-binomial
#'     distribution to characterize forest health. Canadian journal of forest
#'     research 25, 462–469.
#------------------------------------------------------------------------------#
"dogwood_anthracnose"

#------------------------------------------------------------------------------#
#' Ray blight disease of pyrethrum.
#'
#' An assessment of the incidence of foliar symptoms due to ray blight disease
#' of pyrethrum in 62 quadrats.
#'
#' @format A data frame with 62 rows and 2 variables:
#' \tabular{rll}{
#'     [, 1] \tab d \tab Number of diseased samples in the quadrat (i.e., 0, 1, 2, 3, 4, 5, or 6). \cr
#'     [, 2] \tab n \tab Sampling unit size. Number of samples assessed per quadrat, n = 6. \cr
#' }
#'
#' @source Pethybridge SJ, Esker P, Hay F, Wilson C, Nutter FW. 2005.
#'     Spatiotemporal description of epidemics caused by Phoma ligulicola in
#'     Tasmanian pyrethrum fields. \emph{Phytopathology} 95, 648–658.
#'     \href{http://dx.doi.org/10.1094/PHYTO-95-0648}{doi:10.1094/PHYTO-95-0648}
#------------------------------------------------------------------------------#
"pyrethrum_ray_blight"

#------------------------------------------------------------------------------#
#' Offspring survival of rats experiencing different diets
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
#'     [, 2]   \tab d   \tab Pups weaned. \cr
#'     [, 3]   \tab n   \tab Pups alive at 4 days. \cr
#' }
#'
#' @source Weil CS. 1970. Selection of the valid number of sampling units and a
#'     consideration of their combination in toxicological studies involving
#'     reproduction, teratogenesis or carcinogenesis. Food and Cosmetics
#'     Toxicology 8: 177-182.
#------------------------------------------------------------------------------#
"offspring_survival"

#------------------------------------------------------------------------------#
#' Aphid counts
#'
#' TODO.
#'
#' @format A data frame with 63 rows and 3 variables:
#' \tabular{lll}{
#'     [, 1:2] \tab x,y \tab Spatial coordonates. \cr
#'     [, 3]   \tab d   \tab Counts. \cr
#' }
#'
#' @source Perry JN, Winder L, Holland JM, Alston RD. 1999. Red-blue plots for
#'     detecting clusters in count data. Ecology Letters 2, 106-13.
#'     \href{http://dx.doi.org/10.1046/j.1461-0248.1999.22057.x}{10.1046/j.1461-0248.1999.22057.x}
#------------------------------------------------------------------------------#
"aphids"

#------------------------------------------------------------------------------#
#' Arthropods counts
#'
#' TODO.
#'
#' @format A data frame with 378 rows and 4 variables:
#' \tabular{lll}{
#'     [, 1]   \tab set \tab Grid id. \cr
#'     [, 2:3] \tab x,y \tab Spatial coordonates. \cr
#'     [, 4]   \tab d   \tab Counts. \cr
#' }
#'
#' @source Holland JM, Winder L, Perry JN. 1999. Arthropod prey of farmland
#'     birds: Their spatial distribution within a sprayed field with and without
#'     buffer zones. Aspects of Applied Biology 54: 53-60.
#------------------------------------------------------------------------------#
"arthropods"

#------------------------------------------------------------------------------#
#' Codling moth counts
#'
#' TODO.
#'
#' @format A data frame with 30 rows and 3 variables:
#' \tabular{lll}{
#'     [, 1:2] \tab x,y \tab Spatial coordonates. \cr
#'     [, 3]   \tab d   \tab Counts. \cr
#' }
#'
#' @source Lavigne C, Ricci B, Franck P, Senoussi R. 2010. Spatial analyses of
#'     ecological count data: A density map comparison approach. Basic and
#'     Applied Ecology 11: 734-42.
#'     \href{http://dx.doi.org/10.1016/j.baae.2010.08.011}{10.1016/j.baae.2010.08.011}
#------------------------------------------------------------------------------#
"codling_moths"





