#------------------------------------------------------------------------------#
# Useful packages
#------------------------------------------------------------------------------#
library(magrittr)
library(tidyverse)
library(devtools)
#options(stringsAsFactors = FALSE) # To use "standard" data frame, do noy use that.

#------------------------------------------------------------------------------#
# Useful functions
#------------------------------------------------------------------------------#
wide2long <- function(data) {
    data$x <- as.integer(row.names(data))
    res    <- data %>% gather(y, value, -x)
    # substring(...) because R gives 'V#' as default column names
    res$y  <- as.integer(substring(res$y, 2))
    res
}

listOfValuesSup0 <- function(data, vals)
    lapply(vals, function(val) filter(data, value == val))

generateGrid <- function(data, times) {
    res <- expand.grid(x = unique(data$x),
                       y = unique(data$y),
                       t = times,
                       KEEP.OUT.ATTRS = FALSE)
    res   <- res %>% arrange(x, y, t)
    res$r <- 0
    res$n <- 0
    res
}

fillDis <- function(data, listVals) {
    lapply(seq_len(length(listVals)), function(id) {
        # ida: row-ids in data where (x, y) are also present in listVals[[id]]
        ida <- which(
            apply(select(data, x, y), 1, paste0, collapse = "_") %in%
                apply(select(listVals[[id]], x, y), 1, paste0, collapse = "_"))
        # idb: row-ids in data where t >= id
        idb <- which(data$t >= id)
        data[intersect(ida, idb), ]$r <<- 1 # Here, 1 means diseased (0: healthy)
    })
    data
}

#------------------------------------------------------------------------------#
# Dataset: tomato_tswv (ex. dataCochran1936)
#------------------------------------------------------------------------------#
# Plants recorded as diseased at time:
# - 1: 18 December 1929
# - 2: 31 December 1929
tomato_tswv_1929_1plot <- read.csv("data-raw/tomato_tswv_1929_1plot.csv", header = TRUE)
#write.csv(tomato_tswv_1928_4plots, file = "data-raw/tomato-tswv-1928-4plots.csv", quote = FALSE, row.names = FALSE)

#------------------------------------------------------------------------------#
### Export tomato_tswv_4plots (ex. dataBald197)
# 0 = Healthy plants
# Below: diseased plants
# 1 = 6 November 1928
# 2 = 14 November 1928
# 3 = 21 November 1928
# 4 = 28-29 November 1928
# 5 = 5 December 1928
# 6 = 12 December 1928
tomato_tswv_1928_4plots <- read.csv("data-raw/tomato_tswv_1928_4plots.csv", header = TRUE)

tomato_tswv <- list("1928_4plots" = tomato_tswv_1928_4plots,
                    "1929_1plot"  = tomato_tswv_1929_1plot)
use_data(tomato_tswv)



#------------------------------------------------------------------------------#
# Dataset: citrus_ctv (ex. dataGottwald1996)
#------------------------------------------------------------------------------#
citrus_ctv <- list("IVI3&4"      = read.csv("data-raw/citrus_ctv_IVIA3and4.csv", header = TRUE),
                   "IVI6&7"      = read.csv("data-raw/citrus_ctv_IVIA6and7.csv", header = TRUE),
                   "El-Realengo" = read.csv("data-raw/citrus_ctv_El-Realengo.csv", header = TRUE))
use_data(citrus_ctv)

#------------------------------------------------------------------------------#
### Export tobacco_viruses (ex. dataMadden1987)
tobacco_viruses <- read.csv("data-raw/tobacco_viruses.csv", header = TRUE)
use_data(tobacco_viruses)

#------------------------------------------------------------------------------#
### Export pyrethrum_ray_blight (ex. dataPethybridge2005)
pyrethrum_ray_blight <- read.csv("data-raw/pyrethrum_ray_blight.csv", header = TRUE)
use_data(pyrethrum_ray_blight)

#------------------------------------------------------------------------------#
### Export onion_bacterial_blight (ex. dataRoumagnac2004)
# 1 = new diseased at first date
# 2 = new diseased at second data
onion_bacterial_blight <- read.csv("data-raw/onion_bacterial_blight.csv", header = TRUE)
use_data(onion_bacterial_blight)

#------------------------------------------------------------------------------#
### Export dataSkellam1948
#dataSkellam1948 <- read.csv("data-raw/Skellam1948.csv", header = TRUE)
#use_data(dataSkellam1948)

#------------------------------------------------------------------------------#
### Export offspring_survival (ex. dataWilliams1975)
# Interesting : n not the same size everywhere
offspring_survival <- read.csv("data-raw/offspring_survival.csv", header = TRUE)
use_data(offspring_survival)

#------------------------------------------------------------------------------#
### Export simulated_epidemics (ex. dataXuMadden2004)
# pattern = clumped, regular or random
# mu =  median spore dispersal parameter
simulated_epidemics <- read.csv("data-raw/simulated_epidemics.csv", header = TRUE)
use_data(simulated_epidemics)

#------------------------------------------------------------------------------#
### Export dogwood_anthracnose (ex. dataZarnoch1995)
dogwood_anthracnose <- read.csv("data-raw/dogwood_anthracnose.csv", header = TRUE)
use_data(dogwood_anthracnose)









#------------------------------------------------------------------------------#
# Export aphid counts (Perry et al., 1999)
#------------------------------------------------------------------------------#
aphids <- read.csv("data-raw/aphids.csv", header = TRUE)

with(aphids, plot(x, y, pch = NA))
with(aphids, text(x, y, labels = r))

devtools::use_data(aphids)

#------------------------------------------------------------------------------#
# Export arthropods counts (Holland et al., 1999)
#------------------------------------------------------------------------------#
arthropods <- read.csv("data-raw/arthropods.csv", header = TRUE)
arthropods_splitted <- split(arthropods, arthropods$set)

opar <- par()
par(mar = c(2,2,2,2)) # We do not see x and y labels anymore.
layout(matrix(1:6, nrow = 3, ncol = 2, byrow = TRUE))
invisible(lapply(arthropods_splitted, function(set) {
    with(set, plot(x, y, pch = NA))
    with(set, text(x, y, labels = r))
}))
par(opar)
layout(1)

devtools::use_data(arthropods)

#------------------------------------------------------------------------------#
# Export codling moth counts (Lavigne et al., 2010)
#------------------------------------------------------------------------------#
codling_moths <- read.csv("data-raw/codling_moths.csv", header = TRUE)

coef <- 80 # Coef to convert arbitrary units into meters (m = coef * u.a.)
codling_moths[, 1:2] %<>% multiply_by(coef)

with(codling_moths, plot(x, y, pch = NA))
with(codling_moths, text(x, y, labels = r))

devtools::use_data(codling_moths)


#------------------------------------------------------------------------------#

# reorg <- function(my_data, year) {
#
#     my_data$x  <- round((my_data$x * (50.4 / 50)) / 2.1) + 1
#     my_data$y  <- round((my_data$y * (90   / 90)) / 1.8) + 1
#     my_data$r  <- 1
#
#     tmp <- expand.grid(x = min(my_data$x):max(my_data$x),
#                        y = min(my_data$y):max(my_data$y))
#
#     my_data <- left_join(tmp, my_data)
#     my_data[is.na(my_data$r), ]$r <- 0
#     my_data$n  <- 1
#     my_data$xm <- (my_data$x - 1) * 2.1
#     my_data$ym <- (my_data$y - 1) * 1.8
#     my_data$t  <- year
#
#     my_data <- my_data %>%
#         select(x, y, xm, ym, t, r, n) %>%
#         arrange(x, y, t)
#
#     my_data
# }
#
# year_1996 <- c("data-raw/Pethybridge_apple_A_1996.csv",
#                "data-raw/Pethybridge_apple_B_1996.csv",
#                "data-raw/Pethybridge_apple_C_1996.csv")
# year_1997 <- c("data-raw/Pethybridge_apple_A_1997_subsequent.csv",
#                "data-raw/Pethybridge_apple_B_1997_subsequent.csv",
#                "data-raw/Pethybridge_apple_C_1997_subsequent.csv")
# viruses   <- c("HpLV", "HpMV", "ApMV")
#
# for (i in 1:3) {
#     my_data_1 <- read.csv2(year_1996[i], header = TRUE)
#     my_data_2 <- read.csv2(year_1997[i], header = TRUE)
#     my_data_2 <- bind_rows(my_data_1, my_data_2)
#     my_data_3 <- bind_rows(reorg(my_data_1, 1996),
#                            reorg(my_data_2, 1997))
#
#     #with(my_data_3 %>% filter(r == 1, t == 1997),
#     #     plot(xm, ym, pch = 19, col = rgb(0,0,0,1), yaxt = "n"))
#     #axis(2, at = seq(0, 90, by = 10), las = 2)
#
#     write.csv(my_data_3, file = paste0("data-raw/hop_", i, "_", viruses[i], ".csv"),
#               quote = FALSE, row.names = FALSE)
# }

hop_viruses <- list("HpLV" = read.csv("data-raw/hop_1_HpLV.csv", header = TRUE),
                    "HpMV" = read.csv("data-raw/hop_2_HpMV.csv", header = TRUE),
                    "ApMV" = read.csv("data-raw/hop_3_ApMV.csv", header = TRUE))
use_data(hop_viruses)

#------------------------------------------------------------------------------#

