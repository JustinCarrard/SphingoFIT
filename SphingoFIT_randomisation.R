#==========================================================================================
# Randomization for SpingoFIT
# Author: Denis Infanger
#==========================================================================================
#------------------------------------------------------------------------------------------
# Load packages
#------------------------------------------------------------------------------------------

library(ggplot2)
library(randomizeR)

#------------------------------------------------------------------------------------------
# Set up parameters for randomization
#------------------------------------------------------------------------------------------

# Total number of participans

N <- 100

# Number of groups

K <- 2

# Name of the groups

groups <- c("HIIT", "noHIIT")

# Number of participants in each stratum

N_per_group <- N/4

#------------------------------------------------------------------------------------------
# Randomized permuted block randomization with block sizes 2, 4, 6
# Ranodmization is stratified by age (40-49, 50-60) and sex (male, female)
#------------------------------------------------------------------------------------------

# Block sizes

blk_size <- c(2, 4, 6)

# Define randomization procedure for each stratum

params_male_40_49 <- rpbrPar(N_per_group, rb = blk_size, K = K, groups = groups)
params_male_50_60 <- rpbrPar(N_per_group, rb = blk_size, K = K, groups = groups)
params_female_40_49 <- rpbrPar(N_per_group, rb = blk_size, K = K, groups = groups)
params_female_50_60 <- rpbrPar(N_per_group, rb = blk_size, K = K, groups = groups)

# Generate sequences

mySeq_male_40_49 <- genSeq(params_male_40_49, seed = 919871)
mySeq_male_50_60 <- genSeq(params_male_50_60, seed = 654983511)
mySeq_female_40_49 <- genSeq(params_female_40_49, seed = 5698587)
mySeq_female_50_60 <- genSeq(params_female_50_60, seed = 354923831)

# Get randomized lists

list_male_40_49 <- getRandList(mySeq_male_40_49)
list_male_50_60 <- getRandList(mySeq_male_50_60)
list_female_40_49 <- getRandList(mySeq_female_40_49)
list_female_50_60 <- getRandList(mySeq_female_50_60)

# Combine sublists into the final list

final_list <- data.frame(
  group = c(list_male_40_49, list_male_50_60, list_female_40_49, list_female_50_60)
  , sex = rep(c("m", "f"), each = N/2)
  , ages = rep(c("40-49", "50-60"), each = N/4, times = K)
)

