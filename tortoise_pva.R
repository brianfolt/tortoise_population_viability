################################################################################
##### Population viability analysis for gopher tortoises
##### What demographic rates are necessary to yield population viability?
##### B Folt, February 2024

## Set up workspace
# Clear memory
rm(list = ls())

## Load packages
# devtools::install_github("BruceKendall/mpmtools")
library(mpmtools)
library(popbio)


##### Population model -----

## Create a matrix-population model for a female-only population

# What is the approximate age of sexual maturity?
ma <- 18

# Juvenile ages
juv <- 1:(ma-1)

# What is the productivity of an adult female?
# A product of breeding probability (bp), fecundity (f), nest survival (s_n),
#   egg viability (ve), probability female (pf), and hatchling survival (s_h).
bp <- 0.97
f <- 8
s_n <- 0.35
ve <- 0.85
pf <- 0.5
s_h <- 0.13

# What is the survival rate of juvenile and adult females?
s_j <- 0.75
s_a <- 0.96


# Create a demography schedule, with juvenile and mature age classes
# The model has a pre-breeding census
demog_sched <- data.frame(x = 1:ma,
                          sx = rep(NA, length(1:ma)),
                          mx = rep(NA, length(1:ma)))

# Specify productivity of juvenile and adult females
demog_sched[1:(ma-1), "mx"] <- 0
demog_sched[ma, "mx"] <- bp * f * s_n * ve * pf *s_h

# Specify juvenile and adult survival rates
demog_sched[1:(ma-1), "sx"] <- s_j
demog_sched[ma, "sx"] <- s_a

# Construct a Leslie matrix from this demography schedule
A1 <- make_Leslie_matrix(demog_sched)

# Calculate the asymptotic growth rate of the population governed by this 
#   demography schedule:
lambda1(A1)



### Specify alternative conditions
bp <- 0.97
f <- 8
s_n <- 0.35
ve <- 0.85
pf <- 0.5
s_h <- 0.25 #increased from 0.13 (meta-analysis) to 0.25
s_j <- 0.85 #increased from 0.75 (apparent survival) to 0.85
s_a <- 0.98 #increased from 0.96 (apparent survival) to 0.98

# Specify productivity of juvenile and adult females
demog_sched[ma, "mx"] <- bp * f * s_n * ve * pf *s_h

# Specify juvenile and adult survival rates
demog_sched[1:(ma-1), "sx"] <- s_j
demog_sched[ma, "sx"] <- s_a

# Construct a Leslie matrix from this demography schedule
A <- make_Leslie_matrix(demog_sched)

# Calculate the asymptotic growth rate of the population governed by this 
#   demography schedule:
lambda1(A)



##### Population projection -----

## Specify features of simulation
nreps <- 10
nyears <- 50
N <- matrix(0, nyears, nreps)  # initializes totalpop matrix to store trajectories

## Initial population size and structure
# What is the stable-stage distribution (SSD) of the matrix
ssd <- stable.stage(A)

# How many individuals in the population?
n <- 50

# Spread individuals across the SSD
n_ssd <- ssd * n

# Use poisson draws to randomly populate numbers per stage for each simulation
#   replicate
n_i <- matrix(NA, nreps, dim(A)[1])
for (i in 1:nreps){
   n_i[i,]  <- rpois(length(n_ssd), n_ssd)
}
# rowSums(n_i) # usually between 40-60 females to start

## Partition matrix for simulation
x <- splitA(A)
x_T <- x$T
x_F <- x$F

## Run the simulation!
for (j in 1:nreps){ # for each replicate
  n <- n_i[j,] # specify initial population size, randomly drawn above
  for (i in 1:nyears){ # for each year
    n <- multiresultm(n, x_T, x_F)
    N[i,j] <- sum(n)
  } #year
} #rep
matplot(N, type = 'l', log="y",
        xlab = 'Time (years)', ylab = 'Total population')

# This projection works well at tracking total population size,
# but does not track population structure through time.




##### Project and track population structure -----

## Specify features of simulation
nreps <- 10
nyears <- 50
nstages <- dim(A)[1]

# Matrices to save population size and structure
N_tot <- matrix(0, nyears, nreps)
N_stages <- array(0, c(nyears, nreps, nstages))
# array to save population structure across years, reps, and ages

## Initial population size and structure
# What is the stable-stage distribution (SSD) of the matrix
ssd <- stable.stage(A)

# How many individuals in the population?
n <- 50

# Spread individuals across the SSD
n_ssd <- ssd * n

# Use poisson draws to randomly populate numbers per stage for each simulation
#   replicate
n_i <- matrix(NA, nreps, dim(A)[1])
for (i in 1:nreps){
  n_i[i,]  <- rpois(length(n_ssd), n_ssd)
}
# rowSums(n_i) # usually between 40-60 females to start

## Partition matrix for simulation
x <- splitA(A)
x_T <- x$T
x_F <- x$F

## Run the simulation!
for (j in 1:nreps){ # for each replicate
  n_i_j <- n_i[j,] # specify initial population size, randomly drawn above
  for (i in 1:nyears){ # for each year
    # If it's the first year, specify n_i_j; else,
    #   specify N_stage from previous yeaar
    if(i ==1){
      n_sim <- n_i_j
    } else {
        n_sim <- N_stages[i-1,j,]
        }
    N_stages[i,j,] <- multiresultm(n_sim, x_T, x_F) # project pop
    N_tot[i,j] <- sum(N_stages[i,j,]) # save N
  } #year
} #rep
matplot(N_tot, type = 'l', log="y",
        xlab = 'Time (years)', ylab = 'Total population')

# This works well, and we can track age-specific abundances