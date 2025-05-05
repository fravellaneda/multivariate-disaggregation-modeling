#Simulations

library(spatstat)
library(ggplot2)
library(viridis)
library(mvtnorm)
library(dplyr)
library(sf)
library(terra)
require(gridExtra)
library(INLA)
library(SpatialEpi)

set.seed(1234)

# Region of interest
boundary <- st_as_sf(as.polygonal(square(1)))

# Simulations parameters
alpha <- c(0.1, 0.05, -0.1) 
m.var <- c(0.25,0.2, 0.15)
range <- c(0.1,0.2,0.1)

beta <- c(0.3, 0.1, -0.3) 
e.sd <- c(0.1, 0.2, 0.15)

# Simulation grid
resolution<- 48
x <- seq(0+(1/(2*resolution)),1-(1/(2*resolution)), 1/resolution)
dp_pre <- as.matrix(expand.grid(x, x))

# Initial values
theta.ini <- (c(log(range), log((1/m.var)**2), log(1 / e.sd^2),
                beta))[c(9,10,11,1,5,2,6,3,7,12,13,14)]
theta.ini.e = theta.ini + rnorm(length(theta.ini), 0, 0.1)

# Number of areas in each case
areas <- c(36,64,144,36,64,144)

for (j in 1:6){
  n_areas <- areas[j]
  for (seed in 1:100){
    print(paste0("Group: ",j," Iteration:",seed))
    #set.seed(seed)
    
    source("create_dataset.R")
    
    source("disaggregation_model.R")
  }
}

