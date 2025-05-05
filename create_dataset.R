#Simulations

library(spatstat)
library(dplyr)
library(sf)
library(INLA)

source("aux_fun.R")

fun_gen <- function(x,y) {n_areas*(3/2)*((x)**2+(y-0)**2)}

# Create latent proceses z1,z2 and z3

# Create the directory if it doesn't exist
if (!dir.exists("./results/sims/")) {
  dir.create("./results/sims/", recursive = TRUE)
}

if (j==1){
  z1 <- book.rMatern(1, dp_pre, range = range[1],
                     sigma = sqrt(m.var[1]))
  z2 <- book.rMatern(1, dp_pre, range = range[2],
                     sigma = sqrt(m.var[2]))
  z3 <- book.rMatern(1, dp_pre, range = range[3],
                     sigma = sqrt(m.var[3]))
  save(z1,z2,z3,
       file = paste0("./results/sims/sim_",seed,".RData"))
} else{
  load(paste0("./results/sims/sim_",seed,".RData"))
}

# Use the latent variables to create the response variables acording with
# The corregionalizated model

pre.l.1 <- alpha[1] + z1 
pre.l.2 <- alpha[2] + beta[1] * z1 + z2
pre.l.3 <- alpha[3] + beta[2] * z1 + beta[3] * z2 + z3

dpdf <- as.data.frame(dp_pre)
names(dpdf) <- c("x","y")

dpdf$pre.l.1 <- pre.l.1
dpdf$pre.l.2 <- pre.l.2
dpdf$pre.l.3 <- pre.l.3

dpdf_sf <- st_as_sf(dpdf,coords = c("x","y"))
dpdf_sf$id.1 <- 1:nrow(dpdf)

# Create the mesh for the INLA estimation

dp <- dp_pre
bnd <- inla.nonconvex.hull(dp)

mesh <- inla.mesh.2d(
  boundary = bnd,
  max.edge = c(0.04)
)

mesh_sf <- st_as_sf(data.frame(x=mesh$loc[,1],y=mesh$loc[,2]),
                    coords = c("x","y"))

mesh_df <- as.data.frame(mesh$loc[,1:2])
names(mesh_df) <- c("x","y")

areas.0 <- 1

while (areas.0>0){
  if (j>3){
    voronoi <- st_make_grid(boundary, cellsize = 1/sqrt(n_areas))
    areas.0 <- 0
  } else{
    p_gen <- rpoispp(fun_gen,win=square(1))
    p_gen_sf <- st_as_sf(data.frame(x=p_gen$x,y=p_gen$y),coords = c("x","y"))
    voronoi <- terra::voronoi(x = terra::vect(p_gen_sf),bnd = boundary)
    
    vor <- st_as_sf(voronoi)
    
    vor$counts.mesh <- lengths(st_contains(vor,mesh_sf))
    vor$counts.pred <- lengths(st_contains(vor,dpdf_sf))
    
    areas.0 <- nrow(vor[vor$counts.mesh==0,])+
      nrow(vor[vor$counts.mesh==0,])
  }
}

voronoi_sf <- st_as_sf(voronoi)
voronoi_sf$id <- 1:nrow(voronoi_sf)
voronoi_int <- st_intersection(voronoi_sf,dpdf_sf)
values <- as.data.frame(voronoi_int)
values_agg_pre <- aggregate(values[,c("pre.l.1","pre.l.2","pre.l.3")],
                            by=list(values$id), mean)
names(values_agg_pre)[1] <- "id"
values_agg <- data.frame(id=values_agg_pre$id)

for (f in 1:3){
  values_agg[,paste0("error",f)] <- rnorm(
    nrow(values_agg), 0, e.sd[f])
  values_agg[,paste0("pre.linear.",f)] <- values_agg_pre[,paste0("pre.l.",f)]
  values_agg[,paste0("linear",f)] <- values_agg_pre[,paste0("pre.l.",f)] +
    values_agg[,paste0("error",f)]
}

map <- merge(voronoi_sf,values_agg,on="id",all.x=TRUE)

spde <- inla.spde2.pcmatern(
  alpha = 2,
  mesh = mesh,
  prior.range = c(0.6, 0.99), # P(range < 0.5) = 0.01
  prior.sigma = c(1, 0.01)) # P(sigma > 1) = 0.01

#spde <- inla.spde2.matern(mesh = mesh, alpha = 2, constr = TRUE)

# Number of mesh nodes
nv <- mesh$n

# number of groups
n.pp <- 3

#number of areas
n.ar <- nrow(map)

l.1 <- matrix(NA, nrow = n.ar, ncol = n.pp)
l.1[, 1] <- map$linear1

l.2 <- matrix(NA, nrow = n.ar, ncol = n.pp)
l.2[, 2] <- map$linear2

l.3 <- matrix(NA, nrow = n.ar, ncol = n.pp)
l.3[, 3] <- map$linear3

# Projection matrix

t1 <-Sys.time()

A.m <- matrix(0, n.ar, nv)
inter <- as.list(st_contains(map[,c()],mesh_sf))

for (i in 1:n.ar){
  idx <- inter[[i]]
  A.m[i,idx] <- 1/length(idx)
}
t2 <-Sys.time()
print("Time projection matrix construction ")
print(t2-t1)
