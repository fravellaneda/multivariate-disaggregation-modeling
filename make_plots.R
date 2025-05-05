# to make plots of the results after modeling

library(ggplot2)
library(latex2exp)
require(gridExtra)
library(sf)
library(viridis)
library(zeallot)

total <- 100

gen.plots <- function(s =sim,seed = sample(1:total,1)){
  dpdf <- read.csv(paste0("./results/sim.results.",s,"/pred_",seed,".csv"))
  map <- st_read(paste0("./results/sim.results.",s,"/maps/map_",seed,".shp"),
                 crs="")
  
  plots <- list()
  f <- 1
  
  for (i in 1:3){
    
    mini <- min(dpdf[,c(paste0("mean.l",i),paste0("pre.l.",i))])
    maxi <- max(dpdf[,c(paste0("mean.l",i),paste0("pre.l.",i))])
    
    d1 <- dpdf[,c("x","y",paste0("mean.l",i))]
    names(d1) <- c("x","y","value")
    
    plots[[f]] <- ggplot(map) +
      geom_tile(data = d1,aes(x = x, y = y, fill = value)) +
      geom_sf(fill=NA,color="black") +
      ggtitle("Multivariate Disaggregation Modeling") +
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_fill_viridis(TeX(paste0("$\\w_{",i,"}(\\textbf{s})$")),
                         limits=c(mini,maxi)) 
    
    f <- f+1
    
    d2 <- dpdf[,c("x","y",paste0("pre.l.",i))]
    names(d2) <- c("x","y","value")
    plots[[f]] <- ggplot() + 
      geom_tile(data = d2,aes(x = x, y = y, fill = value)) +
      geom_sf(fill=NA,color="gray") +
      ggtitle("True values") +
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_fill_viridis(TeX(paste0("$\\w_{",i,"}(\\textbf{s})$")),
                         limits=c(mini,maxi))
    
    
    f <- f+1
  }
  return(list(plots,seed))
}

c(plots, seed) %<-% gen.plots(s=3)

grid.arrange(plots[[1]],plots[[3]],plots[[5]],
             plots[[2]],plots[[4]],plots[[6]],nrow=2,
             top=paste0("Seed ",seed," Scenario ",s))
