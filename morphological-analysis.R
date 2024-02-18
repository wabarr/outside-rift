library(sf)
library(geometry)
library(dplyr)
library(ggplot2)
library(gridExtra)
theme_set(theme_bw(12))

## set your working directory to the directory where the script file and data files live

##read data files
baboons <- read_sf("./baboons-PCs.gpkg")
guenons <- read_sf("./guenons-PCs.gpkg")

baboons$whichrift[which(is.na(baboons$whichrift))] <- "outside"
guenons$whichrift[which(is.na(guenons$whichrift))] <- "outside"


## read PC summary files

baboon_PC_summary <- read.table("./baboon-PC-variance-explained.csv", header=T, sep = ",")
guenon_PC_summary <- read.table("./guenons-PC-variance-explained.csv", header=T, sep = ",")



fullBaboonHull <- baboons %>% 
  select(PC1:PC3) %>% 
  st_drop_geometry() %>% 
  convhulln(options="FA")

riftBaboonHull <- baboons %>% 
  filter(whichrift=="northern") %>%
  select(PC1:PC3) %>% 
  st_drop_geometry() %>% 
  convhulln(options="FA")

fullGuenonHull <- guenons %>% 
  select(PC1:PC3) %>% 
  st_drop_geometry() %>% 
  convhulln(options="FA")

riftGuenonHull <- guenons %>% 
  filter(whichrift=="northern") %>%
  select(PC1:PC3) %>% 
  st_drop_geometry() %>% 
  convhulln(options="FA")

#we want the vol paramater which per documentation is 'the generalised volume of the hull. This is volume of a 3D hull or the area of a 2D hull'
riftBaboonHull$vol/ fullBaboonHull$vol* 100
riftGuenonHull$vol/ fullGuenonHull$vol* 100

plotcolz <- c("#1f78b4", "#b2df8a")

StatChull <- ggproto("StatChull", Stat,
                     compute_group = function(data, scales) {
                       data[chull(data$x, data$y), , drop = FALSE]
                     },
                     
                     required_aes = c("x", "y")
)

stat_chull <- function(mapping = NULL, data = NULL, geom = "polygon",
                       position = "identity", na.rm = FALSE, show.legend = NA, 
                       inherit.aes = TRUE, ...) {
  layer(
    stat = StatChull, data = data, mapping = mapping, geom = geom, 
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}

guenonPC1PC2 <- ggplot(guenons, aes(x=PC1, y=PC2)) + 
  geom_point(aes(color=whichrift)) + 
  stat_chull(aes(color=whichrift), fill="transparent", linewidth=1.3) + 
  labs(x=sprintf("PC1 (%0.1f%% of variance)", guenon_PC_summary[1,2]*100),
       y=sprintf("PC2 (%0.1f%% of variance)", guenon_PC_summary[2,2]*100)
       ) + 
  scale_color_manual(values=plotcolz) + 
  theme(legend.position = "none")

guenonPC2PC3 <- ggplot(guenons, aes(x=PC2, y=PC3)) + 
  geom_point(aes(color=whichrift)) + 
  stat_chull(aes(color=whichrift), fill="transparent", linewidth=1.3) + 
  labs(x=sprintf("PC2 (%0.1f%% of variance)", guenon_PC_summary[2,2]*100),
       y=sprintf("PC3 (%0.1f%% of variance)", guenon_PC_summary[3,2]*100)
  ) + 
  scale_color_manual(values=plotcolz) + 
  theme(legend.position = "none")


  

baboonPC1PC2 <- ggplot(baboons, aes(x=PC1, y=PC2)) + 
  geom_point(aes(color=whichrift)) + 
  stat_chull(aes(color=whichrift), fill="transparent", linewidth=1.3) + 
  labs(x=sprintf("PC1 (%0.1f%% of variance)", baboon_PC_summary[1,2]*100),
       y=sprintf("PC2 (%0.1f%% of variance)", baboon_PC_summary[2,2]*100)
  ) + 
  scale_color_manual(values=plotcolz) + 
  theme(legend.position = "none")


baboonPC2PC3 <- ggplot(baboons, aes(x=PC2, y=PC3)) + 
  geom_point(aes(color=whichrift)) + 
  stat_chull(aes(color=whichrift), fill="transparent", linewidth=1.3) + 
  labs(x=sprintf("PC2 (%0.1f%% of variance)", baboon_PC_summary[2,2]*100),
       y=sprintf("PC3 (%0.1f%% of variance)", baboon_PC_summary[3,2]*100)
  ) + 
  scale_color_manual(values=plotcolz) + 
  theme(legend.position = "none")




## read in random rift polygons 
randomRifts <- st_read("./randomRifts.gpkg")


computeHullCountSpecimens <- function(randomRift, datfram) {
  ## for the spatial join to work to see which primates fall in the random rift
  ## you must make the polygon a sf feature, with at least one field, which will be joined to datframe 
  randomRift$col_for_finding_specimens_in_polygon <- 1
  joined <- st_join(datfram, randomRift, join=st_within)
  nspecimens <- sum(!is.na(joined$col_for_finding_specimens_in_polygon))
  if(nspecimens < 4) return(NA) #4 is the minimum number of points needed by qhull for a 3D hull calculation
  filtered <- filter(joined, !is.na(col_for_finding_specimens_in_polygon))
  convHull <- select(filtered, PC1, PC2, PC3) %>% st_drop_geometry %>% convhulln(options="FA")

  return(cbind(randomRift,
    nspecimens=nspecimens, 
    hullvol=convHull$vol))
}


## calculate hulls for baboons in random rifts
baboonRandomResults <- vector("list", nrow(randomRifts))
for(i in 1:nrow(randomRifts)) {
  baboonRandomResults[[i]] <- computeHullCountSpecimens(randomRifts[i,], baboons)
}
baboonRandomResults <- baboonRandomResults[-which(is.na(baboonRandomResults))]
baboonRandomResults <- do.call(rbind, baboonRandomResults)

## baboons in actual rift
baboonsInActualRift <- data.frame(
  hullvol=riftBaboonHull$vol,
  nspecimens=nrow(filter(baboons, whichrift=="northern"))
)

modBaboonRandom <- lm(hullvol^(1/3)~log(nspecimens), data=baboonRandomResults)
baboonRandomResults$samplingEfficiency <- resid(modBaboonRandom)
summary(modBaboonRandom)

predictedBaboonVal <- predict(modBaboonRandom, newdata = baboonsInActualRift)

baboonRiftEfficiencyPercentile <- mean(resid(modBaboonRandom) < baboonsInActualRift$hullvol^(1/3) - predictedBaboonVal) * 100
baboonRiftEfficiencyPercentile

baboonResidPlot <- ggplot(data=data.frame(resid=resid(modBaboonRandom)), aes(x=resid)) + 
  geom_histogram(color="white") + 
  geom_vline(color="red", 
             size=1.4,
             xintercept=baboonsInActualRift$hullvol^(1/3) - predict(modBaboonRandom, newdata = baboonsInActualRift)) + 
  annotate(x=-0.0005, y=23, geom="text", label="eastern rift", color="red", size=8) + 
  labs(x="sampling efficiency in random rifts", title=sprintf("Baboons %0.0frd percentile for sampling efficiency", baboonRiftEfficiencyPercentile)) + 
  theme_bw(14) 
baboonResidPlot



## calculate hulls for guenons in random rifts

guenonRandomResults <- vector("list", nrow(randomRifts))
for(i in 1:nrow(randomRifts)) {
  guenonRandomResults[[i]] <- computeHullCountSpecimens(randomRifts[i,], guenons)
}
guenonRandomResults <- guenonRandomResults[-which(is.na(guenonRandomResults))]
guenonRandomResults <- do.call(rbind, guenonRandomResults)

## guenons in actual rift
guenonsInActualRift <- data.frame(
  hullvol=riftGuenonHull$vol,
  nspecimens=nrow(filter(guenons, whichrift=="northern"))
)



modGuenonRandom <- lm(hullvol^(1/3)~log(nspecimens), data=guenonRandomResults)
summary(modGuenonRandom)

predictedGuenonVal <- predict(modGuenonRandom, newdata = guenonsInActualRift)


guenonRiftEfficiencyPercentile <- mean(resid(modGuenonRandom) < guenonsInActualRift$hullvol^(1/3) - predictedGuenonVal) * 100
guenonRiftEfficiencyPercentile

guenonResidPlot <- ggplot(data=data.frame(resid=resid(modGuenonRandom)), aes(x=resid)) + 
  geom_histogram(color="white") + 
  geom_vline(color="red", 
             size=1.4, 
             xintercept=guenonsInActualRift$hullvol^(1/3) - predictedGuenonVal) + 
  labs(x="sampling efficiency in random rifts", title=sprintf("Guenons+ %0.0fst percentile for sampling efficiency", guenonRiftEfficiencyPercentile)) + 
  annotate(x=-0.012, y=28, geom="text", label="eastern rift", color="red", size=8) + 
  theme_bw(14)
guenonResidPlot




sumVars <- function(df) {
  vars <- sapply(df, FUN=function(x){
    if(is(x, "numeric")) var(x)
  })
  return(sum(unlist(vars)))
}

## use sum of variance for first 3 PCs of baboons
rezInRiftBaboon <- filter(baboons, whichrift=="northern") %>% select(PC1:PC3) %>% sumVars()
rezBaboon <- baboons %>% select(PC1:PC3) %>% sumVars()
rezInRiftBaboon/rezBaboon

rezInRiftGuenon <- filter(guenons, whichrift=="northern") %>% select(PC1:PC3) %>% sumVars()
rezGuenon <- guenons %>% select(PC1:PC3) %>% sumVars()
rezInRiftGuenon/rezGuenon



