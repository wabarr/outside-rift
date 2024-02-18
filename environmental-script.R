library(sf)
library(tidyr)
library(dplyr)
library(ggplot2)
library(pander)

blue <- "#033C5A"
buff <- "#AA9868"

theme_set(theme_bw(14))

## set your working directory to the directory where the script file and data files live
mamms_full <- st_read("./IUCN_mammals_presence1_selectedOrders_disolvedBinomial_intersectGREGORY.gpkg")
mamms_intersection <- st_read("./IUCN_mammals_presence1_selectedOrders_disolvedBinomial_intersectGREGORY_intersection.gpkg")

africa <- st_read("./africa_simple_WGS84_geographic.geojson")

rift <- st_read("./great_rift_valley_ofb.json")
rift <- rift[1,]

## find the area, using planar geometry to avoid validation errors
## for some of the polygons
sf_use_s2(FALSE)
fullAreas <- st_area(mamms_full)
intersectingAreas <- st_area(mamms_intersection)
sf_use_s2(TRUE)


## do a comparison between areas computed on the plane versus on the sphere
tryArea <- function(object, i) {
  tryCatch(st_area(object[i,]), error=function(e) return(NA))
}
sphericalFullAreas <- lapply(1:nrow(mamms_full), function(i) tryArea(mamms_full, i))
# all spherical full areas are larger than the planar full areas
diffs <- as.numeric(fullAreas) - as.numeric(sphericalFullAreas)

## median percentage error between planar and spherical os 0.4%
median(diffs / as.numeric(sphericalFullAreas), na.rm=TRUE) * 100

## largest percentage error between planar and spherical is 0.45%
min(diffs / as.numeric(sphericalFullAreas), na.rm=TRUE) * 100

#decision made to go with planar geometry! it's not a big deal
rm(sphericalFullAreas)
mamms_full$planar_area <- as.numeric(fullAreas) / 1000^2
mamms_intersection$planar_area <- as.numeric(intersectingAreas) / 1000^2

areas_full <- select(mamms_full, binomial, planar_area) %>% st_drop_geometry()
areas_intersection <- select(mamms_intersection, binomial, planar_area) %>% st_drop_geometry()

overlaps <- data.frame(overlaps=areas_intersection$planar_area / areas_full$planar_area)
ggplot(overlaps, aes(x=overlaps*100)) + 
  geom_density(fill=blue, alpha=0.7) + 
  labs(title = "Percentage of total range within EARS", x="% overlap") + 
  theme_bw(16) + 
  annotate("text", x=50, y=0.15, label=sprintf("median overlap = %0.2f%%", median(overlaps$overlaps)*100))
ggsave("~/Dropbox/thinking_outside_the_rift/environmental-analyses-Gregory-rift-only/range_percent_overlaps.jpg", width=6, height=6)


bioclim_intersection <- read.table("./GEE-exports/intersection-stats-bioclim-gregory.csv", header=T, sep=",")
#convert bioclim temp variables to degrees celcius
bioclim_intersection$bio01_mean <- bioclim_intersection$bio01_mean/10
bioclim_intersection$range_subset <- "rift"
bioclim_intersection_long <- 
  bioclim_intersection %>%
  left_join(areas_intersection) %>% 
  select(binomial, order_, family, genus, range_subset, contains("_mean")) %>%
  pivot_longer(names_to="metric", cols = contains("_mean"))

bioclim_full <- read.table("./GEE-exports/fullrange-stats-bioclim-gregory.csv", header=T, sep=",")
#convert bioclim temp variables to degrees celcius
bioclim_full$bio01_mean <- bioclim_full$bio01_mean/10
bioclim_full$range_subset <- "full"
bioclim_full_long <- 
  bioclim_full %>%
  left_join(areas_full) %>% 
  select(binomial, order_, family, genus, range_subset, contains("_mean")) %>%
  pivot_longer(names_to="metric", cols = contains("_mean"))

copernicus_intersection <- read.table("./GEE-exports/intersection-stats-copernicus-gregory.csv", header=T, sep=",")
copernicus_intersection$range_subset <- "rift"
copernicus_intersection_long <-
  copernicus_intersection %>% 
  left_join(areas_intersection) %>% 
  select(binomial, order_, family, genus,contains("coverfraction_mean"), range_subset, planar_area) %>%
  pivot_longer(names_to="metric", cols = contains("coverfraction_mean"))

copernicus_full <- read.table("./GEE-exports/fullrange-stats-copernicus-gregory.csv", header=T, sep=",")
copernicus_full$range_subset <- "full"
copernicus_full_long <-
  copernicus_full %>% 
  left_join(areas_full) %>% 
  select(binomial, order_, family, genus,contains("coverfraction_mean"), range_subset, planar_area) %>%
  pivot_longer(names_to="metric", cols = contains("coverfraction_mean"))

copernicus <- do.call(rbind, list(copernicus_full_long, copernicus_intersection_long))
copernicus$metric <- factor(copernicus$metric)
levels(copernicus$metric) <- gsub(" mean", " cover (%)", gsub(".coverfraction_", " ", levels(copernicus$metric)))


bioclim <- do.call(rbind, list(bioclim_full_long, bioclim_intersection_long))
bioclim$metric <- factor(bioclim$metric)
levels(bioclim$metric) <- c("annual mean temperature (°C)", "annual precipitation (mm)", "precip. wettest month (mm)", "precip. driest month (mm)", "precipitation seasonality (CV)" )


copernicus$range_subset <- ordered(copernicus$range_subset, levels=c("rift", "full"))
bioclim$range_subset <- ordered(bioclim$range_subset, levels=c("rift", "full"))

legend_labs <- c("species range in rift", "full species range")

mods_copernicus <- 
  group_by(copernicus, metric) %>%
  do(mod=aov(value~range_subset + Error(binomial), data=.)) 
names(mods_copernicus$mod) <- mods_copernicus$metric

mods_bioclim <- 
  group_by(bioclim, metric) %>%
  do(mod=aov(value~range_subset + Error(binomial), data=.))
names(mods_bioclim$mod) <- mods_bioclim$metric

all_combined <- rbind(select(bioclim, range_subset, metric, value), select(copernicus, range_subset, metric, value))

ggplot(data=all_combined, aes(x=value)) + 
  geom_density(aes(fill=range_subset, color=range_subset), alpha=0.7) + 
  facet_wrap(~metric, scale="free", nrow = 3) + 
  labs(fill="", color="") +
  scale_fill_manual(values=c(blue, buff), labels=legend_labs, name="") + 
  scale_color_manual(values=c(blue, buff),labels=legend_labs, name="") + 
  theme_bw(11) + 
  theme(legend.position="bottom", axis.text = element_text(size=9))
