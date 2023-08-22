### Genomic Sample Database Work ####

library(geosphere)
library(sp)
library(mapdata)
library(maps)
library(ggmaps)

AllPos = read.delim("~/Documents/Postdoc Research/Genomic Sample Bank/Allpositives.txt")

#### 1. Plot all unique sample locations across CA ####

# First, keep only unique collection sites (i.e. unique lat/long)
# 1249 total samples, 238 unique locations
uniap = AllPos[!duplicated(AllPos[ , c("latitude", "longitude")]), ] 
uniaps = uniap
coordinates(uniaps) = ~longitude+latitude

# Obtain base CA county map
counties <- map_data("county")
ca_county <- counties %>% filter(region == "california")

# Plot 
ggplot(data = ca_county) + theme_bw() + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = alpha("#B09C85FF", 0.5), color = "white") + 
  coord_quickmap(xlim = c(-123, -118),  ylim = c(34, 39)) + 
  guides(fill = FALSE) +
  geom_point(data = uniap, mapping = aes(x = longitude, y = latitude),
             shape = 4, size = 3)


#### 2. Calculate euclidean distance between all sampled locations ####

# Function to calculate distance between 2 points
distm(c(uniap$longitude[1], uniap$latitude[1]), c(uniap$longitude[2], uniap$latitude[2]), fun = distHaversine)

# Now calculate for all possible pairs
distances = rep(NA, 238 * 238)

for (i in 1:238) {
  for (j in 1:238) {
  distances[(i + (238*j - 238))] = distm(c(uniap$longitude[i], uniap$latitude[i]), c(uniap$longitude[j], uniap$latitude[j]), fun = distHaversine)
}
}



