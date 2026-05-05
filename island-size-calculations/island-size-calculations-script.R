library(sf)
library(sp)
library(terra)
library(PleistoDist)

#load shapefile
islands <- st_read("E:/Jenna_Maps/Todiramphus_specific.shp")

#define Equi7 Oceania projection
prj_oceania <- st_crs("+proj=aeqd +lat_0=-19.5 +lon_0=131.5 +x_0=6988408.536 +y_0=7654884.537 +datum=WGS84 +units=m +no_defs +type=crs")

#project islands to Equal Earth Asia-Pacific
islands_projected <- st_transform(islands,crs="EPSG:8859")

islands_projected

#load taxa names as list
taxa <- c("alb_alb","alb_orii","alb_owston","albonotatu","aus_aus","aus_dam","aus_inter","aus_odites","aus_tring","chl_azelus","chl_chlori","chl_chloro","chl_collar","chl_davis","chl_occip","chl_terao","cinnamo","colonus","diops","enigma","farquhari","funebris","gertrudae","godeffroyi","lazuli","leucopy","pelewensis","recurviros","reichen","ruficollar","sac_amoe","sac_aneity","sac_brachy","sac_erroma","sac_eximi","sac_juliae","sac_mala","sac_manuae","sac_marin","sac_melano","sac_ornatu","sac_pavuvu","sac_pealei","sac_regina","sac_sacer","sac_santo","sac_solomo","sac_soror","sac_tannen","sac_torres","sac_tucopi","sac_utupua","sac_vicina","sac_viti","sanc_cana","sanc_macmi","sanc_norf","sanc_vagan","sauro_adm","sauro_ana","tris_alber","tris_benne","tris_matt","tris_novae","tris_nusae","tris_stres","tris_tris","tutus_aitu","tutus_mauk","tutus_tutu","ven_venera","ven_young","winch_alf","winch_mind","winch_nesy","winch_nigo","winch_winc")

#create data frame
output <- tibble(.rows = 77)
output$taxon <- taxa
ta <- c()
ma <- c()
cha <- c()
centroids <- data.frame()
ch <- data.frame()


#run for loop that iterates through taxa names, subsets islands, and calculates area
for (i in taxa) {
  #subset relevant features
  islandssubset <- subset(islands_projected,islands_projected[[i]] == "1")

  #calculate centroid
  centroids <- rbind(centroids,st_centroid(islandssubset))

  #calculate total area
  total_area <- sum(st_area(islandssubset))
  ta <- c(ta,total_area)

  #calculate mean area
  mean_area <- mean(st_area(islandssubset))
  ma <- c(ma,mean_area)

  #calculate minimum convex hull area
  conv_hull <- st_as_sf(st_convex_hull(st_union(islandssubset)))
  conv_hull$name = i
  ch <- rbind(ch,conv_hull)
  convex_hull_area <- st_area(conv_hull)
  cha <- c(cha,convex_hull_area)
}

output$totalarea <- ta
output$meanarea <- ma
output$convexhull <- cha

#write output file with taxon areas
write.csv(output,"E:/Jenna_Maps/Todiramphus_areas.csv")

#write output shapefile with convex hull polygons
st_write(ch,"E:/Jenna_Maps/convex_hulls.shp")

#project convex hull shapefile to equidistant projection
convexhulls_projected <- st_transform(ch,crs=prj_oceania)

#extract centroids from convex hull polygons
convexhulls_centroids <- st_centroid(convexhulls_projected)

#calculate inter-convex hull distances
convexhulls_distances <- st_distance(convexhulls_centroids)

#write distance matrix to file
write.csv(convexhulls_distances,"E:/Jenna_Maps/centroid_distances.csv")
