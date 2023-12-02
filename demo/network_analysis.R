#Demonstration case II: Analyze an Interchange (a small network)

#Load required packages
library(osmdata) # This package is for loading OSM data.
library(tidyverse) # this package includes ggplot2 and dplyr for data visualization and manipulation.
library(sf) # simple features (sf) is a standard way to encode spatial vector data.
library(leaflet) # Ths package is used for plotting maps.
library(elevatr)   #This package is for requesting the USGS elevation data
library(forecast)  #This package is used for outlier identification in vertical alignment
library(strucchange) # This package is used for segmenting roads in the vertical alignment analysis

#load function files
source('fun_get_routes.R') #identify unique routes
source('fun_horizontal_analysis.R') #horizontal alignment analysis
source('fun_get_elevation.R') #get elevation for an sf object
source('fun_vertical_analysis.R') #vertical analysis

#Define bounding box (minlng,maxlng,minlat,maxlat) of AOI.
#An alternative is to use getbb() function such as bbox <- getbb("Tucker Georgia")

cb<-c(32.40536704789187, -84.92423455993193) #center of a interchange near Columbus, GA

bbox <-
  matrix(
    c(cb[2] - 0.009, cb[2] + 0.009, cb[1] - 0.009, cb[1] + 0.009), #size of the bbox
    byrow = TRUE,
    nrow = 2,
    ncol = 2,
    dimnames = list(c('x', 'y'), c('min', 'max'))
  ) 

cropbox <- array(bbox,
                 dimnames = list(c('xmin', 'ymin', 'xmax', 'ymax')))

#Load the motorway (controlled access one-way roads) and motorway link (ramps) shape data within the bbox.
#Reference for keys and values: https://wiki.openstreetmap.org/wiki/Key:highway
mtwy <- bbox %>% opq() %>%
  add_osm_feature(key = "highway", value = c("motorway")) %>%
  osmdata_sf()

mtwy_line <- mtwy$osm_lines %>% st_crop(cropbox)
mtwy_pt <- mtwy$osm_points %>% st_crop(cropbox)

cat(length(st_geometry(mtwy_line)), " motorway segments loaded from OSM.") #Number of road segments
cat(length(st_geometry(mtwy_pt)), " motorway nodes loaded from OSM.") #Number of road segments

mtlk <- bbox %>% opq() %>%
  add_osm_feature(key = "highway", value = c("motorway_link")) %>%
  osmdata_sf()

mtlk_line <- mtlk$osm_lines %>% st_crop(cropbox)
mtlk_pt <- mtlk$osm_points %>% st_crop(cropbox)

cat(length(st_geometry(mtlk_line)), " motorway ramp segments loaded from OSM.") #Number of road segments
cat(length(st_geometry(mtlk_pt)), " motorway ramp nodes loaded from OSM.") #Number of road segments


#Plot raw data using leaflet for a visual inspection.
#Reference for map options: https://rstudio.github.io/leaflet/markers.html
leaflet() %>%
  addTiles() %>% #There are other interesting base maps (https://rstudio.github.io/leaflet/basemaps.html)
  addPolylines(data=mtwy_line, color = 'blue', weight = 4) %>%
  addPolylines(data=mtlk_line, color = 'green', weight = 4) %>%
  addRectangles(lng1=cropbox[1], lat1=cropbox[2],lng2=cropbox[3], lat2=cropbox[4],
                weight=2,color='purple',fillColor = "transparent") %>%
  addScaleBar(position='bottomleft')

#Gall find_routes function 
sfc_mtwy <- get_routes(mtwy_line) #Identify routes from motorway
sfc_mtlk <- get_routes(mtlk_line) #Identify routes from ramps
routes<-rbind(sfc_mtwy,sfc_mtlk)  #Merge the two sfcs together into one sfc object


################ Horizontal Analysis ##################
#Call the horizontal analysis function for each route in the network
#This function returns horizontal tangents and curves in an sfc object 
#All results are merged into one sfc
#Track time

start.time <- Sys.time()
hr <- horizontal_analysis(routes[1, ])
for (i in 2:nrow(routes)) {
  hr <- rbind(hr, horizontal_analysis(routes[i, ]))
}
end.time <- Sys.time()
end.time - start.time

#Seperate the horizontal curves and tangents
hc <- hr %>% filter(nmpd >=0) #h curves
ht <- hr %>% filter(is.na(nmpd)) #h tangents


#filtering for quality control
hc <- hc[!duplicated(hc$geom), ] #Remove duplicates
hc <- hc%>%filter(nmpd <=0.05 & length%>%as.numeric()>100) #filter based on fit quality and length
hc <- hc[st_covered_by(hc)%>%lengths()==1,] #remove shorter segments that are covered by a longer segment

ht <- ht[!duplicated(ht$geom),] #Remove duplicates
ht <- ht%>%filter(length%>%as.numeric()>100) #filter based on fit quality and length
ht <- ht[st_covered_by(ht)%>%lengths()==1,] #remove shorter segments that are covered by a longer segment


#plot curves
leaflet() %>%
  addProviderTiles(providers$CartoDB.Positron) %>% 
  addPolylines(data=hc, color='blue', weight = 3) %>%
  addPolylines(data=ht, color='black', weight = 3) %>%
  addCircles(data=hc,lng = ~lngc, lat = ~latc, radius = ~rad, color = "red",  
             weight = 1,fillColor = "transparent") %>%
  addScaleBar(position='bottomright')


################ Vertical Analysis ##################

#Get elevation data to transform routes from xy to xyz format, recommended resolution is 10m
#Speed depends on 3DEP server status
routes_xyz<-get_elevation(routes[1,],10)
for (i in 2:nrow(routes)){
  cat("download elevation for ,", i, " of ", nrow(routes))
  routes_xyz<-rbind(routes_xyz,get_elevation(routes[i,],10)) 
} 

#Call the vertical_analysis function.
#This function returns vertical tangents and curves in an sfc object 
#All results are merged into one sfc object named vr
#Speed depends on the number and lengths of the routes.
#Track time
start.time <- Sys.time() #start time
vr <- vertical_analysis(routes_xyz[1,]) 
for (i in 2:nrow(routes)){
  cat("vertical analysis,", i, " of ", nrow(routes))
  vr <- rbind(vr,vertical_analysis(routes_xyz[i,]))
}
end.time <- Sys.time()
end.time - start.time

#Apply Curve and tangent classification criteria
vr <- vr %>% mutate(type=ifelse(r2>=0.5,"vc",
                            ifelse(r2<=0.3 & std<=1,"vt","u")))

vr<-vr[!duplicated(vr2$geom),] #Remove duplicates. This is very unlikely

#Filter based on length and duplicates
vrc<-vr[vr$type=='vc',]
vrc<-vrc[vrc$length%>%as.numeric()>100,] 
vrc<-vrc[st_covered_by(vrc)%>%lengths()==1,] #remove shorter segments that are covered by a longer segment

vrt<-vr[vr$type=='vt',]
vrt<-vrt[vrt$length%>%as.numeric()>100,]
vrt<-vrt[st_covered_by(vrt)%>%lengths()==1,] #remove shorter segments that are covered by a longer segment

#Draw map of curves and tangents
leaflet() %>%
  addProviderTiles(providers$CartoDB.Positron) %>% 
  addPolylines(data=vrc%>%st_zm(), color='red', weight = 3) %>%
  addPolylines(data=vrt%>%st_zm(), color='black', weight = 3) %>%
  addScaleBar(position='bottomright')

#Analyze a specific route in the network
exp_route<-routes[4,]
plot(exp_route)
exp_routes_xyz<-get_elevation(routes[4,],10)
exp_vr <- vertical_analysis(exp_routes_xyz)
exp_vr <- exp_vr %>% mutate(type=ifelse(r2>=0.5,"vc",
                            ifelse(r2<=0.3 & std<=1,"vt","u")))
coords <- exp_routes_xyz %>% st_coordinates()
len_tot <- st_length(routes_xyz) %>% as.numeric()
df <- data.frame(dist=(c(1:nrow(coords))-1)/(nrow(coords)-1)*len_tot, elev=coords[,3])


plot(df,xlab="Distance (m)",ylab="Elevation (m)",type="l",col = 'blue',
     cex.lab = 1.5, cex.axis = 1.5,, ylim=c(70,120))
abline(v = df[exp_vr[-1,]$from,]$dist,lty = 'dashed')
for (i in 1:nrow(exp_vr)){
  if (exp_vr[i,]$type=='vc') {lines(x=c(df[exp_vr[i,]$from,]$dist, df[exp_vr[i,]$to,]$dist),y=c(77,77),col='red',lwd = 4)}
  if (exp_vr[i,]$type=='vt') {
    lines(x=c(df[exp_vr[i,]$from,]$dist, df[exp_vr[i,]$to,]$dist),y=c(75,75),col='black',lwd = 4)
    text(x=df[exp_vr[i,]$from,]$dist+100,y=222,paste('G = ', round(exp_vr[i,]$avg,2),'%'),srt=90)}
}
legend(200,120,lty = 1,lwd = 4,col = c('red','black'),legend=c("Vertical curve", "Vertical tangent"),cex = 1.5)

