#Demonstration case I: Analyze an 11-km single route

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
bbox <- matrix(
  c(-85.18, -85.08, 34.898, 34.95),
  byrow = TRUE,
  nrow = 2,
  ncol = 2,
  dimnames = list(c('x', 'y'), c('min', 'max'))
) 

#Load the motorway (controlled access one-way roads) shape data within the bbox.
#Reference for keys and values: https://wiki.openstreetmap.org/wiki/Key:highway
roads <- bbox %>% opq() %>%
  add_osm_feature(key = "highway", value = c("motorway")) %>%
  osmdata_sf()

#Crop the segments outside the bbox
#and then get unique routes from the raw data
cropbox <- array(bbox,
                dimnames = list(c('xmin', 'ymin', 'xmax', 'ymax')))

routes <- roads$osm_lines %>% 
  st_crop(cropbox) %>%
  get_routes()

cat(nrow(roads$osm_lines), " road segments loaded from OSM.") #Number of road segments
cat(nrow(routes), " unique routes identified.") #Number of road segments

#Plot raw data using leaflet for a visual inspection.
#Reference for map options: https://rstudio.github.io/leaflet/markers.html

leaflet() %>%
  addTiles() %>% #There are other interesting basemaps (https://rstudio.github.io/leaflet/basemaps.html)
  addPolylines(data=routes, color = 'black', weight = 2) %>%
  addCircleMarkers(data=routes%>%st_geometry()%>%st_cast('POINT'), color = 'blue',weight = 2,radius = 3) %>%
  addRectangles(lng1=cropbox[1], lat1=cropbox[2],lng2=cropbox[3], lat2=cropbox[4],
                weight=2,color='purple',fillColor = 'transparent') %>%
  addScaleBar(position='bottomleft')


#choose one route to analyze. Here we choose the south bound.
route <- routes[2, ]

leaflet() %>%
  addTiles() %>% #There are other interesting base maps (https://rstudio.github.io/leaflet/basemaps.html)
  addPolylines(data=route, color = 'black', weight = 2) %>%
  addCircleMarkers(data=route%>%st_geometry()%>%st_cast("POINT"), color = 'blue',weight = 2,radius = 3) %>%
  addRectangles(lng1=cropbox[1], lat1=cropbox[2],lng2=cropbox[3], lat2=cropbox[4],
                weight=2,color='purple',fillColor = "transparent") %>%
  addScaleBar(position='bottomleft')

################ Horizontal Analysis ##################
#Call the horizontal analysis function
#This function returns horizontal tangents and curves in an sfc object 
hr<-horizontal_analysis(route)

#Seperate the horizontal curves and tangents
hc <- hr %>% filter(nmpd < 0.05) #h curves
ht <- hr %>% filter(is.na(nmpd)) #h tangents

#Plot curves with best fit circles
leaflet() %>%
  addProviderTiles(providers$CartoDB.Positron) %>% 
  addPolylines(data=hc, color='blue', weight = 8,
               label = ~paste("r =",round(rad,0), " m"),
               labelOptions = labelOptions(noHide = T,
                                             textOnly=T, textsize = "16px",
                                             offset = c(20, -10))) %>%
  addPolylines(data=ht, color='black', weight = 2) %>%
  addCircles(data=hc,lng = ~lngc, lat = ~latc, radius = ~rad, color = "red",  
             weight = 1,fillColor = "transparent") %>%
  addScaleBar(position='bottomleft')

#Plot tangents
leaflet() %>%
  addProviderTiles(providers$CartoDB.Positron) %>% 
  addPolylines(data=hc, color='blue', weight = 2) %>%
  addPolylines(data=ht, color='black', weight = 8) %>%
  addScaleBar(position='bottomleft')


################ Vertical Analysis ##################

#Get elevation data to transform routes from xy to xyz format, recommended resolution is 10m
#Speed depends on 3DEP server status
routes_xyz <- get_elevation(route, 10)

#Call the vertical_analysis function.
#Speed depends on the length of the route.
vr <- vertical_analysis(routes_xyz) 

#Apply Curve and tangent classification criteria
vr <- vr %>% mutate(type=ifelse(r2>=0.5,"vc",
                            ifelse(r2<=0.3 & std<=1,"vt","u")))

#Draw map of vertical curves and tangents
leaflet() %>%
  addProviderTiles(providers$CartoDB.Positron) %>% 
  addPolylines(data=vr[vr$type=="vc",]%>%st_zm(), color='red', weight = 3) %>%
  addPolylines(data=vr[vr$type=="vt",]%>%st_zm(), color='black', weight = 3) %>%
  addScaleBar(position='bottomleft')

#draw vertical profile vs distance for a specific route
coords <- routes_xyz %>% st_coordinates()
len_tot <- st_length(routes_xyz) %>% as.numeric()
df <- data.frame(dist = (c(1:nrow(coords)) - 1) / (nrow(coords) - 1) * len_tot,
                 elev = coords[, 3])

plot(df,xlab="Distance (m)",ylab="Elevation (m)",type="l",col = 'blue',
     cex.lab = 1.5, cex.axis = 1.5,
     ylim=c(200,300))
abline(v = df[vr[-1,]$from,]$dist,lty = 'dashed')

for (i in 1:nrow(vr)){
  if (vr[i,]$type=='vc') {lines(x=c(df[vr[i,]$from,]$dist, df[vr[i,]$to,]$dist),y=c(207,207),col='red',lwd = 4)}
  if (vr[i,]$type=='vt') {
    lines(x=c(df[vr[i,]$from,]$dist, df[vr[i,]$to,]$dist),y=c(205,205),col='black',lwd = 4)
    text(x=df[vr[i,]$from,]$dist+100,y=222,paste('G = ', round(vr[i,]$avg,2),'%'),srt=90)}
}
  legend(6000,300,lty = 1,lwd = 4,col = c('red','black'),legend=c("Vertical curve", "Vertical tangent"),cex = 1.5)

  
