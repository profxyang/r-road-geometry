#This function requests elevation of routes from USGS at regular intervals
#The number of sample points is determined by the target resolution (tres)
#The input is a linestring of routes with XY coordinates
#The output is a linestring  of routes with XYZ coordinates

get_elevation <- function(route, tres) {
  
  #the total length of the route should be evaluated before tranformation
  tot_length <- st_length(route) %>% as.numeric()
  
  #first transform the route to XY crs for the st_line_sample function to work.
  route <- route %>% st_transform(3857)
  
  #a collection of points where elevations are to be evaluated.
  vector <- seq(from = 0, to = 1, by = 1 / (tot_length / tres))
  
  sample_pts <- st_line_sample(route, sample = vector)
  
  #convert the multipoint feature "sample_pts" to a simple feature collection of points
  pts <- sample_pts %>%
    st_transform(4326) %>%
    st_coordinates() %>%
    data.frame() %>%
    st_as_sf(coords = c("X", "Y"), crs = 4326)
  
  #get elevation from sample points pts. It will take some time.
  usgs_elev <-
    get_elev_point(location = pts, prj = 4326,src = "epqs")
  
  cat("Requested", nrow(usgs_elev), "elevation records from USGS.", 
      sum(is.na(usgs_elev)), "records missing.\n")
  
  usgs_elev$elevation <-
    imputeTS::na_interpolation(usgs_elev$elevation)
  
  #Interpolate missing data. If no missing data, this does nothing.
  
  #I have problem using st_zm function to add the Z coordinate. Here is a workaround.
  
  xyz <-
    cbind(pts %>% st_geometry() %>% st_coordinates(), 'Z' = usgs_elev$elevation) #XYZ matrix
  
  geometry <- st_linestring(xyz, dim = "XYZ") %>% st_sfc(crs = 4326)
  
  route2 <- st_set_geometry(route, geometry)
  
  return(route2)
}
