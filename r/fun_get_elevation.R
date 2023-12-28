# This function requests elevation of routes from USGS at regular intervals
# The number of sample points is determined by the target resolution (tres)
# The input is a linestring with XY coordinates
# The output is a linestring with XYZ coordinates
library(elevatr)

get_elevation <- function(route,tres) { 
  #Determine elevation sample points
  tot_length <- st_length(route) %>% as.numeric()
  vector <- seq(from = 0, to = 1, by = 1 / (tot_length / tres))
  pts <- st_line_sample(route %>% st_transform(3857), sample = vector) %>%
    st_transform(4326) %>%
    st_cast("POINT")
    
  #get elevation from sample points'pts'. Be patient.
  usgs_elev <- get_elev_point(pts %>% st_as_sf(), prj = 4326, src = "epqs")
  cat(
    "Requested", nrow(usgs_elev), "elevation records from USGS.",
    sum(is.na(usgs_elev)), "records missing.\n"
  )
  usgs_elev$elevation <- imputeTS::na_interpolation(usgs_elev$elevation)
  
  #build the XYZ linestring  
  route3d <- cbind(pts %>% st_geometry() %>% st_coordinates(),
                   "Z" = usgs_elev$elevation) %>%
              st_linestring(dim = "XYZ") %>%
              st_sfc() %>% 
              st_set_crs(4326)
  return(route3d)
}
