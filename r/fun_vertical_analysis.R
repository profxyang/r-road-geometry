# This function performs outlier treatment and segmentation
# The input is the sf object which contains the elevation data of the route
# The function returns an sfc object with information of all segments .
# This function requires several libraries: "strucchange" "forecast" "sf"

vertical_analysis <- function(sf) {

  #convert sf to dataframe. This dataframe is for vertical profile
  df<-sf %>%
    st_transform(4326) %>%
    st_coordinates() %>% 
    data.frame()
  
  #add a distance (dist) column
  len_tot <- st_length(sf) %>% as.numeric()
  df$dist <- (c(1:nrow(df))-1)/(nrow(df)-1)*len_tot #distance from the first vertex of route
  
  #calculate and report actual data resolution
  res <- len_tot/(nrow(df)-1) #resolution in m
  
  #Find and replace outliers (bridges) based on both elevation and vertical grade G profiles
  otlrs_z <- tsoutliers(df$Z)
  g <- diff(df$Z)/res*100 # calculate G profile
  otlrs_g <- tsoutliers(g)
  otlrs <- c(otlrs_z$index-1,otlrs_g$index) %>% unique() %>% sort()
  g[otlrs] <- NA
  
  #plot outliers on vertical profile
  plot(df$dist,df$Z,type="o",col = 'blue', xlab="Distance (m)",  #Plot elevation profile
       ylab="Elevation (m)",cex.lab = 1.5, cex.axis = 1.5)
  points(df[otlrs+1,]$dist,df[otlrs+1,]$Z, cex = 2, col = 'red')
  
  #Data treatment
  g<-imputeTS::na_interpolation(g)
  
  # Segmentation for vertical alignment 
  df2 <- data.frame(x = df$dist[-1], y = g) #create a data frame of distance vs G. 
  bp<-breakpoints(y~x, data=df2, h=ceiling(100/res)-1)$breakpoints
  
  #Output segmentation result and plot
  if (sum(is.na(bp))>0){
    cat('1 segments detected!')
  } else {cat(length(bp)+1, 'segments detected!')}
  
  abline(v=df$dist[-1][bp],lty = 'dashed') #Show segmentation boundaries on the plot
  plot(df$dist[-1],g,type="o",col = 'blue', xlab="Distance (m)",  #Plot elevation profile
       ylab="Vertical Grade (%)",cex.lab = 1.5, cex.axis = 1.5)
  abline(v=df$dist[-1][bp],lty = 'dashed') #Show segmentation boundaries on the plot
  
  #Prepare the result data frame
  #The data frame starts with the from and to of the segments (in G profile)
  if (sum(is.na(bp))>0){
    result<-data.frame(from = 1,to=length(g))
  } else {result <- data.frame(from = c(1,bp), to = c(bp,length(g)+1))}
  
  
  #regression coefficients a and b, r2, std, avg
  result$a<-mapply(function (a,b){lm(y ~ x, data = df2[a:b,])$coefficients[[2]]}, result$from, result$to-1)
  result$b<-mapply(function (a,b){lm(y ~ x, data = df2[a:b,])$coefficients[[1]]}, result$from, result$to-1)
  result$r2<-mapply(function (a,b){summary(lm(y ~ x, data = df2[a:b,]))$r.square}, result$from, result$to-1)
  result$std<-mapply(function (a,b){sd(df2[a:b,]$y)}, result$from, result$to-1)
  result$avg<-mapply(function (a,b){mean(df2[a:b,]$y)}, result$from, result$to-1)
            
  #G1 and G2 are calculated from vertical curve segments
  result <- result %>% mutate(G1=df2$x[from]*a+b) %>% 
    mutate(G2=df2$x[to-1]*a+b) %>% 
    mutate(DG = abs(G2-G1)) %>%
    mutate(length = (to-from)*res)
  
  #Generate line geometries for the segments
  gen_geom <- function (a,b){
    points<-sf%>%st_cast("POINT")
    points[a:b] %>% st_coordinates() %>% st_linestring() %>% st_geometry()
  }
  
  result$geom <- mapply(gen_geom, result$from, result$to) #generate linestring geometry
  
  return(result %>% st_as_sf() %>% st_set_crs(4326))
  # The output dataframe "result" is a simple features collection (sfc) object which contains information 
  # a, b, r2, avg, std, G1, G2, DeltaG, Length
  # and the geometry of each segment.
  # The classification of vertical tangents and curves are done outside of this function so that users
        # can try different classification criteria without re-run the segmentation process.

}
