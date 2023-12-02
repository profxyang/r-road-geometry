# This function performs outlier treatment and segmentation
# The input is the sf object which contains the elevation data of the route
# The function returns an sfc object with information of all segments .
# This function requires several libraries: "strucchange" "forecase" "sf"

vertical_analysis <- function(sf) { #sf<-a   ,    for testing

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
  cat('total length =', len_tot, 'm, elevation data extracted at', res,"m resolution.  ") 
  
  #Find and replace outliers (bridges) based on vertical grade G profile
  G <- diff(df$Z)/res*100 # calculate G provile
  otlrs <- tsoutliers(G) #Find outliers of G based on the forecast package
  
  #plot outliers on vertical profile
  plot(df$dist,df$Z,type="l",col = 'blue', xlab="Distance (m)",  #Plot elevation profile
       ylab="Elevation (m)",cex.lab = 1.5, cex.axis = 1.5)
  points(df[otlrs$index+1,]$dist,df[otlrs$index+1,]$Z, col = 'red') #Plot outliers identified.
  
  #The default outlier replacement algorithm is tsclean(). It is not good enough
  whole<-1:length(G)
  normal<-whole[-otlrs$index]
  cleanG<-G
  for (i in otlrs$index){
    leftindex<-mean(normal[normal<i]%>%tail(3))
    leftvalue<-mean(G[normal[normal<i]%>%tail(3)])
    rightindex<-mean(normal[normal>i]%>%head(3))
    rightvalue<-mean(G[normal[normal>i]%>%head(3)])
    cleanG[i]<-approx(c(leftindex,rightindex),c(leftvalue,rightvalue),i)$y
  }
  
  
  # Segmentation for vertical alignment 
  df2 <- data.frame(x = df$dist[-1], y = cleanG) #create a data frame of distance vs G. 
  
  #For long routes with more than 1000 data points, the "split-lapply-combine" method can save significant time.
  #The code below split the df dataframe into 1000 row sub data frames and apply segmt function and then merge results
  #Skip these lines if you want to run the segmentation algorithm on the route at once. 
  
  
  segmt <- function (dfs) {
      bps <- breakpoints(y~x, data = dfs, h=ceiling(100/res)) # Minimum length of segment is 100m, or 110m?
      dfs[bps$breakpoints,]%>%row.names%>%as.numeric()
  }

  if (len_tot>1000) {
  bp <- df2 %>% split(rep(1:ceiling(nrow(df)/1000),each=1000)[1:nrow(df2)]) %>% 
      lapply(segmt) %>% unlist() %>% array()
  } else{
  # If you want to run the segmentation algorithm on the route at once, un-comment the following line.
  bp<-breakpoints(y~x, data=df2, h=10)$breakpoints
  }
  
  #Output segmentation result and plot
  if (sum(is.na(bp))>0){
    cat('1 segments detected!')
  } else {cat(length(bp)+1, 'segments detected!')}
  
  abline(v=df$dist[-1][bp],lty = 'dashed') #Show segmentation boundaries on the plot
  plot(df2,xlab="Distance (m)",ylab="Vertical Grade (%)",cex.lab = 1.5, cex.axis = 1.5)
  abline(v=df$dist[-1][bp],lty = 'dashed') #Show segmentation boundaries on the plot
  
  #Prepare the result data frame
  #The data frame starts with the from and to of the segments (in G profile)
  if (sum(is.na(bp))>0){
    result<-data.frame(from = 1,to=length(G))
  } else {result <- data.frame(from = c(1,bp), to = c(bp,length(G)))}
  
  
  #regression coefficients a and b, r2, std, avg
  result$a<-mapply(function (a,b){lm(y ~ x, data = df2[a:b,])$coefficients[[2]]}, result$from, result$to)
  result$b<-mapply(function (a,b){lm(y ~ x, data = df2[a:b,])$coefficients[[1]]}, result$from, result$to)
  result$r2<-mapply(function (a,b){summary(lm(y ~ x, data = df2[a:b,]))$r.square}, result$from, result$to)
  result$std<-mapply(function (a,b){sd(df2[a:b,]$y)}, result$from, result$to)
  result$avg<-mapply(function (a,b){mean(df2[a:b,]$y)}, result$from, result$to)
            
  #G1 and G2 are calculated from vertical curve segments
  result <- result %>% mutate(G1=df2$x[from]*a+b) %>% 
    mutate(G2=df2$x[to]*a+b) %>% 
    mutate(DG = abs(G2-G1)) %>%
    mutate(length = (to-from)*res)
  
  #Generate line geometries for the segments
  gen_geom <- function (a,b){
    points<-sf%>%st_cast("POINT")
    points[a:b,]$geometry %>% st_coordinates() %>% st_linestring() %>% st_geometry()
  }
  
  result$geom <- mapply(gen_geom, result$from, result$to) #generate linestring geometry
  
  return(result %>% st_as_sf() %>% st_set_crs(4326))
  # The output dataframe "result" is a simple features collection (sfc) object which contains information 
  # a, b, r2, avg, std, G1, G2, DeltaG, Length
  # and the geometry of each segment.
  # The classification of vertical tangents and curves are done outside of this function so that users
        # can try different classification criteria without re-run the segmentation process.

}
