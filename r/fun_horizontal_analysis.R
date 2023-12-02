# This function identifies horizontal curves based on the angle (alpha) between route legs
# The input is a sfc, which is a long string.
# The output is a sfc object with a collection of linestring shapes.
# This function requires libraries "sf" and "tidyverse"

horizontal_analysis <- function(route){
  #Convert route to a table of X, Y using reference system of 3857 (google map)
  sf <- st_as_sf(route) %>% st_set_crs(4326) # Convert route to sf object

  #generate a table (data frame) of vertices coordinates (X-Y)
  tb <- sf %>%
    st_transform(3857) %>% 
    st_coordinates() %>%
    as.data.frame()
  
  #calculate the scale factor k, for measuring the true distance
  tb$k <- with (sf %>% st_coordinates() %>% as.data.frame(), cos(Y/180*3.14159))
  
  #detect and label vertex type : 0-tangent, 1-left, -1-right
  tb <- tb %>% 
    mutate(difxab = X - lag(X)) %>%
    mutate(difyab = Y - lag(Y)) %>%
    mutate(difxbc = lead(X) - X) %>%
    mutate(difybc = lead(Y) - Y) %>%
    mutate(lab = sqrt(difxab^2+difyab^2))%>%
    mutate(lbc = sqrt(difxbc^2+difybc^2))%>%
    mutate(lac = sqrt((lead(X)-lag(X))^2+(lead(Y)-lag(Y))^2))%>%
    mutate(rad2 = lab*lbc*lac/sqrt((lab+lbc+lac)*(-lab+lbc+lac)*(lab-lbc+lac)*(lab+lbc-lac))*k)%>%
    #The above formula find the radius of the road based on the leg lengths lab, lbc, and lac
    #https://www.mathopenref.com/trianglecircumcircle.html
    mutate(turn = difybc/difxbc-difyab/difxab) %>% #note that difxbc = 0 or difxab = 0 will not be a problem.
    mutate(type = if_else(rad2>5000, 0, sign(turn))) # R>5000m is used as a criterion for tangent sections 
  
  #merge very short "irregular" sections (with only one vertices) inside a curve
  # 1,0,1, is changed to 1,1,1,
  #-1,0,-1, is changed to -1,-1,-1
  #Also replace 1,-1,1 and -1,1,-1 with -1 and 1
  
  rolling_length <- rle(tb$type)
  
  rolling_length$values <- with(rolling_length, 
                                ifelse (is.na(values),values,
                                        ifelse(is.na(lag(values)) | is.na(lead(values)), values,
                                        ifelse(lag(values) == -1 & values == 0 & lead(values) == -1 & lengths ==1, -1, 
                                        ifelse(lag(values) == 1 & values == 0 & lead(values) == 1 & lengths ==1, 1, 
                                        ifelse(lag(values) == 1 & values == -1 & lead(values) == 1 & lengths ==1, 1,
                                        ifelse(lag(values) == -1 & values == 1 & lead(values) == -1 & lengths ==1, -1,
                                                values)))))))
  
  tb$type_mod <- inverse.rle(rolling_length)
  #remove very short curves, should I?
  # 0,1,0, is changed to 0,0,0,
  # 0,-1,0, is changed to 0,0,0,
  
  rolling_length <- rle(tb$type_mod)
  
  rolling_length$values <- with(rolling_length, 
                                ifelse (is.na(values),values,
                                        ifelse(lag(values) == 0 & values == 1 & lead(values) == 0 & lengths ==1, 0, 
                                               ifelse(lag(values) == 0 & values == -1 & lead(values) == 0 & lengths ==1, 0, 
                                                      values)))) 
 
   # regenerate the new vertex type column in tb
  tb$type_mod <- inverse.rle(rolling_length)
  
  
  # generate the segment list
  slist <- with(rle(tb$type_mod), data.frame(values, lengths, cumsum = cumsum(lengths))) %>%
    mutate(sec_start=lag(cumsum)+1) %>%       # calculate starting point
    mutate(sec_end=cumsum) %>%       # and ending point numbers of each segment
    drop_na()                       #  drop the first and last vertices
  
  # create a function to make geometry for each segment
  gen_geom <- function (a,b){
    cbind(tb[a:b,1],tb[a:b,2]) %>%
      st_linestring() %>% 
      st_geometry()
  }
  
  slist$geom <- mapply(gen_geom, slist$sec_start, slist$sec_end) #generate linestring geometry for each segment
  
  #create a sf object with geometry of curves. 
  #curves are determined from the segment list, with at least three curve vertices.
  curves <- slist[slist$values!=0 & slist$length>=3,] %>% st_as_sf() %>% st_set_crs(3857)

  ################################Calculate Radius and Center###########################
  #Calculate the geometric information (R and center) of the curve
  #curve fitting function - Modified Least Square Method
  #outputs are center coordinates, radius, standard deviation of radius, cov of radius
  
  MLS <- function (geom){
    ctb <- as.data.frame(st_coordinates(geom))[,1:2] #calculation table
    ctb <- ctb %>% mutate (X2=X^2) %>% mutate (Y2=Y^2)
    n <- length(ctb$X)
    A <- n*(n-1)*cov(ctb$X,ctb$X)
    B <- n*(n-1)*cov(ctb$X,ctb$Y)
    C <- n*(n-1)*cov(ctb$Y,ctb$Y)
    D <- 0.5*n*(n-1)*(cov(ctb$X,ctb$Y2)+cov(ctb$X,ctb$X2))
    E <- 0.5*n*(n-1)*(cov(ctb$Y,ctb$X2)+cov(ctb$Y,ctb$Y2))
    
    am <- (D*C-B*E)/(A*C-B^2)
    bm <- (A*E-B*D)/(A*C-B^2)
    r <- mean(sqrt((ctb$X-am)^2+(ctb$Y-bm)^2)) #At this point, scale factor is not applied to this radius
    sd <- sd(sqrt((ctb$X-am)^2+(ctb$Y-bm)^2))
    return(c(am,bm,r,sd,sd/r))
  }
  
  circles <- mapply(MLS,curves$geom)%>%t()%>%as.data.frame()
  colnames(circles) <- c('xc','yc','rad','sd','cov')
  curves<-cbind(curves,circles)     #Join the circle fitting results into the sf object.
  
  #Remove curves with radius > 5000m. There shouldn't be, but just in case.
  curves<- curves[curves$rad<=6100,] 
  
  ################################Final Touch###########################
  #Sometimes the point before PC or after PT also should be included int the circle.
  #To get the best estimate of the PC, PT and curve length, a final touch is needed.
  #Final Touch (FT) function. If the point before PC or after PT is less than 1m from the current circle, 
  #then include them in the curve.
  
  FT <- function (cur_start,cur_end, cur_xc,cur_yc,cur_r){
    xpt <- tb[cur_start-1,1]
    ypt <- tb[cur_start-1,2]
    start <- cur_start- ifelse(abs(sqrt((xpt-cur_xc)^2+(ypt-cur_yc)^2)-cur_r)<=1,1,0)
    xpt <- tb[cur_end+1,1]
    ypt <- tb[cur_end+1,2]
    end <- cur_end + ifelse(abs(sqrt((xpt-cur_xc)^2+(ypt-cur_yc)^2)-cur_r)<=1,1,0)
    cbind(start,end)
  }
  
  #Repeat the function for every row of data frame "curves".
  b <- mapply(FT, curves$sec_start,curves$sec_end, curves$xc,curves$yc,curves$rad)
  
  curves$sec_start2 <- b[1,] #the new curve segment start point
  curves$sec_end2 <- b[2,] #the new curve segment end point
  
  
  
  hc <- with(curves, data.frame(sec_start2, sec_end2))
  
  # rebuild geometry with modified starting and ending points
  hc$geom <- mapply(gen_geom, hc$sec_start2, hc$sec_end2) 
  
  hc <- hc %>% st_as_sf() %>% st_set_crs(3857) #convert to sf
  
  #recall MLS function for the new geometry
  a <- mapply(MLS,hc$geom)%>%t()
  
  hc$xc <-a[,1]
  hc$yc <-a[,2]
  hc$rad <-a[,3]
  hc$sd <-a[,4]
  hc$cov <-a[,5]
  
  hc<-hc%>%st_transform(4326)
  
  # apply scale factor to correct the length of r, 
  # the correction factor is calculated based on the center of the curve
  
  c<-with(hc, data.frame(xc,yc)) %>% 
    st_as_sf(coords = c('xc', 'yc'),crs=3857) %>%
    st_transform(4326) %>% 
    st_coordinates()
  hc$k <- cos(c[,2]/180*3.14159)
  hc$latc <- c[,2]
  hc$lngc <- c[,1]
  hc$rad <- hc$rad*hc$k
  
  hc$length <- st_length(hc)
  
  ################################Tangent Analysis###########################
  
  #build the tangent geometry, assuming all non-curve segments are tangents first.
  ht <- with(curves, data.frame(sec_start2, sec_end2)) %>% 
    mutate(t_start=ifelse(is.na(lag(sec_end2)),1,lag(sec_end2))) %>%
    mutate(t_end=sec_start2) 
  
  if (tail(ht$sec_end2,1)==nrow(tb)) {
    ht<-ht[c('t_start','t_end')]
  } else{
    ht<-ht[c('t_start','t_end')]
    ht[nrow(ht)+1,] = c(tail(curves$sec_end2,1),tail(nrow(tb)))
  }
  
  if (ht$t_end[1]==1){ht<-ht[-1,]} # remove the first tangent if the route starts with a curve
  
  
  ht$geom <- mapply(gen_geom, ht$t_start, ht$t_end)
  
  
  ht <- ht %>% st_as_sf() %>% st_set_crs(3857)
  
  #Measure tangents, calculate r2 and rmse. this should be done in the XY coordinate system
  
  rmse <- 0
  r2<-0
  maxe<-0
  
  for (i in 1:nrow(ht)){
    df<-ht[i,]%>%st_geometry()%>%st_coordinates()%>%data.frame()
    rmse[i]<-sqrt(mean((lm(Y~X,df[,1:2])$residuals)^2))
    r2[i]<-summary(lm(Y~X,df[,1:2]))$r.squared
    maxe[i]<-max(abs(lm(Y~X,df[,1:2])$residuals))
  }
  
  ht$rmse<-rmse
  ht$r2<-r2
  ht$maxe<-maxe
  ht<- ht %>% st_transform(4326)
  ht$length <- st_length(ht)
  
  return(list(hc[c('sec_start2','sec_end2','length','latc','lngc','rad','sd','cov','geom')],
              ht[c('t_start','t_end','length','r2','rmse','maxe','geom')]))
  
}
