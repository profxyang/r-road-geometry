# This function identifies horizontal curves based on the angle (alpha) between route legs
# The input is a sfc, which is a long string.
# The output is a sfc object with a collection of linestring shapes.
# This function requires libraries "sf" and "tidyverse"

horizontal_analysis <- function(route){ #route<-routes[4,]
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
    mutate(rad = lab*lbc*lac/sqrt((lab+lbc+lac)*(-lab+lbc+lac)*(lab-lbc+lac)*(lab+lbc-lac))*k)%>%
    #The above formula find the radius of the road based on the leg lengths lab, lbc, and lac
    #https://www.mathopenref.com/trianglecircumcircle.html
    mutate(turn = difybc/difxbc-difyab/difxab) %>% #note that difxbc = 0 or difxab = 0 will not be a problem.
    mutate(type = if_else(rad>5000, 0, sign(turn))) # R>5000m is used as a criterion for tangent sections
  
  #fix remove short irragular segments
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
 
   # generate the segment list
  slist <- with(rle(tb$type_mod), data.frame(values, lengths, cumsum = cumsum(lengths))) %>%
    mutate(sec_start=lag(cumsum)+1) %>%       # calculate starting point
    mutate(sec_end=cumsum) %>%       # and ending point numbers of each segment
    drop_na()                       #  drop the first and last vertices
  # create a function to make geometry for each segment
  # the function takes the begining vertex index a and the ending vertex index b
  # this geometry will be in 3857 reference system since tb is in XY coordinates
  # However a st_set_crs() function is needed for the sfc after the geometry is generated.
  gen_geom <- function (a,b){
    cbind(tb[a:b,1],tb[a:b,2]) %>%
      st_linestring() %>%
      st_geometry()
  }

  slist$geometry <- mapply(gen_geom, slist$sec_start, slist$sec_end) #generate linestring geometry for each segment

  #create a sf object with geometry of curves.
  #curves are determined from the segment list, with at least three curve vertices.
  curves <- slist[slist$values!=0 & slist$length>=3,] %>% st_as_sf() %>% st_set_crs(3857)
  #curve analysis is only necessary when there is at least one curve
  #otherwise the whole route should be sent to tangent analysis.
  if (nrow(curves) >= 1) {
    ################################Calculate Radius and Center###########################
    #Calculate the geometric information (R and center) of the curve
    #curve fitting function - Modified Least Square Method
    #outputs are center coordinates, radius, standard deviation of radius, cov of radius
    MLS <- function (geom) {
      ctb <- as.data.frame(st_coordinates(geom))[, 1:2] #calculation table
      ctb <- ctb %>% mutate (X2 = X ^ 2) %>% mutate (Y2 = Y ^ 2)
      n <- length(ctb$X)
      A <- n * (n - 1) * cov(ctb$X, ctb$X)
      B <- n * (n - 1) * cov(ctb$X, ctb$Y)
      C <- n * (n - 1) * cov(ctb$Y, ctb$Y)
      D <- 0.5 * n * (n - 1) * (cov(ctb$X, ctb$Y2) + cov(ctb$X, ctb$X2))
      E <- 0.5 * n * (n - 1) * (cov(ctb$Y, ctb$X2) + cov(ctb$Y, ctb$Y2))
      am <- (D * C - B * E) / (A * C - B ^ 2)
      bm <- (A * E - B * D) / (A * C - B ^ 2)
      r <-
        mean(sqrt((ctb$X - am) ^ 2 + (ctb$Y - bm) ^ 2)) #At this point, scale factor is not applied to this radius
      mpd <- mean(abs(sqrt((ctb$X - am) ^ 2 + (ctb$Y - bm) ^ 2) - r))
      return(c(am, bm, r, mpd, mpd / r))
    }
    
    #apply circle calculation function to all geometries in the sfc
    circles <- mapply(MLS, curves$geometry) %>% t() %>% as.data.frame()
    colnames(circles) <- c('xc', 'yc', 'rad', 'mpd', 'nmpd')
    curves <- cbind(curves, circles)     #Join the circle fitting results into the sf object.
    
    ################################Final Touch###########################
    #Sometimes the point before PC or after PT also should be included int the circle.
    #To get the best estimate of the PC, PT and curve length, a final touch is needed.
    #Final Touch (FT) function. If the point before PC or after PT is less than 1m from the current circle,
    #then include them in the curve.
    FT <- function (cur_start,
                    cur_end,
                    cur_xc,
                    cur_yc,
                    cur_r) {
      xpt <- tb[cur_start - 1, 1]
      ypt <- tb[cur_start - 1, 2]
      start <-
        cur_start - ifelse(abs(sqrt((xpt - cur_xc) ^ 2 + (ypt - cur_yc) ^ 2) -
                                 cur_r) <= 1, 1, 0)
      xpt <- tb[cur_end + 1, 1]
      ypt <- tb[cur_end + 1, 2]
      end <-
        cur_end + ifelse(abs(sqrt((xpt - cur_xc) ^ 2 + (ypt - cur_yc) ^ 2) - cur_r) <=
                           1, 1, 0)
      cbind(start, end)
    }
    
    #Repeat the function for every row of data frame "curves".
    b <- mapply(FT,
             curves$sec_start,
             curves$sec_end,
             curves$xc,
             curves$yc,
             curves$rad)
    
    curves$sec_start2 <- b[1, ] #the new curve segment start point
    curves$sec_end2 <- b[2, ] #the new curve segment end point
    
    #start to rebuild the curve sf
    hc <- with(curves, data.frame(sec_start2, sec_end2))
    # rebuild geometry with modified starting and ending points
    hc$geometry <- mapply(gen_geom, hc$sec_start2, hc$sec_end2)
    hc <- hc %>% st_as_sf() %>% st_set_crs(3857) #until now the 3857 crs has been assigned.
    #*But Is this necessary?
    #since the coordinates are still in X-Y format.
    #recall MLS function for the new geometry
    a <- mapply(MLS, hc$geom) %>% t()
    circles <- mapply(MLS, hc$geometry) %>% t() %>% as.data.frame()
    colnames(circles) <- c('xc', 'yc', 'rad', 'mpd', 'nmpd')
    hc <- cbind(hc, circles)
    hc <- hc %>% st_transform(4326) # this transformation is necessary to calculate k.
    # apply scale factor to correct the length of r, mpd, and nmpd,
    # the correction factor is calculated based on the center of the curve
    c <- with(hc, data.frame(xc, yc)) %>%
      st_as_sf(coords = c('xc', 'yc'), crs = 3857) %>%
      st_transform(4326) %>%   #This is to transform xc, yc coordinates to 4326.
      st_coordinates()
    hc$k <- cos(c[, 2] / 180 * 3.14159)
    hc$latc <- c[, 2]
    hc$lngc <- c[, 1]
    hc$rad <- hc$rad * hc$k
    hc$mpd <- hc$mpd * hc$k
    hc$nmpd <- hc$nmpd * hc$k
    hc$length <- st_length(hc)
    names(hc)[1] <- "from"
    names(hc)[2] <- "to"
    
    #at the end build the initial tangent geometry
    #assuming all non-curve segments are tangents first.
    ht <- with(hc, data.frame(from, to)) %>%
      mutate(t_start=ifelse(is.na(lag(to)),1,lag(to))) %>%
      mutate(t_end=from)
    if (tail(ht$to,1)==nrow(tb)) {
      ht<-ht[c('t_start','t_end')]
    } else{
      ht<-ht[c('t_start','t_end')]
      ht[nrow(ht)+1,] = c(tail(hc$to,1),tail(nrow(tb)))
    }
    if (ht$t_end[1]==1){ht<-ht[-1,]} # remove the first tangent if the route starts with a curve
    
    ht<-ht%>%filter(t_start<t_end) # this is to ensure proper geometry
    
    
    ht$geometry <- mapply(gen_geom, ht$t_start, ht$t_end)
    ht <- ht %>% st_as_sf() %>% st_set_crs(3857)
    names(ht)[1] <- "from"
    names(ht)[2] <- "to"
    #initial tangent geometry built
    # now sfc (ht) is ready for tangent analysis
    #create a coordinate check table.
    #This is for finding out the index of vertex in the original route.
    tab<-data.frame('X'=tb$X,'Y'=tb$Y,"N"=1:nrow(tb))
    
    tangent_analysis<-function(sf){
      pts<-st_simplify(sf,dTolerance = 5)%>% #Douglas-Peucker algorithm
        st_coordinates()%>%
        as.data.frame()
      df<-data.frame("fromX"=pts[-nrow(pts),1],
                     "fromY"=pts[-nrow(pts),2],
                     "toX"=pts[-1,1],
                     "toY"=pts[-1,2])
      from<-merge(x=df,y=tab, by.x=c("fromX","fromY"), by.y=c("X","Y"))
      to<-merge(x=df,y=tab, by.x=c("toX","toY"), by.y=c("X","Y"))
      return(cbind(from$N,to$N))
    }
    ht <- do.call(rbind.data.frame, mapply(tangent_analysis, ht$geometry,SIMPLIFY = F))
    colnames(ht)<-c("from","to")
    ht$geometry <- mapply(gen_geom, ht$from, ht$to) # build the tangent geometry
    ht <- ht %>% st_as_sf() %>% st_set_crs(3857) %>% st_transform(4326)
    ht <- ht %>% mutate(xc = NA) %>%
      mutate(yc = NA) %>%
      mutate(rad = NA) %>%
      mutate(mpd = NA) %>%
      mutate(nmpd = NA) %>%
      mutate(latc = NA) %>%
      mutate(lngc = NA) %>%
      mutate(k = NA)
    ht$length <- st_length(ht)
    for (i in 1:nrow(ht)){ #calculate mpd
      headlat<-head(ht$geometry[i]%>%st_coordinates(),1)[,2]
      ht$k[i]<- cos(headlat / 180 * 3.14159)
      ht$mpd[i]<-mean(st_distance(ht[i,]$geometry%>%st_cast("POINT"),st_simplify(ht[i,],dTolerance = 5)))*ht$k[i]
    }
    #combine h curves and h tangents
    hr<-rbind(hc,ht)
  } else { #apply tangent analysis to the whole route
    tab<-data.frame('X'=tb$X,'Y'=tb$Y,"N"=1:nrow(tb))
    tangent_analysis<-function(sf){
      pts<-st_simplify(sf,dTolerance = 5)%>% #Douglas-Peucker algorithm
        st_coordinates()%>%
        as.data.frame()
      df<-data.frame("fromX"=pts[-nrow(pts),1],
                     "fromY"=pts[-nrow(pts),2],
                     "toX"=pts[-1,1],
                     "toY"=pts[-1,2])
      from<-merge(x=df,y=tab, by.x=c("fromX","fromY"), by.y=c("X","Y"))
      to<-merge(x=df,y=tab, by.x=c("toX","toY"), by.y=c("X","Y"))
      return(cbind(from$N,to$N))
    }
    ht<-tangent_analysis(route%>%st_transform(3857))%>%as.data.frame()
    colnames(ht)<-c("from","to")
    ht$geometry <- mapply(gen_geom, ht$from, ht$to)
    ht <- ht %>% st_as_sf() %>% st_set_crs(3857)
    ht <- ht %>% mutate(xc = NA) %>%
      mutate(yc = NA) %>%
      mutate(rad = NA) %>%
      mutate(mpd = NA) %>%
      mutate(nmpd = NA) %>%
      mutate(latc = NA) %>%
      mutate(lngc = NA) %>%
      mutate(k = NA) %>%
      st_transform(4326)
    ht$length <- st_length(ht)
    for (i in 1:nrow(ht)){
      ht$mpd[i]<-mean(st_distance(ht[i,]$geometry%>%st_cast("POINT"),st_simplify(ht[i,],dTolerance = 5)))
    }
    #pass tanget segments to results.
    hr<-ht
  }
  return(hr)

}
