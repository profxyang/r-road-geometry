#get unique routes from an sfc object of road geometry
#the road geometry sfc can be extracted from OpenStreetMap
#the output of the function is an sfc with shapes of unique routes
#call this function for motorway and motorway links separately. 
#this function requires the stnetworks package.
library(sfnetworks)

get_routes <- function(sfc){ #sfc<-mtlk_line
  net = as_sfnetwork(sfc,directed=T)
  
  # find all the starting points (never showed up in the "to" column in the net)
  sfc_edges<-st_as_sf(net, "edges")
  start_pts<-setdiff(sfc_edges$from,sfc_edges$to)
  
  paths = st_network_paths(net, from = start_pts[1], type = "all_simple" )
  
  for (i in start_pts[-1]){
    pathi = st_network_paths(net, from = i, type = "all_simple" )
    paths = rbind(paths,pathi)
  }
  
  #remove paths that are part of another path, resulting in a list of unique paths
  pick<-rep(1,nrow(paths))
  
  for (i in 1:nrow(paths)){
    a<-paths[i,]%>%pull(node_paths)%>%unlist() #the node path to be examined
    blist<-paths[-i,] #the node paths to be compared with a
    for (j in 1:nrow(blist)){
      b<-blist[j,]%>%pull(node_paths)%>%unlist()
      if (all(a %in% b)){pick[i]=0}
    }
  }
  
  uniq_paths<-paths[pick==1,]
  
  #build the edge path from the node path (uniq_paths). union edges.
  
  connect_path<-function(arry){
    df_npath<- data.frame(from=arry,to=lead(arry))%>%drop_na()
    ls<-apply(df_npath,1,function(x) sfc_edges%>%filter(from==x[1] & to==x[2]))
    return(do.call(rbind,ls)%>%st_union()) 
    #This union function sometimes returns multilinestring instead of a linestring
    #The reason is unknown.
    #perhaps is because some edge only has two nodes?
    
  }
  
  ls_uniq_path<-uniq_paths%>%pull(node_paths) #get the node path from the unique edge path
  
  ls_route<-lapply(ls_uniq_path,connect_path)
  
  routes <- st_sf(id = 1:nrow(uniq_paths), geometry = do.call(rbind,ls_route), crs = 4326)
  
  cat(nrow(routes),"routes identified.")
  
  #merge multilinestring routes to linestring
  clean_geom<-routes%>%st_geometry()
  for (i in 1:nrow(routes)){ 
    #this loop removes duplicate coordinates *next to each other* and convert multilinestring to linestring
    a<-routes[i,]%>%st_coordinates()
    clean_geom[i]<-a[diff(a[,1])!=0 & diff(a[,2])!=0,1:2]%>%st_linestring(dim = "XY")
  }
  #replace the geometry with all linestrings
  routes<-st_set_geometry(routes,clean_geom)
  
  return(routes)
}
