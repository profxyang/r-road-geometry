#get unique routes from an sfc object of road geometry
#the road geometry sfc can be extracted from OpenStreetMap
#the output of the function is an sfc with shapes of unique routes
#call this function for motorway and motorway links separately. 
#this function requires the stnetworks package.
library(sfnetworks)
library(tidygraph)
get_routes <- function(sfc){
  
  
  net <- as_sfnetwork(sfc$geometry)
  smooth_net<-convert(net,to_spatial_smooth)
  clean_net <- convert(smooth_net, to_spatial_subdivision)
  
  sfc_edges <- st_as_sf(clean_net, "edges")
  start_pts <- setdiff(sfc_edges$from, sfc_edges$to)
  end_pts <- setdiff(sfc_edges$to, sfc_edges$from)
  
  
  n_paths = st_network_paths(clean_net, from = start_pts[1], to=end_pts,type = "all_simple" )
  for (i in start_pts[-1]) {
    pathi <- st_network_paths(clean_net, from = i, to = end_pts, type = "all_simple")
    n_paths <- rbind(n_paths, pathi)
  }
  
  find_epath <- function(array) {
    ep <- st_network_paths(clean_net, from = head(array, 1), to = tail(array, 1)) %>%
      pull(edge_paths) %>%
      unlist()
    return(sfc_edges[ep, ] %>% st_union() %>% st_sf())
  }
  
  ls <- lapply(n_paths %>% pull(node_paths), find_epath)
  
  e_paths <- do.call(rbind, ls)
  return(e_paths)
}
