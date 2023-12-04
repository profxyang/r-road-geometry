#function file of get_routes()
#get unique routes from an sfc object of road geometry
#the road geometry sfc can be extracted from OpenStreetMap
#the output of the function is an sfc with shapes of unique routes
#call this function for motorway and motorway links separately. 
#this function requires the stnetworks package.
library(sfnetworks)
library(tidygraph)
library(igraph)

get_routes <- function(sfc){ # sfc<-mtlk_line[,'lanes']
  
  
  # Step 1 convert the routes sfc to a sf network object
  sfc<-sfc[,'lanes'] # keep the data simple
  net <- as_sfnetwork(sfc,directed=T) 
  
  edge_colors = function(x) rep(sf.colors(12, categorical = TRUE)[-2], 2)[c(1:ecount(x))]
  plot(st_geometry(net, "edges"), col = edge_colors(net), lwd = 2)
  plot(st_geometry(net, "nodes"), pch = 20, cex = 1, add = TRUE)
  title(main='Original Network')
  
  #Step 2 processing the network
  #2.1 simplify - remove loop edges and multiple edges
  simp_net <- net %>%
    activate("edges") %>%
    filter(!edge_is_multiple()) %>%
    filter(!edge_is_loop())
  plot(st_geometry(sm_simp_net, "edges"), col = edge_colors(net), lwd = 2,title='simple')
  plot(st_geometry(sm_simp_net, "nodes"), pch = 20, cex = 1, add = TRUE)
  title(main='Simplified Network')
  #2.2 smooth - remove psudeo nodes (middle nodes)
  sm_simp_net<-convert(simp_net,to_spatial_smooth)
  plot(st_geometry(sm_simp_net, "edges"), col = edge_colors(net), lwd = 2)
  plot(st_geometry(sm_simp_net, "nodes"), pch = 20, cex = 1, add = TRUE)
  title(main='Smoothed Simplified Network')
  
  # find all the starting points (never showed up in the "to" column in the net)
  # and all the end points (never showed up in the "from" column in the net)
  sfc_edges<-st_as_sf(sm_simp_net, "edges")
  start_pts<-setdiff(sfc_edges$from,sfc_edges$to) 
  end_pts<-setdiff(sfc_edges$to,sfc_edges$from) 
  
  # find all simple paths form all "start_pts" to all "end_pts"
  # for some reason this only returns node paths
  n_paths = st_network_paths(sm_simp_net, from = start_pts[1], to=end_pts,type = "all_simple" )
  for (i in start_pts[-1]){
    pathi = st_network_paths(sm_simp_net, from = i, to=end_pts, type = "all_simple" )
    n_paths = rbind(n_paths,pathi)
  }
  
  
  #build the edge path from the node path (uniq_paths).
  find_epath<-function(array){
    ep<-st_network_paths(sm_simp_net, from=head(array,1),to=tail(array,1)) %>% pull(edge_paths) %>% unlist()
    return(sfc_edges[ep,] %>% st_union() %>% st_sf())
  }
  
  ls<-lapply(n_paths %>% pull(node_paths),find_epath)
  
  e_paths<-do.call(rbind,ls)
  
  cat(nrow(e_paths), " unique paths identified")
  
  return(e_paths)
}
