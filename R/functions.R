################################################################################

###
###                               Functions
###

################################################################################


#' Generate example dataset
#'
#' Pre-built function with no arguements to return example trajectories. Each trajectory
#' contains a location (ward), and a corresponding time (date) on that location. 
#'
#' @return Dataframe in long format with each row contiang the unique location-time
#' positions of each trajectory


example_trajectories = function(){
  
  t_vec = seq(as.Date('2010-01-01'),as.Date('2010-10-26'),by = 1)
  
  # ------------------
  # Cluster 1
  
  patient_1 = data.frame("patient.ID" = rep("patient_1",15),
                         "location" = c(rep("ward_1",6),rep("ward_2",9)),
                         "t" = t_vec[1:15])
  
  patient_2 = data.frame("patient.ID" = rep("patient_2",40),
                         "location" = c(rep("ward_2",15),rep("ward_3",5),rep("ward_4",20)),
                         "t" = t_vec[1:40])
  
  patient_3 = data.frame("patient.ID" = rep("patient_3",45),
                         "location" = c(rep("ward_1",30),rep("ward_2",15)),
                         "t" = t_vec[14:58])
  
  patient_4 = data.frame("patient.ID" = rep("patient_4",30),
                         "location" = c(rep("ward_4",23),rep("ward_3",7)),
                         "t" = t_vec[41:70])
  
  patient_5 = data.frame("patient.ID" = rep("patient_5",20),
                         "location" = c(rep("ward_3",10),rep("ward_5",10)),
                         "t" = t_vec[66:85])
  
  patient_6 = data.frame("patient.ID" = rep("patient_6",53),
                         "location" = c(rep("ward_1",10),rep("ward_5",43)),
                         "t" = t_vec[20:72])
  
  
  # ------------------
  # Cluster 2
  
  
  patient_7 = data.frame("patient.ID" = rep("patient_7",53),
                         "location" = c(rep("ward_7",10),rep("ward_8",43)),
                         "t" = t_vec[90:142])
  
  patient_8 = data.frame("patient.ID" = rep("patient_8",53),
                         "location" = c(rep("ward_7",5),rep("ward_6",48)),
                         "t" = t_vec[90:142])
  
  
  # ------------------
  # Cluster 3
  
  patient_9 = data.frame("patient.ID" = rep("patient_9",53),
                         "location" = c(rep("ward_10",10),rep("ward_9",43)),
                         "t" = t_vec[110:162])
  
  patient_10 = data.frame("patient.ID" = rep("patient_10",53),
                          "location" = c(rep("ward_12",53)),
                          "t" = t_vec[115:167])
  
  patient_11 = data.frame("patient.ID" = rep("patient_11",53),
                          "location" = c(rep("ward_10",43),rep("ward_11",10)),
                          "t" = t_vec[70:122])
  
  patient_12 = data.frame("patient.ID" = rep("patient_12",53),
                          "location" = c(rep("ward_9",20),rep("ward_10",33)),
                          "t" = t_vec[100:152])
  
  # ------------------
  # Combine trajectories
  
  trajectories = rbind(patient_1,patient_2,patient_3,patient_4,patient_5,patient_6,
                       patient_7,patient_8,patient_9,patient_10,patient_11,patient_12)
  
  return(trajectories)
  
}

#' Plot pathways function
#'
#' Build ggplot visualisation for time components of trajectories
#' 
#' @param traj.l List: list of trajectory dataframes
#' 
#' @return GGplot visulisation


plot_trajectories = function(traj.l){
  
  plot.df = bind_rows(lapply(traj.l, function(x){
    x = x[,c(1,3)]
    return(data.frame("patient.ID"=x[1,1],"min"=min(x[,2]),"max"=max(x[,2]), stringsAsFactors=FALSE))
  }))
  
  plot.df$patient.ID = factor(plot.df$patient.ID,levels = plot.df$patient.ID[order(plot.df$min)] )
  
  # Plot
  p = ggplot(plot.df) +
    geom_segment( aes(x=patient.ID, xend=patient.ID, y=min, yend=max), color="black") +
    geom_point( aes(x=patient.ID, y=min), color=rgb(0.2,0.7,0.1,0.5), size=2 ) +
    geom_point( aes(x=patient.ID, y=max), color=rgb(0.7,0.2,0.1,0.5), size=2 ) +
    coord_flip()+
    theme_bw() +
    theme(
      legend.position = "none",
    ) +
    xlab("Patient") +
    ylab("time")
  
  return(p)
}


#' Effective distances from background movement
#'
#' Return shortest paths between locations in terms of effective distance given 
#' a background movement. Effective distance (as described in https://science.sciencemag.org/content/342/6164/1337),
#' 
  #' Routes of disease spread are often dominated by a set of most probable trajectories. 
  #' These probable trajectories can be derived using the effective distance given a
#' flux connectivity matrix $P$ from a systems mobility patterns. In $P$, elements 
#' $P_{ij}$ confer to the fraction of departing individuals leaving node (representing 
#' hospital wards in this study) $v_i$ and arriving at node $v_j$, such that $1 < P_{ij} \leq  1$. 
#' Given $P$ define effective distance $d_ij$ between as, 
#' 
#' d_{ij}= (1 - \log P_{ij} ) \leq 1.
#' 
#' In $d_ij$ a low proportion of movement from $ v_i \rightarrow v_j$ is corresponds 
#' to a large effective distance. Explained in greater detail in https://science.sciencemag.org/content/342/6164/1337, 
#' the logarithm captures effective distances are additive along multi-step pathways. 
#' Hence, deriving the most probable trajectories of disease spread from $ v_i \rightarrow v_j$ 
#' is given by the directed path $\delta_{ij} = \{v_i \rightarrow,...,\rightarrow v_j\}$ 
#' with the smallest total effective distance $\lambda ( \Gamma )$:
#' 
#' \hat{\delta_{ij}} =  \min_{ \lambda } \lambda \;( \delta_{ij} ).
#' 
#' Typically, $\hat{\delta_{ij}}  \neq \hat{\delta_{hi}}$, since patient movement 
#' patterns vary depending on procedures a given ward offers, and their starting 
#' locations. Also note, for paths a shortest path $\hat{\delta_{ij}}$ may not 
#' exist (i.e for wards which are purely entry or exit points in the hospital).
#' 
#' @param edges dataframe: directed and weighted edge dataframe capturing bacgkround 
#' movement between source nodes and target nodes.
#' 
#' @return effective distance matrix, with elements the shortest path effective 
#' distance weight between nodes i and j

eff_dist = function(edges){
  colnames(edges)[3] = "weight"
  
  # Construct adjancey matrix of systems mobility patterns
  g <- graph.data.frame(edges)
  adj = get.adjacency(g, sparse = FALSE, attr='weight')
  
  # Convert to transistion (weights fraction of node volume) matrix P by 
  # normalising by indegree 
  degree =  function(x) rowSums(x) 
  P = adj/degree(adj)
  P[which(P == "NaN")] = 0
  
  # Convert to effective distance matrix
  d_mat = 1 - log(P)
  d_mat[which(d_mat > 1000)] = 0
  d_mat[which(d_mat == 0)] = NA
  diag(d_mat) <- 0
  
  # Get shortest path matrix
  D = allShortestPaths(d_mat)$length
  row.names(D) = row.names(d_mat)
  colnames(D) = colnames(d_mat)
  
  return(D)
}


#' Spatial temporal proximity kernal (matrix computation)
#'
#' Computes total spatial temporal proximities using the kernal. For every pair 
#' of patients, we define spatio-temporal proximity between ward-time locations 
#' $l_{i}$ and $l_{j}$ with the kernel:
#' 
#' \kappa(l_{i},l_{j})  = e^{-\delta_{ij}-\beta \tau_{ij}},
#' 
#' where $\tau_{ij}=\begin{vmatrix}t_i - t_j\end{vmatrix}$, parameter $\beta$ 
#' represents a propagation speed, and $\delta_{ij}$ denotes the shortest-path 
#' distance (the most probable pathway for disease propagation) between wards 
#' $v_i$ and $v_j$ across the background movement network $G$. 
#' Note the equation reaches a maximum of one when 
#' $l_{ia} = l_{jb}$ (exact overlaps), and decays to zero as spatial-temporal 
#' proximity becomes more distant. We then measure overall similarity between 
#' trajectories $T_m$ and $T_n$ by summing over pairwise proximity measures 
#' between $l_i \in T_n$ and $l_j \in T_m$:
#' 
#' \mathcal{S}(T_n,T_m) =  \sum_{l_i \in T_n} \, \sum_{l_j \in T_m} \kappa(l_{i},l_{j}).
#' 
#' 
#' Full description found in XXXX LINK PAPER XXXXX
#' 
#' @param beta double: paramter regulating the effect of time in spatial-temporal
#' proximity. beta can also be intepreted as speed of proporgation across the 
#' background mobility data.
#' @param space_d matrix: pairwise space (network) distances between locations of
#' two trajectories.
#' @param time_d matrix: pairwise time distances between locations of
#' two trajectories.
#' 
#' @return sum of spatial-temporal proximities

kernal_dist = function(space_d,time_d,beta=0){
  sp_prox = exp(-1*abs(space_d) - abs(time_d)*beta)
  return(sum(sp_prox))
}


#' Compute space and time distance matricies
#'
#' Given two trajectories of length m and n respectively, form 2 m x n matricies 
#' corresponding to trajectory locations. One matrix dist_spatial, confers to 
#' distance across a background movement graph, and another matrix dist_temporal 
#' is the temporal distance between trajectory locations.
#' 
#' @param traj_a dataframe: space-time locations of trajectory a.
#' @param traj_a dataframe: space-time locations of trajectory b.
#' @param D matrix: effective distance matrix, with elements the shortest path effective 
#' distance weight between nodes i and j
#' 
#' @return list of ditance matricies

comp_dist =  function(traj_a,traj_b,D){
  
  # Check both df's have greater than 0 rows
  if(!(nrow(traj_a)>0 & nrow(traj_b)>0)){
    print("Missing trajectory data")
    return(list("spatial"=NA,"temporal"=NA))
  }
  
  # Matrix to save distances
  dist_spatial = matrix(NA, nrow = nrow(traj_a), ncol = nrow(traj_b))
  dist_temporal = matrix(NA, nrow = nrow(traj_a), ncol = nrow(traj_b))
  
  # Pairwise distance comparison
  for (i in 1:nrow(traj_a)) {
    for (j in 1:nrow(traj_b)) {
      l_i = traj_a$location[i]
      l_j = traj_b$location[j]
      
      # Check overlap on same location
      if(l_i %in% colnames(D) & l_j %in% colnames(D)){
        
        # Time distances
        t_i = traj_a$t[i]
        t_j = traj_b$t[j]
        dist_temporal[i,j] = abs(as.numeric(t_i - t_j))
        
        # Select elements representing bidirectional shortest path distances 
        sp_ij = D[which(colnames(D) == l_i),which(colnames(D) == l_j)]
        sp_ji = D[which(colnames(D) == l_j),which(colnames(D) == l_i)]
        
        # True/False condition of which path to choose - depenedent of location time ordering
        dist_spatial[i,j] = ifelse((t_i>t_j),sp_ji,ifelse((t_i>t_j),sp_ij,min(sp_ij,sp_ji)))
        
      }
    }
  }
  return(list("spatial"=dist_spatial,"temporal"=dist_temporal))
}


#' Compute full space and time between all trajectories
#'
#' Given the total set of trajectories, compute proxmity between all combinations.
#' 
#' @param traj.l list: list of all trajectories
#' @param D matrix: effective distance matrix, with elements the shortest path effective 
#' distance weight between nodes i and j
#' @param beta double: paramter regulating the effect of time in spatial-temporal
#' proximity. beta can also be intepreted as speed of proporgation across the 
#' background mobility data. To pass when for the spatial temporal proximity kernal (matrix computation).
#' 
#' @return weighted edge list of similarities between trajectories.


getSpatialTempProx = function(traj.l,D,beta=0.6){
  
  # User update
  print("Computing proximities")
  
  # Set all combinations of patient trajectories - duplicated combinations should not be included here
  combs = combinations(n=length(names(traj.l)),r=2,v=names(traj.l),repeats.allowed=F)
  
  #Progress bar
  total <- nrow(combs)
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  
  # Lapply over combinations to compute a time and space distance measure
  dist.list = lapply(1:total, function(x,comb.df,patient_tr,D,pb,beta){
    setTxtProgressBar(pb, x)
    
    # Select trajectories and columns
    traj_i = patient_tr[[which(names(patient_tr) == comb.df[x,1])]]
    traj_j = patient_tr[[which(names(patient_tr) == comb.df[x,2])]]
    idx_ward = which(colnames(traj_i) == "location")
    idx_time = which(colnames(traj_i) == "t")
    
    # Compute space and time distances matrices 
    dists = comp_dist(traj_i[,c(idx_ward,idx_time)],
                      traj_j[,c(idx_ward,idx_time)],
                      D)
    
    # Compute kernal distance
    kernal_sim = kernal_dist(dists$spatial,dists$temporal,beta)
    
    return(list("traj_i" = comb.df[x,1],
                "traj_j" = comb.df[x,2],
                "kernal_sim" = kernal_sim
    ))
    
  },comb.df = combs,traj.l,D = D,pb=pb,beta)
  
  
  
  # Lapply and bind rows to output edges with weights
  edges  = bind_rows(lapply(dist.list, function(x){
    if(!is.na(x)){
      return(data.frame("source" = x$traj_i,
                        "target"=x$traj_j,
                        "kernal_sim"= x$kernal_sim))
    }else{
      return(data.frame("source" = NA,
                        "target"=NA,
                        "kernal_sim"= NA))
    }
  }))
  
  edges$source = as.character(edges$source)
  edges$target = as.character(edges$target)
  
  return(edges)
}


#' Implementation of of Continuous k-Nearest Neighbors (CKNN) in R
#'
  #' See https://arxiv.org/pdf/1606.02353.pdf. CKNN connect points x,y 
  #' 
  #' if d(x,y) < lambda * ( d(x,x_k) * d(y,y_k) )^0.5
  #' 
  #' Where d(x,y) is the distance between points x and y. And x_k, and y_k, are the 
  #' k-th closest neighbors of points x and y respectively.
#' 
#' @param k integer: number of neighbors
#' @param lambda double: positive parameter regulating density.
#' @param edges dataframe: weighted undirected edge list for similarity between 
#' trajectories.
#' 
#' @return dataframe of edges with CKNN selected edges only.

cknneighbors_graph = function(k = 2,lambda = 1,edges){
  
  # ------------------
  # Internal function to check if cknn edge present between nodes
  
  check_cknn_edge = function(lambda,k,node_i,node_j,node_edges,all_edges){
    
    #Distance from i to j
    d_ij = edges_2$w[intersect(which(edges_2$source == node_i),
                               which(edges_2$target == node_j))]
    
    #Distances from i(j) to the k-th nearest neighbor of i(j) 
    d_i_ik = node_edges[[node_i]]$w[k]
    d_j_jk = node_edges[[node_j]]$w[k]
    
    # If true edge between nodes
    return(ifelse(d_ij < lambda * sqrt(d_i_ik*d_j_jk),1,0))
  }
  
  # ------------------
  # Step 1: Pre-process edge data
  
  # Add edge in both directions
  edges_2 = rbind(edges,data.frame("source" = edges$target,
                                   "target" = edges$source,
                                   "kernal_sim"= edges$kernal_sim))
  # Similarity to distance
  edges_2$w = max(edges$kernal_sim)-edges$kernal_sim 
  edges_2 = edges_2[,-3]
  
  
  # ------------------
  # Step 2: Ordered list of each nodes 1-step neighbors
  
  nodes = unique(c((edges_2$source),(edges_2$target)))
  node_edges = lapply(nodes,function(x,nodes){
    e = filter(edges_2, source == x)
    e = e[order(e$w),]
    rownames(e) <- NULL
    return(e)
  },nodes)
  names(node_edges) = nodes
  
  
  # ------------------
  # Step 3: check over edges
  
  edges$cknn_edge = 0
  for (i in 1:nrow(edges)) {
    edges$cknn_edge[i] = check_cknn_edge(lambda,k,edges$source[i],
                                         edges$target[i],
                                         node_edges,edges_2)
  }
  
  # ------------------
  # Step 4: Final edges
  
  ret.edges = filter(edges, cknn_edge == 1)
  
  return(ret.edges[,1:3])
}

#' Preprocess data for nodes and edges
#' 
#' @param trajectories dataframe: full trajectory dataframe (long format)
#' @param edges dataframe: weighted undirected edge list for similarity between 
#' trajectories.
#' 
#' @return list: list of nodes and edges for visualisation.

preproNet = function(trajectories,edges){
  nodes =  data.frame("id" = 1:length(unique(trajectories$patient.ID)),
                      "label" = unique(trajectories$patient.ID),
                      "group" = "B")
  
  edges.net = left_join(left_join(edges,nodes[,1:2],by= c("source" = "label")),
                        nodes[,1:2],by= c("target" = "label"))
  
  edges.net = edges.net[,c(4,5,3)]
  colnames(edges.net) = c("from","to","value")
  
  return(list("nodes"=nodes,"edges"=edges.net))
}
