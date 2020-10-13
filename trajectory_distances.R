## ---------------------------
##
## Script name: Compute Trajectory Similarity
##
## Purpose of script: 
##
## Author: Mr. Ashleigh C. Myall
##
## Date Created: 04-06-2020
##
## Copyright (c) Ashleigh C. Myall 2020
## Email: a.myall19@imperial.ac.uk
##
## ---------------------------
##
## Notes: 
##
## ---------------------------

################################################################################

### Load Libs

################################################################################


library(sna)
library(tsna)
library(ndtv)
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plyr)
library(readr)
library(igraph)
library(e1071)
library(gtools)


################################################################################

###
###                               Functions
###

################################################################################

# ------------------------------------------------------------------------------

## Function to measure distance between Trajectories and return mxn distance matrix

get_dist_list =  function(traj_a,traj_b,sp_m){
  
  # Check both df's have greater than 0 rows
  if(!(nrow(traj_a)>0 & nrow(traj_b)>0)){
    #print("Empty Trajectory")
    return(list("spatial"=NA,"temporal"=NA))
  }
  
  #dist_m = matrix(0, nrow = nrow(traj_a), ncol = nrow(traj_b))
  dist_spatial = matrix(0, nrow = nrow(traj_a), ncol = nrow(traj_b))
  dist_temporal = matrix(0, nrow = nrow(traj_a), ncol = nrow(traj_b))
  
  for (i in 1:nrow(traj_a)) {
    for (j in 1:nrow(traj_b)) {
      w_i = traj_a$raw_name[i]
      w_j = traj_b$raw_name[j]
      
      sp_indx_w_i = which(colnames(sp_m) == w_i)
      sp_indx_w_j = which(colnames(sp_m) == w_j)
      
      sp_ij = sp_m[sp_indx_w_i,sp_indx_w_j]
      sp_ji = sp_m[sp_indx_w_j,sp_indx_w_i]
      
      t_i = traj_a$WardDate[i]
      t_j = traj_b$WardDate[j]
      
      t_dif = abs(as.numeric(t_i - t_j))
      
      dist_spatial[i,j] = min(sp_ij,sp_ji)
      dist_temporal[i,j] = t_dif 
    }
  }
  #print("Dist computed")
  return(list("spatial"=dist_spatial,"temporal"=dist_temporal))
}

# ------------------------------------------------------------------------------

SP_Dist = function(sp_ij,sp_ji,t_dif,beta=0.5){
  return(min(sp_ij,sp_ji) + beta * t_dif)
}

# ------------------------------------------------------------------------------

SP_Dist_exp = function(space_d,time_d ,R=5,beta=7){
  
  sd_m = apply(space_d, 2, function(x,R){R^(-x)},R=R)
  
  t_m = apply(time_d, 2, function(x,beta){exp(-(x/beta))},beta=beta)
  
  D = sd_m * t_m
  
  return(sum(D))
}

# ------------------------------------------------------------------------------

## Space-time dist 2 from Mauricio

# y=e^{-space -beta*t}

# beta >= 0 : param
# space_d: distance matrix
# time-d: distance matrix

SP_Dist_exp_m = function(space_d,time_d,beta=0){
  D = exp(-1*abs(space_d) - abs(time_d)*beta)
  return(sum(D))
}


# ------------------------------------------------------------------------------

## Function to construct knn mst graph by taking the union if the two algo's

get_knn_mst_adj = function(m,k=3){
  neigh = knn_dist_mat(m,k)
  
  # Convert edge lists into new Adjancey matrix 
  edge_df = rbind(data.frame("Node1" = 1:nrow(neigh),
                             "Node2" = neigh[,1]))
  
  # For each Kth neighbor add the edge connections to df
  for (i in 1:k) {
    if(i > 1){
      edge_df = rbind(edge_df,
                      data.frame("Node1" = 1:nrow(neigh),
                                 "Node2" = neigh[,i]))
    }
  }
  
  K_adj = as.matrix(get.adjacency(graph.data.frame(edge_df)))
  colnames(K_adj) = colnames(m)
  row.names(K_adj) = row.names(m)
  
  # Do same for MST 
  mst.g <- graph.adjacency(as.matrix(m), mode="undirected", weighted=TRUE)
  mst <- minimum.spanning.tree(mst.g)
  
  # Take Union of two Adj
  mst_adj = as.matrix(get.adjacency(mst))
  knn_mst_adj = K_adj + mst_adj 
  knn_mst_adj[knn_mst_adj > 1] <- 1 # replace duplicated
  diag(knn_mst_adj) <- 0
  
  return(knn_mst_adj)
}

# ------------------------------------------------------------------------------

## Do KNN over pre computed distance matrix

# Take distrance matrix, k, and return matrix of K nearest neighbors

knn_dist_mat <- function(dist_mat,k){
  n = dim(dist_mat)[1]
  dist_mat = dist_mat + diag(x = 1, n) #replace diagnal elements with high num
  nn = matrix(0,n,k) # n x k
  for (i in 1:n)
    nn[i,] = k.nearest.neighbors(i, dist_mat, k = k)
  return(nn)
}

# ------------------------------------------------------------------------------

## Function to Get Effective Distance Shortest Path Matrix from edge list 

eff_dist_sp = function(edges){
  colnames(edges)[3] = "weight"
  
  # Shortest path matrix
  g <- graph.data.frame(edges)
  adj = get.adjacency(g, sparse = FALSE, attr='weight')
  
  # Convert to transistion matrix m by nromalising by indegree 
  degree_centralty <- function(x) rowSums(x) 
  m = adj/degree_centralty(adj)
  m[which(m == "NaN")] = 0
  
  # Convert to effective distance matrix
  D = 1 - log(m)
  D[which(D > 1000)] = 0
  D[which(D == 0)] = NA
  diag(D) <- 0
  
  sp_m = allShortestPaths(D)$length
  row.names(sp_m) = row.names(D)
  colnames(sp_m) = colnames(D)
  
  # Compute matrix of concat stirngs for wards on shortest path seperateed by ";"
  sp_ward_m = sp_m
  z = allShortestPaths(D)
  
  for (i in 1:length(row.names(D))) {
    for (j in 1:length(colnames(sp_m))) {
      if(!is.na(sp_ward_m[i,j])){
        sp_ward_m[i,j] = paste(colnames(D)[extractPath(z,i,j)],collapse = ";")
      }
    }
  }
  
  return(list("sp_m" =  sp_m, "sp_ward_m" = sp_ward_m))
  #return(sp_m)
}


################################################################################

##########
##########                       Main Work Flow 
##########

################################################################################


# ------------------------------------------------------------------------------

##
##  1. Load in Data
##

## load RDS

# CPE or Covid
CPE = T

if(CPE){
  patient_tr = readRDS("./Input Data/CPE/patient.list.clean.RDS")
}else{
  patient_tr = readRDS("./Input Data/Covid/covid_patient_traj.RDS")
  # Filter for 14 days prior to positive sample date (incubation period)
  patient_tr = lapply(patient_tr, function(x){
    x = x[order(x$WardDate),]
    x = x[max(1,(nrow(x))-14):nrow(x),]
  })
}

# Background Movement Graph + attributes
#background_CPE_g <- read_csv("Input Data/background_CPE_g.csv")
#background_CPE_g <- read_csv("./Input Data/Background Graphs/background_all_patients.csv")
background_CPE_g <- read_csv("./Input Data/Background Graphs/background_CPE_and_background_patients.csv")
#background_CPE_g <- read_csv("./Input Data/Background Graphs/background_CPE_patients.csv")


ward_attributes <- read_csv("Input Data/ward_attributes.csv")

# ------------------------------------------------------------------------------

##
##  2. Preprocess
##

## 2.1 Background Graph

sp_m = eff_dist_sp(background_CPE_g)$sp_m
sp_ward_m = eff_dist_sp(background_CPE_g)$sp_ward_m

## 2.2 Correct Mapping for Patient Wards

patient_tr = lapply(patient_tr, function(x){
  x= left_join(x, ward_attributes[,c(1,4)], by = c("ward"))
  x= na.omit(x)
  return(x)
})



# ------------------------------------------------------------------------------

##
##  2.3 Search for intersections
##


# Set all combinations of patient trajectories - duplicated combinations should not be included here
combs = combinations(n=length(names(patient_tr)),r=2,v=names(patient_tr),repeats.allowed=F)

#Progress bar
total <- nrow(combs)
pb <- txtProgressBar(min = 0, max = total, style = 3)

# Lapply over combinations to compute a time and space distance measure
dists.list = lapply(1:total, function(x,comb.df,patient_tr,sp_m,pb){
  setTxtProgressBar(pb, x)
  
  traj_i = patient_tr[[which(names(patient_tr) == comb.df[x,1])]]
  traj_j = patient_tr[[which(names(patient_tr) == comb.df[x,2])]]
  
  idx_ward = which(colnames(traj_i) == "raw_name")
  idx_time = which(colnames(traj_i) == "WardDate")
  
  dists = get_dist_list(traj_i[,c(idx_ward,idx_time)],traj_j[,c(idx_ward,idx_time)],sp_m=sp_m)
  
  return(list("traj_i" = comb.df[x,1],
              "traj_j" = comb.df[x,2],
              "space_d" = dists$spatial,
              "time_d" = dists$temporal))
  
},comb.df = combs,patient_tr,sp_m = sp_m,pb=pb)


# Save distance list
#saveRDS(dists.list,"./Output Data/CPE/distance_list_G.RDS")
saveRDS(dists.list,"./Output Data/CPE/distance_list_D.RDS")
#saveRDS(dists.list,"./Output Data/Covid/distance_list_G.RDS")
#saveRDS(dists.list,"./Output Data/Covid/distance_list_D.RDS")


## Lapply over distances to compute space-time distance metric
st_dist_l = lapply(dists.list, function(x){
  if(is.na(x$space_d) | is.na(x$time_d)){
    return(NA)
  }
  return(list("traj_i" = x$traj_i,
              "traj_j" = x$traj_j,
              "st_d" = SP_Dist_exp_m(x$space_d,x$time_d,beta=0.1)))
})


## Lapply and bind rows to output edges with weights
edges  = bind_rows(lapply(st_dist_l, function(x){
  if(!is.na(x)){
    return(data.frame("source_raw" = x$traj_i,"target_raw"=x$traj_j,"weight"= x$st_d))
  }else{
    return(data.frame("source_raw" = NA,"target_raw"=NA,"weight"= NA))
  }
}))

edges = na.omit(edges)
#edges = filter(edges, weight > 1)





# ------------------------------------------------------------------------------

##
##  2.3 Run Batch
##



folder = "./Output Data/CPE/"

uniq_id = "background_CPE_and_background_patients"

b_params = seq(0,2,0.05)

#mat.clus_coef =  matrix(0, nrow = length(b_params), ncol = length(R_params))
#colnames(mat.clus_coef) = R_params
#row.names(mat.clus_coef) = b_params

#mat.edges =  matrix(0, nrow = length(b_params), ncol = length(R_params))
#colnames(mat.edges) = R_params
#row.names(mat.edges) = b_params

dir.create(paste(folder,"Param_explore_",uniq_id,sep = ""))


for (b in b_params) {
  dir.create(paste(folder,"Param_explore_",uniq_id,"/Run_b",b,sep = ""))
  
  # Compute 
  # Lapply over distances to compute space-time distance metric
  ## Lapply over distances to compute space-time distance metric
  st_dist_l = lapply(dists.list, function(x){
    if(is.na(x$space_d) | is.na(x$time_d)){
      return(NA)
    }
    return(list("traj_i" = x$traj_i,
                "traj_j" = x$traj_j,
                "st_d" = SP_Dist_exp_m(x$space_d,x$time_d,beta=b)))
  })
  
  
  # Lapply and bind rows to output edges with weights
  edges  = bind_rows(lapply(st_dist_l, function(x){
    if(!is.na(x)){
      return(data.frame("source_raw" = x$traj_i,"target_raw"=x$traj_j,"weight"= x$st_d))
    }else{
      return(data.frame("source_raw" = NA,"target_raw"=NA,"weight"= NA))
    }
  }))
  
  edges = na.omit(edges)
  
  # joion exact overlaps..
  edges = left_join(edges, exact_intersects)
  
  # Save to file
  write.csv(edges,paste(folder,"Param_explore_",uniq_id,"/Run_b",b,"/edges.csv",sep = ""))
  
}










## ggplot output from clustering coef

plot.df = as.data.frame(t(mat.clus_coef))
plot.df$R_param = row.names(plot.df)

library(reshape2)

# Specify id.vars: the variables to keep but not split apart on
plot.df = melt(plot.df, id.vars=c("R_param"))
colnames(plot.df)[2:3] = c("Beta","value")

p1 = ggplot(plot.df,aes(x= R_param, y = value,group=Beta, color=Beta))+
  geom_line()+
  geom_point()+
  #scale_x_discrete(breaks = c(1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6))+
  theme_bw()+
  labs(title = "Clustering Coef over Param Range (CPE D")



plot.df = as.data.frame(t(mat.edges))
plot.df$R_param = row.names(plot.df)

library(reshape2)

# Specify id.vars: the variables to keep but not split apart on
plot.df = melt(plot.df, id.vars=c("R_param"))
colnames(plot.df)[2:3] = c("Beta","value")

p2 = ggplot(plot.df,aes(x= R_param, y = value,group=Beta, color=Beta))+
  geom_line()+
  geom_point()+
  #scale_x_discrete(breaks = c(1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6))+
  theme_bw()+
  labs(title = "Edges over Param Range (CPE D)")


library(ggpubr)
ggarrange(p1,p2,
          ncol = 1, nrow = 2)