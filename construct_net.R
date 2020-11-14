## ---------------------------
##
## Script name: construct_net.R
##
## Purpose of script: Construct edge list from patient pathways (with example)
##
## Author: Mr. Ashleigh C. Myall
##
## Date Created: 14-11-2020
##
## Copyright (c) Ashleigh C. Myall 2020
## Email: a.myall19@imperial.ac.uk
##
## ---------------------------
##
## Notes: This script computes an edge list based on the pathways of patient movement
##        history. Here we demonstrate a small example which can handle either daily,
##        or hourly data.
##
## ---------------------------

################################################################################

### Load Libs

################################################################################


library(dplyr)
library(tidyr)
library(gtools)
library(ggplot2)
library(visNetwork)


################################################################################

###
###                               Functions
###

################################################################################

# ------------------------------------------------------------------------------

## Generate an example dataset (either hourly or daily times)

example_paths = function(gran = "hour"){
  
  # use either hour or days to measure time
  if(gran == "hour"){
    t_vec = seq(as.POSIXct("2010-02-22 23:00:00"), as.POSIXct("2010-02-28 08:00:00"), by="hour")
  }else{
    t_vec = seq(as.Date('2010-01-01'),as.Date('2010-04-26'),by = 1)
  }
  
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
  
  paths = rbind(patient_1,patient_2,patient_3,patient_4,patient_5,patient_6)
  
  return(paths)
  
}

# ------------------------------------------------------------------------------

## Plot pathways function

plot_paths = function(path.l){
  plot.df = bind_rows(lapply(path.l, function(x){
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


# ------------------------------------------------------------------------------

## Space-time distance function

# Given two trajectories, form an m x n matrix corresponding to trajectory locations.
# Then check which wards are identical, and enter the time difference between being on identical
# wards. Return the matrix, with 0 being exact overlaps, and NA's being different wards.

traj_dist =  function(traj_a,traj_b,gran = "days"){
  
  # Check both df's have greater than 0 rows
  if(!(nrow(traj_a)>0 | nrow(traj_b)>0)){
    print("Empty trajectory")
    return(NA)
  }
  
  dist_m = matrix(NA, nrow = nrow(traj_a), ncol = nrow(traj_b))
  
  for (i in 1:nrow(traj_a)) {
    for (j in 1:nrow(traj_b)) {
      w_i = traj_a$location[i]
      w_j = traj_b$location[j]
      if(w_i == w_j){
        #dist_m[i,j] = abs(as.numeric(difftime(traj_a$t[i] - traj_b$t[j])))
        dist_m[i,j] = abs(as.numeric(difftime(traj_a$t[i], traj_b$t[j], units = gran)))
      }
    }
  }
  
  return(dist_m)
}

# ------------------------------------------------------------------------------

## Get edge weights from space-time distance matrix

getWeightFromMat = function(m,t=0){
  m[which(m > t)] = NA #replace those outside of time frame
  if(t == 0){
    m[which(m == t)] = 1
    m[is.na(m)] = 0
    return(sum(m))
  }else if(t > 0){
    #Revere weights - higher weights = closer in time
    m = 1 + t - m 
    m[is.na(m)] = 0
    return(sum(m))
  }else{
    warning("Bad t value")
  }
}


# ------------------------------------------------------------------------------

## Preprocess data for nodes and edges

preproNet = function(paths,edges){
  nodes =  data.frame("id" = 1:length(unique(paths$patient.ID)),
                      "label" = unique(paths$patient.ID))
  
  edges.net = left_join(left_join(edges,nodes[,1:2],by= c("source" = "label")),
                        nodes[,1:2],by= c("target" = "label"))
  
  edges.net = edges.net[,c(4,5,3)]
  colnames(edges.net) = c("from","to","value")
  
  return(list("nodes"=nodes,"edges"=edges.net))
}

################################################################################

##########
##########                       Main Work Flow 
##########

################################################################################


# ------------------------------------------------------------------------------

##
##  0. Pre-process data
##


## 0.1 Load in data here / get exampled data (hourly "hour" data or daily "day" data)
gran = "hour"
#gran = "day"
paths = example_paths(gran)


## 0.2 Split long data into a list, with each list element a single patient dataframe
path.l <- split(paths , f = paths$patient.ID)



# ------------------------------------------------------------------------------

##
##  1. Visualise pathway
##

## 1.1 Visual the pathways with y capturing time and x capturing patients
p = plot_paths(path.l)
plot(p)


# ------------------------------------------------------------------------------

##
##  2. Compute list of distances matricies between trajectories
##


## 2.1 Set all combinations of patient trajectories
combs = combinations(n=length(names(path.l)),r=2,v=names(path.l),repeats.allowed=F)

## 2.2 Lapply over combinations to compute a time and space distance measure

#Progress bar
total <- nrow(combs)
pb <- txtProgressBar(min = 0, max = total, style = 3)

dist.list = lapply(1:total, function(x,comb.df,patient_tr,pb,gran){
  setTxtProgressBar(pb, x)
  traj_i = patient_tr[[which(names(patient_tr) == comb.df[x,1])]]
  traj_j = patient_tr[[which(names(patient_tr) == comb.df[x,2])]]
  return(list("traj_i"=comb.df[x,1],
              "traj_j"=comb.df[x,2],
              "dist_m"=traj_dist(traj_i,traj_j,gran)))
},comb.df = combs,path.l,pb=pb,gran = "hours")



# ------------------------------------------------------------------------------

##
##  3. Extract (weighted) edges
##

## 3.1 Here controls how close two patients can be in time for an edge be between them,
##     have a look how the netwrok below changes when you vary this paramter

t_diff = 0 # <-- Exact overlaps only
#t_diff = 1 # <-- +/-1 time unit
#t_diff = 7 # <-- +/-7 time unit

## 3.2 Compute edge weights

edges = bind_rows(lapply(dist.list, function(x,t){
  return(data.frame("source" = x$traj_i,"target" = x$traj_j,
                    "weight" = getWeightFromMat(x$dist_m,t)))
},t=t_diff))

## 3.3 Filter edge weights for those without weights (w > 0)
edges = filter(edges,weight > 0)

## 3.4 Save edges
#write.csv(edges,"edges.csv",row.names = F)


# ------------------------------------------------------------------------------

##
##  4. Visualisation of network
##

## 4.1 Preprocess data
netDat = preproNet(paths,edges)

## 4.2 Visualise network
visNetwork(nodes = netDat$nodes,edges = netDat$edges,  
           height = "500px",width = "100%")%>% 
  visOptions(highlightNearest = TRUE)

