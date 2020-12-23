## ---------------------------
##
## Script name: main.R
##
## Purpose of script: Construct a weighted sparsified graph from proxmity betweeen trajectories, measured across a background movement network.
##
## Author: Mr. Ashleigh C. Myall
##
## Date Created: 21-12-2020
##
## Copyright (c) Ashleigh C. Myall 2020
## Email: a.myall19@imperial.ac.uk
##
## ---------------------------
##
## Notes: Tdo: 1) link our preprint to function, 2) add node shapes as people for visualisation, 3) put functions in seperate file
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
library(readr)
library(igraph)
library(e1071)

source("./R/functions.R")


################################################################################

##########
##########                       Main Work Flow 
##########

################################################################################


##
##  0. Pre-process data
##


## 0.1 Load in data here / get exampled data
trajectories = example_trajectories()

## 0.2 Split long data into a list, with each list element a single patient dataframe
traj.l <- split(trajectories , f = trajectories$patient.ID)

## 0.3 Convert background movement graph to effective distances matrix D, containing 
# shortest directed path between nodes i -> j.
D = eff_dist(read_csv("data/background_movement.csv")) 



##
##  1. Compute proximity between all trajectories
##

edges = getSpatialTempProx(traj.l,   # list of trajectories
                 D,                  # efffective distance matrix 
                 beta = 0.6)         # paramter for speed of propergation



##
##  2. Sparse graph construction with CKNN
##

edges_cknn = cknneighbors_graph(k=3,             # parameter for k-nearest neighbors
                                #lambda = 1,     # data point density
                                edges = edges)   # fully connected edges


##
##  3. Visualisation of network
##

## 3.1 Preprocess data
netDat = preproNet(trajectories,edges_cknn)


## 3.2 Visualise network
visNetwork(nodes = netDat$nodes,edges = netDat$edges)%>% 
  visOptions(highlightNearest = TRUE)%>%
  visEdges(shadow = TRUE,
           color = list(color = "lightblue", highlight = "red"))%>% 
  visGroups(groupname = "B", shape = "icon", 
            icon = list(code = "f007", color = "red")) %>%
  addFontAwesome()

