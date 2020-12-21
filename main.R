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


# ------------------------------------------------------------------------------

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


# ------------------------------------------------------------------------------

##
##  1. Visualise pathway
##

## 1.1 Visual the pathways with y capturing time and x capturing patients
p = plot_trajectories(traj.l)
plot(p)


# ------------------------------------------------------------------------------

##
##  2. Compute list of distances matricies between trajectories
##

edges = getSpatialTempProx(traj.l,
                 D,
                 beta = 0.6)



# ------------------------------------------------------------------------------

##
##  3. Sparse graph construction with CKNN
##

edges_cknn = cknneighbors_graph(k=3,
                                lambda = 1,
                                edges = edges)


# ------------------------------------------------------------------------------

##
##  4. Visualisation of network
##

## 4.1 Preprocess data
netDat = preproNet(trajectories,edges_cknn)

## 4.2 Visualise network
visNetwork(nodes = netDat$nodes,edges = netDat$edges,  
           height = "500px",width = "100%")%>% 
  visOptions(highlightNearest = TRUE)

