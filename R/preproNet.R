#' @title Preprocess data for creating nodes and edges
#'
#' @description A function creating nodes and edges from input trajectory data.
#'
#' @param trajectories dataframe: full trajectory dataframe (long format)
#' @param edges dataframe: weighted undirected edge list for similarity between
#' trajectories.
#'
#' @return list: A list of nodes and edges for visualisation.
#' @export
#' @author Ashleigh C. Myall (\email{a.myall19@@imperial.ac.uk})
#  Copyright (C) 2020-2021 Ashleigh C. Myall
#  Licensed under the GNU General Public License v3.0

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
