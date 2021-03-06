% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/comp_dist.R
\name{comp_dist}
\alias{comp_dist}
\title{Compute space and time distance matricies}
\usage{
comp_dist(traj_a, traj_b, D)
}
\arguments{
\item{traj_a}{dataframe: space-time locations of trajectory b.}

\item{D}{matrix: effective distance matrix, with elements the shortest path effective
distance weight between nodes i and j}
}
\value{
A list of ditance matricies
}
\description{
Given two trajectories of length m and n respectively, form 2 m x n matricies
corresponding to trajectory locations. One matrix dist_spatial, confers to
distance across a background movement graph, and another matrix dist_temporal
is the temporal distance between trajectory locations.
}
\author{
Ashleigh C. Myall (\email{a.myall19@imperial.ac.uk})
}
