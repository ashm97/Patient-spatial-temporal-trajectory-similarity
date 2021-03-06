% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preproNet.R
\name{preproNet}
\alias{preproNet}
\title{Preprocess data for creating nodes and edges}
\usage{
preproNet(trajectories, edges)
}
\arguments{
\item{trajectories}{dataframe: full trajectory dataframe (long format)}

\item{edges}{dataframe: weighted undirected edge list for similarity between
trajectories.}
}
\value{
list: A list of nodes and edges for visualisation.
}
\description{
A function creating nodes and edges from input trajectory data.
}
\author{
Ashleigh C. Myall (\email{a.myall19@imperial.ac.uk})
}
