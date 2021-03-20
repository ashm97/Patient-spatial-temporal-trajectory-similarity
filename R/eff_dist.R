#' @title Effective distances from background movement
#'
#' @description Return shortest paths between locations in terms of effective distance given
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
#' @export
#' @author Ashleigh C. Myall (\email{a.myall19@@imperial.ac.uk})
#  Copyright (C) 2020-2021 Ashleigh C. Myall
#  Licensed under the GNU General Public License v3.0

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
