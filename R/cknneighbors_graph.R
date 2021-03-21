#' @title Implementation of of Continuous k-Nearest Neighbors (CKNN) in R
#'
#' @description See https://arxiv.org/pdf/1606.02353.pdf. CKNN connect points x,y
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
#' @export
#' @author Ashleigh C. Myall (\email{a.myall19@@imperial.ac.uk})
#  Copyright (C) 2020-2021 Ashleigh C. Myall
#  Licensed under the GNU General Public License v3.0

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
