#' @title Compute space and time distance matricies
#'
#' @description Given two trajectories of length m and n respectively, form 2 m x n matricies
#' corresponding to trajectory locations. One matrix dist_spatial, confers to
#' distance across a background movement graph, and another matrix dist_temporal
#' is the temporal distance between trajectory locations.
#'
#' @param traj_a dataframe: space-time locations of trajectory a.
#' @param traj_a dataframe: space-time locations of trajectory b.
#' @param D matrix: effective distance matrix, with elements the shortest path effective
#' distance weight between nodes i and j
#'
#' @return A list of ditance matricies
#' @export
#' @author Ashleigh C. Myall (\email{a.myall19@@imperial.ac.uk})
#  Copyright (C) 2020-2021 Ashleigh C. Myall
#  Licensed under the GNU General Public License v3.0

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
