#' @title Compute full space and time between all trajectories
#'
#' @description Given the total set of trajectories, compute proxmity between all combinations.
#'
#' @param traj.l list: list of all trajectories
#' @param D matrix: effective distance matrix, with elements the shortest path effective
#' distance weight between nodes i and j
#' @param beta double: paramter regulating the effect of time in spatial-temporal
#' proximity. beta can also be intepreted as speed of proporgation across the
#' background mobility data. To pass when for the spatial temporal proximity kernal (matrix computation).
#'
#' @return weighted edge list of similarities between trajectories.
#' @export
#' @author Ashleigh C. Myall (\email{a.myall19@@imperial.ac.uk})
#  Copyright (C) 2020-2021 Ashleigh C. Myall
#  Licensed under the GNU General Public License v3.0

getSpatialTempProx = function(traj.l,D,beta=0.6){

    # User update
    print("Computing proximities")

    # Set all combinations of patient trajectories - duplicated combinations should not be included here
    combs = combinations(n=length(names(traj.l)),r=2,v=names(traj.l),repeats.allowed=F)

    #Progress bar
    total <- nrow(combs)
    pb <- txtProgressBar(min = 0, max = total, style = 3)

    # Lapply over combinations to compute a time and space distance measure
    dist.list = lapply(1:total, function(x,comb.df,patient_tr,D,pb,beta){
        setTxtProgressBar(pb, x)

        # Select trajectories and columns
        traj_i = patient_tr[[which(names(patient_tr) == comb.df[x,1])]]
        traj_j = patient_tr[[which(names(patient_tr) == comb.df[x,2])]]
        idx_ward = which(colnames(traj_i) == "location")
        idx_time = which(colnames(traj_i) == "t")

        # Compute space and time distances matrices
        dists = comp_dist(traj_i[,c(idx_ward,idx_time)],
                          traj_j[,c(idx_ward,idx_time)],
                          D)

        # Compute kernal distance
        kernal_sim = kernal_dist(dists$spatial,dists$temporal,beta)

        return(list("traj_i" = comb.df[x,1],
                    "traj_j" = comb.df[x,2],
                    "kernal_sim" = kernal_sim
        ))

    },comb.df = combs,traj.l,D = D,pb=pb,beta)



    # Lapply and bind rows to output edges with weights
    edges  = bind_rows(lapply(dist.list, function(x){
        if(!is.na(x)){
            return(data.frame("source" = x$traj_i,
                              "target"=x$traj_j,
                              "kernal_sim"= x$kernal_sim))
        }else{
            return(data.frame("source" = NA,
                              "target"=NA,
                              "kernal_sim"= NA))
        }
    }))

    edges$source = as.character(edges$source)
    edges$target = as.character(edges$target)

    return(edges)
}
