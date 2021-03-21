#' @title Spatial temporal proximity kernal (matrix computation)
#'
#' @description Computes total spatial temporal proximities using the kernal. For every pair
#' of patients, we define spatio-temporal proximity between ward-time locations
#' $l_{i}$ and $l_{j}$ with the kernel:
#'
#' \kappa(l_{i},l_{j})  = e^{-\delta_{ij}-\beta \tau_{ij}},
#'
#' where $\tau_{ij}=\begin{vmatrix}t_i - t_j\end{vmatrix}$, parameter $\beta$
#' represents a propagation speed, and $\delta_{ij}$ denotes the shortest-path
#' distance (the most probable pathway for disease propagation) between wards
#' $v_i$ and $v_j$ across the background movement network $G$.
#' Note the equation reaches a maximum of one when
#' $l_{ia} = l_{jb}$ (exact overlaps), and decays to zero as spatial-temporal
#' proximity becomes more distant. We then measure overall similarity between
#' trajectories $T_m$ and $T_n$ by summing over pairwise proximity measures
#' between $l_i \in T_n$ and $l_j \in T_m$:
#'
#' \mathcal{S}(T_n,T_m) =  \sum_{l_i \in T_n} \, \sum_{l_j \in T_m} \kappa(l_{i},l_{j}).
#'
#'
#' Full description found in XXXX LINK PAPER XXXXX
#'
#' @param beta double: paramter regulating the effect of time in spatial-temporal
#' proximity. beta can also be intepreted as speed of proporgation across the
#' background mobility data.
#' @param space_d matrix: pairwise space (network) distances between locations of
#' two trajectories.
#' @param time_d matrix: pairwise time distances between locations of
#' two trajectories.
#'
#' @return sum of spatial-temporal proximities
#' @export
#' @author Ashleigh C. Myall (\email{a.myall19@@imperial.ac.uk})
#  Copyright (C) 2020-2021 Ashleigh C. Myall
#  Licensed under the GNU General Public License v3.0

kernal_dist = function(space_d,time_d,beta=0){
    sp_prox = exp(-1*abs(space_d) - abs(time_d)*beta)
    return(sum(sp_prox))
}
