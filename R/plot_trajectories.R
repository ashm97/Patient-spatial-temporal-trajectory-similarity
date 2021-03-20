#' @title Plot pathways function
#'
#' @description Build ggplot visualisation for time components of trajectories
#'
#' @param traj.l List: list of trajectory dataframes
#'
#' @return GGplot visulisation
#'
#' @export
#'
#' @author Ashleigh C. Myall (\email{a.myall19@@imperial.ac.uk})
#
#  Copyright (C) 2020-2021 Ashleigh C. Myall
#  Licensed under the GNU General Public License v3.0


plot_trajectories = function(traj.l){

    plot.df = bind_rows(lapply(traj.l, function(x){
        x = x[,c(1,3)]
        return(data.frame("patient.ID"=x[1,1],"min"=min(x[,2]),"max"=max(x[,2]), stringsAsFactors=FALSE))
    }))

    plot.df$patient.ID = factor(plot.df$patient.ID,levels = plot.df$patient.ID[order(plot.df$min)] )

    # Plot
    p = ggplot(plot.df) +
        geom_segment( aes(x=patient.ID, xend=patient.ID, y=min, yend=max), color="black") +
        geom_point( aes(x=patient.ID, y=min), color=rgb(0.2,0.7,0.1,0.5), size=2 ) +
        geom_point( aes(x=patient.ID, y=max), color=rgb(0.7,0.2,0.1,0.5), size=2 ) +
        coord_flip()+
        theme_bw() +
        theme(
            legend.position = "none",
        ) +
        xlab("Patient") +
        ylab("time")

    return(p)
}
