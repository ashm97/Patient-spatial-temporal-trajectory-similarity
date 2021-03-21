#' @title Generate an example dataset
#'
#' @description A pre-built function with no arguements to return example trajectories. Each trajectory
#' contains a location (ward), and a corresponding time (date) on that location.
#'
#' @return Dataframe in long format with each row contiang the unique location-time
#' positions of each trajectory
#'
#' @export
#' @author Ashleigh C. Myall (\email{a.myall19@@imperial.ac.uk})
#  Copyright (C) 2020-2021 Ashleigh C. Myall
#  Licensed under the GNU General Public License v3.0

example_trajectories = function(){

    t_vec = seq(as.Date('2010-01-01'),as.Date('2010-10-26'),by = 1)

    # ------------------
    # Cluster 1

    patient_1 = data.frame("patient.ID" = rep("patient_1",15),
                           "location" = c(rep("ward_1",6),rep("ward_2",9)),
                           "t" = t_vec[1:15])

    patient_2 = data.frame("patient.ID" = rep("patient_2",40),
                           "location" = c(rep("ward_2",15),rep("ward_3",5),rep("ward_4",20)),
                           "t" = t_vec[1:40])

    patient_3 = data.frame("patient.ID" = rep("patient_3",45),
                           "location" = c(rep("ward_1",30),rep("ward_2",15)),
                           "t" = t_vec[14:58])

    patient_4 = data.frame("patient.ID" = rep("patient_4",30),
                           "location" = c(rep("ward_4",23),rep("ward_3",7)),
                           "t" = t_vec[41:70])

    patient_5 = data.frame("patient.ID" = rep("patient_5",20),
                           "location" = c(rep("ward_3",10),rep("ward_5",10)),
                           "t" = t_vec[66:85])

    patient_6 = data.frame("patient.ID" = rep("patient_6",53),
                           "location" = c(rep("ward_1",10),rep("ward_5",43)),
                           "t" = t_vec[20:72])


    # ------------------
    # Cluster 2


    patient_7 = data.frame("patient.ID" = rep("patient_7",53),
                           "location" = c(rep("ward_7",10),rep("ward_8",43)),
                           "t" = t_vec[90:142])

    patient_8 = data.frame("patient.ID" = rep("patient_8",53),
                           "location" = c(rep("ward_7",5),rep("ward_6",48)),
                           "t" = t_vec[90:142])


    # ------------------
    # Cluster 3

    patient_9 = data.frame("patient.ID" = rep("patient_9",53),
                           "location" = c(rep("ward_10",10),rep("ward_9",43)),
                           "t" = t_vec[110:162])

    patient_10 = data.frame("patient.ID" = rep("patient_10",53),
                            "location" = c(rep("ward_12",53)),
                            "t" = t_vec[115:167])

    patient_11 = data.frame("patient.ID" = rep("patient_11",53),
                            "location" = c(rep("ward_10",43),rep("ward_11",10)),
                            "t" = t_vec[70:122])

    patient_12 = data.frame("patient.ID" = rep("patient_12",53),
                            "location" = c(rep("ward_9",20),rep("ward_10",33)),
                            "t" = t_vec[100:152])

    # ------------------
    # Combine trajectories

    trajectories = rbind(patient_1,patient_2,patient_3,patient_4,patient_5,patient_6,
                         patient_7,patient_8,patient_9,patient_10,patient_11,patient_12)

    return(trajectories)

}
