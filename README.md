# Revealing missed transmission in hospital outbreaks of disease with a novel contact recovery mode

## Overview

We introduce a novel contact model which recovers undetected transmission of healthcare-associated infections (HAIs). Outbreaks of HAIs are both burdensome and extremely common. Contact tracing based on direct contacts is often used in outbreaks HAIs to prevent further spread. However, missing and indirect contacts pose severely limit contact tracing and result in misleading conclusions. 

 IMAGE 1

*The limitations missing data presents for outbreak investigations. Traditional contact tracing based on direct contact, a widely used tool in disease outbreak investigations, can be hampered by missing data. Both non-observable contacts between infected patients, or unknown infected patients serving as links between infected patients can result in misleading and limited results. Moreover, transmission facilitated through environmental surfaces, or by staff can all form an abundance of indirect and non-observable links.*

Here, we propose a graph-based model to mitigate these problems by measuring proximity between network-constrained temporal trajectories across background movement patterns. Our model naturally captures already known, but also missed contact likely to have resulted in transmission. 

IMAGE 2

*Model overview. Firstly, patient movement histories are captured as a set of network trajectories $\begin{Bmatrix} T_1,T_2,T_3,...,T_N \end{Bmatrix}$ passing through nodes of a background movement graph $D$ (Panel A). Typically, only direct contact is considered when determine transmission from infected patient movement, however, this inherently missed indirect contact which can also be a source of disease transmission. Our model recovers both indirect and direct contact by measuring points between trajectories in terms of the spatial proximity $\delta_{ij}$ (network-wise) and temporal proximity $\tau_{ij}$ (Panel B). This spatial-temporal proximity allows us to quantifies how close patients have into a mathematical object (a patient trajectory similarity matrix) $S$ (Panel C). Given disease portability of disease transmission increases the longer patients coincide, we recover the contact leading to transmission by forming a contact graph $\hat{S}$ by looking for strongest contact patterns (weighted by proximity) in $S$ (Panel D).

In our preprint (link) we showcased our models ability to capture transmission routes in an outbreak of HAI between 116 hospital patients, and demonstrated its real-time deployment on a cohort of 863 hospital patients. Using a semi-supervised learning framework, and bio-markers obtained from Whole Genome Sequencing, we showed that our model reveals missing patient interactions that improve characterisation of disease transmission.

If used please cite: *preprint link*


### Example data


```R
trajectories = example_trajectories()
```


### Visualising trajectories


### Computing proximity between trajectories


### Graph construction


### Visualising final contact network


## References
