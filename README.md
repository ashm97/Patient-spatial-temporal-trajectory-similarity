# Revealing missed transmission in hospital outbreaks of disease with a novel contact recovery mode

## Paper overview

We introduce a novel contact model which recovers undetected transmission of healthcare-associated infections (HAIs). Outbreaks of HAIs are both burdensome and extremely common. Contact tracing based on direct contacts is often used in outbreaks HAIs to prevent further spread. However, missing and indirect contacts pose severely limit contact tracing and result in misleading conclusions. 

![alt text](https://raw.githubusercontent.com/ashm97/Patient-spatial-temporal-trajectory-similarity/main/images/fig_intro_missing_data.png | width=100) 

*The limitations missing data presents for outbreak investigations. Traditional contact tracing based on direct contact, a widely used tool in disease outbreak investigations, can be hampered by missing data. Both non-observable contacts between infected patients, or unknown infected patients serving as links between infected patients can result in misleading and limited results. Moreover, transmission facilitated through environmental surfaces, or by staff can all form an abundance of indirect and non-observable links.*

Here, we propose a graph-based model to mitigate these problems by measuring proximity between network-constrained temporal trajectories across background movement patterns. Our model naturally captures already known, but also missed contact likely to have resulted in transmission. 

![](https://raw.githubusercontent.com/ashm97/Patient-spatial-temporal-trajectory-similarity/main/images/methodology_explained.png) 

*Model overview. Firstly, patient movement histories are captured as a set of network trajectories $\begin{Bmatrix} T_1,T_2,T_3,...,T_N \end{Bmatrix}$ passing through nodes of a background movement graph $D$ (Panel A). Typically, only direct contact is considered when determine transmission from infected patient movement, however, this inherently missed indirect contact which can also be a source of disease transmission. Our model recovers both indirect and direct contact by measuring points between trajectories in terms of the spatial proximity $\delta_{ij}$ (network-wise) and temporal proximity $\tau_{ij}$ (Panel B). This spatial-temporal proximity allows us to quantifies how close patients have into a mathematical object (a patient trajectory similarity matrix) $S$ (Panel C). Given disease portability of disease transmission increases the longer patients coincide, we recover the contact leading to transmission by forming a contact graph $\hat{S}$ by looking for strongest contact patterns (weighted by proximity) in $S$ (Panel D).*

In our preprint (link) we showcased our models ability to capture transmission routes in an outbreak of HAI between 116 hospital patients, and demonstrated its real-time deployment on a cohort of 863 hospital patients. Using a semi-supervised learning framework, and bio-markers obtained from Whole Genome Sequencing, we showed that our model reveals missing patient interactions that improve characterisation of disease transmission.

If used please cite: *preprint link*


## Repo overview

This repo provides an implementable example of the model proposed in ref.  `R/` is the folder for scripts that contain R functions. All functions are documented with [roxygen2](https://roxygen2.r-lib.org/) syntax.

### Example data

There are 2 datasets required for our indirect contact model. Firstly, is the background mobility patterns over which indirect contact is measured. Secondly are the trajectories of infected individuals under investigation.

Background mobility patterns capture how diseases will spread since individuals and their person-person are a primary vector in transmission. Routes of disease spread are often dominated by a set of most probable trajectories (Brockmann and Helbing 2013). However, as we previously showed heterogeneous structure existed amongst movements patterns (Myall et al. 2020). Probable trajectories then must be derived based on background mobility which captures general movement patterns. Specifically, we capture these most probable trajectories `D` which contains the lengths of shortest paths from nodes v_i -> v_j in terms of effective distance (Brockmann and Helbing 2013). Here we can call the function `eff_dist()` which computes matrix `D` when provided a dataframe containing the weighted directed edges of background mobility patterns.

```R
D = eff_dist(read_csv("data/background_movement.csv")) 
```

The dataset `data/background_movement.csv` is a simple example with 12 locations (wards) which the subsequent `trajectories` are recorded over:

![alt text](https://raw.githubusercontent.com/ashm97/Patient-spatial-temporal-trajectory-similarity/main/images/background_movement_example.png) 


The function `example_trajectories()` generates a long dataframe with 12 pre-set trajectories. Each trajectory is a series of locations (rows in the dataframe) with both a location (ward or network nove v_i) and a time component. For later computation, we also split the `trajectories` dataframe into a list of dataframes `traj.l` based on the `trajectories$patient.ID` column.


```R
trajectories = example_trajectories()

traj.l <- split(trajectories , f = trajectories$patient.ID)
```

### Visualising trajectories

We provide a function `plot_trajectories()` to visualise the time that `trajectories` were recorded in `data/background_movement.csv`.

```R
plot_trajectories(traj.l)
```

*Trajectory summary. Each row refers to one individual (patient), and x-axis records time. The line plot shows the beggining of each trajectory (green marker), and the finishing time (red marker).*

### Computing proximity between trajectories

The function `getSpatialTempProx()` computes total spatial temporal proximities using our kernal function and returns a weighted undirected dataframe of `edges`. For every pair of patients, we define spatio-temporal proximity between ward-time locations $l_{i}$ and $l_{j}$ with the kernel:

\kappa(l_{i},l_{j})  = e^{-\delta_{ij}-\beta \tau_{ij}},

where $\tau_{ij}=\begin{vmatrix}t_i - t_j\end{vmatrix}$, parameter $\beta$ represents a propagation speed, and $\delta_{ij}$ denotes the shortest-path 
distance across `D` (the most probable pathway for disease propagation) between wards. We then measure overall similarity between trajectories $T_m$ and $T_n$ by summing over pairwise proximity measures between $l_i \in T_n$ and $l_j \in T_m$:

\mathcal{S}(T_n,T_m) =  \sum_{l_i \in T_n} \, \sum_{l_j \in T_m} \kappa(l_{i},l_{j}).

Full description found in *preprint link*

```R
edges = getSpatialTempProx(traj.l,   # list of trajectories
                 D,                  # efffective distance matrix 
                 beta = 0.6)         # paramter for speed of propergation
```

### Graph construction with Continuous *k*-nearest neighbors (Cknn)

Incorporating graph structure between data points can aid classification through the emphasis of strong relationships. Hence, to reveal stronger contacts capturing transmission between patients, we remove weak connections by sparsifying the corresponding fully connected graph made from the `edges` in `getSpatialTempProx()`. 

Although several alternative graph construction methods exist, we focus on Continuous *k*-nearest neighbors (Cknn) an extension to *k*-nearest neighbors (Berry and Timothy 2016). Implemented in the function `cknneighbors_graph`, it takes provide the fully connected `edges`, the paramter `k` for number of nearest neighbors to include, and an optinal paramter `lambda` (set by defaut to `lambda=1` if not provided), and returns a sparsified edges `edges_cknn` according to the Cknn algorithm: Cknn connect points x,y if d(x,y) < lambda * ( d(x,x_k) * d(y,y_k) )^0.5. Where d(x,y) is the distance between points x and y. And x_k, and y_k, are the k-th closest neighbors of points x and y respectively.

```R
edges_cknn = cknneighbors_graph(k=3,             # parameter for k-nearest neighbors
                                #lambda = 1,     # data point density
                                edges = edges)   # fully connected edges
```


### Visualising final contact network

```R
# Preprocess data
netDat = preproNet(trajectories,edges_cknn)

# Visualise network
visNetwork(nodes = netDat$nodes,edges = netDat$edges,  
           height = "500px",width = "100%")%>% 
  visOptions(highlightNearest = TRUE)
```

## References

<div id="refs" class="references">

<div id="ref-brockmann_2013">

Brockmann, Dirk, and Helbing, Dirk. 2013. “The hidden geometry of complex, network-driven contagion phenomena.” science 342.6164 (2013): 1337-1342. https://science.sciencemag.org/content/342/6164/1337

</div>

</div>

<div id="ref-myall_2020">
Myall, Ashleigh C and Peach, Robert L and Weiße, Andrea Y and Davies, Frances and Mookerjee, Siddharth and Holmes, Alison and Barahona, Mauricio. 2020. Network memory in the movement of hospital patients carrying drug-resistant bacteria. arXiv preprint. https://arxiv.org/abs/2009.14480v2

</div>

</div>

<div id="ref-myall_2020">
Berry, Tyrus, and Timothy Sauer. 2016. Consistent manifold representation for topological data analysis. Foundations of Data Science. http://aimsciences.org//article/id/2556e6c9-b4b9-455a-9d9e-886ef0cd166f

</div>
