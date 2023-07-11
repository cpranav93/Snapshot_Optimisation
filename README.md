# Snapshot_Optimisation
A computationally efficient flow reconstruction technique, exploiting homogeneity, to recreate 3D instantaneous turbulent velocity fields from snapshots of 2D planar fields. This methodology, termed as ‘snapshot optimisation’ or SO, can help provide 3D data-sets from 2D data restricted by the limitations of experimental measurement techniques. For more details refer to the paper: 'Fast 3D flow reconstructions from 2D cross-plane observations' (accessible here: https://hal.science/hal-02012453/document)

Algorithms:
SO_Algorithm.f90: Classical snapshot optimisation algorithm for reconstructing 3D data-sets from 2D snapshots

Avg_SO_Algorithm.f90: Algorithm for the Averaging SO method 

POD_4_SO.f90: Reduced order formulation of the SO algorithm using proper orthogonal decomposition (POD)
