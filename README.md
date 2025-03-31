# FC-SpatialHeterogeneity
Spatial Heterogeneity and Subtypes of Functional Connectivity Development in Youth

This repository contains the code and scripts used for our study investigating [Spatial Heterogeneity and Subtypes of Functional Connectivity (FC) Development in Youth](<https://doi.org/10.1101/2025.01.24.634828>).

Traditional brain age metrics provide a global overview, potentially masking important local variations. To address this, we introduce a regional brain development index for assessing spatially fine-grained FC maturation. ![FC-Figure](https://github.com/ML-DataAnalytics/FC-SpatialHeterogeneity/blob/main/Spatial%20Heterogeneity%20and%20Subtypes%20of%20Functional%20Connectivity%20Development%20in%20Youth.jpg) 

Using fMRI data from the Philadelphia Neurodevelopmental Cohort (PNC) and replicating in the Human Connectome Project Development (HCP-D) cohort, we: 1) Characterized spatial patterns of FC development across the cortex; 2) Identified distinct subtypes of individuals based on their unique spatial FC development patterns using clustering; 3) Demonstrated that subtypes reflecting advanced development along the sensorimotor-association axis show superior cognitive performance; and 4) Linked these developmental patterns to gene expression profiles enriched for neural differentiation, synaptogenesis, and myelination. This work highlights the importance of spatial heterogeneity in understanding neurocognitive development. 

The code provided here allows for the calculation of the regional development index, replication of the clustering analysis, and exploration of the links between FC development patterns, cognition, and gene expression.

## Brief descriptions of the code and scripts:

1. **regionwise_age_prediction**: scripts for region-wise age prediction based on the functional connectivity (FC) profile of each cortical region.

2. **regionwise_fc_dev_subtyping**: scripts for computing regional brain development (RBD) index maps and subtyping.

3. **gene_association_analysis**: scripts for association analysis between RBD patterns and gene data.

4. **replication_hcpd**: scripts for replication analysis of the RBD patterns obtained.


## Reference:

    @article{li2025spatial,
      title = {Spatial Heterogeneity and Subtypes of Functional Connectivity Development in Youth},
      author = {Li, Hongming and Cui, Zaixu and Cieslak, Matthew and Salo, Taylor and Moore, Tyler M and Gur, Raquel E and Gur, Ruben C and Shinohara, Russell T and Oathes, Desmond J and Davatzikos, Christos and Satterthwaite, Theodore D and Fan, Yong},
      journal = {bioRxiv},
      pages = {2025--01},
      year = {2025},
      publisher = {Cold Spring Harbor Laboratory},
      doi = {https://doi.org/10.1101/2025.01.24.634828}
    }
