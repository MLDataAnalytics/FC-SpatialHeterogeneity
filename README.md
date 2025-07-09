# FC-SpatialHeterogeneity: Spatial Heterogeneity and Subtypes of Functional Connectivity Development in Youth

[![DOI](https://img.shields.io/badge/DOI-10.1101%2F2025.01.24.634828-blue?style=for-the-badge&logo=biorxiv)](https://doi.org/10.1101/2025.01.24.634828)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) 

This repository contains the official code and scripts used for our study investigating **Spatial Heterogeneity and Subtypes of Functional Connectivity (FC) Development in Youth**. Our work challenges traditional global brain age metrics by introducing a spatially fine-grained approach to understanding brain maturation.

---

## ✨ Overview

While global brain age metrics offer a general understanding of brain development, they often mask critical local variations. Our research addresses this limitation by introducing a **regional brain development (RBD) index**, which provides a spatially fine-grained assessment of functional connectivity (FC) maturation across the cortex.

![FC-SpatialHeterogeneity Visualization](https://github.com/ML-DataAnalytics/FC-SpatialHeterogeneity/blob/main/Spatial%20Heterogeneity%20and%20Subtypes%20of%20Functional%20Connectivity%20Development%20in%20Youth.jpg)

Using functional Magnetic Resonance Imaging (fMRI) data from the **Philadelphia Neurodevelopmental Cohort (PNC)** for discovery and replicating our findings in the **Human Connectome Project Development (HCP-D) cohort**, we:

1.  **Characterized spatial patterns** of FC development across the cortex.
2.  **Identified distinct subtypes** of individuals based on their unique spatial FC development patterns using advanced clustering techniques.
3.  **Demonstrated superior cognitive performance** in subtypes reflecting advanced development along the sensorimotor-association axis.
4.  **Linked these developmental patterns to gene expression profiles** that are significantly enriched for neural differentiation, synaptogenesis, and myelination.

This work underscores the critical importance of spatial heterogeneity in gaining a deeper understanding of neurocognitive development.

The code provided in this repository enables users to:
* Calculate the regional development index.
* Replicate the clustering analysis to identify FC development subtypes.
* Explore the intricate links between FC development patterns, cognitive performance, and gene expression profiles.

---

## 📁 Repository Structure & Code Descriptions

This repository is organized into distinct modules, each containing scripts for specific parts of our analysis:

1.  **`regionwise_age_prediction/`**:
    Contains scripts for performing region-wise age prediction based on the functional connectivity (FC) profile of each cortical region. This forms the basis for the regional brain development index.

2.  **`regionwise_fc_dev_subtyping/`**:
    Includes scripts dedicated to computing the Regional Brain Development (RBD) index maps and subsequently performing the clustering analysis to identify developmental subtypes.

3.  **`gene_association_analysis/`**:
    Houses the scripts for conducting association analyses between the derived RBD patterns and relevant gene expression data.

4.  **`replication_hcpd/`**:
    Provides scripts specifically designed for the replication analysis of the RBD patterns and subtypes using the independent Human Connectome Project Development (HCP-D) cohort.

---

## 🚀 Getting Started

Detailed instructions for setting up the environment, preparing datasets, and running the code will be provided soon. Please check back for updates!

---

## 📝 Citation

If you find our work useful or inspiring in your research, please cite our study:

```bibtex
@article{li2025spatial,
  title = {Spatial Heterogeneity and Subtypes of Functional Connectivity Development in Youth},
  author = {Li, Hongming and Cui, Zaixu and Cieslak, Matthew and Salo, Taylor and Moore, Tyler M and Gur, Raquel E and Gur, Ruben C and Shinohara, Russell T and Oathes, Desmond J and Davatzikos, Christos and Satterthwaite, Theodore D and Fan, Yong},
  journal = {bioRxiv},
  pages = {2025--01},
  year = {2025},
  publisher = {Cold Spring Harbor Laboratory},
  doi = {[https://doi.org/10.1101/2025.01.24.634828](https://doi.org/10.1101/2025.01.24.634828)}
}
```

---

## 🙏 Acknowledgment

This project has been generously supported in part by the National Institutes of Health (NIH) through grants **U24NS130411** and **R01EB022573**. We are grateful for their support in making this research possible.
