# Description of all scripts used for `semla` analyses
  
General:  

* `custom_colors.R` – Definition of colors. Called within most scripts.  


Analysis:  

* `radial_dist_lung_fibrosis.qmd` – Extract radial distances from regions of interest (human and mouse data), later used in STUtility analysis.  
* `hs_IPF_alveolar_DEG_rdist_fibrosis.qmd` – Analysis of DEGs with regards to distance to fibrotic borders (Ext Data Fig. 2d).  
* `hsF14-C0_spatial_quantification.qmd` – Quantification of hsNMF-F14-C0 cluster spots in relation to areas of fibroblastic foci (Ext Data Fig. 3 d,e).
* `mmF14-C0_spatial_quantification.qmd` – Quantification of mmNMF_d21-F14-C0 cluster spots in relation to areas of fibrosis (Ext Data Fig. 4 f,g).  
* `hsNMF_visualizations.qmd` – Plot selected factors (Suppl Fig. 1a).  


Misc: 

* `hs_prepare_STable_DEA_distance_cor.R ` – Prepare output from `hs_IPF_alveolar_DEG_rdist_fibrosis.qmd` into a Supplementary Table (S. Table 3).  