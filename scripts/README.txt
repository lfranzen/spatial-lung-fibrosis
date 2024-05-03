# Description of all scripts used for `STUtility` analyses
  

General scripts:  

* `custom_colors.R` – Definition of colors. Called within most scripts.  
* `custom_functions.R` – Newly defined functions used in some other scripts.  
* `prep_suppl_data.R` – Generation of tables on NMF data for the supplement.  


Human data analysis ("hs_visium_*.R"):  

* `hs_visium_preprocessing_A.R` – Main processing workflow for human HC/IPF lung Visium data.  
* `hs_visium_nmf.R` – Main NMF analysis, deconvolving all data into 30 factors. (Ext Data Fig 3a)  
* `hs_visium_pseudobulk.R` – Generation of pseudo-bulk data from Visium, by grouping and pooling data in different ways. (Fig 1d)  
* `hs_visium_cell2location_proc_data.R` – Process and analyze cell2location results. (Fig 1f)  
* `hs_visium_nichenet.R` – NicheNet (CCC) analysis for F14hi-C0 and neighboring clusters. (Fig 3 g,h)  
* `hs_visium_radial_distance.R` – Analysis and visualisation of genes/cells correlating with radial distance from F14hi-C0. Distance vectors were extracted using `semla` prior to this analysis. (Fig 3c,d)  
* `hs_visium_preprocessing_B.R` – Alternative processing workflow for human HC/IPF Visium data, where each donor is processed separately.  
* `hs_visium_B_nmf.R` – Alternative NMF analysis performed within each donor.  
* `hs_visium_manus_figs.R` – Collection of scripts for generating various QC and manuscript figures. (Fig 1 c,e,f,g; Fig 2 a,b,c,d,e,f; Fig 3 a,b,e)  
* `IPA_heatmaps.R` – Heatmap visualisation of IPA results. (Fig 3f)

Mouse data analysis ("mm_visium_*.R"):  

* `mm_visium_preprocessing.R` – Main processing workflow for human mouse lung Visium data.  
* `mm_visium_nmf.R` – Main NMF analysis, deconvolving all data, the day 7 data, and the day 21 data separately. (Fig 5b) (Ext Data Fig 4e)  
* `mm_visium_pseudobulk.R` – Generation of pseudo-bulk data from Visium, by grouping and pooling data in different ways.  
* `mm_visium_cell2location_proc_data.R` – Process and analyze cell2location results. (Fig 4d-g)  
* `mm_visium_radial_distance.R` – Analysis and visualisation of genes/cells correlating with radial distance from F14hi-C0. Distance vectors were extracted using `semla` prior to this analysis. (Fig 5f)  
* `mm_visium_plots.R` – Small collection of scripts for generating some QC and manuscript figures. (Ext Data Fig 4 b,c)  


Translational data analysis ("hs_mm_visium_*.R"):  

* `hs_mm_visium_comparison_gene_conversion.R` – Generation of a combined Hs-Mm gene conversion table using `orthogene`.  
* `hs_mm_visium_comparison_DEA.R` – Collection of all Differential Expression Analyses performed on Hs and Mm pseudo-bulk data. (Fig 4b)  
* `hs_mm_visium_comparison_NMF.R` – Jaccard similarity calculations between top contributing NMF factor genes. (Fig 4h)  
* `hs_mm_visium_comparison_AbBa.R` – Translational comparison between Hs-F14-C0 and Mm-F14-C0 ("aberrant basaloid" clusters). (Fig 5 a,c; Ext Data Fig 4e)  
* `hs_mm_visium_comparison_AT2_trans_AT1_trajectory.R` – Slingshot trajectory analysis of AT2-AT1 intermediates. (Fig 5 g-h)  
* `hs_mm_visium_comparison_Mac.R` – Comparison of macrophage NMF factors between human IPF and mouse bleo data. (Ext Data Fig 5)  
* `hs_mm_visium_comparison_TLS.R` – Comparison of NMF factors containing tertiary lymphoid structure-like (or lymphocyte aggregate) signatures. (Fig 6 a-c)  
* `IPA_cnet_plots.R` - Network visualisation of IPA results. (Fig 5d)  