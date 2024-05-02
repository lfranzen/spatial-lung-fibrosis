# DATA FOLDER README


## Visium data (human / mouse)

Download Visium Space Ranger output data from Mendeley and place within their respective folders, and into the `visium/` subfolder:
  
* human: Mendeley DOI 10.17632/nkbjsbmng6.1 or ArrayExpress acc S-BSST1410, 'hs_visium_spaceranger_output.zip'  
* mouse: Mendeley DOI 10.17632/5y6gd44x33.1 or ArrayExpress acc S-BSST1409, 'mm_visium_spaceranger_output.zip'  


## Cell2location results data

Download the cell2location results, unzip, and place in their respective folder (and into the `sc_deconvolution*/` subfolder):

* human, Habermann-based: Mendeley DOI 10.17632/nkbjsbmng6.1 or ArrayExpress acc S-BSST1410, 'cell2location_habermann2020.zip'  
* mouse, Strunz-based: Mendeley DOI 10.17632/5y6gd44x33.1 or ArrayExpress acc S-BSST1409, 'cell2location_strunz2020.zip'  


## Misc

The 'misc/' folder contains additional data files needed or desired for the downstream analysis. For instance; gene annotation tables from BioMart (used for Visium data filtering)), human-mouse orthologue gene conversion (fetched using 'orthogene' R package), and cell type descriptions of the data used for cell2location deconvolution.