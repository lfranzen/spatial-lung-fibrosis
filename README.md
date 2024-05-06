# *Mapping spatially resolved transcriptomes in human and mouse pulmonary fibrosis*

Lovisa Franzén<sup>§</sup>, Martina Olsson Lindvall<sup>§</sup>, Michael Hühn, Victoria Ptasinski, Laura Setyo, Benjamin P Keith, Astrid Collin, Steven Oag, Thomas Volckaert, Annika Borde, Joakim Lundeberg, Julia Lindgren, Graham Belfield, Sonya Jackson, Anna Ollerstam, Marianna Stamou<sup>$</sup>, Patrik L Ståhl<sup>$</sup>, Jorrit J Hornberg

§ These authors contributed equally to the work  
$ Corresponding authors

*bioRxiv* DOI:[https://doi.org/10.1101/2023.12.21.572330]


## Description

This repository contains the R code used to produce all the analyses and figures presented in the article *Mapping spatially resolved transcriptomes in human and mouse pulmonary fibrosis* (Franzén & Olsson Lindvall et al.).

The underlying data used for the analyses have been deposited to ArrayExpress (human data, accession [S-BSST1410](https://www.ebi.ac.uk/biostudies/studies/S-BSST1410); mouse data, accession [S-BSST1410](https://www.ebi.ac.uk/biostudies/studies/S-BSST1410). Space Ranger output found within the zipped files in folders named "V*****-***-*1". To generate these files, raw FastQ files from the NovaSeq sequencing were processed with the Space Ranger pipeline (v. 1.2.2, 10x Genomics), where the reads were mapped to the GRCh38 (human) or mm10 (mouse) reference genome. Manual spot alignment was performed in the Loupe Browser (v. 6, 10x Genomics) software.

Cell type mapping results were obtained using the cell2location (v. 0.1) method, integrating the Space Ranger output data with annotated single cell RNA-seq data produced from human IPF lung, published by [Habermann et al. (2020)](https://doi.org/10.1126/sciadv.aba1972) (GEO accession: GSE135893), or mouse bleomycin-injured lungs, published by [Strunz et al. (2020)](https://doi.org/10.1038/s41467-020-17358-3) (GEO accession: GSE141259).

Seurat/STUtility object was generated from the Space Ranger output files, using the R packages [STUtility](https://github.com/jbergenstrahle/STUtility) (v. 1.1.1) and Seurat (v. 4.1.1) in R (v. 4.0.5) or using [*semla*](https://github.com/ludvigla/semla) (v. 1.1.6) and Seurat (v. 4.3.0.1) in R (v. 4.2.3).

The data processing workflows have been illustrated [here](https://github.com/lfranzen/spatial-lung-fibrosis/blob/master/doc/analysis_workflow_schematic_HsMm.pdf).


<p align="center"><img src="/doc/graphical_abstract.png" height="400" width=400"></p>


## Content

* `bin/`: Installation files for R packages `NNLM`  
* `data/`: Follow instructions in the `README.txt` file placed within this folder to download and populate these folders with the relevant input data  
  * `human/`  
    * `visium/` 
    * `sc_deconvolution_habermann/`  
  * `mouse/`  
    * `visium/` 
    * `sc_deconvolution_strunz/` 
  * `misc`: BioMart gene annotation tables, Orthogene mouse-human conversion tables, and single cell annotation groups  
* `scripts/`: All main R scripts used for processing and analyzing the Visium data  
* `results/`: Directory for saving analysis output objects and figures  
* `doc/`: Contains schematic workflow overview for the analyses   
* `semla_analysis/`: All files used for analyses performed using *semla*    
  * `bin/`: Installation files for R packages `spatstat` and `NNLM`  
  * `data/`: Visium object metadata used for input in semla analyses, together with Visium Seurat objects (needs to be downloaded separately)  
  * `scripts/`: All scripts used for analyzing the Visium data with *semla*


## Get started

If you intend to replicate some of the analyses presented here, it is recommended to do the following once you have forked/downloaded the content of this repo:

1. Move the `semla_analysis` folder and content to another location (this analysis path has its own R project file and package requirements)  
2. Download the spaceranger output data and/or STUtility/Seurat objects and place them in their correct locations (see `data/README.txt`)  
3. Install [`renv`](https://rstudio.github.io/renv/articles/renv.html) and install the necessary packages from the `renv.lock` file by using `renv::restore()`
4. Manually install the packages within the `bin/` folder(s) and [STUtility v. 1.1.1](https://github.com/jbergenstrahle/STUtility/releases/tag/1.1.1)  
5. Read the `scripts/README.txt` for information about all scripts and their content


For details on the analysis steps, please see the overview schematic flowcharts below:

<p align="center"><img src="/doc/analysis_workflow_schematic_Hs.png" width=400"></p>

<p align="center"><img src="/doc/analysis_workflow_schematic_Mm.png" width=400"></p>


## Contact

For questions related to this repo and its content, please contact Lovisa Franzén (lovisa.franzen@scilifelab.se)


