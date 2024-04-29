# *Mapping spatially resolved transcriptomes in human and mouse pulmonary fibrosis*

Lovisa Franzén<sup>§</sup>, Martina Olsson Lindvall<sup>§</sup>, Michael Hühn, Victoria Ptasinski, Laura Setyo, Benjamin P Keith, Astrid Collin, Steven Oag, Thomas Volckaert, Annika Borde, Joakim Lundeberg, Julia Lindgren, Graham Belfield, Sonya Jackson, Anna Ollerstam, Marianna Stamou<sup>$</sup>, Patrik L Ståhl<sup>$</sup>, Jorrit J Hornberg

§ These authors contributed equally to the work  
$ Corresponding authors

*bioRxiv* DOI:[https://doi.org/10.1101/2023.12.21.572330]

## Description

This repository contains the R code used to produce all the analyses and figures presented in the article *Mapping spatially resolved transcriptomes in human and mouse pulmonary fibrosis* (Franzén & Olsson Lindvall et al.).


<p align="center"><img src="/doc/graphical_abstract.png" height="400" width=400"></p>

## Content
* `bin/`: Installation files for R packages `devtools`, `STUtility`, and `NNLM`  
* `data/`: Follow instructions in README.txt file placed here to download and populate these folders with the relevant input data  
  * `visium/`  
    * `human/`  
    * `mouse/`  
  * `misc`  
* `scripts/`: All main R scripts used for processing and analysing the Visium data  
* `results/`: Directory for saving analysis output objects and figures  
* `doc/`: Contains schematic workflow overview for the analyses   
* `semla_analysis/`: All files used for analyses performed using *semla*    
  * `bin/`  
  * `data/`  
  * `scripts/`
  * `results/`  

## Contact

For questions related to this repo and its content, please contact Lovisa Franzén (lovisa.franzen@scilifelab.se)


