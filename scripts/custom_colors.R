#' [custom_colors.R]
#'
#' Collection of color palettes and scales
#' 

#### libs ####
require(RColorBrewer)
require(scico)
require(viridis)
require(ggsci)  # https://cran.r-project.org/web/packages/ggsci/vignettes/ggsci.html

#### color scales ####
col_scale_bupu <- RColorBrewer::brewer.pal(9, "BuPu")
col_scale_purd <- RColorBrewer::brewer.pal(9, "PuRd")
col_scale_spec <- RColorBrewer::brewer.pal(9, "Spectral")
col_scale_or <- RColorBrewer::brewer.pal(9, "Oranges")
col_scale_acton <- scico::scico(10, direction = -1, palette = "acton")
col_scale_batlow <- scico::scico(10, direction = -1, palette = "batlowW")
col_scale_bilbao <- scico::scico(10, palette = "bilbao")
col_scale_tokyo <- scico::scico(10, palette = "tokyo")
col_scale_rocket <- viridis::rocket(10, direction = -1)
col_scale_mako <- viridis::mako(13, direction = -1)[1:10]
col_scale_viridis <- viridis::viridis(10, direction = -1)
col_scale_lajolla <- hcl.colors(n=10, palette = "Lajolla")
col_scale_earth <- hcl.colors(n = 10, palette = "Earth")
col_scale_temps <- hcl.colors(n = 10, palette = "Temps", rev = T)
col_scale_or_blu <- c(hcl.colors(n=14, palette = "Lajolla", rev = T)[2:13], hcl.colors(n=7, palette = "BluYl", rev = T)[2:7])


col_scale_div_custom <- c(rev(col_scale_mako[1:5]), col_scale_acton[1:8])
col_scale_div_custom2 <- c(rev(col_scale_mako[1:8]), "white", col_scale_acton[1:8])
col_scale_div_expr <- c(rev(col_scale_mako[1:8]), col_scale_lajolla[1:8])

#### categorical ####
cols_d3_20 <- ggsci::pal_d3("category20")(20)
cols_cat_custom <- c("#5DA6C6", "#3873B9", "#215189",
                              "#E391BA", "#8D71A6",
                              "#F4D8A8", "#EDC962", "#E7B922", "#DB8712", "#CA6623", "#AB4A3D",
                              "#B1E4C2", "#78D6AE", "#49C1AD")

#### data specific palettes ####
if (SPECIES == "human") {
  message("Reading colors specific for human data set")
  
  cols_donor <- setNames(c(col_scale_mako[c(2,4,6,8)], col_scale_acton[c(2,4,6,8)]), 
                         nm = c(paste0("HC_", 1:4), paste0("IPF_", 1:4)))
  
  cols_cond <- setNames(cols_donor[c(2,7)], nm = c("control", "IPF"))
  
  cols_grade <- setNames(c(cols_donor[2], cols_donor[6:8]), 
                         nm = as.character(0:3))

  cols_annotation <- setNames(object = c("#1F77B4FF",
                                         "#9EDAE5FF",
                                         "#FF9896FF",
                                         "#FF7F0EFF",
                                         "#D62728FF",
                                         "#C5B0D5FF"),
                              nm = c("Normal Alveolar and Other",
                                     "Large Airway",
                                     "Blood Vessel",
                                     "Diseased",
                                     "Inflammation",
                                     "Suspect Fibrosis"))
  cols_annotation2 <- setNames(object = c("#5DA6C6",
                                          "#215189",
                                          "#E391BA",
                                          "#EDC962",
                                          "#AB4A3D",
                                          "#8D71A6"),
                                          nm = c("Normal Alveolar and Other",
                                                 "Large Airway",
                                                 "Blood Vessel",
                                                 "Diseased",
                                                 "Inflammation",
                                                 "Suspect Fibrosis"))
  
} else if (SPECIES == "mouse") {
  message("Reading colors specific for mouse data set")
  
  cols_group <- setNames(RColorBrewer::brewer.pal(4, "Spectral")[c(1,4, 2,3)], 
                         nm = c("d21_BLM", "d21_CTRL","d7_BLM", "d7_CTRL"))

  cols_cond <- setNames(c("#215288", "#7AB3C1"), nm = c("bleomycin", "control"))
  
  cols_annotation <- setNames(object = c("#1F77B4FF",
                                         "#9EDAE5FF",
                                         "#FF9896FF",
                                         "#D62728FF",
                                         "#D62728FF",
                                         "#FF7F0EFF",
                                         "#C5B0D5FF"),
                              nm = c("Normal Alveolar and Other",
                                     "Large Airway",
                                     "Blood Vessel",
                                     "Inflammation",
                                     "Inflammation_d7",
                                     "Inflammation_d21",
                                     "Suspect Fibrosis/Fibroplasia"))
} else if (SPECIES == "translational") {
  message("Reading colors specific for translational analyses")
  cols_annotation <- setNames(object = c("#1F77B4FF",
                                         "#9EDAE5FF",
                                         "#FF9896FF",
                                         "#D62728FF",
                                         "#D62728FF",
                                         "#FF7F0EFF",
                                         "#FF7F0EFF",
                                         "#C5B0D5FF",
                                         "#C5B0D5FF"),
                                         nm = c("Normal Alveolar and Other",
                                                "Large Airway",
                                                "Blood Vessel",
                                                "Inflammation",
                                                "Inflammation_d7",
                                                "Inflammation_d21",
                                                "Diseased",
                                                "Suspect Fibrosis/Fibroplasia",
                                                "Suspect Fibrosis"
                                                ))
  
  cols_cond <- setNames(c("#9EDAE5", "#B26694", "#215288"), # #1F77B4FF
                        nm = c("control", "IPF", "bleomycin"))
  cols_species <- setNames(c("#B26694", "#215288"),
                           nm = c("human", "mouse"))
  cols_donor <- setNames(c(col_scale_mako[c(2,4,6,8)], col_scale_acton[c(2,4,6,8)]), 
                         nm = c(paste0("HC_", 1:4), paste0("IPF_", 1:4)))
  
} else (
  message("Please specify SPECIES (human or mouse) before sourcing")
)

# scales::show_col(cols_d3_20)
