# install pkgs
reqPkg = c("data.table", "Matrix", "hdf5r", "reticulate", "ggplot2", 
           "gridExtra", "glue", "readr", "RColorBrewer", "R.utils", "Seurat")
newPkg = reqPkg[!(reqPkg %in% installed.packages()[,"Package"])]
if(length(newPkg)){install.packages(newPkg)}

reqPkg = c("shiny", "shinyhelper", "data.table", "Matrix", "DT", "hdf5r", 
           "reticulate", "ggplot2", "gridExtra", "magrittr", "ggdendro")
newPkg = reqPkg[!(reqPkg %in% installed.packages()[,"Package"])]
if(length(newPkg)){install.packages(newPkg)}

devtools::install_github("SGDDNB/ShinyCell")


#
setwd('Documents/snRNA/MS/mouse/ten_samples/')
#mice data with six samples
inpFile = 'results/Mouse_harmony_rmBatchbetweenSamples.h5ad'

library(ShinyCell)
scConf = createConfig(inpFile)
makeShinyApp(inpFile, scConf, 
             default.gene1 = "Clic6", default.gene2 = "Folr1",
             default.multigene = c("Clic6","Folr1","Krt18","Htr2c","Car12", "Igfbp2", "Ttr"),
             shiny.dir = "shinyAppH5ad_MiceWithDbl/", 
             shiny.title = "ten mice samples in ChP") 
