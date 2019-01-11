## Prerequisities:

bash, GCC, G++
yum install pandoc mysql-devel libcurl-devel libxml2-devel

conda create --name pf r=3.4.3
conda acticvate pf

install.packages("XML")
install.packages("rmarkdown")
apt-get install libxml2-dev
install.packages("data.table", dependencies=TRUE)
install.packages("ggplot2")
install.packages("dplyr")
install.packages("ggrepel")
install.packages("VennDiagram")
install.packages("htmlTable")
install.packages("tictoc")

 - packages to add:
 - install.packages("rmarkdown")
 - apt-get install pandoc
 - apt-get install libcurl-dev
 - apt-get install libxml2-dev
 - install.packages("data.table", dependencies=TRUE) 
 - install.packages("ggplot2")
 - install.packages("dplyr")
 - install.packages("ggrepel")
 - install.packages("VennDiagram")
 - install.packages("htmlTable")
 - install.packages("kableExtra")

 - source("https://bioconductor.org/biocLite.R")
 - biocLite("VariantAnnotation")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("VariantAnnotation", version = "3.8")
