---
title: "Intra- and inter-individual variation masks inter-species variation in the microbiota of sympatric *Peromyscus* populations"
author: "Nielson T. Baxter & Patrick D. Schloss"
date: ""
output: html_document
    keep_md: yes

---

# Intra- and inter-individual variation masks inter-species variation in the microbiota of sympatric _Peromyscus_ populations  
Nielson T. Baxter, Judy J. Wan, Alyxandria M. Schubert, Matthew L. Jenior, Philip Myers, and Patrick D. Schloss

## Introduction
This is a digital notebook to accompany the paper titled, "Intra- and inter-individual variation masks inter-species variation in the microbiota of sympatric _Peromyscus_ populations" that will be published in *Applied & Environmental Microbiology*. It was written in [R markdown](http://rmarkdown.rstudio.com) and converted to html using the R knitr package. This enables us to embed the results of our analyses directly into the text to allow for a completely reproducible data analysis pipeline. A [github repository is available](https://github.com/schlossLab/wild_mice) where you can pull down your own version of the notebook to modify our analysis or adapt it to your analysis. 

This document was generated using [mothur v.1.34](http://www.mothur.org/wiki) and [R](http://www.r-project.org):


```r
deps = c("randomForest", "maps", "vegan", "date", "knitr");

for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  }
  library(dep, verbose=FALSE, character.only=TRUE)
}
```

```
## randomForest 4.6-10
## Type rfNews() to see new features/changes/bug fixes.
## Loading required package: permute
## Loading required package: lattice
## This is vegan 2.0-10
```

```r
sessionInfo()
```

```
## R version 3.0.1 (2013-05-16)
## Platform: x86_64-redhat-linux-gnu (64-bit)
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=C                 LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] date_1.2-34         vegan_2.0-10        lattice_0.20-29    
## [4] permute_0.8-3       maps_2.3-9          randomForest_4.6-10
## [7] knitr_1.6.21       
## 
## loaded via a namespace (and not attached):
## [1] evaluate_0.5.5 formatR_1.0    grid_3.0.1     stringr_0.6.2 
## [5] tools_3.0.1
```

```r
perm = 1e4
```

To convert the R markdown to html (or any other format), you'll also need to install the [SRA tools package](http://www.ncbi.nlm.nih.gov/books/NBK158900/) and have a live internet connection. You can make the conversion using the [knitr package](http://yihui.name/knitr/) from within R:


```r
knit2html("wild_mice.Rmd")
```

I have used the following knitr settings to compile this document:





























































