## cds2-synthetic-lethal-notebooks

CDS2 is a synthetic lethal target in mesenchymal-like cancers (software notebooks)

## System requirements

For .rmd files:
Rstudio (2023.03.0+386)
Package Rmarkdown version 2.25

For .ijm files:
FIJI 2.14.0
Java 1.8.0
ImageJ 1.54f

For code in word document:
MAGeCK v0.5.9.5 was run in terminal

Provided code was rerun by a second independent researcher. Mac was always used. Running code should be manageable on an average computer. Code for 1d and ext1c may take half an hour to run. Loading large files when running code is fastest if those files are saved locally.

## Installation guide

Packages used in the Studio code can be installed using "install.packages("<INSERT PACKAGE NAME>")". Each Rmarkdown file runs "library("<INSERT PACKAGE NAME>")" for each package used.

Installing R, Rstudio, Java, FIJI, ImageJ and Mageck was previously described:
R:
https://cran.r-project.org/
Rstudio:
https://rstudio-education.github.io/hopr/starting.html
Java:
https://www.java.com/en/
FIJI & ImageJ:
https://imagej.net/software/fiji/downloads
Mageck:
https://sourceforge.net/p/mageck/wiki/demo/#the-first-demo-starting-from-read-count-tables

Install time for these programs should be manageable on an average computer.

## Demo

Rmarkdown files are provided to allow reproduction of every complex computational analysis in the paper with explanation of the code by chunk (see "1-software-notebooks"). The code is named by panel. Input files are provided, or can be downloaded from repositories using "Instructions for obtaining files.rmd" (see "2-files-to-run-software"). If data was plotted using Graphpad Prism and Adobe Illustrator the plot data was generated by the Rmarkdown code. In some cases data was (also) plotted using ggplot. In some cases we added rudimentary plots in Rmarkdown and the actual plots were made in Prism. Plot data for every panel is provided as Source Data in the manuscript.

## Instructions for use

The software notebooks are designed to allow step by step reproduction of the plot data for every complex bioinformatics analysis and to allow running tweaked code for also new questions. They are intended to be opened as .rmd and run chunk by chunk using the play buttons in the upper right corner of each chunk. Intermediate results can appear below a chunk after running it. The Rmarkdown files are not intended to be knitted. Which chunk of code provides what functionality is regularly explained above each chunk or by comments within a chunk.

The code can be tweaked for other research questions. For example, it is possible to switch a gene name, or change input data to cover another dataset (like your own data) or data subset (for example covering a different cancer type). 

The code was designed to be an extension of the manuscript and duplication was avoided when possible. Therefore, we advise also studying the figure panel, figure legend, methods section and manuscript text for context on what code for a panel is doing. All software and related files are provided here on GitHub. All other relevant information is provided in the manuscript or on repositories to allow larger file sizes. Unprocessed proteomic data is provided under accession PXD045833 on the PRIDE ProteomeXchange repository. Other unprocessed data for the manuscript is provided on Figshare.
