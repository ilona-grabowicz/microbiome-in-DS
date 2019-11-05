# microbiome-in-DS
This repository contains the compressed files with code and input files allowing to reproduce the results of the analysis of the microbiome data. Those data were obtained to determine the microbial content of the faecal microbiome samples from Wild Type and Down Syndrome model mice fed high-fat diet.
## What is in the package?
- script to run seed-kraken
- the .csv files with microbial species counts generated by seed-kraken tool
- python code to process the .csv files
- R code to analyse the data and visualise results (analysis.R)
## Basic usage
The script:
```bash run_all.sh```
will run the data pre-processing, data analysis performed in R and generate figures used in the manuscript.

## Detailed usage
Run: 
```analysis.R``` 
in RStudio to investigate all the steps of the analysis.
