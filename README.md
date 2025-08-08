# Individual differences of neurophysiological traits across the lifespan  

## Introduction

Brain-fingerprinting is a novel technique used to explore the inter-individual differences in brain activity. Previous research demonstrates that we are capable of accurately differentiating individuals in a cohort based solely on the resting-stage spectral brain activity using magnetoencephalography (MEG). How these fingerprints are change across the lifespan remain to be explored.

## Description of the project

Here are the different steps we did :

- 1-We preprocessed data taken from the CamCAN and SickKids cohorts. We then derived neural power spectra for each parcel of the Destrieux atlas for every individual from split-half recordings. This resulted in 2 brain fingerprints to differentiate individuals from one another. 
- 2-This code applies the brain-fingerprinting approach to each dataset separately to compute i) differentiation accuracies and ii) differentiability. Here we explore age effects and show that participants become more differentiable from their neurophysiological traits across neurodevelopment.
- 3-We compute salient features of neurophysiological traits across age groups using a sliding window approach. We then assess the alignment of these salient features to the first functional gradient (i.e., unimodal-to-transmodal)
- 4-We assess the alignment between salient features for neurophysiological traits and cortical patterns of gene expression derived from our [previous work](https://github.com/Epideixx/Fingerprints_Twins) 
- 5-We assess the multivariate alignment between cortical patterns of gene expression and salient features for participant differentiation across the lifespan to verify whether genetic signatures for participant differentiation remain stable across neurodevelopment. 
- 6-Empty room pseudo neurophysiological trait analysis to ensure the robustness of our participant differentiation results. 
- 8-Control analysis to ensure the robustness of the 3rd order polynomial effects reported in Figures 2 and 3. Here we randomly permuted the association between participants age and their neurophysiological trait.

An explanation of the gene expression data obtained from [AHBA](https://human.brain-map.org) can be found at [Hansen et al 2021](https://github.com/netneurolab/hansen_genescognition)


## Manuscript and Citation

This work on [bioRxiv](https://www.biorxiv.org/content/10.1101/2024.11.27.624077v1.abstract). Please cite da Silva Castanheira et al., 2024. If you have any questions please contact the authors of the paper.


