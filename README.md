## Mapping locations of cortical single astrocyte RNA-seq

This repository contains the code needed to reproduce the mapping of astrocytes from scRNA-seq data to cortical layers in the somatosensory cortex in mice using smFISH reference.    

Astrocyte layers in the mammalian cerebral cortex revealed by a single-cell in situ transcriptomic map   

Omer Ali Bayraktar, Theresa Bartels, Staffan Holmqvist, Vitalii Kleshchevnikov, Araks Martirosyan, Damon Polioudakis, Lucile Ben Haim, Adam M.H. Young, Mykhailo Batiuk, Kirti Prakash, Alexander Brown, Kenny Roberts, Mercedes F. Paredes, Riki Kawaguchi, John Stockley, Khalida Sabeur, Sandra M. Chang, Eric Huang, Peter Hutchinson, Erik M. Ullian, Martin Hemberg, Giovanni Coppola, Matthew G. Holt, Daniel H. Geschwind & David H. Rowitch     
Nature Neuroscience, 2019

This is not meant to be a standalone software but the R and Matlab scripts should allow reproducing spatial reconstruction figures from our paper.   
Matlab scripts are adapted from Halpern et al, Nature 2017 and were kindly shared by Shalev Itzkovitz.   

Main scripts to note:

1. Holt_astrocytes_spatial_map_Shalev_p56_16markers_n_total180k.Rmd - contains code needed to reproduce spatial reconstruction figures as well as additional diagnostic plots.   
2. reviewer_comments.Rmd - contains code to reproduce figures that rely only on single cell data.  
3. spatial_mapping_code/Zonation_main_with_cell_groups_leave_out.m - master script to run leave-one-out cross-validation analysis.  
4. spatial_mapping_code/Zonation_main_with_cell_groups_p56_new.m - main script that loads data and defines the model. All other scripts define functions used by this script.  
