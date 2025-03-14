# OSTPRE-brain-MRI

This repository contains code used for the analyses in the study
["Real-world brain imaging in a population-based cohort enables accurate markers for dementia"](https://doi.org/10.1101/2025.02.13.25322207).

The study describes how real-world data (RWD) on brain MRI scans on a large cohort of
aged women (*OSTPRE*) can be used and analyzed to extract several volume parameters.
We analyse how these volume parameters associate with cognitive status and investigate
are the detected associations similar than detected in separately collected research
data (*ADNI*).

The repository is organized as follows:

Folder [utils](utils) contains scripts for initial processing.
- [*brain_dicom2nifti.R*](utils/brain_dicom2nifti.R) converts original DICOMs to NIfTI-format and collects relevant metadata tags
- [*sort.py*](utils/sort.py) detects the T1-weighted scans using the extracted metadata tags
- [*create_bids_folders.py*](utils/create_bids_folders.py) (re)organises and renames NIfTI data to BIDS-compatible structure

Folder [image_analysis_matlab](image_analysis_matlab) contains the image analysis scripts done in Matlab and
folder [slurm](slurm) contains the slurm scripts used to run the image analysis scripts on the HPC server.

Folder [analysis](analysis) contains the [R-script](analysis/brain_imagedementia.R) that performs the statistical analyses
and plots the Figures reported in the study.

