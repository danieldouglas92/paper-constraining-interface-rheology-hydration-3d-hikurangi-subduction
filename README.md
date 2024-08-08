This repository is assosciated with the publication

*Constraining Solid Dynamics, Interface Rheology, and Slab Hydration in the Hikurangi Subduction Zone Using 3D Fully Dynamic Models*

by

Douglas, D.
Naliboff, J.
Fraters, M. R. T.
Dannberg, J.
Eberhart-Philips, D.
Ellis, S.

which is currently in review.

# Software
Our numerical simulations were run using the open source geodynamics software ASPECT ([https://aspect.geodynamics.org/](https://aspect.geodynamics.org/)) and the initial conditions software the Geodynamic World Builder ([https://gwb.readthedocs.io/en/latest/](https://gwb.readthedocs.io/en/latest/)). Specifically, we utilize the developmental versions of ASPECT (2.5.0-pre) and the Geodynamic World Builder (0.6.0), with copies of the the source code for both software packages included in this repository.

# Overview
All scripts required to reproduce the results and the figures in Douglas et al. (2024), "Constraining Solid Dynamics, Interface Rheology, and Slab Hydration in the Hikurangi Subduction Zone Using 3D Fully Dynamic Models", submitted to G3. With the exception of the Paraview state files (*.pvsm files), and the two scripts used to create Figures 5 and 6 (Convergence_Comparison.ipynb and additional_tests.ipynb), all scripts can be run without additional modifications if they are run from the location that they were placed in this repository. The two scripts used to create Figures 5 and 6 require the solution files for every model, which exceeds the storage limits of this github repository. The *.pvsm files can be run by specifying the path to the solution files within this repository upon opening the state files. The summary of each directory is outlined below.

## ASPECT_models
Contains solution files for three ASPECT models: the model that most accurately reproduced Pacific-Australian convergence (with a 5e20 Pa s interface viscosity and dry rheology), along with two LPO calculations assuming two hydration conditions, dry (0 H/Si ppm) and wet (1000 H/Si ppm).

## ASPECT_parameter_files
Includes the input files needed to reproduce the ASPECT models, along with the world builder file used to set the initial conditions.

## ASPECT_source
Contains the source code for the version of ASPECT used to run the models.

## WorldBuilder_source
Contains the source code for the version of WorldBuilder used for defining the initial conditions in the models.

## data
Contains relevant data files required by the various python scripts which generate the published figures.

## main_figures
Contains the Python notebooks and ParaView state files used to generate the published figures, along with the figures themselves.

## postprocessing_scripts
Contains various python scripts which are used by the python notebooks in main_figures.

## slab2_worldbuilder_generation
Contains python scripts for generating WorldBuilder files, and a python notebook which generates the WorldBuilder file used in the ASPECT models.

## supplemental_figures
Contains python notebooks which generates the supplemental figures, as well as the supplemental figures.
