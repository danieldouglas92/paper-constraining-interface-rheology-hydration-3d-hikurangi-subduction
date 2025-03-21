This repository is assosciated with the publication

*Constraining Solid Dynamics, Interface Rheology, and Slab Hydration in the Hikurangi Subduction Zone Using 3D Fully Dynamic Models*

by

Douglas, D.,
Naliboff, J.,
Fraters, M. R. T.,
Dannberg, J.,
Eberhart-Philips, D.,
Ellis, S.,

which is currently in review.

# Software
Our numerical simulations were run using the open source geodynamics software ASPECT ([https://aspect.geodynamics.org/](https://aspect.geodynamics.org/)) and the initial conditions software the Geodynamic World Builder ([https://gwb.readthedocs.io/en/latest/](https://gwb.readthedocs.io/en/latest/)). Specifically, we utilize the developmental versions of ASPECT (2.5.0-pre) and the Geodynamic World Builder (0.6.0), with copies of the the source code for both software packages included in this repository.

# Overview
All scripts required to reproduce the results and the figures in Douglas et al. (2024), "Constraining Solid Dynamics, Interface Rheology, and Slab Hydration in the Hikurangi Subduction Zone Using 3D Fully Dynamic Models", submitted to G3. All scripts must be modified to properly include the path to the model output on your local machine. The *.pvsm files can be run by specifying the path to the solution files on your local machine upon opening the state files. The summary of each directory is outlined below.

## ASPECT_parameter_files
Includes the input files needed to reproduce the ASPECT models.

## GWB_files
Includes the world builder file used to setup the initial conditions.

## ASPECTv2.5.0_source
Contains the source code for the ASPECT version 2.5.0 used to run all of the models, with the exception of the heterogeneously hydrated models.


## ASPECTv3.1.0_source
Contains the source code for the ASPECT version 3.1.0 used only for the heterogeneously hydrated models.

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
