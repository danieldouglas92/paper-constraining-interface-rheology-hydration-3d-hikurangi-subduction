This repository contains all scripts required to reproduce the results and the figures in Douglas et al. (2024), "Constraining Solid Dynamics, Interface Rheology, and Slab Hydration in the Hikurangi Subduction Zone Using 3D Fully Dynamic Models", submitted to G3. With the exception of the Paraview state files (*.pvsm files), and the two scripts used to create Figures 5 and 6 (Convergence_Comparison.ipynb and additional_tests.ipynb), all scripts can be run without additional modifications if they are run from the location that they were placed in this repository. The two scripts used to create Figures 5 and 6 require the solution files for every model, which exceeds the storage limits of this github repository. The *.pvsm files can be run by specifying the path to the solution files within this repository upon opening the state files. The summary of each directory is outlined below.

ASPECT_models: Contains the solution files for 3 of the ASPECT models, the model which best reproduced Pacific-Australian convergence (5e20 Pa s interface viscosity, dry rheology), and the two calculations of LPO (dry, 0 H/Si ppm and wet, 1000 H/Si ppm).

ASPECT_parameter_files: Contains the input files required to reproduce the ASPECT models, as well as the world builder file used for the initial conditions.

ASPECT_source: Contains the source code for the version of ASPECT used to run the models.

WorldBuilder_source: Contains the source code for the version of WorldBuilder used for the ASPECT models.

data: Contains relevant data files required by the various python scripts which generate the published figures, or the WorldBuilder file.

main_figures: Contains the python notebooks and paraview state files used to produce the published figures, as well as the published figures.

postprocessing_scripts: Contains various python scripts which are used by the python notebooks in main_figures.

slab2_worldbuilder_generation: Contains python scripts and a python notebook which generates the WorldBuilder file.

supplemental_figures: Contains python notebooks which generate the supplemental figures, as well as the supplemental figures.
