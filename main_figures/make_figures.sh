#!/bin/sh
convert dry_LPO_cross_sections.png wet_LPO_cross_sections.png +smush 25 total_LPO_depth_cross_sections.png
convert schematic_diagram.png 3d_slab_underview.png depth_tectonics.png -smush 25 schematic_figure.png
convert hikurangi_density.png middle_density.png kermadec_density.png -smush 25 density_cross_sections.png
convert 1km_depth_visc.png 150km_depth_visc.png vel_visc_colourbar.png -smush 25 depth_viscosity_cross_sections.png
convert spherical_mesh.png plateau_interface.png temperature_colourbar.png -smush 25 mesh_interface.png

convert convergence_rates_figures/vel_dry_1e20.png convergence_rates_figures/vel_dry_5e20.png convergence_rates_figures/vel_dry_1e21.png +smush 10 convergence_rates_figures/top_row.png
convert convergence_rates_figures/vel_wet_upper_mantle_5e20.png convergence_rates_figures/vel_wet_slab_5e20.png convergence_rates_figures/vel_wet_crust_5e20.png +smush 10 convergence_rates_figures/bot_row.png
convert convergence_rates_figures/top_row.png convergence_rates_figures/bot_row.png -smush 150 total_velocity.png