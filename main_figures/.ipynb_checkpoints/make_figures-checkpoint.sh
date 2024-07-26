#!/bin/sh
convert dry_LPO_cross_sections.png wet_LPO_cross_sections.png +smush 25 total_LPO_depth_cross_sections.png
convert 10.0-melting-regions.png 50.0-melting-regions.png +smush 25 melting-regions.png

convert convergence_rates_figures/vel_dry_1e20.png convergence_rates_figures/vel_dry_5e20.png convergence_rates_figures/vel_dry_1e21.png +smush 10 convergence_rates_figures/top_row.png
convert convergence_rates_figures/vel_wet_upper_mantle_5e20.png convergence_rates_figures/vel_wet_slab_5e20.png convergence_rates_figures/vel_wet_crust_5e20.png +smush 10 convergence_rates_figures/bot_row.png
convert convergence_rates_figures/top_row.png convergence_rates_figures/bot_row.png -smush 150 total_velocity.png