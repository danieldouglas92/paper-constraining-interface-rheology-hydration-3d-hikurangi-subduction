convert vel_dry_1e20.png vel_dry_5e20.png vel_dry_1e21.png +smush 10 top_row.png
convert vel_wet_upper_mantle_5e20.png vel_wet_slab_5e20.png vel_wet_crust_5e20.png +smush 10 bot_row.png
convert top_row.png bot_row.png -smush 150 total.png
