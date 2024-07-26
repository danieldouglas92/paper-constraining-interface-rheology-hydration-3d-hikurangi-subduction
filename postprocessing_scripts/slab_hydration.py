import numpy as np
from scipy.spatial import KDTree
import parallel_curves

def water_concentration_func(LR_polyarr, c_sat_polyarr, Td_polyarr, pressure, temperature, lithology):
    """ 
    This function takes polynomials used to define the phase diagram for hydration of a bulk rock composition from Tian et al. 2019, as well as
    a temperature and pressure, and outputs the partition coefficient of the rock at the specified P-T conditions.
    
    lithology: A string that specifies which rock composition from Tian et al., 2019 is being used. Can either be: "sediment", "MORB", "gabbro", or "peridotite".
    LR_polyarr: An array with the polynomial coefficients for LR. Highest order coefficient first.
    c_sat_polyarr: An array with the polynomial coefficients for c_sat. Highest order coefficient first.
    Td_polyarr: An array with the polynomial coefficients for Td. Highest order coefficient first.
    pressure: The pressure, in GPa, of the given lithology
    temperature: The temperature, in K, of the given lithology
    """
    
    # The polynomials defined in Tian et al., 2019 break down above certain pressures. If the input pressure is above this pressure,
    # set it to the highest pressure where the polynomials do NOT breakdown
    if lithology == "gabbro":
        if pressure > 26:
            pressure = 26
            
    if lithology == "MORB":
        if pressure > 16:
            pressure = 16
            
    if lithology == "peridotite":
        if pressure > 10:
            pressure = 10

    # Initialize the variables
    inv_pressure = 1/pressure
    ln_LR_val = 0
    ln_c_sat_val = 0
    Td_val = 0

    # Calculate the values for LR, c_sat, and Td
    for i in range(len(LR_polyarr)):
        ln_LR_val += LR_polyarr[i] * (inv_pressure**(len(LR_polyarr) - 1 - i))

    for j in range(len(c_sat_polyarr)):
        if lithology == "sediment":
            ln_c_sat_val += c_sat_polyarr[j] * (np.log10(pressure)**(len(c_sat_polyarr) - 1 - j))
            
        else:
            ln_c_sat_val += c_sat_polyarr[j] * (pressure**(len(c_sat_polyarr) - 1 - j))

    for k in range(len(Td_polyarr)):
        Td_val += Td_polyarr[k] * (pressure**(len(Td_polyarr) - 1 - k))

    LR_val = np.exp(ln_LR_val)
    c_sat_val = np.exp(ln_c_sat_val)

    # Calculate the partition coefficient (amount of water that can be stored in the rock)
    partition_coeff = c_sat_val * np.exp(LR_val * (1/temperature - 1/Td_val))
    return partition_coeff

def dehydration_func(LR_polyarr, c_sat_polyarr, Td_polyarr, lithology_density, parallel_slab_pressure, parallel_slab_temperature, initial_water_content, \
                     hydration_container, dehydration_container, parallel_slab_surface_profiles, composition_depth_range, line_spacing, sediment):
    """
    This function uses the water_concentration_func to first compute the water content of the slab at each Pressure-Temperature point in the slab.
    The function assumes a maximum initial_water_content, and then assigns the hydration state of the slab so that it is equal to or less than the
    initial value. An input set of surface slab parallel profiles are specified, and the points in these profiles are used to track the water at even
    intervals of depth. 
    
    After the hydration state of the slab is computed. Each parallel slab profile is iterated along, tracking the change in dehydration from point to point.
    This defines the amount of dehydration.
    
    LR_polyarr                    : The polynomial used to describe the effective enthalpy change for the specified material composition from Tian et al. 2019
    c_sat_polyarr                 : The polynomial used to describe the saturated mass fraction of water for the specified material composition from Tian et al 2019
    Td_polyarr                    : The polynomial used to describe the onset temperature of devolatilization for water for the specified material composition from Tian et al. 2019
    parallel_slab_pressure        : The pressure in the slab along each point in the parallel profiles
    parallel_slab_temperature     : The temperature in the slab along each point in the parallel profiles
    initial_water_content         : the upper bound on water content
    parallel_slab_surface_profiles: array containing the cartesian points of the profiles parallel to the slab surface
    """
    for i in range(len(parallel_slab_surface_profiles)):
        if (i * line_spacing < np.max(composition_depth_range)) & (i * line_spacing >= np.min(composition_depth_range)):
            hydration_profile = np.zeros(len(parallel_slab_surface_profiles[i]))
            dehydration_profile = np.zeros(len(parallel_slab_surface_profiles[i]))

            for k in range(len(parallel_slab_surface_profiles[i])):
                point_temp = parallel_slab_temperature[i][k]
                point_pressure = parallel_slab_pressure[i][k]

                if point_pressure < 1:
                    point_pressure = 1
                point_water = water_concentration_func(LR_polyarr, c_sat_polyarr, Td_polyarr, point_pressure / 1e9, point_temp, sediment) / 100
                
                if k == 0:
                    dehydration_profile[k] = 0
                    if point_water > initial_water_content:
                        hydration_profile[k] = initial_water_content
                    else:
                        hydration_profile[k] = point_water

                else:
                    if point_water <= hydration_profile[k - 1]:
                        hydration_profile[k] = point_water
                    else:
                        hydration_profile[k] = hydration_profile[k - 1]
                    dehydration_profile[k] = np.abs(hydration_profile[k] - hydration_profile[k - 1]) # change in wt % water

            hydration_container[i] = hydration_profile
            dehydration_container[i] = dehydration_profile * lithology_density # kilograms of water

    return hydration_container, dehydration_container

def vertical_surface_water_flux(parallel_slab_surface_profiles, parallel_slab_surface_depths, dehydration_amount, line_spacing):
    """
    This function takes a series of profiles drawn parallel to the surface of the slab, as well as the amount of dehydration
    calculated along each parallel profiles in the function dehydration_func. The true slab surface is iterated along from trench
    to the leading edge of the slab, and the values of dehydration from all points in parallel profiles below the current point are
    summed, giving an effective water flux out of the surface of the slab. The value is then converted to H/Si ppm.
    """
    
    surface_slab_water_flux = np.zeros(len(parallel_slab_surface_profiles[0]), dtype=object)
    # slab_top_dx = abs(parallel_slab_surface_profiles[0][0] - parallel_slab_surface_profiles[0][1])    
    for j in range(len(surface_slab_water_flux)):
        depth_water_profile = np.zeros(len(parallel_slab_surface_profiles))

        for i in range(len(depth_water_profile)):
            # This checks for the points with the same x-position in all profiles below the surface, extracting a vertical profile
            subsurface_vertical_index = np.argmin(np.abs(parallel_slab_surface_profiles[i] - parallel_slab_surface_profiles[0][j]))

            # This if statement checks the distance of the point in the profile below the current profile
            if np.abs(parallel_slab_surface_profiles[i][subsurface_vertical_index] - parallel_slab_surface_profiles[0][j]) >= line_spacing:
                depth_water_profile[i] = 0

            else:
                depth_water_profile[i] = dehydration_amount[i][subsurface_vertical_index]
                
                if i == 1:
                    dxdz = np.abs(parallel_slab_surface_depths[0][j] - parallel_slab_surface_depths[i][subsurface_vertical_index])

        surface_slab_water_flux[j] = np.trapz(depth_water_profile, dx=dxdz)

    return surface_slab_water_flux


def surface_water_flux(parallel_slab_surface_profiles, parallel_slab_surface_depths, dehydration_amount, line_spacing):
    """
    This function takes a series of profiles drawn parallel to the surface of the slab, as well as the amount of dehydration
    calculated along each parallel profiles in the function dehydration_func. The values of dehydration are summed up at
    all parallel profiles perpendicular to the slab surface to get a surface flux. This value is then converted into H/Si ppm.
    """
    surface_slab_water_flux = np.zeros(len(parallel_slab_surface_profiles[0]), dtype=object)

    for j in range(len(surface_slab_water_flux)):
        depth_water_profile = np.zeros(len(parallel_slab_surface_profiles))

        for i in range(len(depth_water_profile)):
            depth_water_profile[i] = dehydration_amount[i][j]
        surface_slab_water_flux[j] = np.trapz(depth_water_profile, dx=line_spacing)
        
        Si_molar_wt = 28.085 # g/mol
        SiO2_molar_wt = 28.085 + 15.999*2
        SiO2_wt_percent = 48.64
        H2O_molar_wt = 1.008*2 + 15.999 # g/mol
        gabbro_wt = 63.7602958 # g/mol
        
        if j <= len(parallel_slab_surface_profiles[0]):
            dxdz = np.sqrt(line_spacing**2 + (parallel_slab_surface_depths[0][j] - parallel_slab_surface_depths[0][j - 1])**2)
        else:
            dxdz = line_spacing
        
        Si_wt = SiO2_wt_percent * SiO2_molar_wt * (Si_molar_wt/SiO2_molar_wt)
        H_wt = surface_slab_water_flux[j] * H2O_molar_wt * (1.008*2/H2O_molar_wt)
        surface_slab_water_flux[j] = H_wt/Si_wt * 1e6 / dxdz
    
    return surface_slab_water_flux