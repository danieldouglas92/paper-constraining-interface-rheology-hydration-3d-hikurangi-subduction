import numpy as np
import matplotlib.pyplot as plt
import os
import pygmt
from scipy.spatial import KDTree

# Use Magali's documentation to calculate spherical velocities:
# Magali defines variables as lamba = latitude deg (Positive for North of 0, negative for South of 0)
#                             phi = longitude deg (Positive for East of 0, negative for West of 0)
#                             omega = angular rotation deg/Myr
# Angles that are entered in sines/cosines must be entered in radians

def MORVEL_surface_velocity(euler_pole_lat, euler_pole_lon, euler_rotation, point_lat, point_lon):
    # Here we calculate the magnitude of the velocity at a point relative to an euler pole
    # First we calculate the distance from the euler pole to the point of interest using:
    # arccos[sin(point_lat) * sin(euler_pole_lat) + cos(point_lat) * cos(euler_pole_lat) * cos(euler_pole_lon - point_lon)]
    # This gives the distance in radians
    point_to_pole_distance_radians = np.arccos(np.sin(np.deg2rad(point_lat)) * np.sin(np.deg2rad(euler_pole_lat)) + \
                                               np.cos(np.deg2rad(point_lat)) * np.cos(np.deg2rad(euler_pole_lat)) * np.cos(np.deg2rad(euler_pole_lon - point_lon)))

    # Now we calculate the angle between the great circle connecting the point of interest to the euler pole,
    # and the point of interest to the North pole
    angle_north_pole_to_euler_pole = np.rad2deg( np.arcsin( (np.cos(np.deg2rad(euler_pole_lat)) * np.sin(np.deg2rad(euler_pole_lon - point_lon))) / \
                                                           np.sin(point_to_pole_distance_radians) ) )

    azimuth = 90 + angle_north_pole_to_euler_pole
    velocity =  np.deg2rad(euler_rotation) * 6371e3 * np.sin(point_to_pole_distance_radians)
    return velocity, azimuth



def ASPECT_surface_velocity(x, y, z, data, point_lats, point_lons, depth, padding):
    
    # Create arrays to store the ouput (the closest x, y, z and velocity from the ASPECT models)
    ASPECT_x = np.zeros(len(point_lats))
    ASPECT_y = np.zeros(len(point_lats))
    ASPECT_z = np.zeros(len(point_lats))
    ASPECT_v = np.zeros(len(point_lats), dtype=object)
    
    # First extract a depth slice from the ASPECT models. This is difficult to do in Cartesian, so calculate
    # ASPECT radius and use that to extract a depth slice
    r = np.sqrt(x**2 + y**2 + z**2)
    depth_index = np.where( (r <= 6371e3 - depth + padding) & (r >= 6317e3 - depth - padding) )
    
    x_slice    = x[depth_index]
    y_slice    = y[depth_index]
    z_slice    = z[depth_index]
    data_slice = data[depth_index]
    
    # Now convert the lat/lon coordinates of the point to global Cartesian for direct comparison to ASPECT data
    point_r = 6371e3 - depth
    point_x = point_r * np.sin(np.deg2rad(90 - point_lats)) * np.cos(np.deg2rad(point_lons))
    point_y = point_r * np.sin(np.deg2rad(90 - point_lats)) * np.sin(np.deg2rad(point_lons))
    point_z = point_r * np.cos(np.deg2rad(90 - point_lats))
    
    # Create a KDTree with the slice through the ASPECT models
    slice_KDTree = KDTree(np.c_[x_slice, y_slice, z_slice])
    
    # Loop through all lat/lon coordinates and use the KDTree to find the closest point
    for p in range(len(point_x)):
        dd, tree_index = slice_KDTree.query([point_x[p], point_y[p], point_z[p]], k=1)
        ASPECT_x[p] = x_slice[tree_index]
        ASPECT_y[p] = y_slice[tree_index]
        ASPECT_z[p] = z_slice[tree_index]
        ASPECT_v[p] = data_slice[tree_index]
    
    # Lastly, convert back to spherical for direct comparison in geographic comparison
    theta = 90 - np.rad2deg( np.arccos( ASPECT_z / (np.sqrt(ASPECT_x**2 + ASPECT_y**2 + ASPECT_z**2)) ) )
    phi =  np.sign(ASPECT_y) * np.rad2deg(np.arccos( ASPECT_x / np.sqrt(ASPECT_x**2 + ASPECT_y**2) ))
    # phi[np.where(phi < 0)] = phi[np.where(phi < 0)] + 360
    # phi[np.where(phi == 0)] = 180
    
    return theta, phi, ASPECT_v

def velocity_converter(theta, phi, surface_velocity):
    spherical_velocity = np.zeros(len(surface_velocity), dtype=object)
    
    for i in range(len(spherical_velocity)):
        theta_rad = np.deg2rad(theta[i])
        phi_rad = np.deg2rad(phi[i])

        i_hat = np.array([np.sin(theta_rad) * np.cos(phi_rad), np.cos(theta_rad) * np.cos(phi_rad), -np.sin(phi_rad)])
        j_hat = np.array([np.sin(theta_rad) * np.sin(phi_rad), np.cos(theta_rad) * np.sin(phi_rad), np.cos(phi_rad)])
        k_hat = np.array([np.cos(theta_rad), -np.sin(theta_rad), phi_rad * 0])
        spherical_velocity[i] = surface_velocity[i][0] * i_hat + surface_velocity[i][1] * j_hat + surface_velocity[i][2] * k_hat

        # r_hat = np.array([np.sin(theta_rad) * np.cos(phi_rad), np.sin(phi_rad) * np.sin(theta_rad), np.cos(theta_rad)])
        # theta_hat = np.array([np.cos(theta_rad) * np.cos(phi_rad), np.cos(theta_rad) * np.sin(phi_rad), -np.sin(theta_rad)])
        # phi_hat = np.array([-np.sin(phi_rad), np.cos(phi_rad), 0])
        # spherical_velocity[i] = np.array([surface_velocity[i][0] * r_hat + surface_velocity[i][1] * theta_hat + surface_velocity[i][2] * phi_hat])
    
    return spherical_velocity
    