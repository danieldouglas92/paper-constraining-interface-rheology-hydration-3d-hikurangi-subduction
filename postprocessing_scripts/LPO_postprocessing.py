import matplotlib.pyplot as plt
import numpy as np
import os
import pygmt
from scipy.spatial import KDTree

def slice_donna_data(data, orientation, constant_value):
    """
    Takes the anisotropy data file provided by Donna Eberhart-Philips that was published
    in her 2009 EPR paper. This function returns a slice through this dataset based on an
    orientation of normal direction to the slice ('x', 'y', or 'z'), and the constant value
    of the plane.
    """
    if orientation == 'x':
        slice_index = np.where(data[:, 2] == constant_value)
    elif orientation == 'y':
        slice_index = np.where(data[:, 3] == constant_value)
    elif orientation == 'z':
        slice_index = np.where(data[:, 4] == constant_value)
    return data[:, 0][slice_index], data[:, 1][slice_index], \
           data[:, 2][slice_index], data[:, 3][slice_index], data[:, 4][slice_index], \
           data[:, 5][slice_index], data[:, 6][slice_index], data[:, 7][slice_index]

def slice_ASPECT(x, y, z, data, orientation, constant_value, padding):
    
    if orientation == 'x':
        slice_index = np.where( (x <= constant_value + padding) & (x >= constant_value - padding) )
    elif orientation == 'y':
        slice_index = np.where( (y <= constant_value + padding) & (y >= constant_value - padding) )
    elif orientation == 'z':
        slice_index = np.where( (z/1e3 <= constant_value + padding) & (z/1e3 >= constant_value - padding) )
        
    slice_x = x[slice_index]
    slice_y = y[slice_index]
    slice_z = z[slice_index]
    slice_data = data[slice_index]
    
    return slice_x, slice_y, slice_z, slice_data

def cartesian_to_spherical(x, y, z):
    """
    Takes an x, y, z and converts it to spherical coordinates. Returns r, theta, phi
    """
    r = np.sqrt(x**2 + y**2 + z**2)
    theta = 90 - np.rad2deg( np.arccos( z / (np.sqrt(x**2 + y**2 + z**2)) ) )
    phi =  np.sign(y) * np.rad2deg(np.arccos( x / np.sqrt(x**2 + y**2) ))
    phi[np.where(phi < 0)] = phi[np.where(phi < 0)] + 360
    phi[np.where(phi == 0)] = 180
    
    return r, theta, phi

def spherical_to_global_cartesian(r, theta, phi):
    """
    Takes spherical coordinates r, theta, and phi and converts to Cartesian coordinates.
    Returns x, y, z
    """
    x = r * np.sin(np.deg2rad(90 - theta)) * np.cos(np.deg2rad(phi))
    y = r * np.sin(np.deg2rad(90 - theta)) * np.sin(np.deg2rad(phi))
    z = r * np.cos(np.deg2rad(90 - theta))
    
    return x, y, z

def spherical_chunk_extractor(x, y, z, field_values, lon_bounds, lat_bounds, r_bounds):
    """
    Takes the x, y, and z coordinates of all points (mesh or particles) from the ASPECT model, as well as the desired field, and
    extracts only points within specified longitude, latitude, and radius bounds. Returns the longitude, latitude, radius,
    and field of the ASPECT model within these bounds. 
    """
    ASPECT_radius, ASPECT_theta, ASPECT_phi = cartesian_to_spherical(x, y, z)
    
    bounded_index = np.where( (ASPECT_radius <= np.max(r_bounds)) & (ASPECT_radius >= np.min(r_bounds)) & \
                              (ASPECT_phi <= np.max(lon_bounds)) & (ASPECT_phi >= np.min(lon_bounds)) & \
                              (ASPECT_theta <= np.max(lat_bounds)) & (ASPECT_theta >= np.min(lat_bounds)) )
    
    return ASPECT_phi[bounded_index], ASPECT_theta[bounded_index], ASPECT_radius[bounded_index], field_values[bounded_index]


def global_cartesian_anisotropy(x_particles, y_particles, z_particles, field_particles, slice_longitude, slice_latitude, slice_depth):
    """
    This function first converts the data provided by Donna Eberhart-Philips and converts it to the same global Cartesian coordinate system
    as the ASPECT model. Then, this function iterates through all of Donna's points in the global coordinate system and finds the closet ASPECT
    particle using a KDTree. Then, return the x,y,z coordinates of both Donna's points and ASPECT particles in the global cartesian coordinate
    system, and the particle field.
    """
    # Set up the bounds based on the slice above to extract particles from the ASPECT model within these bounds
    slice_radius = 6371e3 - slice_depth * 1e3
    lon_bounds = np.array([np.min(slice_longitude), np.max(slice_longitude)])
    lat_bounds = np.array([np.min(slice_latitude), np.max(slice_latitude)])
    r_bounds = np.array([6000e3, 6371e3])

    # Use the spherical chunk extractor to pull data from the entire model domain that are between the lon, lat, and radius bounds
    # of the slice through Donna's data
    lon_chunk, lat_chunk, r_chunk, field_chunk = spherical_chunk_extractor(x_particles, y_particles, z_particles, \
                                                                           field_particles, lon_bounds, lat_bounds, r_bounds)

    # Convert the extract Lon, Lat, R values to Cartesian, this is important for the KDTree, as x, y, and z must
    # all be the same units
    x_chunk, y_chunk, z_chunk = spherical_to_global_cartesian(r_chunk, lat_chunk, lon_chunk)

    # Create the KDTree from the Cartesian values calculated above (ASPECT points)
    chunk_KDTree = KDTree(np.c_[x_chunk, y_chunk, z_chunk])

    # Convert Donna's local coordinate system to the global Coordinate system that ASPECT uses, for use in the
    # KDTree
    anisotropy_global_X = slice_radius * np.sin(np.deg2rad(90 - slice_latitude)) * np.cos(np.deg2rad(slice_longitude))
    anisotropy_global_Y = slice_radius * np.sin(np.deg2rad(90 - slice_latitude)) * np.sin(np.deg2rad(slice_longitude))
    anisotropy_global_Z = slice_radius * np.cos(np.deg2rad(90 - slice_latitude))

    # Initialize arrays for storage
    ASPECT_ANISOTROPY_X = np.empty(0)
    ASPECT_ANISOTROPY_Y = np.empty(0)
    ASPECT_ANISOTROPY_Z = np.empty(0)

    ASPECT_X = np.empty(0)
    ASPECT_Y = np.empty(0)
    ASPECT_Z = np.empty(0)

    # Iterate through all of Donna's data points and find the closest point in the ASPECT models to Donna's data points
    for p in range(len(anisotropy_global_X)):
        dd, tree_index = chunk_KDTree.query([anisotropy_global_X[p], anisotropy_global_Y[p], anisotropy_global_Z[p]], k=1)
        if len(field_chunk.T) == 3:
            ASPECT_ANISOTROPY_X = np.append(ASPECT_ANISOTROPY_X, field_chunk[tree_index][0])
            ASPECT_ANISOTROPY_Y = np.append(ASPECT_ANISOTROPY_Y, field_chunk[tree_index][1])
            ASPECT_ANISOTROPY_Z = np.append(ASPECT_ANISOTROPY_Z, field_chunk[tree_index][2])

            ASPECT_X = np.append(ASPECT_X, x_chunk[tree_index])
            ASPECT_Y = np.append(ASPECT_Y, y_chunk[tree_index])
            ASPECT_Z = np.append(ASPECT_Z, z_chunk[tree_index])
        else:
            ASPECT_ANISOTROPY_X = np.append(ASPECT_ANISOTROPY_X, field_chunk[tree_index])
            ASPECT_ANISOTROPY_Y = np.append(ASPECT_ANISOTROPY_Y, field_chunk[tree_index])
            ASPECT_ANISOTROPY_Z = np.append(ASPECT_ANISOTROPY_Z, field_chunk[tree_index])
            
            ASPECT_X = np.append(ASPECT_X, x_chunk[tree_index])
            ASPECT_Y = np.append(ASPECT_Y, y_chunk[tree_index])
            ASPECT_Z = np.append(ASPECT_Z, z_chunk[tree_index])
        
    return np.array([ASPECT_X, ASPECT_Y, ASPECT_Z]).T, \
           np.array([ASPECT_ANISOTROPY_X, ASPECT_ANISOTROPY_Y, ASPECT_ANISOTROPY_Z]).T, \
           np.array([anisotropy_global_X, anisotropy_global_Y, anisotropy_global_Z]).T

def local_cartesian_coordinates(x_global, y_global, z_global, field_global, azimuth, center, orientation):
    """
    Takes the global cartesian coordinates output by global_cartesian_anisotropy and first converts them to spherical
    coordinates. Then, uses pygmt to project the spherical coordinates onto a local cartesian grid centered around the
    North Island of New Zealand. This projection is specified by an azimuth (clockwise starting from North), and lon,lat
    point describing the center that the azimuth rotates around.
    """
    ASPECT_radius, ASPECT_theta, ASPECT_phi = cartesian_to_spherical(x_global, y_global, z_global)
    
    data_to_project = np.array([ASPECT_phi, ASPECT_theta]).T
    np.savetxt(fname='load_data.txt', X=data_to_project)

    pygmt.project(data='load_data.txt', outfile='test.txt', \
                  center=center, azimuth=azimuth, convention='pq', flat_earth=False, unit=True, \
                  verbose=False)
    
    projected_data = np.loadtxt(fname='test.txt')
    
    if orientation == 'z':
        return projected_data[:, 0], projected_data[:, 1]
    elif orientation == 'y':
        return projected_data[:, 0], 6371e3 - ASPECT_radius.T
    elif orientation == 'x':
        return projected_data[:, 1], 6371e3 - ASPECT_radius.T
    elif orientation == 'xyz':

        data_to_project = np.array([ASPECT_phi, ASPECT_theta, ASPECT_radius, field_global]).T
        np.savetxt(fname='load_data.txt', X=data_to_project)

        pygmt.project(data='load_data.txt', outfile='test.txt', \
                      center=center, azimuth=azimuth, convention='pqz', flat_earth=False, unit=True, \
                      verbose=False)
        
        projected_data = np.loadtxt(fname='test.txt')
        
        return projected_data[:, 0], projected_data[:, 1], 6371e3 - projected_data[:, 2], projected_data[:, 3]
    else:
        return
    
    
    
def spherical_slice(x_solution, y_solution, z_solution, field_solution, tree_unstructured, r_min_max, lon_lat_slice, spacing, slice_type):
    """
    This function takes the mesh output from a spherical ASPECT model and takes a slice through it, and outputs the slice
    in Cartesian coordinates.
        x_solution, y_solution, z_solution - The Cartesian coordinates of the spherical model
        field_solution - The ASPECT solution variable to be sliced
        tree_unstructured - A KDTree generated from the ASPECT mesh points x_solution, y_solution, z_solution
        r_min_max - The radial extent of the slice
        lon_lat_slice - The longitudinal extent (if slice_type == 'Longitude'), latitudinal extent (if slice_type == 'Latitude')
        or both (if slice_type != 'Longitude' | slice_type != 'Latitude')
        slice_type defines whether to slice along a constant longitude, latitude, or neither.
    """
    ASPECT_radius = np.sqrt(x_solution**2 + y_solution**2 + z_solution**2)
    ASPECT_theta = 90 - np.rad2deg( np.arccos( z_solution / (np.sqrt(x_solution**2 + y_solution**2 + z_solution**2)) ) )
    ASPECT_phi =  np.sign(y_solution) * np.rad2deg(np.arccos( x_solution / np.sqrt(x_solution**2 + y_solution**2) ))
    ASPECT_phi[np.where(ASPECT_phi < 0)] = ASPECT_phi[np.where(ASPECT_phi < 0)] + 360
    ASPECT_phi[np.where(ASPECT_phi == 0)] = 180
    
    if slice_type == 'Longitude':
        slice_index = np.where( (ASPECT_phi <= lon_lat_slice[0] + spacing) & (ASPECT_phi >= lon_lat_slice[0] - spacing) & \
                                (ASPECT_radius > np.max(r_min_max) - 500) )
        
    elif slice_type == 'Latitude':
        slice_index = np.where( (ASPECT_theta <= lon_lat_slice[1] + spacing) & (ASPECT_theta >= lon_lat_slice[1] - spacing) & \
                                (ASPECT_radius > np.max(r_min_max) - 500) )
        
    else:
        m, b = np.polyfit(lon_lat_slice[:, 0], lon_lat_slice[:, 1], deg=1)
        lon_array = np.linspace(np.min(ASPECT_phi), np.max(ASPECT_phi), 10000)
        lat_array = m * lon_array + b
        circle_bounded_by_model = np.array([lon_array[np.where( (lat_array <= np.max(ASPECT_theta)) & (lat_array >= np.min(ASPECT_theta)))], \
                                            lat_array[np.where( (lat_array <= np.max(ASPECT_theta)) & (lat_array >= np.min(ASPECT_theta)))]]).T
        
        # fwd_az1, back_az1, dist1 = geodesic.inv(circle_bounded_by_model[0][0], circle_bounded_by_model[0][1], \
        #                                         circle_bounded_by_model[-1][0], circle_bounded_by_model[-1][1])
        
        x_slice = np.max(r_min_max) * np.sin(np.deg2rad(90 - circle_bounded_by_model[:, 1])) * np.cos(np.deg2rad(circle_bounded_by_model[:, 0]))
        y_slice = np.max(r_min_max) * np.sin(np.deg2rad(90 - circle_bounded_by_model[:, 1])) * np.sin(np.deg2rad(circle_bounded_by_model[:, 0]))
        z_slice = np.max(r_min_max) * np.cos(np.deg2rad(90 - circle_bounded_by_model[:, 1]))
        
        slice_index = np.zeros(len(circle_bounded_by_model), dtype=int)
        for p in range(len(circle_bounded_by_model)):
            dd, tree_index = tree_unstructured.query([x_slice[p], y_slice[p], z_slice[p]])
            slice_index[p] = int(tree_index)
       
    radius_vals = np.linspace(np.max(r_min_max), np.min(r_min_max), 100)

    ave_lon = np.average(ASPECT_phi[slice_index])
    ave_lat = np.average(ASPECT_theta[slice_index])
    ave_rad = np.average(radius_vals)

    ave_x = ave_rad * np.sin(np.deg2rad(90 - ave_lat)) * np.cos(np.deg2rad(ave_lon))
    ave_y = ave_rad * np.sin(np.deg2rad(90 - ave_lat)) * np.sin(np.deg2rad(ave_lon))
    ave_z = ave_rad * np.cos(np.deg2rad(90 - ave_lat))

    rotation_angle_around_z_axis = np.arctan(-ave_y / ave_x) 
    rotated_ave_x = ave_x * np.cos(rotation_angle_around_z_axis) - ave_y * np.sin(rotation_angle_around_z_axis)

    rotation_angle_around_y_axis = np.pi - np.arctan(rotated_ave_x / ave_z)

    x_interpolated = np.empty(0)
    y_interpolated = np.empty(0)
    xy_interpolated = np.empty(0)
    z_interpolated = np.empty(0)
    
    if field_solution.ndim == 1:
        field_interpolated = np.empty(0)
        for m in range(len(radius_vals)):

            field_holder = np.zeros(len(ASPECT_phi[slice_index]))
            x_vals = radius_vals[m] * np.sin(np.deg2rad(90 - ASPECT_theta[slice_index])) * np.cos(np.deg2rad(ASPECT_phi[slice_index]))
            y_vals = radius_vals[m] * np.sin(np.deg2rad(90 - ASPECT_theta[slice_index])) * np.sin(np.deg2rad(ASPECT_phi[slice_index]))
            z_vals = radius_vals[m] * np.cos(np.deg2rad(90 - ASPECT_theta[slice_index]))

            x_vals_after_z_rotation = x_vals * np.cos(rotation_angle_around_z_axis) - \
                                      y_vals * np.sin(rotation_angle_around_z_axis)
            y_vals_after_z_rotation = x_vals * np.sin(rotation_angle_around_z_axis) + \
                                      y_vals * np.cos(rotation_angle_around_z_axis)

            x_vals_after_y_rotation = x_vals_after_z_rotation * np.cos(rotation_angle_around_y_axis) + \
                                      z_vals * np.sin(rotation_angle_around_y_axis)
            z_vals_after_y_rotation = -x_vals_after_z_rotation * np.sin(rotation_angle_around_y_axis) + \
                                      z_vals * np.cos(rotation_angle_around_y_axis)
            for i in range(len(x_vals)):               
                dd, tree_index = tree_unstructured.query([x_vals[i], y_vals[i], z_vals[i]], k=1)
                field_holder[i] = field_solution[tree_index]

            x_interpolated = np.concatenate( (x_interpolated, x_vals_after_y_rotation), axis=0)
            y_interpolated = np.concatenate( (y_interpolated, y_vals_after_z_rotation), axis=0)
            xy_interpolated_abs = np.sqrt(x_vals_after_y_rotation**2 + y_vals_after_z_rotation**2)
            xy_interpolated_abs[np.where( (x_vals_after_y_rotation + y_vals_after_z_rotation) < 0 )] *= -1
            xy_interpolated = np.concatenate( (xy_interpolated, xy_interpolated_abs) , axis=0)
            z_interpolated = np.concatenate( (z_interpolated, z_vals_after_y_rotation) , axis=0)
            field_interpolated = np.concatenate( (field_interpolated, field_holder.T))
  

    else:
        field_interpolated = np.empty( (0, len(field_solution)) )
        for m in range(len(radius_vals)):

            field_holder = np.zeros((len(field_solution), len(ASPECT_phi[slice_index])))
            x_vals = radius_vals[m] * np.sin(np.deg2rad(90 - ASPECT_theta[slice_index])) * np.cos(np.deg2rad(ASPECT_phi[slice_index]))
            y_vals = radius_vals[m] * np.sin(np.deg2rad(90 - ASPECT_theta[slice_index])) * np.sin(np.deg2rad(ASPECT_phi[slice_index]))
            z_vals = radius_vals[m] * np.cos(np.deg2rad(90 - ASPECT_theta[slice_index]))

            x_vals_after_z_rotation = x_vals * np.cos(rotation_angle_around_z_axis) - \
                                      y_vals * np.sin(rotation_angle_around_z_axis)
            y_vals_after_z_rotation = x_vals * np.sin(rotation_angle_around_z_axis) + \
                                      y_vals * np.cos(rotation_angle_around_z_axis)

            x_vals_after_y_rotation = x_vals_after_z_rotation * np.cos(rotation_angle_around_y_axis) + \
                                      z_vals * np.sin(rotation_angle_around_y_axis)
            z_vals_after_y_rotation = -x_vals_after_z_rotation * np.sin(rotation_angle_around_y_axis) + \
                                      z_vals * np.cos(rotation_angle_around_y_axis)
            for q in range(len(field_solution)):
                for i in range(len(x_vals)):               
                    dd, tree_index = tree_unstructured.query([x_vals[i], y_vals[i], z_vals[i]], k=1)
                    field_holder[q][i] = field_solution[q][tree_index]

            x_interpolated = np.concatenate( (x_interpolated, x_vals_after_y_rotation), axis=0)
            y_interpolated = np.concatenate( (y_interpolated, y_vals_after_z_rotation), axis=0)
            xy_interpolated_abs = np.sqrt(x_vals_after_y_rotation**2 + y_vals_after_z_rotation**2)
            xy_interpolated_abs[np.where( (x_vals_after_y_rotation + y_vals_after_z_rotation) < 0 )] *= -1
            xy_interpolated = np.concatenate( (xy_interpolated, xy_interpolated_abs) , axis=0)
            z_interpolated = np.concatenate( (z_interpolated, z_vals_after_y_rotation) , axis=0)
            field_interpolated = np.concatenate( (field_interpolated, field_holder.T))

    return x_interpolated, y_interpolated, z_interpolated, xy_interpolated, field_interpolated
    


## LOOK AT THE CHOICE OF INDEXING FOR THE KDTREE
## RIGHT NOW I AM CREATING A KDTREE FROM THE STRUCTURED GRID BUT I THINK
## THIS IS INCORRECT. I SHOULD PROBABLY BE DOING THE KDTREE OVER THE 
## UNSTRUCTURED GRID AND THEN ITERATING THROUGH X,Y CREATED BY NP.MESHGRID
## THE ISSUE IS THAT THE Z VALUES ARE ONLY ON THE UNSTRUCTURED GRID, HOW TO GET
## THEM ONTO THE STRUCTURED GRID?
def bilinear_interpolation(x, y, z, dx, dy, padding):
    x_structured = np.arange(np.min(x), np.max(x) + dx, dx)
    y_structured = np.arange(np.min(y), np.max(y) + dy, dy)
    X, Y = np.meshgrid(x_structured, y_structured)
    x_raveled = np.ravel(X)
    y_raveled = np.ravel(Y)
    
    structured_tree = KDTree(np.c_[x_raveled, y_raveled])
    Z = np.zeros(len(x))
    for i in range(len(x)):
        point = np.array([x[i], y[i]])
        dd, tree_index = structured_tree.query(point, k=1)
        closest_point = np.array([x_raveled[tree_index], y_raveled[tree_index]])

        ###### IN THIS CASE WE ARE EXACTLY ON THE VERTEX OF A MESH AND DO NOT NEED TO INTERPOLATE #####
        if abs(np.sum(closest_point - point)) <= 1e-5:
            Z[i] = z[tree_index]
            break

        ##### IN THIS CASE WE ARE ALONG AN EDGE OF THE MESH WITH A CONSTANT Y VALUE #####
        elif abs(closest_point - point)[1] <= 1e-5:

            ######## IF THE POINT IS CLOSER TO THE LEFT SIDE OF THE CELL ######
            if (closest_point - point)[1] < 0:
                greater_than_index = np.where( (x_raveled > point[0]) & (y_raveled == point[1]) )
                next_point = np.min(x_raveled[greater_than_index])
                next_index = np.where( (x_raveled == next_point) & (y_raveled == point[1]) )[0][0]
                Z[i] = (z[tree_index] * (next_point - point[0]) + z[next_index] * (point[0] - closest_point[0])) / (next_point - closest_point[0])
                break

            ####### IF THE POINT IS CLOSER TO THE RIGHT SIDE OF THE CELL #######
            else:
                less_than_index = np.where( (x_raveled < point[0]) & (y_raveled == point[1]) )
                next_point = np.max(x_raveled[less_than_index])
                next_index = np.where( (x_raveled == next_point) & (y_raveled == point[1]) )[0][0]
                Z[i] = (z[tree_index] * (point[0] - next_point) + z[next_index] * (closest_point[0] - point[0])) / (closest_point[0] - next_point)
                break

        ##### IN THIS CASE WE ARE ALONG AN EDGE OF THE MESH WITH A CONSTANT X VALUE #####
        elif abs(closest_point - point)[0] <= 1e-5:

            ######### IF THE POINT IS CLOSER TO THE BOTTOM SIDE OF THE CELL #######
            if (closest_point - point)[0] < 0:
                greater_than_index = np.where( (y_raveled > point[1]) & (x_raveled == point[0]) )
                next_point = np.max(y_raveled[greater_than_index])
                next_index = np.where( (y_raveled == next_point) & (x_raveled == point[0]) )[0][0]
                Z[i] = (z[tree_index] * (next_point - point[1]) + z[next_index] * (point[1] - closest_point[1])) / (next_point - closest_point[1])
                break

            ######### IF THE POINT IS CLOSER TO THE TOP SIDE OF THE CELL #######
            else:
                less_than_index = np.where( (y_raveled < point[1]) & (x_raveled == point[0]) )
                next_point = np.min(y_raveled[less_than_index])
                next_index = np.where( (y_raveled == next_point) & (x_raveled == point[0]) )[0][0]
                Z[i] = (z[tree_index] * (point[1] - next_point) + z[next_index] * (closest_point[1] - point[1])) / (closest_point[1] - next_point)
                break

        ###### IN THISE CASE WE ARE WITHIN THE INTERIOR OF THE CELL ############   
        else:
            if closest_point[0] - point[0] < 0:

                ####### In this case we are in the bottom left corner of the cell ########
                if closest_point[1] - point[1] < 0:

                    greater_than_x_index = np.where( (y_raveled == closest_point[1]) & (x_raveled > closest_point[0]) )
                    greater_than_y_index = np.where( (x_raveled == closest_point[0]) & (y_raveled > closest_point[1]) )
                    
                    print(point)
                    print(closest_point)
                    
                    next_x_point = np.array([np.min(x_raveled[greater_than_x_index]), closest_point[1]])
                    next_y_point = np.array([closest_point[0], np.min(y_raveled[greater_than_y_index])])
                    next_x_index = np.where( (x_raveled == next_x_point[0]) & (y_raveled == closest_point[1]) )
                    next_y_index = np.where( (y_raveled == next_y_point[1]) & (x_raveled == closest_point[0]) )
                    farthest_point = np.array([next_x_point[0], next_y_point[1]])
                    farthest_index = np.where( (x_raveled == farthest_point[0]) & (y_raveled == farthest_point[1]) )

                ####### In this case we are in the top left corner of the cell ########   
                elif closest_point[1] - point[1] > 0:

                    greater_than_x_index = np.where( (y_raveled == closest_point[1]) & (x_raveled > closest_point[0]) )
                    less_than_y_index    = np.where( (x_raveled == closest_point[0]) & (y_raveled < closest_point[1]) )                        
                    next_x_point = np.array([np.min(x_raveled[greater_than_x_index]), closest_point[1]])

                    print('Point is ' + str(point))
                    print('Closest point is ' + str(closest_point))

                    next_y_point = np.array([closest_point[0], np.max(y_raveled[less_than_y_index])])
                    next_x_index = np.where( (x_raveled == next_x_point[0]) & (y_raveled == closest_point[1]) )
                    next_y_index = np.where( (y_raveled == next_y_point[1]) & (x_raveled == closest_point[0]) )
                    farthest_point = np.array([next_x_point[0], next_y_point[1]])
                    farthest_index = np.where( (x_raveled == farthest_point[0]) & (y_raveled == farthest_point[1]) )

            elif closest_point[0] - point[0] > 0:

                ####### In this case we are in the bottom right corner of the cell ########
                if closest_point[1] - point[1] < 0:

                    less_than_x_index    = np.where( (y_raveled == closest_point[1]) & (x_raveled < closest_point[0]) )
                    greater_than_y_index = np.where( (x_raveled == closest_point[0]) & (y_raveled > closest_point[1]) )
                    next_x_point = np.array([np.max(x_raveled[less_than_x_index]), closest_point[1]])
                    next_y_point = np.array([closest_point[0], np.min(y_raveled[greater_than_y_index])])
                    next_x_index = np.where( (x_raveled == next_x_point[0]) & (y_raveled == closest_point[1]) )
                    next_y_index = np.where( (y_raveled == next_y_point[1]) & (x_raveled == closest_point[0]) )
                    farthest_point = np.array([next_x_point[0], next_y_point[1]])
                    farthest_index = np.where( (x_raveled == farthest_point[0]) & (y_raveled == farthest_point[1]) )

                ####### In this case we are in the top right corner of the cell ########
                elif closest_point[1] - point[1] > 0:

                    less_than_x_index = np.where( (y_raveled == closest_point[1]) & (x_raveled < closest_point[0]) )
                    less_than_y_index = np.where( (x_raveled == closest_point[0]) & (y_raveled < closest_point[1]) )
                    next_x_point = np.array([np.max(x_raveled[less_than_x_index]), closest_point[1]])
                    next_y_point = np.array([closest_point[0], np.max(y_raveled[less_than_y_index])])
                    next_x_index = np.where( (x_raveled == next_x_point[0]) & (y_raveled == closest_point[1]) )
                    next_y_index = np.where( (y_raveled == next_y_point[1]) & (x_raveled == closest_point[0]) )
                    farthest_point = np.array([next_x_point[0], next_y_point[1]])
                    farthest_index = np.where( (x_raveled == farthest_point[0]) & (y_raveled == farthest_point[1]) )

            #### Now we need to calculate the areas of the four boxes that are split by the location of the point is the cell ####
            #### A1 is always assumed to be the smallest area, A4 is always the biggest, and A2 and A3 are intermediate ##########
            A1 = abs(closest_point[0] - point[0]) * abs(closest_point[1] - point[1])
            A2 = abs(next_x_point[0] - point[0]) * abs(next_x_point[1] - point[1])
            A3 = abs(next_y_point[0] - point[0]) * abs(next_y_point[1] - point[1])
            A4 = abs(farthest_point[0] - point[0]) * abs(farthest_point[1] - point[1])
            A_total = abs(closest_point[0] - farthest_point[0]) * abs(closest_point[1] - farthest_point[1])

            z1 = z[tree_index]
            z2 = z[next_x_index]
            z3 = z[next_y_index]
            z4 = z[farthest_index]

            print(next_x_point)
            print(next_y_point)
            print(farthest_point)

            Z[i] = ((z1*A4) + (z2*A3) + (z3*A2) + (z4*A1)) / A_total
    return X, Y, Z