import numpy as np
from scipy.spatial import KDTree

def KDTree_mesh_interpolator(unstructured_x, unstructured_y, unstructured_z, unstructured_data, structured_array_x, structured_array_y, structured_array_z):
    """
    unstructured_coordinates is the 3D positional coordinates of the slice
    unstructured_data is the data field that is being interpolated over the coordinates in unstructured_coordinates
    structured_array_x is the structured array corresponding to the first index of unstructured_coordinates
    structured_array_y is the structured array corresponding to the second index of unstructured_coordinates
    structured_array_z is the structured array corresponding to the third index of unstructured_coordinates
    """
    
    tree_unstructured = KDTree(np.c_[unstructured_x, unstructured_y, unstructured_z])   
    field_interpolated = np.zeros(len(structured_array_x))
        
    for i in range(len(field_interpolated)):
        dd, tree_index = tree_unstructured.query([structured_array_x[i], structured_array_y[i], structured_array_z[i]], k=1)
        field_interpolated[i] = unstructured_data[tree_index]
            
    return field_interpolated


def spherical_slice(x_solution, y_solution, z_solution, field_solution, tree_unstructured, r_min_max, radius_spacing, lon_lat_slice, lon_lat_spacing, spacing, slice_type):
    """
    solution_coordinates are the 3D positional coordinates of the slice
    field_solution is the data field that is being interpolated over the coordinates in unstructured_coordinates
    tree_unstructured is the KD_Tree of the unstructured coordinates
    r_min_max are the radius bounds
    radius_spacing is the number of radius points
    lon_lat_slice is an array of 2 lon lat points that specify the surface projection of the slice
    spacing is the +/- lon/lat range for searching for points within the slice
    slice_type: 'Longitude' for constant longitude, 'Latitude' for constant latitude, or 'Both' for variable lon/lat
    
    returns: the coordinates and field on the slice
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
        lon_array = np.linspace(np.min(ASPECT_phi), np.max(ASPECT_phi), lon_lat_spacing)
        lat_array = m * lon_array + b
        circle_bounded_by_model = np.array([lon_array[np.where( (lat_array <= np.max(ASPECT_theta)) & (lat_array >= np.min(ASPECT_theta)))], \
                                            lat_array[np.where( (lat_array <= np.max(ASPECT_theta)) & (lat_array >= np.min(ASPECT_theta)))]]).T
        
        x_slice = np.max(r_min_max) * np.sin(np.deg2rad(90 - circle_bounded_by_model[:, 1])) * np.cos(np.deg2rad(circle_bounded_by_model[:, 0]))
        y_slice = np.max(r_min_max) * np.sin(np.deg2rad(90 - circle_bounded_by_model[:, 1])) * np.sin(np.deg2rad(circle_bounded_by_model[:, 0]))
        z_slice = np.max(r_min_max) * np.cos(np.deg2rad(90 - circle_bounded_by_model[:, 1]))
        
        slice_surface_index = np.zeros(len(circle_bounded_by_model), dtype=int)
        for p in range(len(circle_bounded_by_model)):
            dd, tree_index = tree_unstructured.query([x_slice[p], y_slice[p], z_slice[p]])
            slice_surface_index[p] = int(tree_index)
    
    radius_vals = np.linspace(np.max(r_min_max), np.min(r_min_max), radius_spacing)
    
    ave_lon = np.average(ASPECT_phi[slice_surface_index])
    ave_lat = np.average(ASPECT_theta[slice_surface_index])
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
    field_interpolated = np.empty( (0, len(field_solution)) )

    for m in range(len(radius_vals)):
        
        field_holder = np.zeros((len(field_solution), len(ASPECT_phi[slice_surface_index])))
        x_vals = radius_vals[m] * np.sin(np.deg2rad(90 - ASPECT_theta[slice_surface_index])) * np.cos(np.deg2rad(ASPECT_phi[slice_surface_index]))
        y_vals = radius_vals[m] * np.sin(np.deg2rad(90 - ASPECT_theta[slice_surface_index])) * np.sin(np.deg2rad(ASPECT_phi[slice_surface_index]))
        z_vals = radius_vals[m] * np.cos(np.deg2rad(90 - ASPECT_theta[slice_surface_index]))

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

        x_interpolated = np.concatenate( (x_interpolated, x_vals_after_y_rotation) , axis=0)
        y_interpolated = np.concatenate( (y_interpolated, y_vals_after_z_rotation) , axis=0)
        xy_interpolated_abs = np.sqrt(x_vals_after_y_rotation**2 + y_vals_after_z_rotation**2)
        xy_interpolated_abs[np.where( (x_vals_after_y_rotation + y_vals_after_z_rotation) < 0 )] *= -1
        xy_interpolated = np.concatenate( (xy_interpolated, xy_interpolated_abs) , axis=0)
        z_interpolated = np.concatenate( (z_interpolated, z_vals_after_y_rotation) , axis=0)
        field_interpolated = np.concatenate( (field_interpolated, field_holder.T))

    return x_interpolated, y_interpolated, z_interpolated, xy_interpolated, field_interpolated