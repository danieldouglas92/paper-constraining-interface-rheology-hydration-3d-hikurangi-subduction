##########################################################################################################################
#################################### FUNCTIONS USED TO MAKE PROFILES ON SLAB2.0 DATA #####################################
##########################################################################################################################

import numpy as np
import geopy
import pygmt
import os
import pyproj
from geopy.distance import geodesic
from geopy.distance import distance
from os import system as sys

def profile_generator(slab_depth_grid, slab_depth_contour, output_dir, depth_of_contour, num_of_profiles, profile_spacing, region_lon, region_lat,\
                      lon_spacing, lat_spacing, max_point_distance, rotation_angle, trench_orientation, use_trench):
    
    '''
    Takes a depth contour and the depth-to-slab-surface data from the slab2.0 database, draws downdip profiles along-strike of the depth contour,
    extracts the depth-to-slab-surface along these profiles, and writes them to files which can then be used to create a worldbuilder file.

    slab_depth_grid    = the slab2.0 file containing lon, lat, and depth to slab surface
    slab_depth_contour = the slab2.0 file containing lon, lat, and depth contours of the slab surface
                         alternatively, this is the contour defining the trench.
    depth_of_contour   = the depth, in km, of the contour to be used to generate the profiles. Note slab2.0 specifies depths as negative
    output_dir         = pathway to the directory where the profiles are written
    num_of_profiles    = the number of profiles to be generated along the slab
    profile_spacing    = the distance between points along the profiles, in km
    region_lon         = the longitudinal extent of the slab
    region_lat         = the latitudinal extent of the slab
    lon_spacing        = the longitudinal spacing around a profile point where slab data will be considered
    lat_spacing        = the latitudinal spacing around a profile point where slab data will be considered
    max_point_distance = the maximum allowed distance from a profile point to a slab data point for it to 
                         be considered
    rotation_angle     = angle in degrees to rotate the profiles by
    trench_orientation = the direction of the trench relative to the OVERRIDING plate. For example, Cascadia would be 'West'.
                         Valid options are 'East', 'West', 'North', 'South'
    use_trench         = True or False, whether to use the trench data from slab2.0 to define the depth contour
    '''

    # Create the directories for storing the profiles (and remove old profiles for clean output)
    sys('mkdir ' + output_dir)
    sys('rm -rf ' + output_dir + '/geographic')
    sys('rm -rf ' + output_dir + '/cartesian')

    sys('mkdir ' + output_dir + '/geographic')
    sys('mkdir ' + output_dir + '/cartesian')
    
    # Define the points along strike of the depth contour (if use_trench = False) or the trench (if use_trench = True) to be 
    # used in the profile_locator() function later.
    if use_trench == False:
        
        # First constrain the depth contour to be within the specified lon/lat bounds
        slab_depth_contour_cut = slab_depth_contour[ np.where( (slab_depth_contour[:, 0] <= np.max(region_lon)) & (slab_depth_contour[:, 0] >= np.min(region_lon)) \
                                                        & (slab_depth_contour[:, 1] <= np.max(region_lat)) & (slab_depth_contour[:, 1] >= np.min(region_lat)) ) ]   

        # Extract the contour at the specified depth
        contour_indices = np.where(slab_depth_contour_cut[:, 2] == depth_of_contour)
        contour_arr = np.flipud(slab_depth_contour_cut[contour_indices])
        
        # Use the length of the specified depth contour and the desired number of profiles to create an interval for along-strike profiles
        contour_points = contour_arr.shape[0]
        contour_interval = np.floor(contour_points / num_of_profiles)
    
    elif use_trench == True:
        
        # First constrain the depth contour to be within the specified lon/lat bounds
        contour_arr = slab_depth_contour[ np.where( (slab_depth_contour[:, 0] <= np.max(region_lon)) & (slab_depth_contour[:, 0] >= np.min(region_lon)) \
                                                           & (slab_depth_contour[:, 1] <= np.max(region_lat)) & (slab_depth_contour[:, 1] >= np.min(region_lat)) ) ]
        
        # Use the length of the depth_contour and the desired number of profiles to create an interval for along-strike profiles
        contour_points = contour_arr.shape[0]
        contour_interval = np.floor(contour_points / num_of_profiles)
    
    # Slab2.0 data may contain NaNs, remove them
    slab_data_no_nans = []
    for j in range(len(slab_depth_grid)):
        if np.isnan(slab_depth_grid[j, 2]) == False:
            slab_data_no_nans.append(slab_depth_grid[j])
    slab_data_no_nans = np.array(slab_data_no_nans)
    
    # Define a max length for profiles for computational efficiency
    max_profile_length = 1000
    profile_within_slab = [0, 0]
    
    # Call profile_locator() and depth_extractor() to finish writing profiles to text files
    for i in np.arange(0, contour_points, contour_interval, dtype=int):
        combine_profile = profile_locator(contour_points, contour_arr, profile_spacing, max_profile_length, trench_orientation, i)
        
        profile_within_slab = depth_extractor(combine_profile, output_dir, lon_spacing, lat_spacing, slab_data_no_nans, region_lon, region_lat, rotation_angle, \
                                              max_point_distance, i)



def profile_locator(contour_points, contour_arr, spacing_of_profiles, max_profile_length, trench_orientation, i):
    
    '''
    This function generates profiles down-dip of the slab at a given point along a specified slab2.0 depth contour
    
    contour_points      = number of points that make up the depth contour array
    contour_arr         = array containing lon, lat, depth of a specified depth contour
    spacing_of_profiles = distance between points on the profiles
    max_profile_length  = the length of the profiles
    i                   = index of the current profile
    '''
    from geopy.distance import geodesic
    geodesic = pyproj.Geod(ellps='WGS84')
  
    if i == (contour_points - 1):
        fwd_az1, back_az1, dist1 = geodesic.inv(contour_arr[i-1, 0], contour_arr[i-1, 1], contour_arr[i, 0], contour_arr[i, 1]) # azimuth 1
        fwd_az2 = fwd_az1
        # The following Point (this condition is for the first point of the profile)
    elif i == 0:
        fwd_az1, back_az1, dist1 = geodesic.inv(contour_arr[i, 0], contour_arr[i, 1], contour_arr[i+1, 0], contour_arr[i + 1, 1]) # azimuth 1
        fwd_az2 = fwd_az1
        # All other points, using the previous and following point
    else:
        fwd_az1, back_az1, dist1 = geodesic.inv(contour_arr[i-1, 0], contour_arr[i-1, 1], contour_arr[i, 0], contour_arr[i, 1]) # azimuth 1
        fwd_az2, back_az2, dist2 = geodesic.inv(contour_arr[i, 0], contour_arr[i, 1], contour_arr[i+1, 0], contour_arr[i+1, 1]) # azimuth 2

        # Average azimuth, which is used to create a great circle profile perpendicular to the depth contour 
    az = (fwd_az1 + fwd_az2)/2.

    # Create two great circle profiles perpendicular (+/-90.) to the strike of the trench
    pygmt.project(center=[contour_arr[i, 0], contour_arr[i, 1]], outfile='positive_profile.xyd', azimuth=az + 90., generate=spacing_of_profiles, length=[0, max_profile_length], unit=True)
    pygmt.project(center=[contour_arr[i, 0], contour_arr[i, 1]], outfile='negative_profile.xyd', azimuth=az - 90., generate=spacing_of_profiles, length=[0, max_profile_length], unit=True)
    
    # Load the profiles to combine them into one array
    pos_profile = np.loadtxt('positive_profile.xyd')
    neg_profile = np.loadtxt('negative_profile.xyd')
    
    # Trench orientation influeces the way you combine the profiles
    if trench_orientation == 'East' or trench_orientation == 'North':
        combine_profile = np.concatenate( (np.flipud(pos_profile), np.delete(neg_profile, 0, axis=0)) ) # delete single duplicating point
    elif trench_orientation == 'West' or trench_orientation == 'South':
        combine_profile = np.concatenate( (np.flipud(neg_profile), np.delete(pos_profile, 0, axis=0)) ) # delete single duplicating point
    
    # Want longitudes to go from 0 - 360, not -180 - 180
    for k in range(len(combine_profile)):
        if combine_profile[:, 0][k] <= 0:
            combine_profile[:, 0][k] = combine_profile[:, 0][k] + 360
            
    sys('rm negative_profile.xyd')
    sys('rm positive_profile.xyd')
    return combine_profile
    

    

def depth_extractor(combine_profile, output_dir, lon_spacing, lat_spacing, slab_data, region_lon, region_lat, rotation_angle, max_point_distance, i):
    
    '''
    This function assigns depths to the profiles generated in the other functions from the slab 2.0 dataset.
    
    combine_profile    = Profile generated in profile_maker.py
    output_dir         = the pathway to the directory where the text files containing the profiles will be written
    lon_spacing        = specifies slab points that will be used within +/- longitude around a given profile point 
    lat_spacing        = specifies slab points that will be used within +/- latitude around a given profile point
    slab_data          = the slab2.0 depth file with nan values removed
    region_lon         = the longitudinal extent of the data you want
    region_lat         = the latitidunal extent of the data you want
    rotation_angle     = rotates the profiles by a given azimuth
    max_distance_point = only slab points that are less than or equal to this value will be considered
    i                  = index of the current profile
    '''
    
    from geopy.distance import distance
    # Calculate the distance between profile and slab point
    # Loop through all profile
    for j in range(len(combine_profile)):
        # Find all slab points within a specified lon/lat of a given profile point
        slab_indices = np.where( (combine_profile[j, 1] <= slab_data[:, 1] + lat_spacing) & (combine_profile[j, 1] >= slab_data[:, 1] - lat_spacing) & \
                                 (combine_profile[j, 0] <= slab_data[:, 0] + lon_spacing) & (combine_profile[j, 0] >= slab_data[:, 0] - lon_spacing))[0]
 
        # If there aren't any slab points, return 0
        if len(slab_indices) == 0:
            combine_profile[j] = 0
            profile_slab_distances = [0]
            
        # Otherwise calculate distance for all slab points
        else:
            profile_slab_distances  = [] #  Variable to store the distance between profile point and slab point
            for k in slab_indices:
                # Calculate the distance between the current profile and slab point
                profile_slab_distances.append(distance((combine_profile[j, 1], combine_profile[j, 0]), (slab_data[k, 1], slab_data[k, 0])).km)
            
            profile_slab_weights = np.zeros(len(profile_slab_distances)) # Variable for storing weights for each slab point
            
            # Iterate through each distance to compute a corresponding weight for later averaging
            for w in range(len(profile_slab_weights)):
                # if distance is very small, assign high weight (avoids potential divide by 0 error)
                if profile_slab_distances[w] <= 1e-10:
                    profile_slab_weights[w] = 1e10
                    
                # Otherwise weight is inverse of slab distance
                else:
                    profile_slab_weights[w] = 1./profile_slab_distances[w]
                    
                    # If the distance to point is greater than the allowed max_distance, discard by setting weight to 0 
                    if profile_slab_weights[w] < (1./ max_point_distance):
                        profile_slab_weights[w] = 0

            # If all weights are 0, set all weights to 1 to avoid averaging error
            if np.sum(profile_slab_weights) == 0.:
                profile_slab_weights[:] = 1
            
            for w in range(len(slab_data[:, 2][slab_indices])):
                if np.isnan(slab_data[:, 2][slab_indices][w]):
                    slab_data[:, 2][slab_indices][w] = 0
            combine_profile[j, 2] = np.average(slab_data[:, 2][slab_indices], weights=profile_slab_weights)
            
            if (np.min(profile_slab_distances) > max_point_distance):
                combine_profile[j, 0] = 0.
                combine_profile[j, 1] = 0.

    # Remove profile points past slab edge
    nonzero             = np.where(combine_profile[:, 0] != 0.)
    profile_within_slab = np.copy(combine_profile[nonzero])

    # Save profiles 
    np.savetxt(output_dir + '/geographic/profile_' + str(i).zfill(5) + '.txt', profile_within_slab, delimiter=' ', fmt='%13.8f')
    
    # For pygmt.project to transform from geographic to cartesian coordinates, it seems like we need to use the method of a center point with an azimuth rotation, and set the flag
    # flat_earth=True to tell the function that we are transforming into a Cartesian plane.
    
    
    pygmt.project(data=output_dir + '/geographic/profile_' + str(i).zfill(5) + '.txt', outfile=output_dir + '/cartesian/profile_' + str(i).zfill(5) + '.txt', \
                  azimuth=rotation_angle, center=np.array([np.average(region_lon), np.average(region_lat)]), flat_earth=False, convention='pqz', unit=True)
    return profile_within_slab


def hikurangi_plateau_depth(plateau_edge_file, profile_directory, rotation_angle, region_lon, region_lat, output_dir, create_plots):
    
    '''
    Takes a file that has a map view projection of the leading edge of the subducted Hikurangi Plateau, and compares the geographic
    coordinates of the profiles to the leading subducted edge. If the profiles are trenchward of the leading edge, then the profile
    is modified to account for thickened oceanic crust, else the profile is unchanged
    
    plateau_edge_file = the pathway to the file containing the leading subducted edge of the Hikurangi Plateau
    profile_directory = the pathway to the directory containing all the profiles created by profile_generator
    rotation_angle    = angle to rotate the profiles during the projection from Spherical to Cartesian
    region_lon        = the longitudinal extents of the region
    region_lat        = the latitudinal extents of the region
    output_dir        = the path to write the modified profiles
    create_plots      = True or False, will plot the modified tracks to confirm this function is working correctly
    '''
    
    index = 0
    
    sys('rm -rf ' + output_dir)
    sys('mkdir ' + output_dir)
    sys('mkdir ' + output_dir + '/geographic')
    sys('mkdir ' + output_dir + '/cartesian')
    for file in np.sort(os.listdir(profile_directory)):
        index += 1
        file_path = os.path.join(profile_directory, file)
        plateau_edge = np.loadtxt(plateau_edge_file)

        sys('gmt select ' + file_path + ' -F' + plateau_edge_file + ' -If -Ef > ' + output_dir + '/plateau_location_along_profile.xyd')
        intersected = np.loadtxt(fname=output_dir + '/plateau_location_along_profile.xyd')
        np.savetxt(fname=output_dir + '/geographic/profile_' + str(index).zfill(5) + '.xyd', X=intersected)

        if create_plots == True:
            import matplotlib.pyplot as plt
            track = np.loadtxt(fname=file_path)
            plt.scatter(track[:, 0], track[:, 1], c = 'k', s = 20)
            plt.scatter(intersected[:, 0], intersected[:, 1], s = 5, c = 'r')

        pygmt.project(data=output_dir + '/geographic/profile_' + str(index).zfill(5) + '.xyd', outfile=output_dir + '/cartesian/profile_' + str(index).zfill(5) + '.xyz', \
                      azimuth=rotation_angle, center=np.array([np.average(region_lon), np.average(region_lat)]), convention='pqz', unit=True)
    return None





# This function keeps profiles from intersecting one another, no longer necessary
# ################## INPUTS ##################
# # truncate_profiles: whether or not you want to truncate profiles to prevent them from intersecting one another
# # profile_within_slab: profile generated by other functions
# # polygon: array containing the polygon used to keep profiles from intersecting one another
# # region_lat: the latitudinal extent of the slab you want
# # region_lon: the longitudinal extent of the slab you want
# # i: profile index
# # trench_orientation: cardinal direction of the trench relative to the OVERRIDING plate (i.e. Cascadia is WEST)

# def polygon_truncation(truncate_profiles, profile_within_slab, left_poly, right_poly, region_lat, region_lon, i, trench_orientation):
    
#     # Two cases, East/West Subducting Slabs and North/South Subducting Slabs
#     # This is for East/West Subducting Slabs
#     if truncate_profiles == True:
#         if trench_orientation == 'East' or trench_orientation == 'West':
#             # Now generate a polygon using the terminal ends of each profile, this polygon grows as the profiles move along the slab
#             if i == 0:
#                 # Special case, save empty file for first profile since no intersections are possible
#                 np.savetxt(fname='polygon.txt', X=left_poly)
#             else:
#                 # Find the terminal points of the previous profile
#                 top_left_ind = np.where( profile_within_slab[:, 0] == np.min(profile_within_slab[:, 0]) )
#                 top_right_ind = np.where( profile_within_slab[:, 0] == np.max(profile_within_slab[:, 0]) )

#                 right_poly.append( [np.max(profile_within_slab[:, 0]), profile_within_slab[:, 1][top_right_ind][0]] )
#                 left_poly.append( [np.min(profile_within_slab[:, 0]), profile_within_slab[:, 1][top_left_ind][0]] )

#                 # Sort the points into the proper format for GMT to recognize the polygons shape
#                 sorted_polygon = []
#                 for j in range(len(right_poly)):
#                     sorted_polygon.append(right_poly[j])

#                 for arr in reversed(left_poly):
#                     sorted_polygon.append(arr)
#                 np.savetxt(fname='polygon.txt', X=sorted_polygon)

#         # This is for North/South Subducting Slabs
#         elif trench_orientation == 'North' or trench_orientation == 'South':
#             if i == 0:
#                 np.savetxt(fname='polygon.txt', X=left_poly)
#             else:
#                 top_right_ind = np.where( profile_within_slab[:, 1] == np.max(profile_within_slab[:, 1]) )
#                 bot_right_ind = np.where( profile_within_slab[:, 1] == np.min(profile_within_slab[:, 1]) )

#                 right_poly.append( [profile_within_slab[:, 0][bot_right_ind], np.min(profile_within_slab[:, 1])] )
#                 left_poly.append( [profile_within_slab[:, 0][top_right_ind], np.max(profile_within_slab[:, 1])] )
                
#                 sorted_polygon = []
#                 for j in range(len(right_poly)):
#                     sorted_polygon.append(right_poly[j])

#                 for arr in reversed(left_poly):
#                     sorted_polygon.append(arr)
#                 np.savetxt(fname='polygon.txt', X=sorted_polygon)

#         # Use GMT to remove any points in the great circle that intersect previous profiles
#         if i > 0:
#             sys('gmt select positive_profile.xyd -Fpolygon.txt -If -Ef > positive_truncated.xyd')
#             sys('gmt select negative_profile.xyd -Fpolygon.txt -If -Ef > negative_truncated.xyd')
#             pos_profile = np.loadtxt('positive_truncated.xyd')
#             neg_profile = np.loadtxt('negative_truncated.xyd')
#         else:
#             pos_profile = np.loadtxt('positive_profile.xyd')
#             neg_profile = np.loadtxt('negative_profile.xyd')
#     else:
#         pos_profile = np.loadtxt('positive_profile.xyd')
#         neg_profile = np.loadtxt('negative_profile.xyd')
#     if trench_orientation == 'East' or trench_orientation == 'North':
#         combine_profile = np.concatenate( (np.flipud(pos_profile), np.delete(neg_profile, 0, axis=0)) ) # delete single duplicating point
#     elif trench_orientation == 'West' or trench_orientation == 'South':
#         combine_profile = np.concatenate( (np.flipud(neg_profile), np.delete(pos_profile, 0, axis=0)) ) # delete single duplicating point
        
#     for k in range(len(combine_profile)):
#         if combine_profile[:, 0][k] <= 0:
#             combine_profile[:, 0][k] = combine_profile[:, 0][k] + 360
#     return combine_profile, left_poly, right_poly