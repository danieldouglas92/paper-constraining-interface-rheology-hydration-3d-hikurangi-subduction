# This script uses Paraview's pvpython to take output from a global mantle
# convection model and extract velocities along user-specified boundaries. 
# Next, the script uses the extracted velocities to create ASCII files which
# can be applied as boundary conditions in regional models to account for
# far-field effects outside of the regional model boundaries. This is achieved
# by specifying the bounds of the regional model, slicing through the 
# global model along planes which coincide with the location of the
# regional model boundaries, and saving the velocities at each point
# into an ASCII file. The regional model is expected to be a 3D spherical
# chunk, and the global model is expected to be a 3D spherical shell.

# The user must specify the following parameters in the script:
# 1. The location of the .pvd file for a global convection model
# 2. The output directory where the ASCII files will be saved
# 3. The refinement level of the global model 
# 4. The radial resolution of the regional model
# 5. The lateral resolution of the regional model
# 6. The maximum and minimum radius, latitude, and  longitude for the regional model.

# This shows an example for applying this script to the global model
# presented in the S2ORTS cookbook, and applying the extracted velocities 
# to a regional model that spans from 20 degrees latitude to 50 degrees
# latitude, 190 degrees longitude to 230 degrees longitude, and from a
# radius of 5070 km to a radius of 6370 km. Before running this script, 
# make sure that you run the S2ORTS cookbook and generate the solution.pvd file.

# This script defines 3 functions
# 1. spherical_to_cartesian, which converts from spherical coordinates to Cartesian
#    coordinates
# 2. cartesian_to_spherical, which converts from Cartesian coordinates to spherical
#    coordinates
# 3. slice_plane_calculator, which determines the normal to the plane that defines
#    the east and west model boundaries in the regional chunk.

# Import packages
import numpy as np
import pandas as pd
import os
from paraview import simple
from paraview import servermanager

def spherical_to_cartesian(radii, latitudes, longitudes, for_slice):
    """
    Converts from spherical coordinates to Cartesian coordinates.
    radii: the radius, in m
    latitudes: the latitude, in degrees ranging from -90 to 90
    longitudes: the longitude, in degrees ranging from 0 to 360
    for_slice: boolean, if false the provided points are directly 
    converted into Cartesian coordinates. If true, radii, latitudes,
    and longitudes must be arrays, and the function returns a uniform 
    structured grid defining the slice.
    Returns x, y, z, in m
    """

    if for_slice:
        cartesian_coordinates = []
        for lat in latitudes:
            for lon in longitudes:
                for r in radii:
                    x = r * np.sin(np.deg2rad(90 - lat)) * np.cos(np.deg2rad(lon))
                    y = r * np.sin(np.deg2rad(90 - lat)) * np.sin(np.deg2rad(lon))
                    z = r * np.cos(np.deg2rad(90 - lat))
    
                    if lon == 0:
                        y = 0
                    elif lon == np.pi/2:
                        x = 0
                    if lat == np.pi/2:
                        z = 0          
                    cartesian_coordinates.append([x, y, z])
    
        return np.array(cartesian_coordinates)
        
    else:
        x = radii * np.sin(np.deg2rad(90 - latitudes)) * np.cos(np.deg2rad(longitudes))
        y = radii * np.sin(np.deg2rad(90 - latitudes)) * np.sin(np.deg2rad(longitudes))
        z = radii * np.cos(np.deg2rad(90 - latitudes))
    
        return np.array([x, y, z])
    
def cartesian_to_spherical(x, y, z):
    """
    Takes an x, y, z point and converts it to spherical coordinates. 
    Returns r in m, latitude in degrees, and longitude, ranging from 0 to 360, in degrees 
    """
    r = np.sqrt(x**2 + y**2 + z**2)
    latitude = 90 - np.rad2deg( np.arccos( z / (np.sqrt(x**2 + y**2 + z**2)) ) )
    longitude =  np.sign(y) * np.rad2deg(np.arccos( x / np.sqrt(x**2 + y**2) ))
    longitude[np.where(longitude < 0)] = longitude[np.where(longitude < 0)] + 360
    longitude[np.where(longitude == 0)] = 180
    
    return r, latitude, longitude

def slice_plane_calculator(boundary_name, radius_bounds, latitude_bounds, longitude_bounds):
    """
    Calculates the normal of a plane which is used for slicing the global models. This is
    achieved by defining three points on either the east or west model boundary using the 
    values provided by radius_bounds, latitude_bounds, and longitude_bounds.
    boundary_name: the name of the model boundary
    radius_bounds: the maximum and minimum radius of the regional models
    latitude_bounds: the maximum and minimum latitude of the regional models
    longitude_bounds: the maximum and minimum longitude of the regional models
    """
    # Define the 3 points on the west or east boundary. If west, we are on the 
    # minimum longitude, and if east we are on the maximum longitude.
    if boundary_name == "west":
        spherical_point_1 = np.array([np.max(radius_bounds), \
                                      np.max(latitude_bounds), \
                                      np.min(longitude_bounds)])
        spherical_point_2 = np.array([np.max(radius_bounds), \
                                      np.min(latitude_bounds), \
                                      np.min(longitude_bounds)])
        spherical_point_3 = np.array([np.min(radius_bounds), \
                                      np.max(latitude_bounds), \
                                      np.min(longitude_bounds)])

    elif boundary_name == "east":
        spherical_point_1 = np.array([np.max(radius_bounds), \
                                     np.max(latitude_bounds), \
                                     np.max(longitude_bounds)])
        spherical_point_2 = np.array([np.max(radius_bounds), \
                                     np.min(latitude_bounds), \
                                     np.max(longitude_bounds)])
        spherical_point_3 = np.array([np.min(radius_bounds), \
                                     np.max(latitude_bounds), \
                                     np.max(longitude_bounds)])

    else:
        raise Exception("Unknown boundary name: " + boundary_name)
    # Convert spherical points to Cartesian
    cartesian_point_1 = spherical_to_cartesian(spherical_point_1[0], \
                                               spherical_point_1[1], \
                                               spherical_point_1[2], \
                                               for_slice=False)
    cartesian_point_2 = spherical_to_cartesian(spherical_point_2[0], \
                                               spherical_point_2[1], \
                                               spherical_point_2[2], \
                                               for_slice=False)
    cartesian_point_3 = spherical_to_cartesian(spherical_point_3[0], \
                                               spherical_point_3[1], \
                                               spherical_point_3[2], \
                                               for_slice=False)
    # Calculate 2 in-plane orthogonal vectors using the 3 Cartesian points
    vector_1_2 = cartesian_point_2 - cartesian_point_1
    vector_1_3 = cartesian_point_3 - cartesian_point_1
    # Taking the cross product yields a vector normal to the model boundary
    normal_vector_to_plane = np.cross(vector_1_2, vector_1_3)
    # Normalize
    unit_normal = normal_vector_to_plane / np.linalg.norm(normal_vector_to_plane)

    return unit_normal    

####################################################################################################################################

""" Usage: This script requires 8 input arguments:
input_directory: solution file for the global model (*.pvd)
output_directory: Where the .txt files for each boundary are saved
refinement_level: The number of mesh refinements in the global model
output_radius_resolution: The radial resolution of the regional slice (in meters)
output_lateral_resolution: The lateral resolution of the regional slice (in degrees)
radius_bounds: Array with the minimum and maximum radius (meters) of the regional model
latitude_bounds: Array with the minimum and maximum latitude (degrees) of the regional model
longitude_bounds: Array with the minimum and maximum longitude (degrees) of the regional model
"""

# Define the input arguments for the S2ORTS cookbook
input_directory = "/Users/danieldouglas/FINAL_SLABS/plate_model/HIGHRES_LITHO/"
radius_bounds = np.array([5870e3, 6370e3])
latitude_bounds = np.array([-30, -45])
longitude_bounds = np.array([170, 190])

model = simple.OpenDataFile(input_directory + "solution.pvd")
model.PointArrays = ["Points", "plastic_strain", "depth", "Crust", "Base_Subducting", "T", "p"]
model.UpdatePipeline()


# Create a calculator filter to determine the latitude in the global models
latitude_calc = simple.PythonCalculator(registrationName="latitude_calc", Input=model)
latitude_calc.ArrayAssociation = "Point Data" # "Point Data" not "Cell Data"
latitude_calc.CopyArrays = True # Copy all of the other variables into the calculator filter
latitude_calc.Expression = "90 - np.arccos( points[:, 2] / (sqrt(points[:, 0]**2 + points[:, 1]**2 + points[:, 2]**2)) ) * 180 / np.pi" # Calculate the latitude
latitude_calc.ArrayName = "latitude"
latitude_calc.UpdatePipeline()

# Create threshold filter
latitude_threshold = simple.Threshold(Input=latitude_calc)
latitude_threshold.Scalars = ("POINTS", "latitude") # Threshold the latitude variable
latitude_threshold.ThresholdMethod = "Between" 

# Threshold on either side of the maximum lat_bound (north) or the minimum lat_bound (south)
# based on the refinement_level of the global models.
latitude_threshold.LowerThreshold = np.min(latitude_bounds)
latitude_threshold.UpperThreshold = np.max(latitude_bounds)

latitude_threshold.UpdatePipeline()

# Create a calculator filter to determine the latitude in the global models
radius_calc = simple.PythonCalculator(registrationName="radius_calc", Input=latitude_threshold)
radius_calc.ArrayAssociation = "Point Data" # "Point Data" not "Cell Data"
radius_calc.CopyArrays = True # Copy all of the other variables into the calculator filter
radius_calc.Expression = "sqrt(points[:, 0]**2 + points[:, 1]**2 + points[:, 2]**2)" # Calculate the latitude
radius_calc.ArrayName = "radius"
radius_calc.UpdatePipeline()

radius_threshold = simple.Threshold(Input=radius_calc)
radius_threshold.Scalars = ("POINTS", "radius") # Threshold the latitude variable
radius_threshold.ThresholdMethod = "Between" 

radius_threshold.LowerThreshold = np.min(radius_bounds)
radius_threshold.UpperThreshold = np.max(radius_bounds)

radius_threshold.UpdatePipeline()


# Create a calculator filter to determine the latitude in the global models
wrong_longitude_calc = simple.PythonCalculator(registrationName="wrong_longitude_calc", Input=radius_threshold)
wrong_longitude_calc.ArrayAssociation = "Point Data" # "Point Data" not "Cell Data"
wrong_longitude_calc.CopyArrays = True # Copy all of the other variables into the calculator filter
wrong_longitude_calc.Expression = "points[:, 1] / abs(points[:, 1]) * np.arccos( points[:, 0] / sqrt(points[:, 0]**2 + points[:, 1]**2 + 1e-10) ) * 180 / np.pi" # Calculate the latitude
wrong_longitude_calc.ArrayName = "wrong_longitude"
wrong_longitude_calc.UpdatePipeline()

longitude_calc = simple.Calculator(registrationName="longitude_calc", Input=wrong_longitude_calc)
# longitude_calc.AttributeType = "Point Data" # "Point Data" not "Cell Data"
# longitude_calc.CopyArrays = True # Copy all of the other variables into the calculator filter
longitude_calc.Function = "if(wrong_longitude <= 0, wrong_longitude + 360, wrong_longitude)" # Calculate the latitude
longitude_calc.ResultArrayName = "longitude"
longitude_calc.UpdatePipeline()

longitude_threshold = simple.Threshold(Input=longitude_calc)
longitude_threshold.Scalars = ("POINTS", "longitude") # Threshold the latitude variable
longitude_threshold.ThresholdMethod = "Between" 

longitude_threshold.LowerThreshold = np.min(longitude_bounds)
longitude_threshold.UpperThreshold = np.max(longitude_bounds)

longitude_threshold.UpdatePipeline()

writer = simple.DataSetWriter(Input=longitude_threshold, FileName=input_directory + "subset.vtk")
writer.UpdatePipeline()
