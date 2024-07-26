############################################# MODEL FEATURES FOR THE WB FILE #############################################
################### Functions which write strings for various model features to the world builder file ###################
##########################################################################################################################

import numpy as np
import os

def model_feature_string(model_name, feature_name, min_depth, max_depth, coordinates, dip_point, is_subducting=False):
    
    '''
    Creates a string which is then written to a world builder file that initializes a world builder feature

    model_name    = string specifying the WB model type you want for the feature 
    feature_name  = user specified string naming the feature (can be anything, for ease of reading the WB file)
    min_depth     = minimum depth extent of the feature
    max_depth     = maximum depth extenet of the feautre
    coordinates   = array specifying the points which create a bounding volume of the feature
    is_subducting = True or False, whether the feature is a subducting plate
    dip_point     = if is_subducting=True, location of the dip point of the subducting slab    
    '''
    
    if is_subducting == False:
        feature_string = '"model":"' + str(model_name) + '", "name":"' + str(feature_name) + '", "min depth":' + str(min_depth) +', "max depth":' + \
                         str(max_depth) + ', "coordinates":' + str(coordinates) + ',\n'
    else:
        feature_string = '"model":"' + str(model_name) + '", "name":"' + str(feature_name) + '", "min depth":' + str(min_depth) + ', "max depth":' + \
                         str(max_depth) + ', "coordinates":' + str(coordinates) + ', "dip point":' + str(dip_point) + ',\n'
    return feature_string


def trench_extractor(profile_directory, xshift, yshift, trench_orientation, coordinate_sys):

    '''
    Extracts the location of the trench from the files containing slab profiles

    profile_directory  = the pathway to the directory containing the profiles across slab2.0 data
    xshift             = Shifts the x coordinate of the trench
    yshift             = Shifts the y coordinate of the trench
    trench_orientation  = direction of the trench relative to the overriding plate
    coordinate_sys     = string specifying whether the coordinate system of the profiles is 'Cartesian' or 'Spherical'
    '''
    
    trench_x = []
    trench_y = []
    if coordinate_sys == 'Cartesian':
        for file in np.sort(os.listdir(profile_directory)):
            profile_file = np.loadtxt(fname=profile_directory + file)
            if trench_orientation == 'East':
                trench_x.append( (np.max(profile_file[:, 0]) + xshift) * 1000)
                trench_y.append( (profile_file[:, 1][np.where(profile_file[:, 0] == np.max(profile_file[:, 0]))][0] + yshift) * 1000)
                
            elif trench_orientation == 'West':
                trench_x.append( (np.min(profile_file[:, 0]) + xshift) * 1000)
                trench_y.append( (profile_file[:, 1][np.where(profile_file[:, 0] == np.min(profile_file[:, 0]))][0] + yshift) * 1000)

            elif trench_orientation == 'South':
                trench_y.append( (np.min(profile_file[:, 1]) + yshift) * 1000)
                trench_x.append( (profile_file[:, 0][np.where(profile_file[:, 1] == np.min(profile_file[:, 1]))][0] + xshift) * 1000)

            elif trench_orientation == 'North':
                trench_y.append( np.max(profile_file[:, 1] + yshift) * 1000)
                trench_x.append( (profile_file[:, 0][np.where(profile_file[:, 1] == np.min(profile_file[:, 1]))][0] + xshift) * 1000)
                
    elif coordinate_sys == 'Spherical':
        for file in np.sort(os.listdir(profile_directory)):
            profile_file = np.loadtxt(fname=profile_directory + file)
            if trench_orientation == 'East':
                trench_x.append( (np.max(profile_file[:, 0]) + xshift))
                trench_y.append( (profile_file[:, 1][np.where(profile_file[:, 0] == np.max(profile_file[:, 0]))][0] + yshift))
                
            elif trench_orientation == 'West':
                trench_x.append( (np.min(profile_file[:, 0]) + xshift))
                trench_y.append( (profile_file[:, 1][np.where(profile_file[:, 0] == np.min(profile_file[:, 0]))][0] + yshift))

            elif trench_orientation == 'South':
                trench_y.append( (np.min(profile_file[:, 1]) + yshift))
                trench_x.append( (profile_file[:, 0][np.where(profile_file[:, 1] == np.min(profile_file[:, 1]))][0] + xshift))

            elif trench_orientation == 'North':
                trench_y.append( np.max(profile_file[:, 1] + yshift))
                trench_x.append( (profile_file[:, 0][np.where(profile_file[:, 1] == np.min(profile_file[:, 1]))][0] + xshift)) 
        
    return trench_x, trench_y

def trench_splitter(trench_x, trench_y, x_bounds, y_bounds, trench_orientation):

    '''
    Takes the boundaries of the model and creates an array of points which splits the surface into two areas:
    An overriding plate bound by the landward side of the trench, and a subducting plate bound by the seaward
    side of the trench. If the trench does not intersect the boundaries of the model, vertical lines are drawn
    from the terminii of the trench to the boundaries of the model which close the regions describing the
    overriding and subducting plates.

    trench_x          = array containing x-points of the trench output by trench_extractor function
    trench_y          = array containing y-points of the trench output by trench_extractor function
    x_bounds          = array with the max and min x-values of the model domain
    y_bounds          = array with the max and min y-values of the model domain
    trench_orientation = string specifying the direction of the trench relative to the overriding plate 
                       (required for extending the 'trench' to the boundaries of the model. For example,
                       an East-West Striking trench will be extended East-West to the Model boundaries)
    '''
    
    if trench_orientation == 'East':
        sub_coords = [[np.max(x_bounds), np.min(y_bounds)]]
        over_coords = [[np.min(x_bounds), np.min(y_bounds)]]
        sub_coords.append([trench_x[0], np.min(y_bounds)])
        over_coords.append([trench_x[0], np.min(y_bounds)])
        for i in range(len(trench_x)):
            sub_coords.append([trench_x[i], trench_y[i]])
            over_coords.append([trench_x[i], trench_y[i]])
        sub_coords.append([trench_x[-1], np.max(y_bounds)])
        over_coords.append([trench_x[-1], np.max(y_bounds)])
        sub_coords.append([np.max(x_bounds), np.max(y_bounds)])
        over_coords.append([np.min(x_bounds), np.max(y_bounds)])
        
    elif trench_orientation == 'West':
        sub_coords = [[np.min(x_bounds), np.min(y_bounds)]]
        over_coords = [[np.max(x_bounds), np.min(y_bounds)]]
        sub_coords.append([trench_x[0], np.min(y_bounds)])
        over_coords.append([trench_x[0], np.min(y_bounds)])
        for i in range(len(trench_x)):
            sub_coords.append([trench_x[i], trench_y[i]])
            over_coords.append([trench_x[i], trench_y[i]])
        sub_coords.append([trench_x[-1], np.max(y_bounds)])
        over_coords.append([trench_x[-1], np.max(y_bounds)])
        sub_coords.append([np.min(x_bounds), np.max(y_bounds)])
        over_coords.append([np.max(x_bounds), np.max(y_bounds)])
    
    elif trench_orientation == 'North':
        sub_coords = [[np.min(x_bounds), np.min(y_bounds)]]
        over_coords = [[np.min(x_bounds), np.max(y_bounds)]]
        sub_coords.append([np.min(x_bounds), trench_y[0]])
        over_coords.append([np.min(x_bounds), trench_y[0]])
        for i in range(len(trench_y)):
            sub_coords.append([trench_x[i], trench_y[i]])
            over_coords.append([trench_x[i], trench_y[i]])
        sub_coords.append([np.max(x_bounds), trench_y[-1]])
        over_coords.append([np.max(x_bounds), trench_y[-1]])
        sub_coords.append([np.max(x_bounds), np.min(y_bounds)])
        over_coords.append([np.max(x_bounds), np.max(y_bounds)])
        
    elif trench_orientation == 'South':
        over_coords = [[np.min(x_bounds), np.max(y_bounds)]]
        sub_coords = [[np.min(x_bounds), np.min(y_bounds)]]
        over_coords.append([np.min(x_bounds), trench_y[0]])
        sub_coords.append([np.min(x_bounds), trench_y[0]])
        for i in range(len(trench_y)):
            over_coords.append([trench_x[i], trench_y[i]])
            sub_coords.append([trench_x[i], trench_y[i]])
        over_coords.append([np.max(x_bounds), trench_y[-1]])
        sub_coords.append([np.max(x_bounds), trench_y[-1]])
        over_coords.append([np.max(x_bounds), np.max(y_bounds)])
        sub_coords.append([np.max(x_bounds), np.min(y_bounds)])
    return sub_coords, over_coords

########################################## TEMPERATURE FEATURES FOR THE WB FILE ##########################################
################ Functions which write strings for various temperature features to the world builder file ################
##########################################################################################################################

def cooling_model(model_name, max_depth, min_depth, bottom_temp, top_temp, spr_vel, ridge_coords, first_or_last='both'):

    '''
    Defines the geotherm in an oceanic plate as the half-space cooling or plate cooling model

    model_name    = 'plate model' or 'half space model', specifiying which cooling model you want
    max_depth     = the maximum depth to which this temperature feature will be used
    min_depth     = the minimum depth to which this temperature feature will be used
    bottom_temp   = the temperature the maximum depth
    top_temp      = the temperature at the minimum depth
    spreading_vel = the velocity of the plate
    ridge_coords  = the coordinates specifying the axis of the spreading ridge
    first_or_last = whether this feature is the first or last temperature feature. 'first' is the first feature, 'last' is the last feature, 'both'
                    is if this is the only feature (i.e. the 'first' and the 'last')
    '''

    if first_or_last == 'first' or first_or_last == 'both':
        string = '"temperature models":['
    else:
        string = ''
        
    if model_name == 'plate model' or model_name == 'half space model':
        string += '{"model":"' + str(model_name) + '", "max depth":' + str(max_depth) + ', "min depth":' + str(min_depth) + ', "top temperature":' + \
                  str(top_temp) + ', "bottom temperature":' + str(bottom_temp) + ', "spreading velocity":' + str(spr_vel) + ', "ridge coordinates":' + \
                  str(ridge_coords) + '}'
    if first_or_last == 'last' or first_or_last == 'both':
        string += '], \n'
    else:
        string += ',\n'
    return string



def linear_model(model_name, max_depth, min_depth, bottom_temp, top_temp, first_or_last='both'):
    
    '''
    Defines the temperature in a continental plate as a linear conductive geotherm

    model_name    = 'linear'
    max_depth     = the maximum depth to which this temperature feature will be used
    min_depth     = the minimum depth to which this temperature feature will be used
    bottom_temp   = the temperature the maximum depth
    top_temp      = the temperature at the minimum depth
    first_or_last = whether this feature is the first or last temperature feature. 'first' is the first feature, 'last' is the last feature, 'both'
                    is if this is the only feature (i.e. the 'first' and the 'last')
    '''
    
    if first_or_last == 'first' or first_or_last == 'both':
        string = '"temperature models":['
    string += '{"model":"' + str(model_name) + '", "max depth":' + str(max_depth) + ', "min depth":' + str(min_depth) + ', "top temperature":' + \
               str(top_temp) + ', "bottom temperature":' + str(bottom_temp) + '}'
    if first_or_last == 'last' or first_or_last == 'both':
        string += '], \n'
    else:
        string += ',\n'
    return string



def uniform_model(model_name, uniform_temp, operation, first_or_last='both'):
    
    '''
    Defines the temperature within the feature as a constant value

    model_name   = 'uniform'
    uniform_temp = the constant temperature value
    operation    = whether uniform_temp overrides the existing temperature field, or adds or subtracts from it. 'replace', 'add', or 'subtract'
    first_or_last = whether this feature is the first or last temperature feature. 'first' is the first feature, 'last' is the last feature, 'both'
                    is if this is the only feature (i.e. the 'first' and the 'last')
    '''
    
    if first_or_last == 'first' or first_or_last == 'both':
            string = '"temperature models":['
    else:
        string = ''
            
    if model_name == 'uniform':
        string += '{"model":"' + str(model_name) + '", "temperature":' + str(uniform_temp) + ', "operation":"' + str(operation) + '"}'   
        
    if first_or_last == 'last' or first_or_last == 'both':
        string += '], \n'
        
    else:
        string += ',\n'
        
    return string
   
   
    
def mass_conserving_model(thermal_model, density, plate_vel, coupling_depth, ridge_coords, taper, max_slab_top, min_slab_top, first_or_last='both'):
    
    '''
    Defines the temperature in a subducting slab based on a half space cooling model (bottom portion of the slab) and a 
    semi-infinite heating model (top portion of the slab).
    
    thermal_model  = either "halfspace model", or "plate model"
    density        = density of the slab
    plate_vel      = velocity of the plate entering the trench
    coupling_depth = the depth where the subducting slab interfaces the overriding plate
    ridge_coords   = the coordinates specifying the axis of the spreading ridge
    max_slab_top   = the distance to the top of the slab from the mid plane of the slab
    min_slab_top   = the distance to the bottom of the slab from the mid plane of the slab
    first_or_last  = whether this feature is the first or last temperature feature. 'first' is the first feature, 'last' is the last feature, 'both'
                     is if this is the only feature (i.e. the 'first' and the 'last')
    '''
    
    if first_or_last == 'first' or first_or_last == 'both':
        string = '"temperature models":['
    else:
        string = ''
        
    string += '{"model":"mass conserving", "reference model name": "' + thermal_model + '", "adiabatic heating":true, "density":' + str(density) + ', "plate velocity":' + str(plate_vel) + ', "coupling depth":' + str(coupling_depth) + \
              ', "ridge coordinates":' + str(ridge_coords) + ', "taper distance":' + str(taper) + \
              ', "max distance slab top":' + str(max_slab_top) + ', "min distance slab top":' + str(min_slab_top) + '}'
    
    if first_or_last == 'last' or first_or_last == 'both':
        string += '], \n'
    else:
        string += ', \n'
        
    return string


########################################## COMPOSITION FEATURES FOR THE WB FILE ##########################################
############################# Functions which write strings segments for composition fields ##############################
##########################################################################################################################



def composition_feature_string(model_name, comp_index, max_depth, min_depth, is_subducting=False, first_or_last='both', operation='replace'):
    
    '''
    Writes the string for a composition feature to the world builder file
    
    model_name    = string specifying the composition model
    comp_index    = the index of the current composition field. Indexing starts at 0
    max_depth     = maximum depth extent of the compositional field
    min_depth     = minimum depth extent of the compositional field
    is_subducting = True or False, whether the current field is a subducting slab
    first_or_last = whether this feature is the first or last temperature feature. 'first' is the first feature, 'last' is the last feature, 'both'
                    is if this is the only feature (i.e. the 'first' and the 'last')
    '''
    if first_or_last == 'first' or first_or_last == 'both':
        string = '"composition models":['
    else:
        string = ''
    
    if is_subducting == False:
        string += '{"model":"' + str(model_name) + '", "compositions":[' + str(comp_index) + '], "max depth":' + str(max_depth) + \
                  ', "operation":"' + str(operation) + '", "min depth":' + str(min_depth) + '}'
    else:
        string += '{"model":"' + str(model_name) + '", "compositions":[' + str(comp_index) + '], "operation":"' + str(operation) + '", ' + \
                  '"max distance slab top":' + str(max_depth) + ', "min distance slab top":' + str(min_depth) + '}'
    if first_or_last == 'last' or first_or_last == 'both':
        string += ']\n'
    else:
        string += ',\n'
    return string

######################################## SUBDUCTION ZONE FEATURES FOR THE WB FILE ########################################
############################# Functions which write strings segments of the subducting slab ##############################
##########################################################################################################################



def segment_string_at_end(length, thickness, angle, top_truncation, coordinate, total_sections, segment_num, current_segment):

    '''
    Writes the string that defines segments of a slab to the world builder file

    When initializing a slab in Worldbuilder, you must specify the same number of segments for each section along strike
    of the slab, which is troublesome since some sections of slab are much longer than others. To get around this,
    the slab is initialized with the same number of segments as the section with the maximum number of segments. 
    Then, sections with less than this number of segments are assigned filler sections with lengths of 0 m until
    they reach the required number of segments

    length          = the length of a given segment
    thickness       = the thickness of a given segment
    angle           = the dip of a given segment
    coordinate      = the trench coordinate index coupled to the current section
    total_sections  = the number of sections making up the slab
    segment_num     = the number of segments making up each section
    current_segment = the current segment index of a section
    '''
 

    string_total = ''
    
    if current_segment <= len(thickness) - 1:
        if current_segment != (segment_num - 1):
            string_total = '{"length":' + str(length[current_segment]) + ', "thickness":' + str(thickness[current_segment]) + \
                            ', "angle":' + str(angle[current_segment]) + ', "top truncation":' + str(top_truncation) + '},\n'
            return string_total

        elif current_segment == (segment_num - 1) and coordinate == total_sections:
            string_total += '{"length":' + str(length[current_segment]) + ', "thickness":' + str(thickness[current_segment]) + \
                             ', "angle":' + str(angle[current_segment]) + ', "top truncation":' + str(top_truncation) + '}]}\n'
            return string_total

        else:
            string_total += '{"length":' + str(length[current_segment]) + ', "thickness":' + str(thickness[current_segment]) + \
                             ', "angle":' + str(angle[current_segment]) + ', "top truncation":' + str(top_truncation) + '}]},\n'
            return string_total
        
    if current_segment > len(thickness) - 1:
        if len(thickness[0]) > 1:
            if current_segment != (segment_num - 1):
                string_total += '{"length":' + str(0.0) + ', "thickness":' + str([np.max(thickness), np.max(thickness)]) + ', "angle":' + str(angle[len(thickness) - 1]) + '},\n'# \
                return string_total
            
            elif current_segment == (segment_num - 1) and coordinate == total_sections:
                string_total += '{"length":' + str(0.0) + ', "thickness":' + str([np.max(thickness), np.max(thickness)]) + ', "angle":' + str(angle[len(thickness) - 1]) + '}]}\n'
                return string_total
                
            else:
                string_total += '{"length":' + str(0.0) + ', "thickness":' + str([np.max(thickness), np.max(thickness)]) + ', "angle":' + str(angle[len(thickness) - 1]) + '}]},\n'
                return string_total
        
        else:
            if current_segment != (segment_num - 1):
                string_total += '{"length":' + str(0.0) + ', "thickness":' + str([np.max(thickness)]) + ', "angle":' + str(angle[len(thickness) - 1]) + '},\n'
                return string_total
            
            elif current_segment == (segment_num - 1) and coordinate == total_sections:
                string_total += '{"length":' + str(0.0) + ', "thickness":' + str([np.max(thickness)]) + ', "angle":' + str(angle[len(thickness) - 1]) + '}]}'
                return string_total
            
            else:
                string_total += '{"length":' + str(0.0) + ', "thickness":' + str([np.max(thickness)]) + ', "angle":' + str(angle[len(thickness) - 1]) + '}]},\n'
                return string_total
            
def segment_string(length, thickness, angle, top_truncation, coordinate, total_sections, segment_num, current_segment):

    '''
    Writes the string that defines segments of a slab to the world builder file

    When initializing a slab in Worldbuilder, you must specify the same number of segments for each section along strike
    of the slab, which is troublesome since some sections of slab are much longer than others. To get around this,
    the slab is initialized with the same number of segments as the section with the maximum number of segments. 
    Then, sections with less than this number of segments are assigned filler sections with lengths of 0 m until
    they reach the required number of segments

    length          = the length of a given segment
    thickness       = the thickness of a given segment
    angle           = the dip of a given segment
    coordinate      = the trench coordinate index coupled to the current section
    total_sections  = the number of sections making up the slab
    segment_num     = the number of segments making up each section
    current_segment = the current segment index of a section
    '''
 
    string_total = ''
    filler_segment_number = abs(len(thickness) - segment_num)
    if current_segment < filler_segment_number:
        if len(thickness[0]) > 1:
            # Here we proceed as above, adding the segment strings until we reach the end of the array
            string_total += '{"length":' + str(0.0) + ', "thickness":' + str([np.min(thickness), np.min(thickness)]) + ', "angle":' + str([0.0, 0.0]) + '},\n'# \
                           # ', "angle":' + str([0,0]) + '},\n'
            return string_total
        
        else:
            # Here we proceed as above, adding the segment strings until we reach the end of the array
            string_total += '{"length":' + str(0.0) + ', "thickness":' + str([np.min(thickness)]) + ', "angle":' + str([0.0, 0.0]) + '},\n'# \
                           # ', "angle":' + str([0,0]) + '},\n'
            return string_total
        
        # Now we have reached the end of the array, but need to fill out the amount of segments to reach the maximum
        # segment number. Add segments with lengths of 1m, thicknesses of 1m, and dips of 1 degree until this is the case.
    else:
        if current_segment != segment_num - 1:
            string_total = '{"length":' + str(length[current_segment - filler_segment_number]) + ', "thickness":' + str(thickness[current_segment - filler_segment_number]) + \
                            ', "angle":' + str(angle[current_segment - filler_segment_number]) + ', "top truncation":' + str(top_truncation) + '},\n'
            return string_total

        elif current_segment == (segment_num - 1) and coordinate == total_sections:
            string_total += '{"length":' + str(length[current_segment - filler_segment_number]) + ', "thickness":' + str(thickness[current_segment - filler_segment_number]) + \
                        ', "angle":' + str(angle[current_segment - filler_segment_number]) + ', "top truncation":' + str(top_truncation) + '}]}\n'
            return string_total

        else:
            string_total += '{"length":' + str(length[current_segment - filler_segment_number]) + ', "thickness":' + str(thickness[current_segment - filler_segment_number]) + \
                        ', "angle":' + str(angle[current_segment - filler_segment_number]) + ', "top truncation":' + str(top_truncation) + '}]},\n'
            return string_total

def segment_section(world_builder_file, profile_directory, xshift, yshift, slab_thickness, top_truncation, coordinate_system, fill_segments):
    
    '''
    Calculates dip, thickness, and length of segments then uses segment_string() to create the string for the world builder file.
    This function has the option to vary the slab thickness between certain depths. If the varialbe slab_thickness is a scalar,
    a constant thickness is assumed for the entire section. To vary thickness, slab_thickness must be an array where each entry is an 
    array with length two. The entries of this array are the thickness of the slab and the depth to where 
    that thickness occurs. For example, to have the slab be 100km thick between 0 <= depth <= 200km, and then
    have the thickness increase to 150km for depth > 200km, the variable slab_thickness would need to be set as:

    slab_thickness = [ [100e3, 200e3], [150e3, 1e10] ]

    Where 1e10 was chosen to be a depth so high it would never be reached. Slab thickness and dip are varied gradually 
    down dip

    world_builder_file = the name of the world builder file
    profile_directory  = directory containing the world builder file
    xshift             = the amount of shift in the x direction, m
    yshift             = the amount of shift in the y direction, m
    slab_thickness     = the thickness of the slab, scaler for uniform thickness, array for variable, km
    '''
  
    # Import packages for spherical coordinate systems
    import geopy
    from geopy.distance import geodesic
    from geopy.distance import distance
    import pyproj
    geodesic = pyproj.Geod(ellps='WGS84')
    
    # Create an array which stores the length of each profile, the longest profile determines the number of
    # segments required for initializing the slab
    length = []
    for file in np.sort(os.listdir(profile_directory)):
        length.append(len(np.loadtxt(fname=profile_directory + file)))
    segment_num = max(length) - 1
    total_sections = len(os.listdir(profile_directory)) - 1
    
    # This loop initializes the slab with the correct number of segments
    world_builder_file.write('"segments":[\n')
    for i in range(segment_num):
        if i != segment_num - 1:
            world_builder_file.write('{"length":0, "thickness":[0.0], "angle":[0]},\n')
        else:
            world_builder_file.write('{"length":0, "thickness":[0.0], "angle":[0]}],\n\n')

    # This loops through all sections and determines the thickness, dip, and length of each segment
    world_builder_file.write('    "sections":[')
    for section_index, file in enumerate(np.sort(os.listdir(profile_directory))):
        filename = os.path.join(profile_directory, file)

        track_x = np.loadtxt(fname=filename, usecols=0) + xshift
        track_y = np.loadtxt(fname=filename, usecols=1) + yshift
        track_z = np.loadtxt(fname=filename, usecols=2)

        slab_length = 0
        segment_length = []
        dip_holder = [0]
        
        # Here we check to see if the thickness is set to vary along the slab by checking if slab_thickness is
        # a scalar or not. Dip is also computed and stored for output to the world builder file here
        if hasattr(slab_thickness, "__len__"):
            thick_holder = []
            for i in range(1, len(track_z)):
                thick_index = np.min(np.where( slab_thickness[section_index][:, 1] > np.abs(track_z[i]) ))
                thick_holder.append(slab_thickness[section_index][:, 0][thick_index])
                
                if coordinate_system == 'Cartesian':
                    slab_cart_proj = np.sqrt((track_x[i] - track_x[i - 1])**2 + (track_y[i] - track_y[i - 1])**2 )
                    slab_length = np.sqrt((track_x[i] - track_x[i - 1])**2 + (track_y[i] - track_y[i - 1])**2 + (track_z[i] - track_z[i - 1])**2)
                    riserun = (track_z[i - 1] - track_z[i]) / (slab_cart_proj)
                    
                elif coordinate_system == 'Spherical':
                    fwd_az, back_az, slab_cart_proj_m = geodesic.inv(track_x[i-1], track_y[i-1], track_x[i], track_y[i])
                    slab_cart_proj = slab_cart_proj_m / 1e3
                    slab_length = np.sqrt(slab_cart_proj**2 + (track_z[i] - track_z[i - 1])**2)
                    riserun = (track_z[i - 1] - track_z[i]) / (slab_cart_proj)
                    
                segment_length.append(slab_length * 1000)
                dip_holder.append(np.rad2deg(np.arctan(abs(riserun))))
                
            segment_thickness = []
            for k in range(len(thick_holder) - 1):
                if k != int(len(thick_holder) - 1):
                    segment_thickness.append([thick_holder[k], thick_holder[k + 1]])
                else:
                    segment_thickness.append([thick_holder[k]])
         
        else:
            segment_thickness = []
            for i in range(1, len(track_z)):
                segment_thickness.append([slab_thickness])
                
                if coordinate_system == 'Cartesian':
                    slab_cart_proj = np.sqrt( (track_x[i] - track_x[i - 1])**2 + (track_y[i] - track_y[i - 1])**2 )
                    slab_length = np.sqrt( (track_x[i] - track_x[i - 1])**2 + (track_y[i] - track_y[i - 1])**2 + (track_z[i] - track_z[i - 1])**2 )
                    
                elif coordinate_system == 'Spherical':
                    fwd_az, back_az, slab_cart_proj_m = geodesic.inv(track_x[i-1], track_y[i-1], track_x[i], track_y[i])
                    slab_cart_proj = slab_cart_proj_m / 1e3
                    slab_length = np.sqrt( slab_cart_proj**2 + (track_z[i] - track_z[i - 1])**2 )
                    
                riserun = (track_z[i - 1] - track_z[i]) / (slab_cart_proj)
                segment_length.append(slab_length * 1000)
                dip_holder.append(np.rad2deg(np.arctan(abs(riserun))))
                
        dips = []
        for k in range(len(dip_holder) - 1):
            if k != int(len(dip_holder) - 1):
                dips.append([dip_holder[k], dip_holder[k + 1]])
            else:
                dips.append([dip_holder[k]])
        
        # Write to the World Builder File
        world_builder_file.write('    {"coordinate":' + str(section_index) + ',\n')
        world_builder_file.write('     "segments":[')
        for current_segment in range(segment_num):
            if fill_segments == 'end':
                world_builder_file.write('    ' + segment_string_at_end(segment_length, segment_thickness, dips, top_truncation, section_index, total_sections, segment_num, current_segment))
            elif fill_segments == 'start':
                world_builder_file.write('    ' + segment_string(segment_length, segment_thickness, dips, top_truncation, section_index, total_sections, segment_num, current_segment))