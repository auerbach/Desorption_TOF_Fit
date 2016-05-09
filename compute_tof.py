# -*- coding: utf-8 -*-
"""
functions to compute the time of flight signal vs time
"""

import numpy as np
from scipy import special
from cutoff import cutoff_function


# -------------------------------------------------------------------------------------------------
#   TOF - compute the signal vs time
# -------------------------------------------------------------------------------------------------
def TOF(time, n_dataset, parms, data, glbl):
    # time of flight signal model. The function takes an array of uncorrected time in seconds and 
    #   returns an array with the corresponing comptute signal
    
    # averaging_type = glbl.averaging_types[n_dataset-1]    
    # prob_curve_type = glbl.functions[n_dataset-1]
    angles = glbl.angles_list
    cutoff_type = glbl.cutoff_types[n_dataset-1]
    mass_molecules = data.mass_molecules
    mass_factor = np.sqrt(mass_molecules[n_dataset - 1] / glbl.massH2)
    y_scale   = parms['y_scale_%i' % n_dataset].value
    baseline = parms['baseline_%i' % n_dataset].value
    # tcutc    = Params['tcutc_%i'     %NDataSet].value
    # tcutw    = Params['tcutw_%i'     %NDataSet].value
    ion_tof = parms['ion_tof_%i' % n_dataset].value * mass_factor
    #Temperature = Params['Temp_%i' %NDataSet].value
    
    # subtract the ion flight time and eliminate singularity that would occur at time = 0
    time = time - ion_tof * 1E-6
    time = np.where(time != 0, time, np.repeat(0.01E-6, len(time)))
    
    cutoff = cutoff_function(parms, data, glbl, n_dataset, time, cutoff_type)
    signal0 = angular_averaging(time, n_dataset, angles, parms, glbl, data)
    signal  = signal0 * cutoff * y_scale + baseline
    
    return signal
    

# -------------------------------------------------------------------------------------------------
#   generate_angle -- generate the angles used for angular averaging.
# -------------------------------------------------------------------------------------------------
#==============================================================================
# def generate_angles(averaging_type, grid_type, \
#                     points_source, points_detector, \
#                     z_source, r_source, \
#                     z_aperture, r_aperture, \
#                     z_detector, length_detector, glbl):
#  
#==============================================================================
def generate_angles(averaging_type, glbl):
    grid_type       = glbl.grid_type
    points_source   = glbl.points_source
    points_detector = glbl.points_detector
    z_source        = glbl.z_source
    r_source        = glbl.r_source
    z_aperture      = glbl.z_aperture
    r_aperture      = glbl.r_aperture
    length_detector = glbl.length_laser
    z_detector      = glbl.z_laser
    
    if averaging_type == 'point_detector':
        angles = np.arange(0., glbl.ang_res + glbl.theta_step, glbl.theta_step)
    elif averaging_type == 'none':
        angles = [0.]
    elif averaging_type == 'line_detector':
        import generate_points
        grid_of_points_source  = generate_points.points_on_the_source(         \
            grid_type, points_source, r_source, z_source)

        grid_of_points_detector = generate_points.points_on_the_detection_line(
            n_points = points_detector, Length = length_detector, ZDetector = z_detector )

        angles = generate_points.ThetaPossibleTrajectories(         \
            grid_source = grid_of_points_source,                          \
            grid_detector= grid_of_points_detector,                        \
            ZAperture = z_aperture,                                      \
            RAperture = r_aperture)


        for i in range(len(angles)):
            angles[i] = np.degrees(angles[i])
        print("Considering ", len(angles),\
            " values of Theta in the angular averaging, minimum: %8.3f" \
              % min(angles), " deg , maximum: %8.3f" % max(angles), " deg.")
    return angles

    

# -------------------------------------------------------------------------------------------------
#   Angular Averaging
# -------------------------------------------------------------------------------------------------
def angular_averaging(time, n_dataset, angles, parms, glbl, data):
    averaging_type = glbl.averaging_types[n_dataset-1]
    prob_curve_type = glbl.functions[n_dataset-1]
    ffr_dist  = parms['ffr_%i' % n_dataset].value * 1E-3
    Temperature = parms['temp_%i' % n_dataset].value

    Signal = 0.  # Initialize Signal to 0
    
    mass = data.mass_molecules[n_dataset - 1]
    if averaging_type == "point_detector":
        # Averaging performed  taking into account different flight time for different angles, 
        # but assuming ionization occurring in one point
        for theta in angles :
        #for theta in [0]:
            time_nonzero = np.where(time != 0, time, np.repeat(0.01E-6, len(time)))
            
            # v = x / t = ( L / cos(theta) ) / t            
            velocity = ffr_dist / (time_nonzero * np.cos(np.radians(theta))) 
            
            # ASSUMING the reaction probability depends on normal energy
            Ekin = (0.5 * mass * velocity**2.) * glbl.eVConst
            
            Enorm = Ekin * np.cos( np.radians(theta) )**glbl.energy_angle_scaling 
            Signal = Signal + (velocity ** 4. * np.exp(-Ekin / (glbl.kb * Temperature)) * \
                               np.cos( np.radians(theta) ) ** 2. * \
                               Prob(Enorm, n_dataset, parms, prob_curve_type)) * \
                              np.sin( np.radians(theta) ) * glbl.theta_step

    elif averaging_type == "none":
        # No angular averaging performed
        time_nonzero = np.where(time != 0, time, np.repeat(0.01E-6, len(time)))
        velocity = ffr_dist / time_nonzero
        Ekin = (0.5 * mass * velocity**2.) * glbl.eVConst
        Enorm = Ekin
        Signal = (velocity ** 4. * np.exp(-Ekin / (glbl.kb * Temperature)) *
                  Prob(Enorm, n_dataset, parms, prob_curve_type))

    elif averaging_type == "line_detector":
        # Averaging along line, accounting for different flight time for different angles
        for theta in angles :
            # v = x / t = ( L / cos(theta) ) / t                        
            velocity = ffr_dist / (time * np.cos(np.radians(theta))) 
            Ekin = (0.5 * mass * velocity**2.) * glbl.eVConst
            Enorm = Ekin * np.cos( np.radians(theta) )**2 # Reaction probability depends on normal energy
            # Here no sin weight, since each theta value has a weight of one
            Signal = Signal + (velocity ** 4. * np.exp(-Ekin / 
            (glbl.kb * Temperature)) * np.cos(np.radians(theta)) ** 2. * 
            Prob(Enorm, n_dataset, parms, prob_curve_type)) * glbl.theta_step

    return Signal


# -------------------------------------------------------------------------------------------------
#   Prob -- reaction probability depending on curve type
# -------------------------------------------------------------------------------------------------
def Prob(En, NDataSet, Params, ProbCurveType):
    # Reaction probability functions. Takes an Energy (eV) returns a reaction probability
    
#==============================================================================
#     A  = Params['A_%i' %NDataSet].value
#     B  = Params['E0_%i' %NDataSet].value
#     # BI = Params['BI_%i' %NDataSet].value
#     C  = Params['W_%i' %NDataSet].value
#     # CI = Params['CI_%i' %NDataSet].value
#     # ni = Params['ni_%i' %NDataSet].value
#==============================================================================

    if ProbCurveType == "erf":
        E0 = Params['e0_%i' %NDataSet].value
        W  = Params['w_%i' %NDataSet].value
        return 1/2. * (1. + special.erf((En - E0)/ W ))
    
#==============================================================================
# 
#     elif fit_function == "GMP":
#         return ( A * np.exp(-np.exp(-(En - B)/C)) )
#     
# 
#     elif fit_function == "LGS":
#         return (A / np.power((1. + ni * np.exp(-(En - B)/C)), (1./ni)) )
#     
# 
#     elif fit_function == "FPC":
#         return A *(np.exp(-np.exp(-(En - B)/C)) / (1. +  np.exp(-(En - BI)/CI)) )
#     
# 
#     elif fit_function == "TPC":
#         return A *(np.exp(-np.exp(-(En - B)/C)) / (1. +  np.exp(-(En - B)/C)) )
#==============================================================================

    elif ProbCurveType.lower().startswith('cal'):
        return 1.
