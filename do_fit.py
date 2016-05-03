# -*- coding: utf-8 -*-
import numpy as np
import lmfit 
from compute_tof import TOF, generate_angles

from compute_tof import TOF, generate_angles
# from Parameters2 import Parameters2

# -------------------------------------------------------------------------------------------------
#   do_fit -- function to perform the fit
# -------------------------------------------------------------------------------------------------
def do_fit(datasets, parms, averaging_type, prob_curve_type, mass_molecules, glbl, data):
    #----------------------------------------------------------------------------------------------    
    # Generate Theta angles employed for angular averaging
    #----------------------------------------------------------------------------------------------
    #glbl.angles_list = generate_angles(
    #                        averaging_type, glbl.grid_type, glbl.points_source,
    #                        glbl.points_detector, glbl.z_source, glbl.r_source, glbl.z_aperture,
    #                        glbl.r_aperture, glbl.z_laser, glbl.length_laser, glbl)
    glbl.angles_list = generate_angles(averaging_type, glbl)
    #----------------------------------------------------------------------------------------------
    # Fit the data
    #----------------------------------------------------------------------------------------------
    # Put data in the form needed by minimize
    X = []
    Y = []
    for DataSet in datasets:
        X = X + list( DataSet[0] )
        Y = Y + list( DataSet[1] )
        
    # Perform the fit
    # result = lmfit.minimize(residual, parms,
    #                   args=(X, Y, datasets, averaging_type, 
    #                         prob_curve_type, mass_molecules, glbl, data))
    result = lmfit.minimize(residual, parms,
                      args=(X, Y, datasets, data, glbl))

    #----------------------------------------------------------------------------------------------
    # The minimizer result object has a Parameters object.  We want to have a Parameters2 object 
    # so that information about whether a parameter is global or not is retained.  
    #
    # The proper way to do this would be to create     #   a new class that inherits from minimize 
    # but uses Parameters2 objects.  For now, we will simply make a new copy the value of global from 
    # the initial Parameters to the results.params object.Even though this is a Parameters object 
    # and has no glbl attribute, it seems we can add an attribute this way
    #----------------------------------------------------------------------------------------------
    for p in result.params:
        result.params[p].glbl = parms[p].glbl

    return result

    
# -------------------------------------------------------------------------------------------------
#   residual function used by minimize in performing the fit
# -------------------------------------------------------------------------------------------------
def residual(parms, x, y, datasets, data, glbl):
    resid = []

    for i in range(len(datasets)):
        NDataSet = i + 1
        resid = resid + list(datasets[i][1] - TOF(datasets[i][0], NDataSet, parms, data, glbl))
    return np.array(resid)