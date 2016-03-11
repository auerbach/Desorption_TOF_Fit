# -*- coding: utf-8 -*-
"""
Functions to compute the time of flight signal vs time
"""

import numpy as np
from scipy import special
from Cutoff import cutoff_function

# -------------------------------------------------------------------------------------------------
#   TOF - compute the signal vs time
# -------------------------------------------------------------------------------------------------
def TOF(Time, NDataSet, Params, data, glbl, AveragingType, ThetaAngles, ProbCurveType,
        cutoff_type, mass_molecules):
    # Time of flight signal model. The function takes an uncorrected time in seconds and returns a signal
    #Dist  = Params['FFRDist_%i'   %NDataSet].value
    mass_factor = np.sqrt(mass_molecules[NDataSet -1] / glbl.massH2)
    Yscale   = Params['Yscale_%i'    %NDataSet].value
    Baseline = Params['Baseline_%i'  %NDataSet].value
    # TCutC    = Params['TCutC_%i'     %NDataSet].value
    # TCutW    = Params['TCutW_%i'     %NDataSet].value
    TimeCorr = Params['IonTOF_%i'    %NDataSet].value * mass_factor
    #Temperature = Params['Temp_%i' %NDataSet].value
    
    # subtract the ion flight time and eliminate singularity that would occur at Time = 0    
    Time = Time - TimeCorr * 1E-6  
    Time = np.where(Time != 0, Time, np.repeat(0.01E-6, len(Time)))    
    
    CutOff = cutoff_function(Params, data, glbl, NDataSet, Time, cutoff_type)
    Signal0 = AngularAveraging(Time, NDataSet, Params, data, AveragingType, ThetaAngles, ProbCurveType, glbl)
    Signal  = Signal0 * CutOff * Yscale + Baseline
    
#==============================================================================
#     for i, sig in enumerate(Signal):
#         if np.isnan(sig) and Time[i] == 0.0:
#             print('NAN found i, ix, Time=', i, ix, Time[i])
#             
#             Signal[i] = 0.
#==============================================================================
        
    return Signal
    

# -------------------------------------------------------------------------------------------------
#   Angular Averaging
# -------------------------------------------------------------------------------------------------
def AngularAveraging(Time, NDataSet, Params, data, 
                     AveragingType, ThetaAngles, ProbCurveType, glbl):

    FFRDist  = Params['FFR_%i'   %NDataSet].value * 1E-3
    Temperature = Params['Temp_%i' %NDataSet].value

    Signal = 0.  # Initialize Signal to 0
    
    mass = data.mass_molecules[NDataSet-1]
    if AveragingType == "PointDetector":
        # Averaging performed  taking into account different flight time for different angles, but assuming ionization occurring in one point
        for Theta in ThetaAngles :
        #for Theta in [0]:
            time_nonzero = np.where(Time != 0, Time, np.repeat(0.01E-6, len(Time)))
            Velocity = FFRDist /(time_nonzero * np.cos( np.radians(Theta) ) ) # v = x / t = ( L / cos(theta) ) / t
            Ekin = (0.5 * mass * Velocity**2.) * glbl.eVConst
            Enorm = Ekin * np.cos( np.radians(Theta) )**2. # Reaction probability depends on normal energy
            Signal = Signal + (Velocity ** 4. * np.exp(-Ekin / (glbl.kb * Temperature)) * \
                               np.cos( np.radians(Theta) ) ** 2. * \
                               Prob(Enorm, NDataSet, Params, ProbCurveType)) *           \
                                np.sin( np.radians(Theta) ) * glbl.ThetaStep

    elif AveragingType == "None":
        # No angular averaging performed
        time_nonzero = np.where(Time != 0, Time, np.repeat(0.01E-6, len(Time)))
        Velocity = FFRDist / time_nonzero
        Ekin = (0.5 * mass * Velocity**2.) * glbl.eVConst
        Enorm = Ekin
        Signal = (Velocity ** 4. * np.exp(-Ekin / (glbl.kb * Temperature)) *
                  Prob(Enorm, NDataSet, Params, ProbCurveType))

    elif AveragingType == "LineDetector":
        # Averaging along line, accounting for different flight time for different angles
        for Theta in ThetaAngles :
                        Velocity = FFRDist/(Time * np.cos( np.radians(Theta) ) ) # v = x / t = ( L / cos(theta) ) / t
                        Ekin = (0.5 * mass * Velocity**2.) * glbl.eVConst
                        Enorm = Ekin * np.cos( np.radians(Theta) )**2 # Reaction probability depends on normal energy
                        # Here no sin weight, since each Theta value has a weight of one
                        Signal = Signal + (Velocity ** 4. * np.exp(-Ekin / (glbl.kb * Temperature)) * np.cos(np.radians(Theta)) ** 2. * Prob(Enorm, NDataSet, Params, ProbCurveType)) * glbl.ThetaStep

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

    if ProbCurveType == "ERF":
        E0 = Params['E0_%i' %NDataSet].value
        W  = Params['W_%i' %NDataSet].value
        return 1/2. * (1. + special.erf((En - E0)/ W ))
    
#==============================================================================
# 
#     elif ProbCurveType == "GMP":
#         return ( A * np.exp(-np.exp(-(En - B)/C)) )
#     
# 
#     elif ProbCurveType == "LGS":
#         return (A / np.power((1. + ni * np.exp(-(En - B)/C)), (1./ni)) )
#     
# 
#     elif ProbCurveType == "FPC":
#         return A *(np.exp(-np.exp(-(En - B)/C)) / (1. +  np.exp(-(En - BI)/CI)) )
#     
# 
#     elif ProbCurveType == "TPC":
#         return A *(np.exp(-np.exp(-(En - B)/C)) / (1. +  np.exp(-(En - B)/C)) )
#==============================================================================

    elif ProbCurveType.lower().startswith('cal'):
        return 1.
