#!/usr/bin/env python
""" 
Desorption_TOF_Fit
Fit data from post permeation desorption TOF
"""
#==============================================================================
# 
# Based on time of Flight fitting program
#   Alessandro Genova, Leiden University, July 2012 
#   Francesco Nattino, Leiden University, June 2015  (vs 7)
#
#==============================================================================

from lmfit import minimize, fit_report
from Parameters2 import Parameters2
import numpy as np
from scipy import special
import shutil
import subprocess
import unicodedata
from glob import glob

from Fit_control import Fit_control
from plot_fit import plot_fit
from cutoff import cutoff_function
from compute_tof import TOF, GenerateThetaAngles
import TOF_fit_global
from Parameters2 import Parameters2

import Data
import write_output

# -------------------------------------------------------------------------------------------------
#   do_fit -- function to perform the fit
# -------------------------------------------------------------------------------------------------
def do_fit(DataSets, Params, AveragingType, ProbCurveType, mass_molecules):
    
    # Generate Theta angles employed for angular averaging
    glbl.angles_list = GenerateThetaAngles(
                AveragingType=AveragingType,
                GridType=glbl.GridType,
                NPointsSource=glbl.NPointsSource,
                NPointsDetector=glbl.NPointsDetector,
                ZSource = glbl.ZSource,
                RSource = glbl.RSource,
                ZAperture = glbl.ZAperture, RAperture = glbl.RAperture,
                ZDetector = glbl.ZLaser, LengthDetector = glbl.LLaser)
                # ZDetector = ZFinal, LengthDetector = 2.*RFinal         \

    # Fit the data
    # Put data in the form needed by minimize
    X = []
    Y = []
    for DataSet in DataSets:
        X = X + list( DataSet[0] )
        Y = Y + list( DataSet[1] )
        
    # Perform non linear least sqaures fit
    result = minimize(residual, Params,
                      args=(X, Y, DataSets, AveragingType, glbl.angles_list,
                            ProbCurveType, mass_molecules))

    #-------------------------------------------------------------------------------------------------------------------
    # The minimizer result object has a Parameters object.  We want to have a Parameters2 object so that information
    #   about whether a parameter is global or not is retained.  The proper way to do this would be to create
    #   a new class that inherits from minimize but uses Parameters2 objects.  For now, we will simply make a new
    #   copy the value of global from the initial Parameters to the results.params object.  Even though this is a
    #   Parameters object and has no glbl attribute, it seems to be allowed to add an attribute this way
    #-------------------------------------------------------------------------------------------------------------------
    for p in result.params:
        result.params[p].glbl = Params[p].glbl

    return result

    
# -------------------------------------------------------------------------------------------------
#   residual function used in fit
# -------------------------------------------------------------------------------------------------
def residual(Params, X, Y, DataSets, AveragingType, ThetaAngles, ProbCurveType, mass_molecules):
    # residual function needed by minimize
    Resid = []

    for i in range( len( DataSets )):
        NDataSet = i + 1
        Resid = Resid + list(DataSets[i][1] \
                - TOF(DataSets[i][0], NDataSet, Params, data, glbl, AveragingType, ThetaAngles,
                      ProbCurveType, cutoff_type, mass_molecules))

    return np.array( Resid )
    


def ProbFromTOFInversion(Time, Signal, NDataSet, Params, AveragingType, ThetaAngles, 
                         ProbCurveType, cutoff_type, mass_molecules):
    # time of flight signal model. The function takes an uncorrected time in seconds and returns a signal
    
    mass_factor = np.sqrt(mass_molecules[NDataSet -1] / glbl.massH2)
    FFRDist  = Params['FFR_%i'      %NDataSet].value * 1E-3
    MaxTOF   = Params['Yscale_%i'   %NDataSet].value
    Baseline = Params['Baseline_%i' %NDataSet].value
    TimeCorr = Params['IonTOF_%i'   %NDataSet].value * mass_factor
    Temperature = Params['Temp_%i'  %NDataSet].value

    # subtract the ion flight time and eliminate singularity that would occur at time = 0
    Time = Time - TimeCorr * 1E-6  
    Time = np.where(Time != 0, Time, np.repeat(0.01E-6, len(Time)))
            
    CutOff = cutoff_function(Params, data, NDataSet, Time, cutoff_type)

    VelocityDistribution = 0.                        # Initialize Signal to 0

    if AveragingType == "PointDetector":    
        # Averaging performed taking into account different flight time for different angles
        for Theta in ThetaAngles:
            # v = x / t = ( L / cos(theta) ) / t
            Velocity = FFRDist /(Time * np.cos(np.radians(Theta))) 
            Ekin = (0.5 * glbl.massAmu * Velocity ** 2.) * glbl.eVConst
            # Enorm = Ekin * np.cos( np.radians(Theta) )**2. # Reaction probability depends on normal energy
            VelocityDistribution = VelocityDistribution                   \
                                   + Velocity**4.                         \
                                   * np.exp(-Ekin / (glbl.kb * Temperature)) \
                                   * np.cos( np.radians(Theta) )**2.      \
                                   * np.sin( np.radians(Theta) ) \
                                     * glbl.ThetaStep

    elif AveragingType == "None":
        # No angular averaging
        Velocity = FFRDist / Time
        Ekin = (0.5 * glbl.massAmu * Velocity ** 2.) * glbl.eVConst
        # Enorm = Ekin
        VelocityDistribution = Velocity**4. * np.exp(-Ekin / (glbl.kb * Temperature))
        
        # kludge to fix divide by zero runtime error        
        for ii in range(len(VelocityDistribution)):
            if VelocityDistribution[ii] == 0:
                VelocityDistribution[ii] = 1E-100
        

    ProbFromTOF = (Signal - Baseline ) / (MaxTOF * VelocityDistribution * CutOff)
    Energy = 0.5 * glbl.massAmu * (FFRDist / Time) ** 2. * glbl.eVConst
    

    return Energy, ProbFromTOF


#==============================================================================
#     TOFInvertedOutFile         = "DataSet_" + str( NDataSet ) + "_" + State + "_" + Label + "_ReacProb_Fit_" + fit_function + "_FromTOFInversion.dat"
#     # Reaction probability curve: minimum energy, maximum energy, deltaE (eV)
#     Emin = 0.0
#     Emax = 2.0
#     DeltaE = 0.01
# 
#     Energy = np.arange( Emin, Emax + DeltaE, DeltaE)
# 
#     # Calculate Reaction probability fitted curve
#     ReactionProbability = Prob(Energy, NDataSet, Params, fit_function)
#     np.savetxt( ReactionProbabilityOutFile, np.column_stack(( Energy, ReactionProbability)))
#     # print("WriteProbOutput: Reaction probability curve written to file ", ReactionProbabilityOutFile)
# 
#     # Reaction probability from TOF inversion; only possible with Point detector or no angular averaging (and normal energy scaling!)
#     if averaging_type != "LineDetector":
#         TOFEnergy, TOFReactionProbability = ProbFromTOFInversion(TOFTime, TOFData, NDataSet, Params, averaging_type, angles_list, fit_function, cutoff_type)
#         np.savetxt( TOFInvertedOutFile, np.column_stack(( TOFEnergy, TOFReactionProbability*Params['Yscale_%i' %NDataSet].value, Prob( TOFEnergy, NDataSet, Params, fit_function)*Params['Yscale_%i' %NDataSet].value )))
#         # print("WriteProbOutput: Inverted TOF written to file ", TOFInvertedOutFile)
#         # print('WriteProbOutput: Exiting')
#    return                                               
#==============================================================================









#============================================================================= 
#============================================================================= 
# ============================  Main PROGRAM  ================================
#============================================================================= 
#============================================================================= 

# import time
# start_time = time.time()

#------------------------------------------------------------------------------
# Create TOF_fit_global, Data, and Parameters, objects
#------------------------------------------------------------------------------
glbl = TOF_fit_global.TOF_fit_global()
data = Data.Data(glbl)
parms = Parameters2()
fit_control = Fit_control()

# define default path to control files and fit output17
pathToFits = 'Fits\\'
# pathToFits = 'd:\\users\\dja\\desktop\\permeation\\All-DJA\\Fits'

#------------------------------------------------------------------------------
# uncomment for testing:
#------------------------------------------------------------------------------
newFitNumber = '034'
newFitNumber = '017'

#------------------------------------------------------------------------------
# begin1 comment out for testing
#------------------------------------------------------------------------------
Get Last fit number from FitNumber.dat
fit_number_file = open(pathToFits + 'FitNumber.dat', 'r+')
fitNumber = '{:03d}'.format(int(fit_number_file.readline()))
oldFitNumber = fitNumber
newFitNumber = fitNumber

while True:
    print('please enter oldfit number: ', '[', oldFitNumber, ']')
    ans = input('?')
    if ans:
        old_n = '{:03d}'.format(int(ans))
    else:
        old_n = oldFitNumber

    if int(old_n) > int(oldFitNumber):
        print('maximum allowed for old fit number is ', oldFitNumber)
    else:
        break

oldFitNumber = old_n

ans = input ('make new command file? [no]')
if ans:
    if ans.upper()[0] == 'Y':
        newFitNumber = '{:03d}'.format(int(fitNumber)+1)
    fit_number_file.seek(0)
    fit_number_file.write(newFitNumber)
else:
    newFitNumber = oldFitNumber

fit_number_file.close()

oldFile = pathToFits + 'Fit' + oldFitNumber + '.tof_in'
newFile = oldFile

if oldFitNumber != newFitNumber:
    newFile = pathToFits + 'Fit' + newFitNumber + '.tof_in'
    shutil.copy2(oldFile, newFile)

subprocess.call(['npp.bat', newFile])

#------------------------------------------------------------------------------
# end of comment out for testing
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Parse the command file
#------------------------------------------------------------------------------
cmdFile = pathToFits + 'Fit' + str(newFitNumber) + '.tof_in'
errors = fit_control.parse_cmd_file(cmdFile, glbl, parms)

if len(errors) > 0:
    print('Errors\n', errors)
    raise SystemExit ('Error in command file' )

signal_filenames     = glbl.signal_filenames
background_filenames = glbl.background_filenames

# initalize lists
fit_datasets  = []
plot_datasets = []
# fit_index_ranges    = []  # this is in Data class

# Read data from input
for i in range(len(signal_filenames)):
    # t_min           = glbl.t_mins[i]
    # t_max           = glbl.t_maxs[i]
    fit_range       = glbl.fit_ranges[i]
    baseline_range  = glbl.baseline_ranges[i]
    fit_function    = glbl.functions[i]
    DataFile        = signal_filenames[i]
    averaging_type  = glbl.averaging_types[i]
    cutoff_type     = glbl.cutoff_types[i]
    
    if background_filenames  != [""]:
        BackgroundFile = background_filenames[i]
    else:
        BackgroundFile = ""
          
    data.read_data(DataFile, BackgroundFile, fit_range, baseline_range)
    
    n_min = data.fit_index_ranges[i][0]
    n_max = data.fit_index_ranges[i][1]
    t_min = data.datasets[i][0][n_min]
    t_max = data.datasets[i][0][n_max]
    
    if data.baselines[i]:
        parms['Baseline_' + str(i+1)].value = data.baselines[i]
        parms['Baseline_' + str(i+1)].vary  = False
    
    time   = data.datasets[i][0]
    Signal = data.datasets[i][1]
    fit_datasets.append([time[n_min:n_max], Signal[n_min:n_max]])
    plot_datasets.append([time, Signal])

#--------------------------------------------------------------------------------------------------
# Fit the data to model 
#--------------------------------------------------------------------------------------------------
fitResult = do_fit(fit_datasets, parms, averaging_type, fit_function, data.mass_molecules)
print(fit_report(fitResult))

result_filename = write_output.write_fit_out(glbl, data, pathToFits, cmdFile, newFitNumber, fitResult, plot_datasets)
# write_output.write_results_excel_file()
plot_fit(pathToFits + result_filename)





    # if not calibration run
#==============================================================================
#     if not glbl.functions[i].lower().startswith('cal'):
#         WriteProbOutput(time, Signal, State, Label, NDataSet, fitResult.params,
#                         averaging_type, angles_list, fit_function)
#==============================================================================
