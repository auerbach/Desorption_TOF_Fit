#!/usr/bin/env python
""" 
Desorption_TOF_Fit
Fit data from post permeation desorption TOF
"""
#==============================================================================
# 
# This program is dervied from time of Flight fitting program
#   Alessandro Genova, Leiden University, July 2012 
#   Francesco Nattino, Leiden University, June 2015  (vs 7)
#
#==============================================================================
import numpy as np
import lmfit
import subprocess
import shutil
# Parameters2 is a class derived from Parameters 
# adds support for a parm.glbl attribute for global fits
from Parameters2 import Parameters2
from Fit_control import Fit_control
# from compute_tof import TOF, generate_angles
from plot_fit import plot_fit
from cutoff import cutoff_function
import TOF_fit_global
from do_fit import do_fit
# from Parameters2 import Parameters2

import Data
import write_output


# -------------------------------------------------------------------------------------------------
#   sticking probability from direct TOF inversion
#--------------------------------------------------------------------------------------------------
def ProbFromTOFInversion(time, signal, n_dataset, parms, averaging_type, theta_angles,
                         prob_curve_type, cutoff_type, mass_molecules):
    # time of flight signal model. The function takes an uncorrected time in seconds and returns a signal
    
    mass_factor = np.sqrt(mass_molecules[n_dataset - 1] / glbl.massH2)
    ffr_dist  = parms['ffr_%i' % n_dataset].value * 1E-3
    max_tof   = parms['y_scale_%i' % n_dataset].value
    baseline = parms['baseline_%i' % n_dataset].value
    TimeCorr = parms['IonTOF_%i' % n_dataset].value * mass_factor
    Temperature = parms['Temp_%i' % n_dataset].value

    # subtract the ion flight time and eliminate singularity that would occur at time = 0
    time = time - TimeCorr * 1E-6  
    time = np.where(time != 0, time, np.repeat(0.01E-6, len(time)))
            
    CutOff = cutoff_function(parms, data, n_dataset, time, cutoff_type)

    velocity_distribution = 0.                        # Initialize Signal to 0

    if averaging_type == "PointDetector":    
        # Averaging performed taking into account different flight time for different angles
        for theta in theta_angles:
            # v = x / t = ( L / cos(theta) ) / t
            velocity = ffr_dist / (time * np.cos(np.radians(theta))) 
            Ekin = (0.5 * glbl.AMU * velocity ** 2.) * glbl.eVConst
            # Enorm = Ekin * np.cos( np.radians(theta) )**2. # Reaction probability depends on normal energy
            velocity_distribution = velocity_distribution                   \
                                   + velocity**4.                         \
                                   * np.exp(-Ekin / (glbl.kb * Temperature)) \
                                   * np.cos( np.radians(theta) )**2.      \
                                   * np.sin( np.radians(theta) ) \
                                     * glbl.theta_step

    elif averaging_type == "None":
        # No angular averaging
        velocity = ffr_dist / time
        Ekin = (0.5 * glbl.AMU * velocity ** 2.) * glbl.eVConst
        # Enorm = Ekin
        velocity_distribution = velocity**4. * np.exp(-Ekin / (glbl.kb * Temperature))
        
        # kludge to fix divide by zero runtime error        
        for ii in range(len(velocity_distribution)):
            if velocity_distribution[ii] == 0:
                velocity_distribution[ii] = 1E-100
        

    ProbFromTOF = (signal - baseline) / (max_tof * velocity_distribution * CutOff)
    Energy = 0.5 * glbl.AMU * (ffr_dist / time) ** 2. * glbl.eVConst
    

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
#         np.savetxt( TOFInvertedOutFile, np.column_stack(( TOFEnergy, TOFReactionProbability*Params['y_scale_%i' %NDataSet].value, Prob( TOFEnergy, NDataSet, Params, fit_function)*Params['y_scale_%i' %NDataSet].value )))
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

# set path to control files and fit output files by reading from the file
# path_to_fits.txt
with open("path_to_fits_and_editor.txt") as file: 
    lines = file.readlines()
    path_to_fits = lines[0].split(':',1)[1].strip()
    editor_cmd   = lines[1].split(':',1)[1].strip()

#------------------------------------------------------------------------------
# for testing:
#------------------------------------------------------------------------------
new_fit_number = '0035'
debug1 = False  # suppress dialog for opening control file
debug2 = False  # suppress writing excel file
#------------------------------------------------------------------------------
# begin1 comment out for testing
#------------------------------------------------------------------------------

if not debug1:
    # Get Last fit number from FitNumber.dat
    fn = path_to_fits + 'FitNumber.dat'
    print('fn =', fn)
    fit_number_file = open(fn, 'r+')
    fit_number = '{:04d}'.format(int(fit_number_file.readline()))
    old_fit_number = fit_number
    new_fit_number = fit_number
    
    while True:
        print('please enter old fit number: ', '[', old_fit_number, ']')
        ans = input('?')
        if ans:
            old_n = '{:04d}'.format(int(ans))
        else:
            old_n = old_fit_number
    
        if int(old_n) > int(old_fit_number):
            print('maximum allowed for old fit number is ', old_fit_number)
        else:
            break
    
    old_fit_number = old_n
    
    ans = input ('make new command file? [no]')
    if ans:
        if ans.upper()[0] == 'Y':
            new_fit_number = '{:04d}'.format(int(fit_number) + 1)
        fit_number_file.seek(0)
        fit_number_file.write(new_fit_number)
    else:
        new_fit_number = old_fit_number
    
    fit_number_file.close()
    
    old_file = path_to_fits + 'Fit' + old_fit_number + '.fit_in'
    new_file = old_file
    
    if old_fit_number != new_fit_number:
        new_file = path_to_fits + 'Fit' + new_fit_number + '.fit_in'
        shutil.copy2(old_file, new_file)
    
    subprocess.run(editor_cmd + ' "' + new_file + '"')
#------------------------------------------------------------------------------
# end of comment out for testing
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Parse the command file
#------------------------------------------------------------------------------
cmd_file = path_to_fits + 'Fit' + str(new_fit_number) + '.fit_in'
errors = fit_control.parse_cmd_file(cmd_file, glbl, parms)

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
    n_dataset       = i + 1
    fit_range       = glbl.fit_ranges[i]
    threshold       = glbl.thresholds[i]
    baseline_range  = glbl.baseline_ranges[i]
    fit_function    = glbl.functions[i]
    data_fn         = signal_filenames[i]
    averaging_type  = glbl.averaging_types[i]
    cutoff_type     = glbl.cutoff_types[i]
    
    if background_filenames:
        background_fn = background_filenames[i]
    else:
        background_fn = ""
          
    data.read_data(n_dataset, data_fn, background_fn, fit_range, baseline_range, threshold = threshold)
    
    n_min = data.fit_index_ranges[i][0]
    n_max = data.fit_index_ranges[i][1]
    t_min = data.datasets[i][0][n_min]
    t_max = data.datasets[i][0][n_max]
    
    if data.baselines[i]:
        parms['baseline_' + str(i+1)].value = data.baselines[i]
        parms['baseline_' + str(i+1)].vary  = False
        
    # check if temperature parameter was set.  If not use temperature from header
    try:
        parms['temp_' + str(i+1)].value
        pass
    except:
        parms.add('temp_' + str(i+1),
                  value = data.temperatures[i],
                  vary = False)
    
    time   = data.datasets[i][0]
    signal = data.datasets[i][1]
    fit_datasets.append([time[n_min:n_max], signal[n_min:n_max]])
    plot_datasets.append([time, signal])

#--------------------------------------------------------------------------------------------------
# Fit the data to model 
#--------------------------------------------------------------------------------------------------
fit_result = do_fit(fit_datasets, parms, averaging_type, fit_function, data.mass_molecules, glbl, data)
print(lmfit.fit_report(fit_result))

result_filename = \
    write_output.write_fit_out(glbl, data, path_to_fits, cmd_file, 
                               new_fit_number, fit_result, plot_datasets)

#==============================================================================
# uncomment when done with testing
#==============================================================================
if not debug2:
    write_output.write_results_excel_file(glbl, data, path_to_fits, new_fit_number, fit_result)

plot_fit(path_to_fits + result_filename)


    # if not calibration run
#==============================================================================
#     if not glbl.functions[i].lower().startswith('cal'):
#         WriteProbOutput(time, Signal, State, Label, NDataSet, fitResult.params,
#                         averaging_type, angles_list, fit_function)
#==============================================================================
