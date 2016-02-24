#!/usr/bin/env python
""" 
Desorption_TOF_Fit
Fit data from post permeation desorption TOF
"""
#==============================================================================
# 
# Based on Time of Flight fitting program
#   Alessandro Genova
#   July 2012, Leiden University
#
#   Francesco Nattino, Leiden University, June 2015  (vs 7)
#
#==============================================================================

# import os
# import sys
import unicodedata
import shutil
import subprocess
import numpy as np
from scipy import special
from lmfit import minimize, fit_report

from ParseCmdFile import parseCmdFile
from PlotFit import PlotFit

# from Parameters2 import Parameter2, Parameters2

import GlobalVariables as glbl
from read_data import Data

def GenerateThetaAngles(AveragingType, GridType,        \
                        NPointsSource, NPointsDetector, \
                        ZSource,    RSource,            \
                        ZAperture, RAperture,           \
                        ZDetector,  LengthDetector):

    if AveragingType == 'PointDetector':
        ThetaAngles = np.arange( 0., glbl.AngRes + glbl.ThetaStep, glbl.ThetaStep )
    elif AveragingType == 'None':
        ThetaAngles = [0.] 
    elif AveragingType == 'LineDetector':
        import GeneratePoints
        GridOfPointsSource  = GeneratePoints.PointsOnTheSource(         \
            GridType = GridType, ZSource = ZSource , RSource = RSource, \
            NPoints = NPointsSource)
        
        GridOfPointsDetector = GeneratePoints.PointsOnTheDetectionLine( \
            ZDetector = ZDetector ,         \
            NPoints = NPointsDetector,      \
            Length = LengthDetector )
               
        ThetaAngles = GeneratePoints.ThetaPossibleTrajectories(         \
            GridSource   = GridOfPointsSource,                          \
            GridDetector = GridOfPointsDetector,                        \
            ZAperture = ZAperture,                                      \
            RAperture = RAperture)
        
            
        for i in range( len( ThetaAngles )):
            ThetaAngles[i] = np.degrees( ThetaAngles[i] )
        print("Considering ", len(ThetaAngles ),\
            " values of Theta in the angular averaging, minimum: %8.3f"\
            %min( ThetaAngles), " deg , maximum: %8.3f" %max( ThetaAngles) ," deg.")
    return ThetaAngles


# -------------------------------------------------------------------------------------------------
#   FitData -- function to perform the fit
# -------------------------------------------------------------------------------------------------
def FitData( DataSets, Params, AveragingType, ThetaAngles, ProbCurveType):
    # Fit the data
    # Give to the datasets a form that "minimize" likes
    X, Y = [], []
    

      
    for DataSet in DataSets:
        X = X + list( DataSet[0] )
        Y = Y + list( DataSet[1] )
        
    # The fit is performed
    result = minimize( Residual, Params, 
                      args=(X, Y, DataSets, AveragingType, ThetaAngles, 
                            ProbCurveType) )
    
    return result

    
# -------------------------------------------------------------------------------------------------
#   Residual -- residual function used in fit
# -------------------------------------------------------------------------------------------------
def Residual( Params, X, Y, DataSets, AveragingType, ThetaAngles, ProbCurveType):
    # Residual function needed by minimize
    Resid = []

    for i in range( len( DataSets )):
        NDataSet = i + 1
        Resid = Resid + list(DataSets[i][1] - TOF(DataSets[i][0], NDataSet, Params,
                                                  AveragingType, ThetaAngles, ProbCurveType))

    return np.array( Resid )
    

# -------------------------------------------------------------------------------------------------
#   TOF - compute the signal vs time
# -------------------------------------------------------------------------------------------------
def TOF(Time, NDataSet, Params, AveragingType, ThetaAngles, ProbCurveType, Debug=False):
    # Time of flight signal model. The function takes an uncorrected time in seconds and returns a signal
    #FFRDist  = Params['FFRDist_%i'   %NDataSet].value
    Yscale   = Params['Yscale_%i'    %NDataSet].value
    Baseline = Params['Baseline_%i'  %NDataSet].value
    TCutC    = Params['TCutC_%i'     %NDataSet].value
    TCutW    = Params['TCutW_%i'     %NDataSet].value
    TimeCorr = Params['IonTOF_%i'    %NDataSet].value
    #Temperature = Params['Temp_%i' %NDataSet].value
    Time   = Time - TimeCorr * 1E-6                            # Correct the time
    CutOff = 0.5 * (1. - np.tanh((Time - TCutC*1E-6) / (TCutW*1E-6)))      # CutOff function
            #used to model the experimental loss of low energy ions
    Signal0 = AngularAveraging( AveragingType, Time, NDataSet, Params, ThetaAngles)       
    Signal  = Signal0 * CutOff * Yscale + Baseline
    
    if Debug :
        print('TOF     : Signal0 = ', Signal0)
        print('        : CutOff  = ', CutOff)
        print('        : Yscale  = ', Yscale) 
        print('        : Signal  = ', Signal)

    return Signal
    

# -------------------------------------------------------------------------------------------------
#   Angular Averaging
# -------------------------------------------------------------------------------------------------
def AngularAveraging( AveragingType, Time, NDataSet, Params, ThetaAngles):

    FFRDist  = Params['FFR_%i'   %NDataSet].value
    Temperature = Params['Temp_%i' %NDataSet].value

    Signal = 0.  # Initialize Signal to 0
    
    mass = data.mass_molecules[NDataSet-1]
    if AveragingType == "PointDetector":
        # Averaging performed  taking into account different flight time for different angles, but assuming ionization occurring in one point
        for Theta in ThetaAngles :
        #for Theta in [0]:
            Velocity = FFRDist /(Time * np.cos( np.radians(Theta) ) ) # v = x / t = ( L / cos(theta) ) / t
            Ekin = (0.5 * mass * Velocity**2.) * glbl.eVConst
            Enorm = Ekin * np.cos( np.radians(Theta) )**2. # Reaction probability depends on normal energy
            Signal = Signal + (Velocity**4. * np.exp( -Ekin / (glbl.kb * Temperature) ) * \
                                np.cos( np.radians(Theta) )**2. *                         \
                                Prob(Enorm, NDataSet, Params, ProbCurveType)) *           \
                                np.sin( np.radians(Theta) ) * glbl.ThetaStep

    elif AveragingType == "None":
        # No angular averaging performed
        Velocity = FFRDist / Time # v = L / t
        Ekin = (0.5 * mass * Velocity**2.) * glbl.eVConst
        Enorm = Ekin
        Signal = (Velocity**4. * np.exp( -Ekin / (glbl.kb * Temperature) ) * 
                  Prob(Enorm, NDataSet, Params, ProbCurveType))

    elif AveragingType == "LineDetector":
        # Averaging along line, accounting for different flight time for different angles
        for Theta in ThetaAngles :
                        Velocity = FFRDist/(Time * np.cos( np.radians(Theta) ) ) # v = x / t = ( L / cos(theta) ) / t
                        Ekin = (0.5 * mass * Velocity**2.) * glbl.eVConst
                        Enorm = Ekin * np.cos( np.radians(Theta) )**2 # Reaction probability depends on normal energy
                        # Here no sin weight, since each Theta value has a weight of one
                        Signal = Signal + (Velocity**4. * np.exp( -Ekin / (glbl.kb * Temperature) ) * np.cos( np.radians(Theta) )**2. * Prob(Enorm, NDataSet, Params, ProbCurveType)) * glbl.ThetaStep

    return Signal


# -------------------------------------------------------------------------------------------------
#   Prob -- reaction probability depending on curve type
# -------------------------------------------------------------------------------------------------
def Prob(En, NDataSet, Params, ProbCurveType):
    # Reaction probability functions. Takes an Energy (eV) returns a reaction probability
    A  = Params['A_%i' %NDataSet].value
    B  = Params['E0_%i' %NDataSet].value
    # BI = Params['BI_%i' %NDataSet].value
    C  = Params['W_%i' %NDataSet].value
    # CI = Params['CI_%i' %NDataSet].value
    # ni = Params['ni_%i' %NDataSet].value
    
    


    if ProbCurveType == "ERF":
        return A/2. * (1. + special.erf((En - B)/ C ))
    
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

    elif ProbCurveType == "Calibration":
        return 1.

def ProbFromTOFInversion(Time, Signal, NDataSet, Params, AveragingType, ThetaAngles, ProbCurveType):
    # Time of flight signal model. The function takes an uncorrected time in seconds and returns a signal
    

    FFRDist  = Params['FFR_%i'      %NDataSet].value
    MaxTOF   = Params['Yscale_%i'   %NDataSet].value
    Baseline = Params['Baseline_%i' %NDataSet].value
    TCutC    = Params['TCutC_%i'    %NDataSet].value
    TCutW    = Params['TCutW_%i'    %NDataSet].value
    TimeCorr = Params['IonTOF_%i'   %NDataSet].value
    Temperature = Params['Temp_%i'  %NDataSet].value

    Time   = Time - TimeCorr * 1E-6                            # Correct the time
    CutOff = 0.5 * (1. - np.tanh((Time - TCutC*1E-6) / (TCutW*1E-6)))    # CutOff function used to model the experimental loss of low energy ions
    

    VelocityDistribution = 0.                        # Initialize Signal to 0

    if AveragingType == "PointDetector":    
        # Averaging performed taking into account different flight time for different angles
        for Theta in ThetaAngles:
            Velocity = FFRDist /(Time * np.cos( np.radians(Theta) ) ) # v = x / t = ( L / cos(theta) ) / t
            Ekin = (0.5 * glbl.massAmu * Velocity**2.) * glbl.eVConst
            # Enorm = Ekin * np.cos( np.radians(Theta) )**2. # Reaction probability depends on normal energy
            VelocityDistribution = VelocityDistribution                   \
                                   + Velocity**4.                         \
                                   * np.exp( -Ekin / (glbl.kb * Temperature) ) \
                                   * np.cos( np.radians(Theta) )**2.      \
                                   * np.sin( np.radians(Theta) )           \
                                   * glbl.ThetaStep

    elif AveragingType == "None":
        # No angular averaging
        Velocity = FFRDist / Time
        Ekin = (0.5 * glbl.massAmu * Velocity**2.) * glbl.eVConst
        # Enorm = Ekin
        VelocityDistribution = Velocity**4. * np.exp( -Ekin / (glbl.kb * Temperature))
        
        # kludge to fix divide by zero runtime error        
        for ii in range(len(VelocityDistribution)):
            if VelocityDistribution[ii] == 0:
                VelocityDistribution[ii] = 1E-100
        

    ProbFromTOF = (Signal - Baseline ) / (MaxTOF * VelocityDistribution * CutOff)
    Energy = 0.5 * glbl.massAmu * ( FFRDist / Time )**2. * glbl.eVConst
    

    return Energy, ProbFromTOF


#==============================================================================
#     TOFInvertedOutFile         = "DataSet_" + str( NDataSet ) + "_" + State + "_" + Label + "_ReacProb_Fit_" + ProbCurveType + "_FromTOFInversion.dat"
#     # Reaction probability curve: minimum energy, maximum energy, deltaE (eV)
#     Emin = 0.0
#     Emax = 2.0
#     DeltaE = 0.01
# 
#     Energy = np.arange( Emin, Emax + DeltaE, DeltaE)
# 
#     # Calculate Reaction probability fitted curve
#     ReactionProbability = Prob(Energy, NDataSet, Params, ProbCurveType)
#     np.savetxt( ReactionProbabilityOutFile, np.column_stack(( Energy, ReactionProbability)))
#     # print("WriteProbOutput: Reaction probability curve written to file ", ReactionProbabilityOutFile)
# 
#     # Reaction probability from TOF inversion; only possible with Point detector or no angular averaging (and normal energy scaling!)
#     if AveragingType != "LineDetector":
#         TOFEnergy, TOFReactionProbability = ProbFromTOFInversion(TOFTime, TOFData, NDataSet, Params, AveragingType, ThetaAngles, ProbCurveType)
#         np.savetxt( TOFInvertedOutFile, np.column_stack(( TOFEnergy, TOFReactionProbability*Params['Yscale_%i' %NDataSet].value, Prob( TOFEnergy, NDataSet, Params, ProbCurveType)*Params['Yscale_%i' %NDataSet].value )))
#         # print("WriteProbOutput: Inverted TOF written to file ", TOFInvertedOutFile)
#         # print('WriteProbOutput: Exiting')
#    return                                               
#==============================================================================




#============================================================================= 
#============================================================================= 
# ============================  Main PROGRAM  ================================
#============================================================================= 
#============================================================================= 

# create a Data object to read and store the data
data = Data()

# define default path to control files and fit output
pathToFits = 'Fits\\'


#==============================================================================
# # uncomment for testing:
# cmdFileName    = 'fit010_test1'
# # cmdFileName    = 'fit010'
# cmdFileTesting = pathToFits + cmdFileName + '.tof_in'
# 
#==============================================================================
# Get Last fit number from FitNumber.dat
fit_number_file = open(pathToFits + 'FitNumber.dat', 'r+')
fitNumber = '{:03d}'.format(int(fit_number_file.readline()))
oldFitNumber = fitNumber
newFitNumber = fitNumber

# comment out for testing
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

oldFile = pathToFits + 'fit' + oldFitNumber + '.tof_in'
newFile = oldFile

if oldFitNumber != newFitNumber:
    newFile = pathToFits + 'fit' + newFitNumber + '.tof_in'
    shutil.copy2(oldFile, newFile)

# subprocess.call(['npp.bat', newFile])
cmdFile = pathToFits + 'fit' + newFitNumber + '.tof_in'

#==============================================================================
# # for testing
# cmdFile = cmdFileTesting
#==============================================================================

#==============================================================================
# # Parse the command file
#==============================================================================
parms, functions, signalFiles, backgroundFiles, errors = parseCmdFile(cmdFile)

if len(errors) > 0:
    print('Errors\n', errors)
    raise SystemExit ('Error in command file' )
    
AveragingType = glbl.AveragingTypes[0]

#==============================================================================
# # use Francesco's names until have time to change
#==============================================================================
DataFiles = signalFiles
BackgroundFiles = backgroundFiles

# initalize lists
states = []         # list if tuples (v, j)
temperatures = []
DataSets = []
PlotDataSets = []
Params = parms

# Read data from input
fit_ranges = []
for i in range( len( DataFiles)):
    Tmin = glbl.Tmins[i]
    Tmax = glbl.Tmaxs[i]
    ProbCurveType = functions[i]
    DataFile = DataFiles[i]
    
    if BackgroundFiles  != [ "" ]:
        BackgroundFile = BackgroundFiles[i]
    else:
        BackgroundFile = BackgroundFiles[0]
          
    data.read_data(DataFile, BackgroundFile, Tmin, Tmax)  
    
    Nmin = data.fit_ranges[i][0]
    Nmax = data.fit_ranges[i][1]    
    Tmin = data.datasets[i][0][Nmin]
    Tmax = data.datasets[i][0][Nmax]
    
    fit_ranges.append([Nmin,Nmax])    
    states.append( data.states[i] )
    temperatures.append(data.temperatures[i])
    
    Time   = data.datasets[i][0]
    Signal = data.datasets[i][1]
    DataSets.append([Time[Nmin:Nmax], Signal[Nmin:Nmax]])
    PlotDataSets.append([Time, Signal])
    pass
    
  

# Generate Theta angles employed for angular averaging
#==============================================================================
# def GenerateThetaAngles(AveragingType, GridType,        \
#                         NPointsSource, NPointsDetector, \
#                         ZSource,    RSource,            \
#                         ZAperture, RAperture,           \
#                         ZDetector,  LengthDetector):
#==============================================================================
ThetaAngles = GenerateThetaAngles(AveragingType=AveragingType, GridType=glbl.GridType,
                                  NPointsSource=glbl.NPointsSource, 
                                  NPointsDetector=glbl.NPointsDetector, 
                                  ZSource = glbl.ZSource, RSource = glbl.RSource,
                                  ZAperture = glbl.ZAperture, RAperture = glbl.RAperture,
                                  ZDetector = glbl.ZLaser, LengthDetector = glbl.LLaser)
#    ZDetector = ZFinal,          LengthDetector = 2.*RFinal         \
  
#==============================================================================
# print()
# print('#-----------------------------')
# print('Angular Averaging Parameters:')
# print('#-----------------------------')
# print('nps =', glbl.NPointsSource, 'npd =', glbl.NPointsDetector)
# print('zs  =', glbl.ZSourceresult_file.write('rs  =', glbl.RSource )
# print('za  =', glbl.ZAperture,     'ra  =', glbl.RAperture)     
# print('zf  =', glbl.ZFinal,        'rf  =', glbl.RFinal)
# print('it mse')
# print('ThetaAngles =', ThetaAngles)
# print()
#==============================================================================

#--------------------------------------------------------------------------------------------------
# Fit the data to model 
#--------------------------------------------------------------------------------------------------
fitResult = FitData( DataSets, Params, AveragingType, ThetaAngles, ProbCurveType)
print(fit_report(fitResult))

state_string = 'v' + str(states[0][0]) + 'j' +str(states[0][1])
result_file_name = 'Fit' + newFitNumber + '_' + state_string + '_' + ProbCurveType + '.fit_out'


#==============================================================================
# # for testing 
# result_file_name = cmdFileName + '_' + state_string + '_' + ProbCurveType + '.fit_out'
#==============================================================================

with open(pathToFits + result_file_name, 'w') as result_file:
    result_file.write('#\n')
    result_file.write('# Fit {} Results\n'.format(newFitNumber))
    result_file.write('#\n')
    result_file.write('#' + 60*'-' + '\n')
    result_file.write('# Control file: {}\n'.format(cmdFile))
    result_file.write('#' + 60*'-' + '\n')
    
    #----------------------------------------------------------------------------------------------
    # write the control file to the result file
    #----------------------------------------------------------------------------------------------
    with open(cmdFile,'r') as cmd:
        cmd_lines = cmd.readlines()

    for line in cmd_lines:
        result_file.write('# ' + line )
    
    result_file.write('\n#' + 60*'-' + '\n')    
    result_file.write('# End Control File\n')
    result_file.write('#' + 60*'-' + '\n')
    result_file.write('\n')
    
    #----------------------------------------------------------------------------------------------
    # write the angular averaging parameters
    #----------------------------------------------------------------------------------------------
    result_file.write('#' + 60*'-' + '\n')
    result_file.write('# Angular Averaging Parameters\n')
    result_file.write('#' + 60*'-' + '\n')
    result_file.write('# NPointsSource   : ' + str(glbl.NPointsSource)   + '\n')
    result_file.write('# NPointsDetector : ' + str(glbl.NPointsDetector) + '\n')
    result_file.write('# ZSource         : ' + str(glbl.ZSource)         + '\n')
    result_file.write('# RSource         : ' + str(glbl.RSource )        + '\n')
    result_file.write('# ZAperture       : ' + str(glbl.ZAperture)       + '\n')
    result_file.write('# RAperture       : ' + str(glbl.RAperture)       + '\n')    
    result_file.write('# ZFinal          : ' + str(glbl.ZFinal)          + '\n')
    result_file.write('# RFinal          : ' + str(glbl.RFinal)          + '\n')
    # result_file.write('# it mse\n')
    
    for i in range(0, len(ThetaAngles), 10):
        result_file.write('# ThetaAngles : ' + str(ThetaAngles[i : i+10]) + '\n')
    result_file.write('#' + 60*'-' + '\n')    
    result_file.write('# End Angular Averaging Parametern')
    result_file.write('#' + 60*'-' + '\n')
    result_file.write('\n')
    
    #----------------------------------------------------------------------------------------------
    # write the fit report to the result file
    #----------------------------------------------------------------------------------------------
    result_file.write('#' + 60*'-' + '\n')
    result_file.write('# Fit Report\n')
    result_file.write('#' + 60*'-' + '\n')    

    report = fit_report(fitResult).split('\n')
    for line in report:
        result_file.write('# ' + line + '\n')
    result_file.write('#' + 60*'-' + '\n')    
    result_file.write('# End Fit Report\n')
    result_file.write('#' + 60*'-' + '\n')
    result_file.write('\n')
    
    #----------------------------------------------------------------------------------------------
    # the number of datasets to plot
    #----------------------------------------------------------------------------------------------
    result_file.write('#' + 60*'-' + '\n')
    result_file.write('# Labels and data for plots\n')
    result_file.write('#' + 60*'-' + '\n')
    result_file.write('# Number of plots : ' + str(len(DataSets)) + '\n')
    result_file.write('#' + 60*'-' + '\n')
    result_file.write('\n')
    
    #----------------------------------------------------------------------------------------------
    # for each data set, write plot title, labels, and data
    #----------------------------------------------------------------------------------------------
    for i in range( len( DataSets )):
        n_dataset = i + 1
        
        result_file.write('#' + 60*'-' + '\n')
        result_file.write('# Plot ' + str(n_dataset) + '\n')
        result_file.write('#' + 60*'-' + '\n')
        
        #------------------------------------------------------------------------------------------
        # write the plot title
        #------------------------------------------------------------------------------------------
        result_file.write('# Title    : Fit {:03d}_{}'.format(int(newFitNumber), n_dataset) + '\n')

        #------------------------------------------------------------------------------------------
        # write the plot label
        #------------------------------------------------------------------------------------------
        result_file.write('# Label    : ' + data.original_signal_names[i] + '\n')
        result_file.write('# Label    : \n')
        # result_file.write('# Label    : ' + '----------------------\n')
        
        avg_type = AveragingType
        if avg_type == 'None' : 
            avg_type = 'No Angular Averaging'
        elif avg_type == 'PointDetector':
            avg_type = 'Point Detector'
        elif avg_type == 'LineDetector':
            avg_type = 'Line Detector' 
        
        # result_file.write('# Label    : Averaging: ' + avg_type + '\n')
        result_file.write('# Label    : Function:  ' + ProbCurveType + '\n')
        result_file.write('# Label    : ' + avg_type + '\n')
        result_file.write('# Label    : \n')         
        parms = fitResult.params     
        parm_list = ['E0', 'W', 'TCutC', 'TCutW']
        for p in parm_list:
            p_ = p + '_' + str(i+1)
            result_file.write('# Label    : {:6s}{:>7.3f}'.format(p, parms[p_].value))
            if not parms[p_].vary:
                result_file.write(' (fixed)\n')
            else:
                stderr = '{:6.3f}'.format(parms[p_].stderr)
                result_file.write(' (' + unicodedata.lookup('plus-minus sign') 
                                  + str(stderr) +')\n')
        
        #------------------------------------------------------------------------------------------
        # write the number of data lines and range of lines fit
        #-----------------------------------------------------------------------------------------
        result_file.write('# Npoints    : ' + str(len(Time)) +'\n')
        result_file.write('# Nmin, Nmax : ' + str(Nmin) + ',' +str(Nmax) + '\n')
        result_file.write('# Tmin, Tmax : ' + str(Tmin * 1E6) + ',' +str(Tmax * 1E6) + '\n')
               
        result_file.write('#' + 60*'-' + '\n')
        
        #------------------------------------------------------------------------------------------
        # write the plot data
        #------------------------------------------------------------------------------------------
        result_file.write('# Time            Signal                Fit\n')
        result_file.write('#' + 60*'-' + '\n')
        result_file.write('# Begin data\n')
        
        Time     =  PlotDataSets[i][0] *1E6  # convert to microseconds
        Signal   =  PlotDataSets[i][1]
        Fit      =  TOF(PlotDataSets[i][0], n_dataset, fitResult.params,
                        AveragingType, ThetaAngles, ProbCurveType) 
        for j in range(len(Time)):
            result_file.write('{:6.2f} {:20.5e} {:20.5e}\n'.format(Time[j], Signal[j], Fit[j]))
        
        result_file.write('# End data')
        result_file.write('#' + 60*'-' + '\n')
        result_file.write('# End Plot ' + str(n_dataset))
            
        
result_file.close()

PlotFit(pathToFits + result_file_name)
    # if DataType != 'Calibration':
#==============================================================================
#     if functions[i] != 'Calibration':
#         WriteProbOutput(Time, Signal, State, Label, NDataSet, fitResult.params, 
#                         AveragingType, ThetaAngles, ProbCurveType)                                                   
#==============================================================================
