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
import shutil
import subprocess
import numpy as np
from scipy import special
from lmfit import minimize, fit_report

from ParseCmdFile import parseCmdFile
from PlotFit import PlotFit
# from Parameters2 import Parameter2, Parameters2

import GlobalVariables as glbl
from read_data import read_data




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
        Resid = Resid + list(  DataSets[i][1] - TOF( DataSets[i][0], NDataSet, Params, AveragingType, ThetaAngles, ProbCurveType) )

    return np.array( Resid )
    

# -------------------------------------------------------------------------------------------------
#   TOF - compute the signal vs time
# -------------------------------------------------------------------------------------------------
def TOF(Time, NDataSet, Params, AveragingType, ThetaAngles, ProbCurveType, Debug=False):
    # Time of flight signal model. The function takes an uncorrected time in seconds and returns a signal
    #FFRDist  = Params['FFRDist_%i'   %NDataSet].value
    MaxTOF   = Params['Yscale_%i'    %NDataSet].value
    Baseline = Params['Baseline_%i'  %NDataSet].value
    TCutC    = Params['TCutC_%i'     %NDataSet].value
    TCutW    = Params['TCutW_%i'     %NDataSet].value
    TimeCorr = Params['IonTOF_%i'    %NDataSet].value
    #Temperature = Params['Temp_%i' %NDataSet].value
    Time   = Time - TimeCorr * 1E-6                            # Correct the time
    CutOff = 0.5 * (1. - np.tanh((Time - TCutC*1E-6) / (TCutW*1E-6)))      # CutOff function
            #used to model the experimental loss of low energy ions
    Signal0 = AngularAveraging( AveragingType, Time, NDataSet, Params, ThetaAngles)       
    Signal = Signal0 * CutOff * MaxTOF + Baseline
    
    if Debug :
        print('TOF     : Signal0 = ', Signal0)
        print('        : CutOff  = ', CutOff)
        print('        : MaxTOF  = ', MaxTOF) 
        print('        : Signal  = ', Signal)

    return Signal
    

# -------------------------------------------------------------------------------------------------
#   Angular Averaging
# -------------------------------------------------------------------------------------------------
def AngularAveraging( AveragingType, Time, NDataSet, Params, ThetaAngles):

    FFRDist  = Params['FFR_%i'   %NDataSet].value
    Temperature = Params['Temp_%i' %NDataSet].value

    Signal = 0.  # Initialize Signal to 0
    

    if AveragingType == "PointDetector":
        # Averaging performed  taking into account different flight time for different angles, but assuming ionization occurring in one point
        for Theta in ThetaAngles :  
            Velocity = FFRDist /(Time * np.cos( np.radians(Theta) ) ) # v = x / t = ( L / cos(theta) ) / t
            Ekin = (0.5 * glbl.MassAmu * Velocity**2.) * glbl.eVConst
            Enorm = Ekin * np.cos( np.radians(Theta) )**2. # Reaction probability depends on normal energy
            Signal = Signal + (Velocity**4. * np.exp( -Ekin / (glbl.kb * Temperature) ) * np.cos( np.radians(Theta) )**2. * Prob(Enorm, NDataSet, Params, ProbCurveType)) * np.sin( np.radians(Theta) ) * glbl.ThetaStep
    elif AveragingType == "None":
        # No angular averaging performed
        Velocity = FFRDist / Time # v = L / t 
        Ekin = (0.5 * glbl.MassAmu * Velocity**2.) * glbl.eVConst
        Enorm = Ekin
        Signal = (Velocity**4. * np.exp( -Ekin / (glbl.kb * Temperature) ) * Prob(Enorm, NDataSet, Params, ProbCurveType))

    elif AveragingType == "LineDetector":
        # Averaging along line, accounting for different flight time for different angles
        for Theta in ThetaAngles :
                        Velocity = FFRDist/(Time * np.cos( np.radians(Theta) ) ) # v = x / t = ( L / cos(theta) ) / t
                        Ekin = (0.5 * glbl.MassAmu * Velocity**2.) * glbl.eVConst
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
            Ekin = (0.5 * glbl.MassAmu * Velocity**2.) * glbl.eVConst
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
        Ekin = (0.5 * glbl.MassAmu * Velocity**2.) * glbl.eVConst
        # Enorm = Ekin
        VelocityDistribution = Velocity**4. * np.exp( -Ekin / (glbl.kb * Temperature))
        
        # kludge to fix divide by zero runtime error        
        for ii in range(len(VelocityDistribution)):
            if VelocityDistribution[ii] == 0:
                VelocityDistribution[ii] = 1E-100
        

    ProbFromTOF = (Signal - Baseline ) / (MaxTOF * VelocityDistribution * CutOff)
    Energy = 0.5 * glbl.MassAmu * ( FFRDist / Time )**2. * glbl.eVConst
    

    return Energy, ProbFromTOF


#============================================================================== 
#   WriteTOFOutput - Write the TOF Data                                                        
#==============================================================================
#==============================================================================
# def WriteTOFOutput( TOFTime, TOFData, State, fit_range, NDataSet, Params, AveragingType, ThetaAngles, ProbCurveType):
#     # Prints the fitted curve. 
#     MinTime =  Params['IonTOF_' + str(NDataSet)].value*1.E-6 + 0.1E-6
#     MaxTime = MinTime + 30.E-6
#     TimeIncr = 0.01E-6
# 
#     Time = np.arange( MinTime, MaxTime + TimeIncr, TimeIncr)
# 
#     # Name Output files
#     result_filename ='fit' +_DataSet_" + str( NDataSet ) + "_" + State + "_" + "_TOF_DataFitted.dat"
#     TOFOutFile                 = Label + "_DataSet_" + str( NDataSet ) + "_" + State + "_" + "_TOF_Fit_"      + ProbCurveType + ".dat"      
# 
#     # First determine the maximum (as average of 5 points) in the experimental signal
#     NMax =  TOFData.argmax()
#     DataMax = ( TOFData[NMax-2] + TOFData[NMax-1] + TOFData[NMax] + TOFData[NMax+1] + TOFData[NMax+2] ) / 5.
# 
#     # Write file with intensity fitted (minus background) - both relative and absolute signal 
#     np.savetxt( FileNoBackground , np.column_stack((np.array( TOFTime ), np.array( TOFData )/DataMax, np.array( TOFData ))) )
#     # print()    
#     # print("WriteTOFOutput: TOF data fitted written to file ", FileNoBackground)       
# 
#     # Calculate TOF fit curve
#     TOFFitted = TOF(Time, NDataSet, Params, AveragingType, ThetaAngles, ProbCurveType)    
# 
#     # Write file with TOF fitted - both relative and absolute signal
#     np.savetxt( TOFOutFile, np.column_stack((np.array( Time ), np.array( TOFFitted )/DataMax, np.array( TOFFitted ))) )
#     # print("WriteTOFOutput: Fitted TOF spectrum written to file ", TOFOutFile)
#     return
#==============================================================================


#==============================================================================
# #--------------------------------------------------------------------------------------------------
# #   write_result                                                      
# #--------------------------------------------------------------------------------------------------
# # write_result(newFitNumber, NDataSet, cmdFile, fitResult, Time, Signal, State, fit_range,
# #                    AveragingType, ThetaAngles, ProbCurveType)
# def write_result(FitNumber, NDataSet, cmdFile, fitresult, Time, Signal, State, fit_range, 
#                  AveragingType, ThetaAngles, ProbCurveType):
# 
#     # FitNumber is a 0 padded string
#     result_file_name = 'Fit' + FitNumber + '_' + str(NDataSet) + '_' + State + '_'   \
#                         + ProbCurveType + '.fit_out'
#     #with open(result_file_name, 'w') as result_file:
#     # fit_curve = TOF(Time, NDataSet, Params, AveragingType, ThetaAngles, ProbCurveType)
#     np.savetxt(result_file_name , np.column_stack((np.array(Time), np.array(Signal))))
#     
#     return 
#==============================================================================
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
# Main PROGRAM
#============================================================================= 

# define default path to control files and default command file name
pathToFits = 'Fits\\'
#cmdFilename = 'fit000.tof_in'


# Get Last fit number from FitNumber.dat, increment and write back
fit_number_file = open(pathToFits + 'FitNumber.dat', 'r+')
fitNumber = '{:03d}'.format(int(fit_number_file.readline()))
oldFitNumber = fitNumber
newFitNumber = fitNumber

#==============================================================================
# while True:
#     print('please enter oldfit number: ', '[', oldFitNumber, ']')
#     ans = input('?')
#     if ans:
#         old_n = '{:03d}'.format(int(ans))
#     else:
#         old_n = oldFitNumber
#     
#     if int(old_n) > int(oldFitNumber):
#         print('maximum allowed for old fit number is ', oldFitNumber)
#     else:
#         break
#     
# oldFitNumber = old_n
# 
# ans = input ('make new command file? [no]')
# if ans:
#     if ans.upper()[0] == 'Y':
#         newFitNumber = '{:03d}'.format(int(fitNumber)+1)
#     fit_number_file.seek(0)
#     fit_number_file.write(newFitNumber)
# else:
#     newFitNumber = oldFitNumber
#==============================================================================

fit_number_file.close()

oldFile = pathToFits + 'fit' + oldFitNumber + '.tof_in'
newFile = oldFile

if oldFitNumber != newFitNumber:
    newFile = pathToFits + 'fit' + newFitNumber + '.tof_in'
    shutil.copy2(oldFile, newFile)

# subprocess.call(['npp.bat', newFile])
cmdFile = pathToFits + 'fit' + newFitNumber + '.tof_in'

# Parse the command file
parms, functions, signalFiles, backgroundFiles, errors = parseCmdFile(cmdFile)

if len(errors) > 0:
    print('Errors\n', errors)
    raise SystemExit ('Error in command file' )
# print('Finished parsing command file')

DataFiles = signalFiles
BackgroundFiles = backgroundFiles
# Label = glbl.Label

angularaveragingList = ['None','PointDetector','LineDetector']
args_angularaveraging = angularaveragingList[0]

# args_input = r'data\referenceline\6008_Au_D2_v1J2_x=35.datv2'
# args_background = r'data\referenceline\6009_Au_D2_v1J2_x=35_offRes.datv2'
args_mintime = 5
args_maxtime = 25

    
# DataFiles = args_input.split(",")
AveragingType = args_angularaveraging

# initialize variables and parametrs
Tmin = glbl.Tmin
Tmax = glbl.Tmax
states = []         # list if tuples (v, j)
DataSets = []
PlotDataSets = []
Params = parms

# Read data from input
fit_ranges = []
for i in range( len( DataFiles)):
    ProbCurveType = functions[0]
    DataFile = DataFiles[i]
    if BackgroundFiles  != [ "" ]:
        BackgroundFile = BackgroundFiles[i]
    else:
        BackgroundFile = BackgroundFiles[0]
        
    state, Temperature, Nmin, Nmax, Time, Signal = read_data(DataFile, BackgroundFile, Tmin, Tmax)    
    states.append( state )
    DataSets.append([Time[Nmin:Nmax], Signal[Nmin:Nmax]])
    PlotDataSets.append([Time, Signal])
    fit_ranges.append([Nmin,Nmax])
  

# Generate Theta angles employed for angular averaging
ThetaAngles = GenerateThetaAngles(                                  \
    AveragingType=AveragingType, GridType=glbl.GridType,                 \
    NPointsSource=glbl.NPointsSource, NPointsDetector=glbl.NPointsDetector,   \
    ZSource = glbl.ZSource,           RSource = glbl.RSource,                 \
    ZAperture = glbl.ZAperture,       RAperture = glbl.RAperture,             \
    ZDetector = glbl.ZLaser,          LengthDetector = glbl.LLaser)
    #    ZDetector = ZFinal,          LengthDetector = 2.*RFinal         \
  
#==============================================================================
# print()
# print('Angular Averaging Parameters:')
# print('nps =', glbl.NPointsSource, 'npd =', glbl.NPointsDetector)
# print('zs  =', glbl.ZSource,       'rs  =', glbl.RSource )
# print('za  =', glbl.ZAperture,     'ra  =', glbl.RAperture)     
# print('zf  =', glbl.ZFinal,        'rf  =', glbl.RFinal)                                                
#==============================================================================

# Fit TOF
fitResult = FitData( DataSets, Params, AveragingType, ThetaAngles, ProbCurveType)
print(fit_report(fitResult))

state_string = 'v' + str(states[0][0]) + 'j' +str(states[0][1])

result_file_name = 'Fit' + newFitNumber + '_' + state_string + '_' + ProbCurveType + '.fit_out'
with open(pathToFits + result_file_name, 'w') as result_file:
    
    result_file.write('#\n')
    result_file.write('# Fit {} Results\n'.format(newFitNumber))
    result_file.write('#\n')
    result_file.write('#' + 40*'-' + '\n')
    result_file.write('# Control file: {}\n'.format(cmdFile))
    result_file.write('#' + 40*'-' + '\n')
    
    with open(cmdFile,'r') as cmd:
        cmd_lines = cmd.readlines()

    for line in cmd_lines:
        result_file.write('# ' + line )
    
    #result_file.write('#\n')
    result_file.write('\n#' + 40*'-' + '\n')    
    result_file.write('# End Control File\n')
    result_file.write('#' + 40*'-' + '\n')
    result_file.write('#\n')
    
    result_file.write('#' + 40*'-' + '\n')    
    result_file.write('# Fit Report\n')
    result_file.write('#' + 40*'-' + '\n')    
    result_file.write('#\n')
    
    report = fit_report(fitResult).split('\n')
    for line in report:
        result_file.write('# ' + line + '\n')
    
    result_file.write('#' + 40*'-' + '\n')    
    result_file.write('# End Fit Report\n')
    result_file.write('#' + 40*'-' + '\n')
    
    for i in range( len( DataSets )):
        # Write output
        NDataSet =  i + 1
        Time     =  PlotDataSets[i][0] * 1E6  # convert to microseconds
        Signal   =  PlotDataSets[i][1]
        state    =  states[i]
        fit_range=  fit_ranges[i]
        
        result_file.write('#\n')
        result_file.write('#' + 40*'-' + '\n')    
        result_file.write('# Begin Data Set {}\n'.format(NDataSet))
        result_file.write('#     States: v ={} J={}\n'.format(states[0][0], states[0][1]))
        result_file.write('#' + 40*'-' + '\n')
        # for j in range(0,len(Time)):
        for j in range(0,20):
            result_file.write('{:6.2f} {:20.5e}\n'.format(Time[j], Signal[j]))
      
#==============================================================================
#     write_result(newFitNumber, NDataSet, cmdFile, fitResult, Time, Signal, State, fit_range,
#                    AveragingType, ThetaAngles, ProbCurveType)
#==============================================================================

    # baseFileName = Label + "_DataSet_" + str( NDataSet ) + "_" + State
    # PlotFit(baseFileName)
    # if DataType != 'Calibration':
#==============================================================================
#     if functions[i] != 'Calibration':
#         WriteProbOutput(Time, Signal, State, Label, NDataSet, fitResult.params, 
#                         AveragingType, ThetaAngles, ProbCurveType)                                                   
#==============================================================================
