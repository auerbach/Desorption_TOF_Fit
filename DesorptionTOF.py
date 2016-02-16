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

import os
import shutil

# import sys
import subprocess
import numpy as np
from scipy import special
from lmfit import minimize, fit_report

from ParseCmdFile import parseCmdFile
from PlotFit import PlotFit
# from Parameters2 import Parameter2, Parameters2

import GlobalVariables as glbl


#------------------------------------------------------------------------------
# read_data  function to read the data and subtract background
#------------------------------------------------------------------------------
def read_data(data_file_name, background_file_name = '', Tmin ='', Tmax ='', Threshold = 0.10):
    # Function to read Datafiles
    # Input: DataFile = name of the file to be read ;
    # BackgroundFile  = name of the file with background; 
    # Tmin, Tmax = minimum (maximum) value of time (in us) to consider ;
    # Threshold = ratio signal/maximum of signal for signal to be considered


    # Check existence of the file
    if not os.path.isfile(DataFile):
        print(data_file_name + " : Data file does not exist! Quitting...")
        raise SystemExit
    
    # Open and read the data file    
    with open(data_file_name, 'r') as data_file:
        lines = data_file.readlines()
    
    DataFormat = lines[glbl.DataFormatLine -1].split()[3]
    
    if DataFormat != '2.1':
        print('File ', data_file_name, ' is not in the right format')
        print('Format = ', DataFormat)
        raise SystemExit
        
    glbl.OriginalFile   = lines[glbl.OriginalFileLine       - 1].split()[3]
    glbl.FileDate       = lines[glbl.FileDateLine           - 1].split()[3]
    glbl.Molecule       = lines[glbl.MoleculeLine           - 1].split()[3]
    Temperature         = float( lines[glbl.TemperatureLine - 1].split()[3])    
    glbl.VibState       = float( lines[glbl.VibStateLine    - 1].split()[3])
    glbl.RotState       = float( lines[glbl.RotStateLine    - 1].split()[3])
    DataCol             = int(lines[glbl.DataColLine        - 1].split()[3])
    DataRow             = int(lines[glbl.DataRowLine        - 1].split()[3])
    
    State = "v" + str( int(glbl.VibState)) + "J" + str( int(glbl.RotState))
    
    if glbl.Molecule  == 'H2':
        glbl.MassAmu = glbl.massH2
    elif glbl.Molecule == 'HD':
        glbl.MassAmu = glbl.massHD
    elif glbl.Molecule == 'D2':
        glbl.MassAmu = glbl.massD2
    else:
        print('Error:  Unknown Molecule ', glbl.Molecule)
        raise SystemExit
    
    State = "v" + str( int(glbl.VibState)) + "J" + str( int(glbl.RotState))

    # If Background file is provided, subtract it from intensities
    if background_file_name :
        if not os.path.isfile(background_file_name):
            print(background_file_name + " : Background file does not exist! Quitting...")
            quit()

        with open(background_file_name, 'r') as background_file:
            BackgroundLines = background_file.readlines()

    # Assuming data last until end of the file
    DataRowEnd = len( lines )

    # Read file
    Time = []
    Data = []
       
    for n in range(DataRow - 1, DataRowEnd):
        # print(n,lines[n])
        T = float( lines[n].split()[0] )*1E-6   # T is ins sec , Gottingen data                                                # is in microsec
        F = float( lines[n].split()[DataCol -1] ) #
        Time.append( T )
        if not BackgroundFile :
            Data.append( F )
        else:
            B = float( BackgroundLines[n].split()[DataCol -1] )
            Data.append( F - B )

    DeltaTime = Time[1] - Time[0]
    Time = np.array(Time)
    Data = np.array(Data)

    #==============================================================================================
    # Select good data: if max or min times are provided use them, otherwise find times where
    #     data is = max of data * Threshold
    # Since the data is noisy, for point n average from n-n_delt to n+n_delt
    #==============================================================================================
    n_delt = 5
    data_max = np.array(Data[100:len(Data)]).max()
        
    # find Nmin = index of first point to use in fitting
    for n in range(100, len( Time )-n_delt):
        if Tmin :
            if float( Tmin ) * 1.E-6 >= Time[n] and float( Tmin )* 1.E-6 < Time[n] + DeltaTime:
                Nmin = n
                break
        else:
            # find time where Data[n] reaches threshold and subtract 1 microseconds
            if np.array(Data[n-n_delt:n + n_delt]).mean() >= data_max * Threshold:
                Nmin = int(n - 1 / (DeltaTime * 1E6))
                break
    
    for n in range( len( Time ) -(n_delt+1), 100, -1 ):
        if Tmax :
            if float( Tmax )* 1.E-6 >= Time[n] and float( Tmax ) * 1.E-6 < Time[n] + DeltaTime:
                Nmax = n
                break
        else:
            # find time where Data[n] falls to threshold and add 1.5 microseconds
            if np.array(Data[n-n_delt:n + n_delt]).mean()  >= data_max * Threshold:
                Nmax = int(n + 1.5 / (DeltaTime * 1E6))
                break
                
    print()
    print('Range of points to use in fitting')                
    print('Nmin=',Nmin,  ' Tmin=',Time[Nmin])
    print('Nmax=',Nmax, ' Tmax=',Time[Nmax])
    print('DataMax=',max( Data[50:len(Data)] ))
    print()

    return State, Temperature, Nmin, Nmax, Time, Data


    
#==============================================================================
#     Params.add('FFRDist_%i'     %(i+1), value=FFRDist,
#                                         min=0.0,
#                                         #min=FFRDist - FFRDist * FFRDistTolerance,
#                                         #max=FFRDist + FFRDist * FFRDistTolerance,
#                                         vary=False)
#     Params.add('TimeCorr_%i'    %(i+1), value=TimeCorr,
#                                         #min=TimeCorr-TimeCorr * TimeCorrTolerance,
#                                         #max=TimeCorr+TimeCorr * TimeCorrTolerance,
#                                         vary=False)
#     Params.add('TCutC_%i'       %(i+1), value=TCutC,
#                                         min=TCutC - TCutC * TCutCTolerance,
#                                         max=TCutC + TCutC * TCutCTolerance,
#                                         vary=False)
#     Params.add('TCutW_%i'       %(i+1), value=TCutW,
#                                         min=TCutW - TCutW * TCutWTolerance,
#                                         max=TCutW + TCutW * TCutWTolerance,
#                                         vary=False)
#     Params.add('Temp_%i' %(i+1), value=Temperature,
#                                        #min=Temperature - TemperatureTolerance * Temperature,
#                                        #max=Temperature + TemperatureTolerance * Temperature,
#                                        vary=False)
#==============================================================================
    

    return Params

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


###############################################################################
def FitData( DataSets, Params, AveragingType, ThetaAngles, ProbCurveType, LabelData):
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
#==============================================================================
#     Params = result.params
#     print('FitData: after minimize')
#     print('Params[E0_1].value = ', Params['E0_1'].value)
#==============================================================================
    return result
#==============================================================================
#     new_stdout =open(Label + '.out', 'w')
#     old_stdout =
#     stdouts=[new_stdout, old_stdout]
#     
#     for stdout in stdouts:
#         sys.stdout=stdout
#==============================================================================
#==============================================================================
#     print()
#     print('Fitdata: ', 25*"#")
#     print(LabelData + " " + ProbCurveType + " curve fitting complete:")
#     print("Success:   ", Result.success)
#     print("chisq:     ", Result.chisqr)
#     print("dof:       ", Result.ndata, Result.nvarys, Result.nfree)
#     print("RedChiSq:  ", Result.redchi)
#     print("Optimized parameters:")
#     for name, par in list(Params.items()):
#         #if par.vary == True:
#             print(name+"  "+str(par.value) + " +/- " + str(par.stderr) +" Vary "+str(par.vary))
#     print('Fitdata: ', 25*"#")
#==============================================================================
    


def Residual( Params, X, Y, DataSets, AveragingType, ThetaAngles, ProbCurveType):
    # Residual function needed by minimize
    Resid = []

    for i in range( len( DataSets )):
        NDataSet = i + 1
        Resid = Resid + list(  DataSets[i][1] - TOF( DataSets[i][0], NDataSet, Params, AveragingType, ThetaAngles, ProbCurveType) )

    return np.array( Resid )
    

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
    
#==============================================================================
# def Debug1(Time=2E-5):
#     # Baseline    = Params['Baseline_1'   ].value
#     TCutC       = Params['TCutC_1'  ].value
#     TCutW       = Params['TCutW_1'  ].value
#     FFRDist     = Params['FFR_1'    ].value
#     MaxTOF      = Params['Yscale_1' ].value
#     TimeCorr    = Params['IonTOF_1' ].value
#     Temperature = Params['Temp_1'   ].value
#     
#     Time2    = Time - TimeCorr *1E-6
#     
#     Signal = TOF(Time,1,Params,'None',0,'Calibration')
#     CutOff = 0.5 * (1. - np.tanh((Time2 - TCutC*1E-6) / (TCutW*1E-6)))
#     print('TOF Sig  : ', Signal)
#     print('MaxTOF   : ', MaxTOF)
#     print('CutOff   : ', CutOff)
#     
#     Velocity = FFRDist / Time2 # v = L / t 
#     Ekin = (0.5 * MassAmu * Velocity**2.) * eVConst
#     print('Velocity : ', Velocity)
#     Signal0 = (Velocity**4. * np.exp( -Ekin / (kb * Temperature) ))
#     print('Signal0  : ', Signal0)
#     print('  *MaxTOF: ', Signal0*MaxTOF)
#     print('  *CutOff: ', Signal0*MaxTOF*CutOff)
#==============================================================================
    
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
def WriteTOFOutput( TOFTime, TOFData, State, Label, NDataSet, Params, AveragingType, ThetaAngles, ProbCurveType):
    # Prints the fitted curve. 
    MinTime =  Params['IonTOF_1'].value*1.E-6 + 0.1E-6
    MaxTime = MinTime + 30.E-6
    MaxTime = 50.E-6
    TimeIncr = 0.01E-6

    Time = np.arange( MinTime, MaxTime + TimeIncr, TimeIncr)

    # Name Output files
    FileNoBackground           = Label + "_DataSet_" + str( NDataSet ) + "_" + State + "_" + "_TOF_DataFitted.dat"
    TOFOutFile                 = Label + "_DataSet_" + str( NDataSet ) + "_" + State + "_" + "_TOF_Fit_"      + ProbCurveType + ".dat"      

    # First determine the maximum (as average of 5 points) in the experimental signal
    NMax =  TOFData.argmax()
    DataMax = ( TOFData[NMax-2] + TOFData[NMax-1] + TOFData[NMax] + TOFData[NMax+1] + TOFData[NMax+2] ) / 5.

    # Write file with intensity fitted (minus background) - both relative and absolute signal 
    np.savetxt( FileNoBackground , np.column_stack((np.array( TOFTime ), np.array( TOFData )/DataMax, np.array( TOFData ))) )
    # print()    
    # print("WriteTOFOutput: TOF data fitted written to file ", FileNoBackground)       

    # Calculate TOF fit curve
    TOFFitted = TOF(Time, NDataSet, Params, AveragingType, ThetaAngles, ProbCurveType)    

    # Write file with TOF fitted - both relative and absolute signal
    np.savetxt( TOFOutFile, np.column_stack((np.array( Time ), np.array( TOFFitted )/DataMax, np.array( TOFFitted ))) )
    # print("WriteTOFOutput: Fitted TOF spectrum written to file ", TOFOutFile)
    return


#============================================================================= 
#   WriteProbOutput - Write the Fitted S0                                                        
#==============================================================================
def WriteProbOutput( TOFTime, TOFData, State, Label, NDataSet, Params, AveragingType, ThetaAngles, ProbCurveType):
        # Name Output files
    ReactionProbabilityOutFile = "DataSet_" + str( NDataSet ) + "_" + State + "_" + Label + "_ReacProb_Fit_" + ProbCurveType + ".dat"
    TOFInvertedOutFile         = "DataSet_" + str( NDataSet ) + "_" + State + "_" + Label + "_ReacProb_Fit_" + ProbCurveType + "_FromTOFInversion.dat"

    # Reaction probability curve: minimum energy, maximum energy, deltaE (eV)
    Emin = 0.0
    Emax = 2.0
    DeltaE = 0.01

    Energy = np.arange( Emin, Emax + DeltaE, DeltaE)

    # Calculate Reaction probability fitted curve
    ReactionProbability = Prob(Energy, NDataSet, Params, ProbCurveType)
    np.savetxt( ReactionProbabilityOutFile, np.column_stack(( Energy, ReactionProbability)))
    # print("WriteProbOutput: Reaction probability curve written to file ", ReactionProbabilityOutFile)

    # Reaction probability from TOF inversion; only possible with Point detector or no angular averaging (and normal energy scaling!)
    if AveragingType != "LineDetector":
        TOFEnergy, TOFReactionProbability = ProbFromTOFInversion(TOFTime, TOFData, NDataSet, Params, AveragingType, ThetaAngles, ProbCurveType)
        np.savetxt( TOFInvertedOutFile, np.column_stack(( TOFEnergy, TOFReactionProbability*Params['Yscale_%i' %NDataSet].value, Prob( TOFEnergy, NDataSet, Params, ProbCurveType)*Params['Yscale_%i' %NDataSet].value )))
        # print("WriteProbOutput: Inverted TOF written to file ", TOFInvertedOutFile)
        # print('WriteProbOutput: Exiting')
    return                                               


#============================================================================= 
# Main PROGRAM
#============================================================================= 

# define default path to control files and default command file name
pathToFits = 'Fits\\'
#cmdFilename = 'fit000.tof_in'


# Get Last fit number from FitNumber.dat, increment and write back
fit_number_file = open(pathToFits + 'FitNumber.dat', 'r+')
oldFitNumber = int(fit_number_file.readline())
oldFitNumber = '{:03d}'.format(oldFitNumber)

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
        newFitNumber = '{:03d}'.format(int(oldFitNumber)+1)
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

subprocess.call(['npp.bat', newFile])
cmdFile = pathToFits + 'fit' + newFitNumber + '.tof_in'

# Parse the command file
parms, functions, signalFiles, backgroundFiles, errors = parseCmdFile(cmdFile)

if len(errors) > 0:
    print('Errors\n', errors)
    raise SystemExit ('Error in command file' )
# print('Finished parsing command file')

DataFiles = signalFiles
BackgroundFiles = backgroundFiles
Label = glbl.Label

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
States = []
DataSets = []
PlotDataSets = []
Params = parms

# Read data from input
for i in range( len( DataFiles)):
    ProbCurveType = functions[0]
    DataFile = DataFiles[i]
    if BackgroundFiles  != [ "" ]:
        BackgroundFile = BackgroundFiles[i]
    else:
        BackgroundFile = BackgroundFiles[0]
        
    State, Temperature, Nmin, Nmax, Time, Signal = read_data(DataFile, BackgroundFile, Tmin, Tmax)    
    States.append( State )
    DataSets.append([Time[Nmin:Nmax], Signal[Nmin:Nmax]])
    PlotDataSets.append([Time, Signal])
  
    # Add default parameters to the dictionary
    # Params.add('Baseline_%i'    %(i+1),  value=0.  , vary=False)    # TOF spectrum Baseline Parameter
  
    


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
fitResult = FitData( DataSets, Params, AveragingType, ThetaAngles, ProbCurveType, Label )
print(fit_report(fitResult))

for i in range( len( DataSets )):
    # Write output
    NDataSet =  i + 1
    Time     =  DataSets[i][0]
    Signal   =  DataSets[i][1]
    State    =  States[i]
      
    WriteTOFOutput( Time, Signal, State, Label, NDataSet, fitResult.params, 
                   AveragingType, ThetaAngles, ProbCurveType)

    baseFileName = Label + "_DataSet_" + str( NDataSet ) + "_" + State
    PlotFit(baseFileName)
    # if DataType != 'Calibration':
    if functions[i] != 'Calibration':
        WriteProbOutput(Time, Signal, State, Label, NDataSet, fitResult.params, 
                        AveragingType, ThetaAngles, ProbCurveType)                                                   
