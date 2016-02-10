#!/usr/bin/env python

###########################################################################
## Time of Flight fitting program
## 
## Written by Alessandro Genova
## July 2012, Leiden University
##
## Fit experimental TOF spectra using several sigmoid functions as a model (and pront the fitted curve to file).
##
## USAGE:
##
## Start the program with the following command and read the provided instructions:
## ./TOF_v5.py --help
##
##
###########################################################################
## Modified by Francesco Nattino
## Last update 06-2014
##
## Changes from vs6
## - various angular averaging included (none and line detector)
## - global fit for calibration implemented
##
###########################################################################

# Importing the extra packages and modules needed

import os
# import sys
# import argparse
import numpy as np
# import scipy as sp
# from scipy.optimize import leastsq
from scipy import special
# from scipy import interpolate
# from lmfit import minimize, Parameters, conf_interval, printfuncs, fit_report
from lmfit import minimize, fit_report

# from lmfit import minimize

from ParseCmdFile import parseCmdFile
from PlotFit import PlotFit
# from Parameters2 import Parameter2, Parameters2

# Physical constants ###################################################

kb          = 8.6173324E-5           # Boltzmann constant in eV/K
eV2J        = 1.602176565E-19           # eV to Joule 
J2eV        = 6.24150934E18           # Joule to eV 
AtomicMass  = 1.660538921E-27           # Atomic mass constant
eVConst     = AtomicMass * J2eV
MassAmu     = 2. * 2.01410178           # Mass of molecule (amu)

# Experimental apparatus constants #####################################

FFRDist           = 29.0E-3        # Distance travelled by the molecule in the fieldfree region as an ion (m)
FFRDistTolerance  = 0.5            # Percentage by which the field free region lenght can vary
TimeCorr          = 4.6            # Time correction (us)
TimeCorrTolerance = 0.8            # (see two lines above )
TCutC             = 28.6           # CutOff function 1st parameter (us) 
TCutCTolerance    = 0.5            # ...
TCutW             = 4.3            # CutOff function 2nd parameter (us)
TCutWTolerance    = 1.0

TemperatureTolerance = 1.         

# Following parameters for point detector
AngRes               = 20.         # Angular resolusion (degrees)
ThetaStep            = 2.          # Theta step in averaging

# Following parameters for point detector
AngRes               = 20.         # Angular resolusion (degrees)
ThetaStep            = 2.          # Theta step in averaging

# Following parameters for line detector
ZDiffWall   = 0.0               # Poistion of inside (detector side)
RDiffWall   = 6.0               # Radius of the hole in differential wall        
ZRef = ZDiffWall                # Use Diff wall as reference point since
                                #   Source and differential wall positions
                                #   might be changed
ZAperture  = ZRef - 3.5        # Position of Aperture in Ta Sheiled   
RAperture   = 1.5               # Radius of aperture in Ta shield
ZSource     = ZAperture - 4.   # Position of Source
RSource     = 0.1               # Radius of the source (Source or knudsen)
    
# Detection volume is determined by length of REMPI laser.
# the actual length is very long, so we should set this parameter to be 
# long enough that it does not limit detection.  
ZLaser      = ZRef + 5.0        # Position of REMP laser beam
LLaser      = 9.3               # Length of REMPI detection volume.  

# In the present code we don't us these parameter.  Instead
# we use an effective detection line based on the acceptance
# of ions at the final field free region, i.e we assume the length of the
# REMPI detection volume is not a limiting factor
       
ZFinal     = ZRef + 34.         # Position of the final grid
RFinal     = 10.0               # Effective acceptance radius for ions at
                                #   at the final grid.  Because of strong
                                #   accleration to extractor take this to
                                #   be equal to the extrator radius

NPointsDetector = 101          # Number of points to consider on the line
    
NPointsSource   = 1             # Number of points to consider on the Source
                                #   If NPointsSource = 1 we treat this as
                                #   point source

GridType = 'Cartesian'          # Generate a cartesian or radial grid on  
                                #   the source. This parameter can have
                                #   values of 'Cartesian' or'Radial'
# Data Format ##########################################################

DataLine         = 34    # Line in DataFile where data start
MassLine         = 1
TemperatureLine  = 2     # Line in DataFile where temperature is reported
                         #  (surface T for desorption experiments,
                         #    nozzle T for Knudsen experiments) 
VibStateLine     = 3
RotStateLine     = 4

########################################################################

def ReadData( DataFile, BackgroundFile = "", Tmin = '', Tmax = '', Threshold = 0.10):
    # Function to read Datafiles
    # Input: DataFile = name of the file to be read ;
    # BackgroundFile  = name of the file with background; 
    # Tmin, Tmax = minimum (maximum) value of time (in us) to consider ;
    # Threshold = ratio signal/maximum of signal for signal to be considered


    # Check existence of the file
    #    (only if Tmin, Tmax are not defined).
    if not os.path.isfile(DataFile):
        print(DataFile + " : Data file does not exist! Quitting...")
        quit()
       
    # Open input file
    lines = open( DataFile ).readlines()
    
    # Read temperature, vib, rot states from head


    Temperature = float( lines[TemperatureLine - 1].split()[3 ] )
    # Mass        = float( lines[MassLine        - 1].split()[3 ] )
    VibState    = float( lines[VibStateLine    - 1].split()[3 ] )
    RotState    = float( lines[RotStateLine    - 1].split()[3 ] )
    # Read rovibrational state (if not in input):
    State = "v" + str( int(VibState)) + "J" + str( int(RotState))

    # If Background file is provided, subtract it from intensities
    if BackgroundFile :
        if not os.path.isfile(BackgroundFile):
            print(BackgroundFile + " : Background file does not exist! Quitting...")
            quit()

        BackgroundLines = open( BackgroundFile ).readlines()

    # Assuming data last until end of the file
    DataLineEnd = len( lines )

    # Read file
    Time = []
    Data = []
       
    for n in range(DataLine - 1, DataLineEnd):
        # print(n,lines[n])
        T = float( lines[n].split()[0] )*1E-6   # T is ins sec , Gottingen data                                                # is in microsec
        F = float( lines[n].split()[1] ) #
        Time.append( T )
        if not BackgroundFile :
            Data.append( F )
        else:
            B = float( BackgroundLines[n].split()[1] )
            Data.append( F - B )

    DeltaTime = Time[1] - Time[0]

    # Select good data: if max or min times are provided use them, otherwise use threshold (ratio signal/maximum signal)
    for n in range( len( Time )):
        if Tmin :
            if float( Tmin ) * 1.E-6 >= Time[n] and float( Tmin )* 1.E-6 < Time[n] + DeltaTime:
                Nin = n
                break
        else:
            #if n != 0 and n != ( len( Time ) -1 ):
            if n >= 100 and n != ( len( Time ) -1 ):
                if ( Data[n] + Data[n-1] + Data[n+1] ) / 3.  >= max( Data[50:len(Data)] ) * Threshold:
                    Nin = n
                    break
    print('Readdata: Nin=',Nin, ' Time[Nin]=',Time[Nin])
    print('Readdata: DataMax=',max( Data[50:len(Data)] ))
    for n in range( len( Time ) -1, -1, -1 ):
        if Tmax :
            if float( Tmax )* 1.E-6 >= Time[n] and float( Tmax ) * 1.E-6 < Time[n] + DeltaTime:
                Nfin = n
                break
        else:
            if n >= 50 and n != ( len( Time ) -1 ):
                if ( Data[n] + Data[n-1] + Data[n+1] ) / 3.  >= max( Data[50:len(Data)] ) * Threshold:
                    Nfin = n
                    break

    # Set arrays for fitting: Time (X) in sec, data (Y) in arbitrary units    
    VarX = np.array( Time[Nin:Nfin] )
    VarY = np.array( Data[Nin:Nfin] )

    return State, Temperature, VarX, VarY

def SetupParametersDictionary( i, Params,  Temperature):
    
    
#==============================================================================
#     # Reaction probability parameters
#     Params.add('A_%i' %(i+1),  value=1.0 , vary=True, max=1. )           # Saturation value
#     Params.add('B_%i' %(i+1),  value=0.6 , vary=False, min=0., max=5.)    # ReacProb Curve E0 Parameter
#     Params.add('BI_%i'%(i+1), value=0.5 , vary=False )                   # ReacProb Curve E0 Parameter (for FPC)
#     Params.add('C_%i' %(i+1),  value=0.16, vary=False, min=0., max=1.)    # ReacProb Curve Width Parameter
#     Params.add('CI_%i'%(i+1), value=0.90, vary=False, min=0.)            # ReacProb Curve Width Parameter (for FPC)
#     Params.add('ni_%i'%(i+1), value=0.05, vary=False, min=1E-7, max=1.)  # ReacProb Curve 'ni' parameter  (for LGS)
#==============================================================================
    

    # Experimental apparatus
    #Params.add('MaxTOF_%i'      %(i+1),  value=1.E-14 , vary=True)    # Max of the TOF spectrum
    Params.add('Baseline_%i'    %(i+1),  value=0.  , vary=False)    # TOF spectrum Baseline Parameter
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
        ThetaAngles = np.arange( 0., AngRes + ThetaStep, ThetaStep )
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
            Ekin = (0.5 * MassAmu * Velocity**2.) * eVConst
            Enorm = Ekin * np.cos( np.radians(Theta) )**2. # Reaction probability depends on normal energy
            Signal = Signal + (Velocity**4. * np.exp( -Ekin / (kb * Temperature) ) * np.cos( np.radians(Theta) )**2. * Prob(Enorm, NDataSet, Params, ProbCurveType)) * np.sin( np.radians(Theta) ) * ThetaStep
    elif AveragingType == "None":
        # No angular averaging performed
        Velocity = FFRDist / Time # v = L / t 
        Ekin = (0.5 * MassAmu * Velocity**2.) * eVConst
        Enorm = Ekin
        Signal = (Velocity**4. * np.exp( -Ekin / (kb * Temperature) ) * Prob(Enorm, NDataSet, Params, ProbCurveType))

    elif AveragingType == "LineDetector":
        # Averaging along line, accounting for different flight time for different angles
        for Theta in ThetaAngles :
                        Velocity = FFRDist/(Time * np.cos( np.radians(Theta) ) ) # v = x / t = ( L / cos(theta) ) / t
                        Ekin = (0.5 * MassAmu * Velocity**2.) * eVConst
                        Enorm = Ekin * np.cos( np.radians(Theta) )**2 # Reaction probability depends on normal energy
                        # Here no sin weight, since each Theta value has a weight of one
                        Signal = Signal + (Velocity**4. * np.exp( -Ekin / (kb * Temperature) ) * np.cos( np.radians(Theta) )**2. * Prob(Enorm, NDataSet, Params, ProbCurveType)) * ThetaStep

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
            Ekin = (0.5 * MassAmu * Velocity**2.) * eVConst
            # Enorm = Ekin * np.cos( np.radians(Theta) )**2. # Reaction probability depends on normal energy
            VelocityDistribution = VelocityDistribution                   \
                                   + Velocity**4.                         \
                                   * np.exp( -Ekin / (kb * Temperature) ) \
                                   * np.cos( np.radians(Theta) )**2.      \
                                   * np.sin( np.radians(Theta) )           \
                                   * ThetaStep

    elif AveragingType == "None":
        # No angular averaging
        Velocity = FFRDist / Time
        Ekin = (0.5 * MassAmu * Velocity**2.) * eVConst
        # Enorm = Ekin
        VelocityDistribution = Velocity**4. * np.exp( -Ekin / (kb * Temperature))
        
        # kludge to fix divide by zero runtime error        
        for ii in range(len(VelocityDistribution)):
            if VelocityDistribution[ii] == 0:
                VelocityDistribution[ii] = 1E-100
        

    ProbFromTOF = (Signal - Baseline ) / (MaxTOF * VelocityDistribution * CutOff)
    Energy = 0.5 * MassAmu * ( FFRDist / Time )**2. * eVConst
    

    return Energy, ProbFromTOF

#==============================================================================
# def SetFreeParameters( Params, DataType, ProbCurveType="ERF", FitTemperature=False, FitCutoffFunction=False, FitCorrectionTime=False, FitBaseline=False ):
#     # Set the parameters free to vary in the fit
#     
# 
#     # Loop over all the parameters
#     for name, param in list(Params.items()):
#         # Scaling constants
#         if name.startswith('MaxTOF_'):
#             param.vary = True
#     
#         if DataType == "FitTOF":
#             # Parameters reaction probability
#             if name.startswith('B_') or name.startswith('C_'):
#                 param.vary = True
#             if ProbCurveType == "LGS":
#                 if name.startswith('ni_'): 
#                     param.vary = True
#             if ProbCurveType == "FPC":
#                 if name.startswith('BI_') or name.startswith('CI_'):
#                     param.vary = True
# 
#         elif DataType == "Calibration":
#             # Constrain parameters of all datasets to be equal
# 
#             # Field free region distance
#             if name.startswith('FFRDist_'):
#                 param.vary = True
#                 if not name.endswith('_1'):       
#                     param.expr = 'FFRDist_1'   
#                 
#             # Correction time
#             if FitCorrectionTime:
#                 if name.startswith('TimeCorr_'): 
#                     param.vary = True
#                     if not name.endswith('_1'):       
#                         param.expr = 'TimeCorr_1'
# 
#             # Cutoff function
#             if FitCutoffFunction :
#                 if name.startswith('TCutC_'): 
#                     param.vary = True
#                     if not name.endswith('_1'):      
#                         param.expr = 'TCutC_1'
#                 if name.startswith('TCutW_'):
#                     param.vary = True
#                     if not name.endswith('_1'):     
#                         param.expr = 'TCutW_1'
# 
#         # Optionally, we can fit baseline or temperature
#         if FitBaseline:
#             if name.startswith('Baseline_'):
#                 param.vary = True
#         if FitTemperature:
#             if name.startswith('Temp_'):
#                 param.vary = True
#
#     return
#==============================================================================
                                       

def LockParam(Params):
    for name, param in list(Params.items()):
        if param.vary == True :
            param.vary = False
    return
                                                       

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
    print()    
    print("WriteTOFOutput: TOF data fitted written to file ", FileNoBackground)       

    # Calculate TOF fit curve
    TOFFitted = TOF(Time, NDataSet, Params, AveragingType, ThetaAngles, ProbCurveType)    

    # Write file with TOF fitted - both relative and absolute signal
    np.savetxt( TOFOutFile, np.column_stack((np.array( Time ), np.array( TOFFitted )/DataMax, np.array( TOFFitted ))) )
    print("WriteTOFOutput: Fitted TOF spectrum written to file ", TOFOutFile)
    return

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
    print("WriteProbOutput: Reaction probability curve written to file ", ReactionProbabilityOutFile)

    # Reaction probability from TOF inversion; only possible with Point detector or no angular averaging (and normal energy scaling!)
    if AveragingType != "LineDetector":
        TOFEnergy, TOFReactionProbability = ProbFromTOFInversion(TOFTime, TOFData, NDataSet, Params, AveragingType, ThetaAngles, ProbCurveType)
        np.savetxt( TOFInvertedOutFile, np.column_stack(( TOFEnergy, TOFReactionProbability*Params['Yscale_%i' %NDataSet].value, Prob( TOFEnergy, NDataSet, Params, ProbCurveType)*Params['Yscale_%i' %NDataSet].value )))
        print("WriteProbOutput: Inverted TOF written to file ", TOFInvertedOutFile)
        print('WriteProbOutput: Exiting')
    return                                               


################################################################################
# Main PROGRAM
################################################################################
filename = 'test2.tof_in'
parms, functions, signalFiles, backgroundFiles, errors = parseCmdFile(filename)


#==============================================================================
# print('\nFuncions: ', functions)    
# print('\nsignalFiles\n', signalFiles)
# print('\nbackgroundFiles\n', backgroundFiles)
# print()
# parms.pretty_print()
#==============================================================================

if len(errors) > 0:
    print('Errors\n', errors)
    raise SystemExit ('Error in command file' )

# print('Finished parsing command file')

#==============================================================================
# francesco command line 
# -t FitTOF -l test1 -f ERF -a None 
# -i 6008_Au_D2_v1J2_x=35.datv2 
# -b 6009_Au_D2_v1J2_x=35_offRes.datv2 
# --mintime 5  --maxtime 25
#==============================================================================
args_task = 'FitTOF'
args_label = 'test2'
args_function = 'ERF'
args_angularaveraging = 'None'
args_input = r'data\referenceline\6008_Au_D2_v1J2_x=35.datv2'
args_background = r'data\referenceline\6009_Au_D2_v1J2_x=35_offRes.datv2'
args_mintime = 5
args_maxtime = 25
args_baseline = False
args_cutoff = False
args_corrtime = False
args_temperature = False

# Translate arguments into variables
DataType = args_task

Label = args_label

if DataType == 'Calibration':
    ProbCurveType = 'Calibration'
else:
    ProbCurveType = args_function
    
DataFiles = args_input.split(",")
AveragingType = args_angularaveraging
BackgroundFiles = args_background.split(",")

if len( BackgroundFiles ) != len( DataFiles ) and  BackgroundFiles  != [ "" ] :
    print("Number of Background files does not match number of Data files. Quit. ")
    quit()

FitBaseline = args_baseline
FitCutoffFunction = args_cutoff
FitCorrectionTime = args_corrtime
FitTemperature = args_temperature
Tmin = args_mintime
Tmax = args_maxtime
    

States = []; DataSets = []
# Params = Parameters()
Params = parms

# Read data from input
for i in range( len( DataFiles)):
    DataFile = DataFiles[i]
    if BackgroundFiles  != [ "" ]:
        BackgroundFile = BackgroundFiles[i]
    else:
        BackgroundFile = BackgroundFiles[0]
        
    State, Temperature, Time, Signal = ReadData( DataFile, BackgroundFile, Tmin, Tmax)
    States.append( State )
    DataSets.append( [ Time, Signal ] )
  
    # Add parameters to the dictionary
    Params = SetupParametersDictionary( i, Params, Temperature )

#==============================================================================
# SetFreeParameters( Params, DataType, ProbCurveType, FitTemperature, 
#                   FitCutoffFunction, FitCorrectionTime, FitBaseline )
#==============================================================================

# Generate Theta angles employed for angular averaging
ThetaAngles = GenerateThetaAngles(                                  \
    AveragingType=AveragingType, GridType=GridType,                 \
    NPointsSource=NPointsSource, NPointsDetector=NPointsDetector,   \
    ZSource = ZSource,           RSource = RSource,                 \
    ZAperture = ZAperture,       RAperture = RAperture,             \
    ZDetector = ZLaser,          LengthDetector = LLaser)
    #    ZDetector = ZFinal,          LengthDetector = 2.*RFinal         \
  
print()
print('Angular Averaging Parameters:')
print('nps =', NPointsSource, 'npd =', NPointsDetector)
print('zs  =', ZSource,       'rs  =', RSource )
print('za  =', ZAperture,     'ra  =', RAperture)     
print('zf  =', ZFinal,        'rf  =', RFinal)                                                

# Fit TOF
fitResult = FitData( DataSets, Params, AveragingType, ThetaAngles, ProbCurveType, Label )

print(fit_report(fitResult))

# LockParam(Params)     

#==============================================================================
# print('\n' + 25*'#')
# print('Main: After  FitData')
# print('Params[E0_1].value = ', Params['E0_1'].value)
# print(lmfit.fit_report(fitResult))
#==============================================================================

                                            
for i in range( len( DataSets )):
    # Write output
    NDataSet =  i + 1
    Time     =  DataSets[i][0]
    Signal   =  DataSets[i][1]
    State    =  States[i]
    
    baseFileName = Label + "_DataSet_" + str( NDataSet ) + "_" + State
    PlotFit(baseFileName)
    
    WriteTOFOutput( Time, Signal, State, Label, NDataSet, fitResult.params, 
                   AveragingType, ThetaAngles, ProbCurveType)

    if DataType != 'Calibration':
        WriteProbOutput(Time, Signal, State, Label, NDataSet, fitResult.params, 
                        AveragingType, ThetaAngles, ProbCurveType)                                                   
