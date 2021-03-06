# Physical constants ###################################################

kb          = 8.6173324E-5           # Boltzmann constant in eV/K
eV2J        = 1.602176565E-19        # eV to Joule 
J2eV        = 6.24150934E18          # Joule to eV 
AtomicMass  = 1.660538921E-27        # AMU in kg
eVConst     = AtomicMass * J2eV
massH       = 1.00782503223
massD       = 2.01410177812
massH2      = 2 * massH
massHD      = massH + massD
massD2      = 2 * massD
# massAmu     = massD2

# Experimental apparatus constants #####################################
#==============================================================================
# FFRDist           = 29.0E-3        # Distance travelled by the molecule in the fieldfree region as an ion (m)
# FFRDistTolerance  = 0.5            # Percentage by which the field free region lenght can vary
# TimeCorr          = 4.6            # time correction (us)
# TimeCorrTolerance = 0.8            # (see two lines above )
#
# TCutC             = 28.6           # CutOff function 1st parameter (us) 
# TCutCTolerance    = 0.5            # ...
# TCutW             = 4.3            # CutOff function 2nd parameter (us)
# 
# TemperatureTolerance = 1.         
# 
# t_min = ''
# t_max = ''
# 
#==============================================================================

#--------------------------------------------------------------------------------------------------
# Angular averaging constants
#--------------------------------------------------------------------------------------------------
# Following parameters for point detector
AngRes               = 20          # Angular resolusion (degrees)
ThetaStep            = 1          # Theta step in averaging

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

# NPointsDetector = 101         # Number of points to consider on the line
NPointsDetector = 11            # Number of points to consider on the line   
# NPointsSource   = 1           # Number of points to consider on the Source
NPointsSource   = 10            # Number of points to consider on the Source
                                #   If NPointsSource = 1 we treat this as
                                #   point source

GridType = 'Cartesian'          # Generate a cartesian or radial grid on  
                                #   the source. This parameter can have
                                #   values of 'Cartesian' or'Radial'
# Data Format ##########################################################


DataFormatLine   =  2;   DataFormatLines  = []
OriginalFileLine =  3;   OriginalFileLines= []
FileDateLine     =  4;   FileDataLines    = []
MoleculeLine     =  5;   MoleculeLines    = []
TemperatureLine  =  6;   TemperatureLines = []     
VibStateLine     =  7;   VibStateLines    = []
RotStateLine     =  8;   RotStateLines    = []
DataColLine      =  9;   DataColLines     = []
DataRowLine      = 10;   DataRowLines     = []
Tmin             = None; Tmins = []
Tmax             = None; Tmaxs = []
Function         = None; Functions = []
AveragingType    = None; AveragingTypes = []
cutoff_type      = None; cutoff_types = []   


#--------------------------------------------------------------------------------------------------
#   Miscelaneous Constants and Variables
#--------------------------------------------------------------------------------------------------
comment_xlsx = None