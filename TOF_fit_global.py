
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 13:42:22 2016

@author: dja
"""

from Parameters2 import Parameters2

class TOF_fit_global(object):
    def __init__(self):
        
        #------------------------------------------------------------------------------------------
        # For testingPhysical Constants
        #------------------------------------------------------------------------------------------
        self.parms = Parameters2()
        self.list1 =[]
        self.float_variable = 0.0
        self.int_variable = 0

        #------------------------------------------------------------------------------------------
        # Physical Constants
        #------------------------------------------------------------------------------------------
        self.kb          = 8.6173324E-5                 # Boltzmann constant in eV/K
        self.eV2J        = 1.602176565E-19              # eV to Joule
        self.J2eV        = 6.24150934E18                # Joule to eV
        self.AtomicMass  = 1.660538921E-27              # AMU in kg
        self.eVConst     = self.AtomicMass * self.J2eV

        self.massH       = 1.00782503223
        self.massD       = 2.01410177812
        self.massH2      = 2 * self.massH
        self.massHD      = self.massH + self.massD
        self.massD2      = 2 * self.massD

        #------------------------------------------------------------------------------------------
        # Angular averaging parameters
        #------------------------------------------------------------------------------------------
        # Parameters for point detector
        self.AngRes    = 20                     # Angular resolusion (degrees)
        self.ThetaStep = 1                      # Theta step in averaging

        # Parameters for line detector
        self.ZDiffWall = 0.0                    # Poistion of inside (detector side)
        self.RDiffWall = 6.0                    # Radius of the hole in differential wall
        self.ZRef = self.ZDiffWall              # Use Diff wall as reference point since
                                                #   Source and differential wall positions
                                                #   might be changed
        self.ZAperture = self.ZRef - 3.5        # Position of Aperture in Ta Sheiled
        self.RAperture = 1.5                    # Radius of aperture in Ta shield
        self.ZSource   = self.ZAperture - 4.    # Position of Source
        self.RSource   = 0.1                    # Radius of the source (Source or knudsen)

        # Detection volume is determined by length of REMPI laser beam.
        # the actual length is very long, so we should set this parameter to be
        # long enough that it does not limit detection.
        self.laser_to_wall = 5.0                # Distance of REMPI laser from differential wall
        self.ZLaser = self.ZRef + self.laser_to_wall # Position of REMPI laser beam
        self.LLaser      = 9.3                  # Length of REMPI detection volume.

        # In the present code we don't us these parameter.  Instead
        # we use an effective detection line based on the acceptance
        # of ions at the final field free region, i.e we assume the length of the
        # REMPI detection volume is not a limiting factor

        self.ZFinal     = self.ZRef + 34.       # Position of the final grid
        self.RFinal     = 10.0                  # Effective acceptance radius for ions at
                                                #   at the final grid.  Because of strong
                                                #   accleration to extractor take this to
                                                #   be equal to the extrator radius

        # self.NPointsDetector = 101            # Number of points to consider on the line
        self.NPointsDetector = 11               # Number of points to consider on the line
        # self.NPointsSource   = 1              # Number of points to consider on the Source
        self. NPointsSource   = 10              # Number of points to consider on the Source
                                                #   If NPointsSource = 1 we treat this as
                                                #   point source

        self.GridType = 'Cartesian'             # Generate a cartesian or radial grid on
                                                #   the source. This parameter can have
                                                #   values of 'Cartesian' or'Radial'

        #------------------------------------------------------------------------------------------
        # Fit Control Variables
        #------------------------------------------------------------------------------------------
        self.AveragingType  = None; self.AveragingTypes  = []
        self.backgroundFile = None; self.backgroundFiles = []
        self.cutoff_type    = None; self.cutoff_types    = []
        self.Function       = None; self.Functions       = []
        self.signalFile     = None; self.signalFiles     = []
        self.Tmin           = None; self.Tmins           = []
        self.Tmax           = None; self.Tmaxs           = []

        #------------------------------------------------------------------------------------------
        # Data Header Format
        #------------------------------------------------------------------------------------------
        self.DataFormatLine   =  2;   self.DataFormatLines  = []
        self.OriginalFileLine =  3;   self.OriginalFileLines= []
        self.FileDateLine     =  4;   self.FileDataLines    = []
        self.MoleculeLine     =  5;   self.MoleculeLines    = []
        self.TemperatureLine  =  6;   self.TemperatureLines = []
        self.VibStateLine     =  7;   self.VibStateLines    = []
        self.RotStateLine     =  8;   self.RotStateLines    = []
        self.DataColLine      =  9;   self.DataColLines     = []
        self.DataRowLine      = 10;   self.DataRowLines     = []

        #------------------------------------------------------------------------------------------
        #   Miscellaneous Variables
        #------------------------------------------------------------------------------------------
        self.comment_xlsx = None
        self.file_label   = None
        self.ThetaAngles  = []