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
        self.kb         = 8.6173324E-5                 # Boltzmann constant in eV/K
        self.eV2J       = 1.602176565E-19              # eV to Joule
        self.J2eV       = 6.24150934E18                # Joule to eV
        self.AMU        = 1.660538921E-27              # AMU in kg
        self.eVConst    = self.AMU * self.J2eV

        self.massH      = 1.00782503223
        self.massD      = 2.01410177812
        self.massH2     = 2 * self.massH
        self.massHD     = self.massH + self.massD
        self.massD2     = 2 * self.massD

        #------------------------------------------------------------------------------------------
        # Angular averaging parameters
        #------------------------------------------------------------------------------------------
        # Parameters for point detector
        self.ang_res    = 20                        # Angular resolusion (degrees)
        self.theta_step = 1                         # Theta step in averaging

        # Parameters for line detector
        self.z_diff_wall = 0.0                      # Poistion of inside (detector side)
        self.r_diff_wall = 6.0                      # Radius of the hole in differential wall
        self.z_grid      = 34.0                     # Distance of final grid from differential wall
        self.z_ref       = self.z_diff_wall         # Use Diff wall as reference point since
                                                    #   Source and differential wall positions
                                                    #   might be changed
        self.aperture_to_shield = 3.5               # Distance of the aperature to Ta shield
        self.z_aperture = self.z_ref - 3.5          # Position of Aperture in Ta shield
        self.r_aperture = 1.5                       # Radius of aperture in Ta shield
        self.z_source   = self.z_aperture - 4.      # Position of Source
        self.r_source   = 0.1                       # Radius of the source (Source or knudsen)

        # Detection volume is determined by length of REMPI laser beam.
        # the actual length is very long, so we should set this parameter to be
        # long enough that it does not limit detection.
        self.laser_to_wall = 5.0                    # Distance of REMPI laser from differential wall
        self.z_laser       = self.z_ref + self.laser_to_wall # Position of REMPI laser beam
        self.length_laser  = 9.3                    # Length of REMPI detection volume.
        # In the present code we don't us these parameter.  Instead
        # we use an effective detection line based on the acceptance
        # of ions at the final field free region, i.e we assume the length of the
        # REMPI detection volume is not a limiting factor

        self.z_final = self.z_ref + self.z_grid       # Position of the final grid
        self.r_final     = 10.0                  # Effective acceptance radius for ions at
                                                #   at the final grid.  Because of strong
                                                #   accleration to extractor take this to
                                                #   be equal to the extrator radius

        self.points_detector = 11               # Number of points to consider on the line
        self.points_source   = 10               # Number of points to consider on the Source
                                                #   If NPointsSource = 1 we treat this as
                                                #   point source

        self.grid_type = 'Cartesian'             # Generate a cartesian or radial grid on
                                                #   the source. This parameter can have
                                                #   values of 'Cartesian' or'Radial'

        #------------------------------------------------------------------------------------------
        # Fit Control Variables and Lists
        #------------------------------------------------------------------------------------------
        self.averaging_type         = None; self.averaging_types      = []  # string
        self.background_filename    = None; self.background_filenames = []  # string
        self.baseline_range         = None; self.baseline_ranges      = []  # tuple of integers (4)
        self.comment                = None; self.comments             = []  # string
        self.cutoff_type            = None; self.cutoff_types         = []  # string
        self.fit_range              = None; self.fit_ranges           = []  # tuple of integers (2)
        self.function               = None; self.functions            = []  # string
        self.n_delt                 = 30  ; self.n_delts              = []  # int
        self.signal_filename        = None; self.signal_filenames     = []  # string
        self.threshold              = 0.05; self.thresholds           = []  # float

        #------------------------------------------------------------------------------------------
        #   Miscellaneous Variables
        #------------------------------------------------------------------------------------------
        self.comment_xlsx          = None; 
        self.file_label            = None
        self.angles_list           = [];   self.angles_lists          = []
        