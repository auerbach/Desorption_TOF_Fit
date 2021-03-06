#! /usr/bin/env python

# Assuming the following geometry for the experiment:
#
#          _ _ _ _            Solid angle acceptance
#         \      /
#          \    /
#           \  /
#            \/
#         --------            Detection line
#
#----------     ------------  Aperture
#
#
#        ____________         Source
#
# Center of the Source, of the aperture and of the detection line are collinear.


import numpy as np

# Generate points on the Source (radius of the Source in mm)
def points_on_the_source(grid_type, n_points, r_source, z_source):
    Grid = []
    if n_points == 1:
        Grid.append( [0., 0., z_source ] )
        return Grid
    
    if grid_type == 'Radial':
        # Center of the circle has weight 0, so we add it to the grid of points
        Grid.append( [0., 0., z_source ] ) 
        #-------------------------------------------------
        # this doesn't make sense check with Francesco
        #-------------------------------------------------
        
        dR = r_source / float( n_points )
        Rmin = dR 
        # We start to add points from dR to z_source
        for r in np.linspace( Rmin, r_source, num = n_points, endpoint=True):
            dThetaAngle = dR / r
            for ThetaAngle in np.arange( 0, 2*np.pi, dThetaAngle):
                XPoint = r * np.cos( ThetaAngle )
                YPoint = r * np.sin( ThetaAngle )
                Grid.append( [XPoint, YPoint, z_source] )
                
    if grid_type == 'Cartesian':
        if ( n_points ) % ( 2 ) == 0.:
            NLPoints = n_points + 1 # Always include the center
        else:
            NLPoints = n_points

        for XPoint in np.linspace( -r_source, r_source, num = NLPoints, endpoint=True):
            for YPoint in np.linspace( -r_source, r_source, num = NLPoints, endpoint=True):
                if (XPoint**2. + YPoint**2.) <= r_source **2. :
                    Grid.append( [XPoint, YPoint, z_source] )
                    #print XPoint, z_source
    return Grid

# Generate points on the detection line (length of the detection line in mm)
def points_on_the_detection_line(n_points, Length, ZDetector):
    Grid = [ ]
    for XPoint in np.linspace( start=-Length /2.,\
                               stop= Length /2.,\
                               num = n_points,\
                               endpoint=True):
        Grid.append( [XPoint, 0., ZDetector] )
    return Grid

def ThetaPossibleTrajectories(grid_source, grid_detector, ZAperture, RAperture):
    
    ThetaAngles = []
    # Generate all possibles trajectories
    for CoordinatesSource in grid_source:
        for CoordinatesDetector in grid_detector:
            # Equation of the line of the form:
            # ( x - q0 ) / m0 = ( y - q1 ) / m1 = ( z - q2 ) / m2
        
            m = [ CoordinatesDetector[0] - CoordinatesSource[0],\
                  CoordinatesDetector[1] - CoordinatesSource[1],\
                  CoordinatesDetector[2] - CoordinatesSource[2] ]
                  
            q = [ CoordinatesSource[0],\
                  CoordinatesSource[1],\
                  CoordinatesSource[2] ]

            # Select trajectories passing through the aperture:
            r2 = ((ZAperture - q[2])/m[2] * m[0] + q[0])**2. +\
                ((ZAperture - q[2])/m[2] * m[1] + q[1])**2.
            print('\ncoord det =', CoordinatesDetector)
            print('coord src =', CoordinatesSource)
            print('sqrt(r2)        =', np.sqrt(r2))
            if (r2 <= RAperture **2.) :   
                # Pass to spherical coordinates (with center in desorbtion point) to determine theta
                r = np.sqrt( m[0]**2. + m[1]**2. + m[2]**2. )
                ThetaAngle = np.arccos( abs( m[2] ) / r )
                ThetaAngles.append( ThetaAngle )
    return ThetaAngles

def histogram( List, DeltaBin ):
    # This function returns an histogram given a list of theta values and a binning width (rad)
    # Only of use for test purposes
    Histogram = [ ]
    for X in np.arange( min( List ), max( List ) , DeltaBin):
        Counter = 0
        for Y in List:
            if Y >= X and Y < ( X + DeltaBin ):
                Counter = Counter + 1
        Histogram.append( [ X , Counter ] )
    return Histogram


if (__name__ == "__main__"):
    # Executing this module prints an histogram of the Theta values produced

    # Following parameters for point detector
    AngRes               = 20.         # Angular resolusion (degrees)
    ThetaStep            = 2.          # Theta step in averaging

    # Following parameters for line detector

    ZDiffWall   = 0.0               # Poistion of inside (detector side)
    RDiffWall   = 6.0               # Radius of the hole in differential wall    
    
    ZRef = ZDiffWall               # Use Diff wall as reference point since
                                    #   Source and differential wall positions
                                    #   might be changed
    
    ZAperture  = ZRef - 3.5        # Position of aperture in Ta Sheiled   
    RAperture   = 1.5               # Radius of aperture in Ta shield

    z_source     = ZAperture - 4.   # Position of Source
    r_source     = 0.1               # Radius of the source (Source or knudsen)
        
    # Detection volume is determined by length of REMPI laser.
    # the actual length is very long, so we should set this parameter to be 
    # long enough that it does not limit detection.  
    
    # In the present code we don't us this parameter.  Instead
    # we use an effective detection line based on the acceptance
    # of ions at the final field free region 
    ZLaser      = ZRef + 5.0        # Position of REMP laser beam
    LLaser      = 9.3               # Length of REMPI detection volume.  
    # the actual length is very long, so we should set this parameter to be 
    # long enough that it does not limit detection.  
    
       
    ZFinal     = ZRef + 34.         # Position of the final grid
    RFinal     = 10.0               # Effective acceptance radius for ions at
                                    #   at the final grid.  Because of strong
                                    #   accleration to extractor take this to
                                    #   be equal to the extrator radius
    
    
    
    NPointsDetector = 11            # Number of points to consider on the line
        
    NPointsSource = 1               # Number of points to consider on the Source
                                    #   If NPointsSource = 1 we treat this as
                                    #   point source
    
    GridType = 'Cartesian'          # Generate a cartesian or radial grid on  
                                    #   the source. This parameter can have
                                    #   values of 'Cartesian' or'Radial'

    AveragingType = 'Line Detector'
        
    if AveragingType == 'No Averaging':
        ThetaAngles = [0.]
        
    elif AveragingType == 'Point Detector':
        ThetaAngles = np.arange( 0., AngRes + ThetaStep, ThetaStep ) 
        
    elif AveragingType == 'Line Detector':
        GridOfPointsSource  = points_on_the_source(grid_type= GridType, \
                                                   z_source = z_source, \
                                                   r_source = r_source, \
                                                   NPoints = NPointsSource)
        print('\nSource Grid')
        for idx, val in enumerate(GridOfPointsSource):
            print(idx, val)            
            
        GridOfPointsDetector = points_on_the_detection_line(\
            NPoints = NPointsDetector,\
            Length  = LLaser,\
            ZDetector = ZLaser)   
        print('\nDetector Grid')
        for idx, val in enumerate(GridOfPointsDetector):
            print(idx, val) 
        
        ThetaAngles = ThetaPossibleTrajectories(
            GridOfPointsSource, GridOfPointsDetector, ZAperture, RAperture)
        print('\nTheta Angles')
        for idx,val in enumerate(ThetaAngles):
            print(idx,np.degrees(val))
        
    for i in range( len( ThetaAngles )):
        ThetaAngles[i] = np.degrees( ThetaAngles[i] )
    print("\nConsidering ", len(ThetaAngles ),\
        " values of Theta in the angular averaging, minimum: %8.3f"\
        %min( ThetaAngles), " deg , maximum: %8.3f" %max( ThetaAngles) ," deg.")
        
    
    Histogram = histogram( ThetaAngles, 1 )
    sum = 0
    for n in Histogram:
        sum = sum + n[1]
        print (n[0], n[1],' #', sum)
    #print ('X=',np.arange( min( angles_list ), max( angles_list ) , 1))
    #print (angles_list)
