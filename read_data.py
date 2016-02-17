#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" 
Read data module for desorption tof fit
"""

import os
import numpy as np
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
    if not os.path.isfile(data_file_name):
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
    
    # State = "v" + str( int(glbl.VibState)) + "J" + str( int(glbl.RotState))
    
    if glbl.Molecule  == 'H2':
        glbl.MassAmu = glbl.massH2
    elif glbl.Molecule == 'HD':
        glbl.MassAmu = glbl.massHD
    elif glbl.Molecule == 'D2':
        glbl.MassAmu = glbl.massD2
    else:
        print('Error:  Unknown Molecule ', glbl.Molecule)
        raise SystemExit
    
    #State = "v" + str( int(glbl.VibState)) + "J" + str( int(glbl.RotState))
    State = (int(glbl.VibState), int(glbl.RotState))
    
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
        if not background_file_name :
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
