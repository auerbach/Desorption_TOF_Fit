#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" 
Read data module for desorption tof fit
"""

import os
import numpy as np
import GlobalVariables as glbl

class Data(object):
    def __init__(self):
        self.run_number         = 0
        
        #initialize lists
        self.background_dates   = []
        self.background_names   = []
        self.datasets           = []
        self.fit_ranges         = []
        self.mass_molecules     = []
        self.molecules          = []
        self.original_signal_names  = []
        self.original_background_names = []        
        self.signal_dates       = []
        self.signal_names       = []        
        self.states             = []
        self.temperatures       = []
        
        


    #------------------------------------------------------------------------------
    # read_data  function to read the data and subtract background
    #------------------------------------------------------------------------------
    def read_data(self, signal_file_name, background_file_name = '', t_min ='', t_max ='', Threshold = 0.10):
        # Function to read Datafiles
        # Input: DataFile = name of the file to be read ;
        # BackgroundFile  = name of the file with background; 
        # t_min, t_max = minimum (maximum) value of time (in us) to consider ;
        # Threshold = ratio signal/maximum of signal for signal to be considered
    
        # increment run_number every time we read a data file
        self.run_number += 1
        
        # Check existence of the file
        if not os.path.isfile(signal_file_name):
            print(signal_file_name + " : signal file does not exist! Quitting...")
            raise SystemExit
        
        # Open and read the signal file    
        with open(signal_file_name, 'r') as signal_file:
            sig_lines = signal_file.readlines()
        
        if background_file_name :
            if not os.path.isfile(background_file_name):
                print(background_file_name + " : Background file does not exist! Quitting...")
                quit()
    
            with open(background_file_name, 'r') as background_file:
                back_lines = background_file.readlines()
                
        data_format = sig_lines[glbl.DataFormatLine -1].split()[3]
        
        if data_format != '2.1':
            print('File ', signal_file_name, ' is not in the right format')
            print('Format = ', data_format)
            raise SystemExit
        
        self.signal_names.append(signal_file_name)
        self.background_names.append(background_file_name)
        self.original_signal_names.append(sig_lines[glbl.OriginalFileLine - 1].split()[3])
        self.molecules.append(sig_lines[glbl.MoleculeLine - 1].split()[3])
        self.temperatures.append(float(sig_lines[glbl.TemperatureLine - 1].split()[3]))
        self.states.append( (int(sig_lines[glbl.VibStateLine - 1].split()[3]),
                             int(sig_lines[glbl.RotStateLine - 1].split()[3])))
                            
        sig_data_col  = int(sig_lines[glbl.DataColLine - 1].split()[3])
        sig_data_row  = int(sig_lines[glbl.DataRowLine - 1].split()[3])
        
        sig_row_end   = len( sig_lines ) #assuming data last until end of file
        back_data_col = int(sig_lines[glbl.DataColLine - 1].split()[3])
        if background_file_name:
            back_data_row = int(sig_lines[glbl.DataRowLine - 1].split()[3])
            back_row_end  = len(back_lines)
        
        
        if background_file_name:        
            if (sig_data_row - sig_row_end) != (back_data_row -back_row_end):
                print('Length of signal and background data are not equal')
                print('  Signal file:    ', signal_file_name)
                print('  Background file:', background_file_name)
                raise SystemExit
                
        # get the mass of the molecule
        molecule = self.molecules[self.run_number - 1]
        if molecule  == 'H2':
            self.mass_molecules.append(glbl.massH2)
        elif molecule == 'HD':
            self.mass_molecules.append(glbl.massHD)
        elif molecule == 'D2':
            self.mass_molecules.append(glbl.massD2)
        else:
            print('Error:  Unknown Molecule ', self.molecules[-1])
            raise SystemExit
        
        time = []
        signal = []
           
        for n_sig in range(sig_data_row - 1, sig_row_end):
            T = float( sig_lines[n_sig].split()[0] )*1E-6   # T is ins sec , Gottingen data                                                # is in microsec
            S = float( sig_lines[n_sig].split()[sig_data_col - 1] ) #
            time.append( T )

            if not background_file_name :
                signal.append(S)
            else:
                n_back = n_sig + back_data_row - sig_data_row
                B = float(back_lines[n_back].split()[back_data_col -1])
                signal.append(S - B)
    
        delta_time = time[1] - time[0]
        self.datasets.append([np.array(time),np.array(signal)])
    
        #==============================================================================================
        # Select good data: 
        #   if max or min times are provided use them, otherwise find times where
        #     data is = max of data * Threshold
        # Since the data is noisy, for point n average from n-n_delt to n+n_delt
        #==============================================================================================
        n_delt = 10
        signal_max = np.array(signal[100:len(signal)]).max()
            
        # find n_min = index of first point to use in fitting
        for n in range(100, len( time )-n_delt):
            if t_min :
                if float(t_min) * 1.E-6 >= time[n] and \
                   float(t_min) * 1.E-6 < time[n] + delta_time:
                    n_min = n
                    break
            else:
                # find time where signal[n] reaches threshold and subtract 1 microseconds
                if np.array(signal[n-n_delt:n + n_delt]).mean() >= signal_max * Threshold:
                    n_min = int(n - 1 / (delta_time * 1E6))
                    break
        
        for n in range( len(time) -(n_delt + 1), 100, -1 ):
            if t_max :
                if float(t_max) * 1.E-6 >= time[n] and \
                   float(t_max) * 1.E-6 <  time[n] + delta_time:
                    n_max = n
                    break
            else:
                # find time where signal[n] falls to threshold and add 1.5 microseconds
                if np.array(signal[n-n_delt:n + n_delt]).mean()  >= signal_max * Threshold:
                    n_max = int(n + 1.5 / (delta_time * 1E6))
                    break
        
        # save the range of signal to use in fitting
        self.fit_ranges.append((n_min, n_max))
            
        print()
        print('Range of points to use in fitting')                
        print('n_min      =',n_min, '     t_min=',time[n_min])
        print('n_max      =',n_max, '     t_max=',time[n_max])
        print('signal_max =',max( signal[50:len(signal)] ))
        print()
    
        #return State, Temperature, n_min, n_max, Time, Data


#==============================================================================
# test class if executed as main
#==============================================================================
if __name__ == '__main__':
    data = Data()
    data.read_data('Data\\2016.02.10\\Au(111)_H2_v1J3.datv2')
    attributes = vars(data)
    for item in attributes:
        if item != 'datasets':
            print('{:25} = {}'.format(item, attributes[item]))
            
    
#==============================================================================
#     datasets = data.datasets
#     print('\ndatasets\n--------\n', datasets)
#     print('\ndatasets[0]\n-----------\n', datasets[0])
#     print('\ndatasets[0][0]\n-------------\n', datasets[0][0])
#     print()
#
#     rows = len(datasets[0][0])
#     print('rows=', rows)
#     for i in range(0,1000,10):
#         print('{:6.2f}    {:10.7f}'.format(datasets[0][0][i]*1E6, datasets[0][1][i]))
#==============================================================================
    
    time   = data.datasets[0][0]
    signal = data.datasets[0][1]
    
    import matplotlib.pyplot as plt

    fig = plt.figure()
    
    file_name = data.original_signal_names[0]
    fig.suptitle(file_name, fontsize=14)
    
    ax = plt.subplot(111)
    ax.set_xlabel('TOF [micro sec]', fontsize=14)
    ax.set_ylabel('Signal', fontsize=14)
    
    plt.xlim(4, 20)
    plt.plot(time * 1E6, signal, 'b.')
    # plt.plot(t2 * 1E6, f1, linewidth = 2.5)
    plt.show()