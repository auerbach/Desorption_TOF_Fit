#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" 
Read data module for desorption tof fit
"""

import os
import numpy as np


class Data(object):
    def __init__(self, glbl):
        self.run_number         = 0
        
        #initialize lists
        self.glbl                       = glbl
        self.background_dates           = []
        self.background_names           = []
        self.baselines                  = []
        self.datasets                   = []
        self.fit_index_ranges           = []
        self.mass_molecules             = []
        self.molecules                  = []
        self.original_signal_names      = []
        self.original_background_names  = []
        self.signal_dates               = []
        self.signal_names               = []
        self.states                     = []
        self.temperatures               = []
        
        
    def is_number(self, a_string):
        try:
            float(a_string)
            return True
        except ValueError:
            return False
        except:
            return False

    #------------------------------------------------------------------------------
    # read_data  function to read the data and subtract background
    #------------------------------------------------------------------------------
    def read_data(self, signal_filename, background_filename ='', fit_range = None,
                  baseline_range = None,  Threshold = 0.05):
        # function to read Datafiles
        # Input: DataFile = name of the file to be read ;
        # BackgroundFile  = name of the file with background; 
        # t_min, t_max = minimum (maximum) value of time (in us) to consider ;
        # Threshold = ratio signal/maximum of signal for signal to be considered
    
        # increment run_number every time we read a data file
        self.run_number += 1
        
        # Check existence of the file
        if not os.path.isfile(signal_filename):
            print(signal_filename + " : signal file does not exist! Quitting...")
            raise SystemExit
        
        # Open and read the signal file    
        with open(signal_filename, 'r') as signal_file:
            sig_lines = signal_file.readlines()
        
        if background_filename :
            if not os.path.isfile(background_filename):
                print(background_filename + " : Background file does not exist! Quitting...")
                raise SystemExit
    
            with open(background_filename, 'r') as background_file:
                back_lines = background_file.readlines()
                
        data_format = sig_lines[self.glbl.DataFormatLine - 1].split()[3]
        
        if data_format != '2.1':
            print('File ', signal_filename, ' is not in the right format')
            print('Format = ', data_format)
            raise SystemExit
        
        self.signal_names.append(signal_filename)
        self.background_names.append(background_filename)
        self.original_signal_names.append(sig_lines[self.glbl.OriginalFileLine - 1].split()[3])
        self.molecules.append(sig_lines[self.glbl.MoleculeLine - 1].split()[3])
        self.temperatures.append(float(sig_lines[self.glbl.TemperatureLine - 1].split()[3]))
        self.states.append((int(sig_lines[self.glbl.VibStateLine - 1].split()[3]),
                             int(sig_lines[self.glbl.RotStateLine - 1].split()[3])))
                            
        sig_data_col  = int(sig_lines[self.glbl.DataColLine - 1].split()[3])
        sig_data_row  = int(sig_lines[self.glbl.DataRowLine - 1].split()[3])
        
        sig_row_end   = len( sig_lines ) #assuming data last until end of file
        back_data_col = int(sig_lines[self.glbl.DataColLine - 1].split()[3])
        if background_filename:
            back_data_row = int(sig_lines[self.glbl.DataRowLine - 1].split()[3])
            back_row_end  = len(back_lines)
        
        
        if background_filename:
            if (sig_data_row - sig_row_end) != (back_data_row -back_row_end):
                print('Length of signal and background data are not equal')
                print('  Signal file:    ', signal_filename)
                print('  Background file:', background_filename)
                raise SystemExit
                
        # get the mass of the molecule
        molecule = self.molecules[self.run_number - 1]
        if molecule  == 'H2':
            self.mass_molecules.append(self.glbl.massH2)
        elif molecule == 'HD':
            self.mass_molecules.append(self.glbl.massHD)
        elif molecule == 'D2':
            self.mass_molecules.append(self.glbl.massD2)
        else:
            print('Error:  Unknown Molecule ', self.molecules[-1])
            raise SystemExit
        
        time = []
        signal = []
           
        for n_sig in range(sig_data_row - 1, sig_row_end):
            T = float( sig_lines[n_sig].split()[0] )*1E-6   # T is ins sec , Gottingen data                                                # is in microsec
            S = float( sig_lines[n_sig].split()[sig_data_col - 1] ) #
            time.append( T )

            if not background_filename :
                signal.append(S)
            else:
                n_back = n_sig + back_data_row - sig_data_row
                B = float(back_lines[n_back].split()[back_data_col -1])
                signal.append(S - B)
    
        delta_time = time[1] - time[0]
        self.datasets.append([np.array(time),np.array(signal)])
    
        #==========================================================================================
        # Select data to fit:
        #   if fit_range is given, use it
        #   otherwise find times where data is = max of data * Threshold
        #
        # Since the data is noisy, for point n average from n-n_delt to n+n_delt
        #==========================================================================================
        t_min = None
        t_max = None
        try:         
            if self.is_number(fit_range[0]):        
                t_min = fit_range[0]
            if self.is_number(fit_range[1]):
                t_max = fit_range[1]
        except:
            pass
        
        n_delt = 30
        signal_max = np.array(signal[100:len(signal)]).max()
            
        # find n_min = index of first point to use in fitting
        # scan through the data and look for point where data exceeds a threshold
        for n in range(n_delt+1, len( time )-(n_delt+1)):

            # ignore data imes less than ion flight time to avoid noise.
            if time[n] < 3.2E-6:      # ion flight time for H2 is 3.2E-6.
                continue

            # if t_min is specified, find the index of the corresponding point in time array
            if self.is_number(t_min):
                if float(t_min) * 1.E-6 >= time[n] and \
                   float(t_min) * 1.E-6 < time[n] + delta_time:
                    n_min = n
                    break

            # otherwise find the index where data exceeds threshold
            else:
                if np.array(signal[n-n_delt:n + n_delt]).mean() >= signal_max * Threshold:
                    n_min = int(n - 1 / (delta_time * 1E6))
                    break

        # repeat for long flight time side of cureve
        for n in range( len(time) -(n_delt + 1), 100, -1 ):
            if self.is_number(t_max) :
                if float(t_max) * 1.E-6 >= time[n] and \
                   float(t_max) * 1.E-6 <  time[n] + delta_time:
                    n_max = n
                    break
            else:
                if np.array(signal[n-n_delt:n + n_delt]).mean()  >= signal_max * Threshold:
                    n_max = int(n + 1.5 / (delta_time * 1E6))
                    break
        
        # save the range of signal to use in fitting
        self.fit_index_ranges.append((n_min, n_max))

        # =========================================================================================
        # if baseline range is supplied - compute the average baseline
        # =========================================================================================
        baseline = None        
        
        if(baseline_range):
            if len(baseline_range) % 2 == 0:
                n_segments = len(baseline_range) / 2

            else:
                print('***** error baseline_range = ', baseline_range)
                print('***** number of segments must be even')
                raise SystemExit('number of baseline_range elements is odd')

            avg_segments = np.array(baseline_range).reshape(n_segments, 2)
            
            first = True        
            for a_baseline_range in avg_segments:
                a_baseline_range = a_baseline_range * 1E-6
                condition = np.logical_and(time >= a_baseline_range[0], time <= a_baseline_range[1])
                if first:
                    conditions = condition
                    first = False
                else: 
                    conditions = np.logical_or(condition, conditions,)
            
            baseline = np.array(signal)[conditions].mean()
        
        self.baselines.append(baseline)

            
        print()
        print('Read data:')
        print('  Range of points to use in fitting')
        print('  n_min      =',n_min, '     t_min=',time[n_min])
        print('  n_max      =',n_max, '     t_max=',time[n_max])
        print('  signal_max =',max( signal[50:len(signal)] ))
        print('  baseline   =', baseline)
        print()
    
        #return State, Temperature, n_min, n_max, time, Data


#==============================================================================
# test class if executed as main
#==============================================================================
if __name__ == '__main__':

    import TOF_fit_global
    glbl= TOF_fit_global.TOF_fit_global()

    data = Data(glbl)
    
    sigfn = 'Data\\2016.02.10\\Au(111)_H2_v1J3.datv2'
    backfn = ''
    fit_range = (7.0, 41.)
    baseline_range = (3., 6., 45., 55.)
    
    data.read_data(sigfn, backfn, fit_range, baseline_range)
    
    attributes = vars(data)
    for item in attributes:
        if item != 'datasets':
            print('{:25} = {}'.format(item, attributes[item]))

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