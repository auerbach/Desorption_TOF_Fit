#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Cutoff function to account for loss of very low energy ions
"""

import numpy as np
from Parameters2 import Parameters2
import GlobalVariables as glbl
from read_data import Data

def cutoff_function(params, data, NDataSet, Time, cutoff_type, debug = False):
    if cutoff_type.lower() == 'tanh':
        TCutC    = params['TCutC_%i'     %NDataSet].value
        TCutW    = params['TCutW_%i'     %NDataSet].value
        CutOff = 0.5 * (1. - np.tanh((Time - TCutC*1E-6) / (TCutW*1E-6)))
    
    elif cutoff_type.lower() == 'exp':
        FFRDist  = params['FFR_%i'     %NDataSet].value * 1E-3
        ECutM    = params['ECutM_%i' %NDataSet].value * 1E-3 # mev --> eV
        ECutS    = params['ECutS_%i' %NDataSet].value * 1E+3 # 1/mev --> 1/eV
        mass = data.mass_molecules[NDataSet-1]
        Velocity = FFRDist / Time 
        Ekin = (0.5 * mass * Velocity**2.) * glbl.eVConst
        
        CutOff = np.zeros(len(Ekin))
        for i in range(len(Ekin)):     
            if Ekin[i] > ECutM:
                CutOff[i] = 1.0 - np.exp(-ECutS * (Ekin[i] - ECutM))
            
    else:
        print('***** Error unknown cutoff type = ', cutoff_type)
        raise SystemExit
        
    return CutOff

if __name__ == '__main__':
    
    import unicodedata
    import matplotlib.pyplot as plt
    
    data = Data()
    params = Parameters2()
    
    data.mass_molecules.append(glbl.massD2)
    
    # data_filename = 'Data\\2016.02.10\\Au(111)_D2_v0J2.datv2'
    # data.read_data(data_filename)
    
    NDataSet = 1    
    debug = True    
    time = np.arange(5,50,.1) * 1E-6
    
    cutoff_type = 'tanh'
    params.add('TCutC_1', 28)
    params.add('TCutW_1', 4)   
    tanh_cut = cutoff_function(params, data, NDataSet, time, cutoff_type )    
    
    cutoff_type = 'exp'
    params.add('ECutM_1', 21.1)
    params.add('ECutS_1', 0.034)
    params.add('FFR_1'  , 31.72861)
    exp_cut = cutoff_function(params, data, NDataSet, time, cutoff_type )    
        
    label=''
    parm_for_plot_label = ['TCutC', 'TCutW', 'ECutM', 'ECutS']
    for p in parm_for_plot_label:
        p_ = p + '_1'
        label += '{:6s}{:>7.3f}\n'.format(p, params[p_].value)
    
    ax = plt.subplot(111)
    ax.set_xlabel('TOF [' + unicodedata.lookup('micro sign') + 's]', fontsize=16)
    ax.set_ylabel('Cutoff Function', fontsize=16)
    
    ax.annotate(label, xy = [35, 0.7], xycoords = 'data', 
                va = 'top', family='monospace', fontsize=12)    
    
    plt.ylim(-0.05, 1.05)
    plt.plot(time*1E6, exp_cut, linewidth = 1.5, label ='exp cutoff')
    plt.plot(time*1E6, tanh_cut, linewidth = 1.5, linestyle ='--', label = 'tanh cutoff')
    plt.legend(loc=1)
    plt.show()
    
    
    