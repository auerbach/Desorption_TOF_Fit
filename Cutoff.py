#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Cutoff function to account for loss of very low energy ions
"""

import numpy as np
import scipy as sp
from Parameters2 import Parameters2

def cutoff_function(params, data, glbl, n_dataset, time, cutoff_type, debug = False):
    if cutoff_type.lower() == 'tanh':
        tcutc    = params['tcutc_%i' % n_dataset].value
        tcutw    = params['tcutw_%i' % n_dataset].value
        CutOff = 0.5 * (1. - np.tanh((time - tcutc * 1E-6) / (tcutw * 1E-6)))
        
    elif cutoff_type.lower() == 'erf':
        tcutc    = params['tcutc_%i' % n_dataset].value
        tcutw    = params['tcutw_%i' % n_dataset].value
        CutOff = 0.5 * (1. - sp.special.erf((time - tcutc * 1E-6) / (tcutw * 1E-6)))
    
    elif cutoff_type.lower() == 'exp':
        ffr_dist  = params['ffr_%i' % n_dataset].value * 1E-3
        ecutm    = params['ecutm_%i' % n_dataset].value * 1E-3 # mev --> eV
        ecuts    = params['ecuts_%i' % n_dataset].value * 1E+3 # 1/mev --> 1/eV
        mass = data.mass_molecules[n_dataset - 1]
        Velocity = ffr_dist / time
        Ekin = (0.5 * mass * Velocity**2.) * glbl.eVConst
        CutOff = (1.0 - np.exp(-ecuts * (Ekin - ecutm)))  * (Ekin > ecutm)
#==============================================================================
#         This code is slow, above line does the same but as vector
#         CutOff = np.zeros(len(Ekin))
#         for i in range(len(Ekin)):     
#             if Ekin[i] > ecutm:
#                 CutOff[i] = 1.0 - np.exp(-ecuts * (Ekin[i] - ecutm))
#==============================================================================
            
    else:
        print('***** Error unknown cutoff type = ', cutoff_type)
        raise SystemExit
        
    return CutOff




#------------------------------------------------------------------------------
#  Main program for testing
#------------------------------------------------------------------------------
if __name__ == '__main__':
    
    import unicodedata
    import matplotlib.pyplot as plt
    import TOF_fit_global
    import Data

    glbl = TOF_fit_global.TOF_fit_global()
    data = Data.Data(glbl)
    params = Parameters2()
    
    data.mass_molecules.append(glbl.massD2)
    
    # data_filename = 'Data\\2016.02.10\\Au(111)_D2_v0J2.datv2'
    # data.read_data(data_filename)
    
    NDataSet = 1    
    debug = True    
    time = np.arange(3,50,.1) * 1E-6
    
    cutoff_type = 'tanh'
    params.add('tcutc_1', 28)
    params.add('tcutw_1', 4)   
    tanh_cut = cutoff_function(params, data, glbl, NDataSet, time, cutoff_type )

    cutoff_type = 'erf'
    params.add('tcutc_1', 28)
    params.add('tcutw_1', 4)   
    erf_cut = cutoff_function(params, data, glbl, NDataSet, time, cutoff_type )
    
    cutoff_type = 'exp'
    params.add('ecutm_1', 20.5)
    params.add('ecuts_1', 0.12)
    params.add('ffr_1'  , 31.72861)
    exp_cut = cutoff_function(params, data, glbl,NDataSet, time, cutoff_type )
        
    label=''
    parm_for_plot_label = ['tcutc', 'tcutw', 'ecutm', 'ecuts']
    for p in parm_for_plot_label:
        p_ = p + '_1'
        label += '{:6s}{:>7.3f}\n'.format(p, params[p_].value)
    
    ax = plt.subplot(111)
    ax.set_xlabel('TOF [' + unicodedata.lookup('micro sign') + 's]', fontsize=16)
    ax.set_ylabel('Cutoff function', fontsize=16)
    
    ax.annotate(label, xy = [35, 0.7], xycoords = 'data', 
                va = 'top', family='monospace', fontsize=12)    
    
    plt.ylim(-0.05, 1.05)
    plt.plot(time*1E6, exp_cut, linewidth = 1.5, label ='exp cutoff')
    plt.plot(time*1E6, tanh_cut, linewidth = 1.5, linestyle ='--', label = 'tanh cutoff')
    plt.plot(time*1E6, erf_cut, linewidth = 1.5, linestyle ='--', label = 'erf cutoff')
    plt.legend(loc=1)
    plt.show()
    
    
    