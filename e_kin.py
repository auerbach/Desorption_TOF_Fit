#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
Function to computer kinetic energy (eV)

    e_kin(parms, data, time, n_dataset)  

    parms:      Parameters2 object with the paramaters of the fit
    data:       Data object with the data (including mass of the molecule (AMU)
    time:       float or nparray (seconds)
    n_dataset:  the number of the data set used in the fit (starts at 1)
    
    return energy (eV)

Created on Wed Mar  2 10:29:48 2016
@author: dja
"""

from Parameters2 import Parameters2
from Data import Data
import global_variables_old as glbl


def e_kin(parms, data, time, n_dataset):

    iontof = parms['IonTOF_%i'    %n_dataset].value
    time = (time - iontof) * 1E-6
    ffr_length = parms['FFR_%i'   %n_dataset].value * 1E-3  # m
    velocity = ffr_length / time                            #
    mass = data.mass_molecules[n_dataset - 1]               # AMU
    
    return (0.5 * mass * velocity**2.) * glbl.eVConst
    
    
if __name__ == '__main__':
    parms = Parameters2()
    data  = Data()
    
    parms.add('FFR_1', 31.72861)        # mm per Sven
    parms.add('IonTOF_1', 4.525)        # micro secconds for D2 per Sven
    data.mass_molecules = [glbl.massD2] # AMU
    
    while True:
        ans = input('time (micro sec) (enter to end)')
        if ans == '': 
            break
        else:
            time = float(ans) 
            print('energy = ', e_kin(parms, data, time, 1))