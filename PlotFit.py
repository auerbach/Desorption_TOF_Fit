#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Plot the data and result of fit """

def PlotFit(filename):

    """ 
    Plot the data and result of fit 
        parameter = base file name
    """
    
    import numpy as np
    import matplotlib.pyplot as plt

    
    file1 = filename + '__TOF_DataFitted.dat'
    file2 = filename + '__TOF_Fit_ERF.dat'
    
    t1, s1, s2 = np.loadtxt(file1, unpack=True)
    t2, f1, f2 = np.loadtxt(file2, unpack=True)
    
    fig = plt.figure()
    fig.suptitle('TOF Fit', fontsize=20)
    ax = plt.subplot(111)
    
    ax.set_xlabel('TOF [micro sec]', fontsize=16)
    # ax.set_ylabel('Sticking Proability', fontsize=16)
    
    plt.xlim(5E-6, 15E-6)
    plt.plot(t1,s1)
    plt.plot(t2,f1)

