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
    
    # time, sig, fit = np.loadtxt(filename, unpack=True)
    
    with open(filename,'r') as fit_out:
        lines = fit_out.readlines()
    
    
    #----------------------------------------------------------------------------------------------
    # get the number of plots to make 
    #----------------------------------------------------------------------------------------------    
    for n_line, line in enumerate(lines):
        if line.startswith('# Number of plots'):
            n_plots = int(line.split(':')[1])
            break
    #----------------------------------------------------------------------------------------------
    # make a plot for each data section
    #----------------------------------------------------------------------------------------------
    for n_plot in range(n_plots):
        
        #------------------------------------------------------------------------------------------
        # get the plot title and label from the result file        
        #------------------------------------------------------------------------------------------
        label = ''        
        for n_line in range(n_line + 1, len(lines)):
            line = lines[n_line]
            
            if line.startswith('# Title'):
                title = line.split(':', 1)[1]
            
            if line.startswith('# Label'):
                label += line.split(':', 1)[1]
            
            if line.startswith('# Begin data'):
                n_points = int(line.split(':')[1])
                break
            
        #------------------------------------------------------------------------------------------
        # get the data to plot
        #------------------------------------------------------------------------------------------
        time = []
        sig  = []
        fit  = []
        for n_line in range(n_line + 1, n_line + 1 + n_points ):
            line = lines[n_line]
            time.append(float(line.split()[0]))
            sig.append (float(line.split()[1]))
            fit.append (float(line.split()[2]))
        time = np.array(time)
        sig  = np.array(sig)            
        fit  = np.array(fit)            
    
        fig = plt.figure(figsize = (6,6), dpi = 200)
        fig = plt.figure()
        fig.suptitle(title, fontsize=16)
        ax = plt.subplot(111)
        ax.set_xlabel('TOF [micro sec]', fontsize=16)
        ax.set_ylabel('Signal', fontsize=16)
    
        ax.annotate(label, xy = [0.5, 0.75, ], xycoords = 'axes fraction', 
                    va = 'top', family='monospace')
    
        print('max data, fit =', sig[100:300].max(), fit[100:300].max())
    
    y1 = sig
    y2 = fit
        
    plt.xlim((3, 15))
    plt.plot(time, y1, 'b.')
    plt.plot(time, y2, 'r', linewidth = 2)
    plt.show()
    
if __name__ == '__main__':
    
    plot_file_name = 'fits\\fit010_test1_v1j3_ERF.fit_out'
    PlotFit(plot_file_name)