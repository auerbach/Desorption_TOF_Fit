#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Plot the data and result of fit """

def PlotFit(filename):


    """ 
    Plot the data and result of fit 
    
    filename should include path/name.ext
    plot will be saved as path/name.pdf
    """
    
    import os
    import numpy as np
    import matplotlib.pyplot as plt
    import unicodedata
    
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
        # get the title, label, information for plot from result file        
        #------------------------------------------------------------------------------------------
        label = ''        
        for n_line in range(n_line + 1, len(lines)):
            line = lines[n_line]
            
            if line.startswith('# Title'):
                title = line.split(':', 1)[1]
            
            if line.startswith('# Label'):
                label += line.split(':', 1)[1]
                
            if line.startswith('# Npoints'):
                n_points = int(line.split(':')[1])
                
            if line.startswith('# Tmin'):
               t_min, t_max = tuple(line.split(':')[1].strip().split(','))
               t_min = float(t_min)
               t_max = float(t_max)
                    
            if line.startswith('# Begin data'):
                break
        #------------------------------------------------------------------------------------------
        # get the number of data points and minimum and maximum times used for the fit       
        #------------------------------------------------------------------------------------------
            
        #------------------------------------------------------------------------------------------
        # get the data to plot
        #------------------------------------------------------------------------------------------
        time   = []
        sig    = []
        fit    = []
        cutoff = []
        for n_line in range(n_line + 1, n_line + 1 + n_points ):
            line = lines[n_line]
            if float(line.split()[0]) < 3.0 :
                continue
            time.append  (float(line.split()[0]))
            sig.append   (float(line.split()[1]))
            fit.append   (float(line.split()[2]))
            cutoff.append(float(line.split()[3]))

        time = np.array(time)
        sig  = np.array(sig)            
        fit  = np.array(fit)
        cutoff = np.array(cutoff) * fit.max()
        
        sig_max = sig.max()
        sig_min = sig.min()
        
    
        fig = plt.figure(figsize = (6,6), dpi = 200)
        fig = plt.figure()
        fig.suptitle(title, fontsize=16)
        ax = plt.subplot(111)
        ax.set_xlabel('TOF [' + unicodedata.lookup('micro sign') + 's]', fontsize=16)
        ax.set_ylabel('Signal', fontsize=16)
    
        ax.annotate(label, xy = [0.55, 0.95, ], xycoords = 'axes fraction', 
                    va = 'top', family='monospace', )    
        ax.annotate('$t_{min}$', xycoords = 'data', xy = (t_min-.1, sig_max * .4), 
                    ha = 'right', va='center', fontsize=14)        
        ax.annotate('$t_{max}$', xycoords = 'data', xy = (t_max+.1, sig_max * .4), 
                    ha = 'left',  va='center', fontsize=14)                
        ax.annotate('', xycoords = 'data', xy = (t_min , 0), xytext =(t_min, sig_max*.5), 
                    arrowprops=dict(linewidth = 1.5, linestyle = '--', arrowstyle = '-')) 
        ax.annotate('', xycoords = 'data', xy = (t_max , 0), xytext =(t_max, sig_max*.5), 
                    arrowprops=dict(linewidth = 1.5, linestyle = '--', arrowstyle = '-'))   
              
        plt.xlim((2, 60))
        plt.ylim((sig_min, sig_max * 1.05))
        plt.plot(time, sig, 'b.')
        plt.plot(time, fit, 'r', linewidth = 2)
        plt.plot(time, cutoff, 'k', linestyle = '--', linewidth = 2)
        
        path_to_file = os.path.dirname(filename)
        base_name    = os.path.basename(filename)
        plot_filename = path_to_file + '\\' + base_name + '.pdf'
        plt.savefig(plot_filename)
        plt.show()
    
    
if __name__ == '__main__':
    
    plot_file_name = 'fits\\fit010_test1_v1j3_ERF.fit_out'
    plot_file_name = 'fits\\fit001_v0j2_ERF.fit_out'
    plot_file_name = 'fits\\fit015_D2_v0j2_ERF.fit_out'
    plot_file_name = 'fits\\fit017_D2_v0j2_Calibration.fit_out'
    PlotFit(plot_file_name)