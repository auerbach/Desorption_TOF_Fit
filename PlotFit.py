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
    from matplotlib.backends.backend_pdf import PdfPages
    import unicodedata

    # get the filename for pdf file; PdfPages appends .pdf
    path_to_file = os.path.dirname(filename)
    fn_no_ext = os.path.basename(filename).rsplit('.',1)[0]  
    plot_filename = path_to_file + '\\' + fn_no_ext

    # create an object for multipage PDF
    pdf_multi = PdfPages(plot_filename + '.pdf')
    
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
        n_start = n_line + 1
        label = ''
        for n_line in range(n_start, len(lines)):
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
               
            if line.startswith('# Baseline'):
                baseline = float(line.split(':')[1])
                    
            # break out of loop when we get to Begin Data
            if line.startswith('# Begin data'):
                break
        #------------------------------------------------------------------------------------------
        # get the number of data points and minimum and maximum times used for the fit       
        #------------------------------------------------------------------------------------------
            
        #------------------------------------------------------------------------------------------
        # get the time signal and fit arrays plot
        #------------------------------------------------------------------------------------------
        time   = []
        sig    = []
        fit    = []
        
        n_start = n_line + 1
        for n_line in range(n_start, n_start + n_points ):
            line = lines[n_line]
            time.append  (float(line.split()[0]))
            sig.append   (float(line.split()[1]))
            fit.append   (float(line.split()[2]))
        
        time = np.array(time)
        sig  = np.array(sig)            
        fit  = np.array(fit)

        #------------------------------------------------------------------------------------------
        # get the cutoff array over the range wher fit > fit.max()
        # use the same range of lines used for the time, signal, and fit arrays
        #  i.e. do not reset n_start
        # do not use n_line as index variable so that n_line remains the last data line
        #------------------------------------------------------------------------------------------
        cutoff = []
        time_cutoff = []
        trigger = False
        for nn in range(n_start, n_start + n_points ):
            line = lines[nn]
            if not trigger and fit[nn - n_start] > 0.5 * fit.max():
                trigger = True
            if trigger and fit[nn - n_start] - baseline < 0.02 * (fit.max() - baseline):
                break
            time_cutoff.append(float(line.split()[0]))
            cutoff.append(float(line.split()[3]))
        
        cutoff = np.array(cutoff) * fit.max()
        time_cutoff = np.array(time_cutoff)
        
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
        ax.annotate('Cutoff', xy = (4., fit.max()*1.03), xycoords ='data', va='bottom')
#==============================================================================
#         ax.annotate('$t_{min}$', xycoords = 'data', xy = (t_min-.1, sig_max * .4), 
#                     ha = 'right', va='center', fontsize=14)        
#         ax.annotate('$t_{max}$', xycoords = 'data', xy = (t_max+.1, sig_max * .4), 
#                     ha = 'left',  va='center', fontsize=14)                
#==============================================================================
        ax.annotate('', xycoords = 'data', xy = (t_min , sig_min), xytext =(t_min, sig_max*.3), 
                    arrowprops=dict(linewidth = 2.5, linestyle = '-', arrowstyle = '-')) 
        ax.annotate('', xycoords = 'data', xy = (t_max , sig_min), xytext =(t_max, sig_max*.3), 
                    arrowprops=dict(linewidth = 2.5, linestyle = '-', arrowstyle = '-'))   

        if t_max <= 20:
            t_max_plot = 30
        else: 
            t_max_plot = 60
            
        plt.xlim((2, t_max_plot))
        plt.ylim((sig_min, sig_max * 1.05))
        plt.plot(time, sig, 'b.')
        plt.plot(time, fit, 'r', linewidth = 2)
        plt.plot(time_cutoff, cutoff, 'g', linestyle = '--', linewidth = 2)
        
#==============================================================================
#         path_to_file = os.path.dirname(filename)
#         base_name    = os.path.basename(filename)
#         plot_filename = path_to_file + '\\' + base_name + '.pdf'
#         plt.savefig(plot_filename)
#==============================================================================
        pdf_multi.savefig()
        
        plt.show()
    
    pdf_multi.close()    
    
    
    
if __name__ == '__main__':
    
    plot_file_name = 'fits\\Fit010_test1_v1j3_ERF.fit_out'
    plot_file_name = 'fits\\Fit001_v0j2_ERF.fit_out'
    plot_file_name = 'fits\\Fit015_D2_v0j2_ERF.fit_out'
    plot_file_name = 'fits\\Fit017_D2_v0j2_CAL.fit_out'
    plot_file_name = 'fits\\Fit018_D2_v0j2_CAL.fit_out'
    plot_file_name = 'fits\\Fit019_D2_v0j2_CAL.fit_out'
    plot_file_name = 'fits\\Fit020_D2_v0j2_CAL.fit_out'
    plot_file_name = 'fits\\Fit021_D2_v0j2_CAL.fit_out'
    plot_file_name = 'fits\\Fit022_D2_v1j2_ERF.fit_out'
    plot_file_name = 'fits\\Fit023_D2_v1j2_ERF.fit_out'
    plot_file_name = 'fits\\Fit024_H2_v1j0_ERF.fit_out'
    plot_file_name = 'fits\\Fit025_H2_v1j0_ERF.fit_out'    
    #plot_file_name = 'fits\\Fit022_d2_v1j2_ERF.fit_out'
    # plot_file_name = 'fits\\Fit023_d2_v1j2_ERF.fit_out'
    
    PlotFit(plot_file_name)