#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Plot the data and result of fit """

def isNumber(s):
    try:
        float(s)
        return True
    except:
        return False
        
def plot_fit(filename):


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

    # get the filename for pdf file
    path_to_file = os.path.dirname(filename)
    fn_no_ext = os.path.basename(filename).rsplit('.',1)[0]  
    plot_filename = path_to_file + '\\' + fn_no_ext

    # create an object for multipage PDF
    pdf_multi = PdfPages(plot_filename + '.pdf')
    
    # read the fit_out file    
    with open(filename,'r') as fit_out:
        lines = fit_out.readlines()
     
    #--------------------------------------------------------------------------
    # get the number of plots to make 
    #--------------------------------------------------------------------------    
    for n_line, line in enumerate(lines):
        if line.startswith('# Number of plots'):
            n_plots = int(line.split(':')[1])
            break
    
    #--------------------------------------------------------------------------
    # list of variable to parse from .fit_out file for plot lables
    # note: t_min and t_max are not included becuase parse routine used here
    #       can't handle tuple time_range
    #--------------------------------------------------------------------------
    variable_list = ['avg_type_label', 
                     'baseline_label', 'baseline',
                     'date', 
                     'e0_label',
                     'ecutm_label', 'ecuts_label', 
                     'ffr_label', 
                     'function', 
                     'ion_tof_label', 
                     'molecule',
                     'n_points',
                     'signal_file',
                     'states', 
                     'tcutc', 'tcutw',
                     'temp_label', 
                     'Title',
                     'w_label',                                                              
                     'y_scale_label']
    
    #--------------------------------------------------------------------------
    # set all variables to None to suppress complaints about undefined variables
    #--------------------------------------------------------------------------
#==============================================================================
#     avg_type_label = None
#     date = None
#     e0_label = None
#     ecutm_label = ecuts = None
#     ion_tof_label = None
#     molecule = None
#     n_points = None
#     states = None
#     tcutc = tcutw = None
#     temp_label = None
#     Title = None
#     w_label = None
#     y_scale_label = None
#==============================================================================
    
    #--------------------------------------------------------------------------
    # make a plot for each data section
    #--------------------------------------------------------------------------
    for n_plot in range(n_plots):
        
        #----------------------------------------------------------------------
        # parse the fit_out file 
        # use exec to set values of variable
        #----------------------------------------------------------------------
        n_start = n_line + 1
        label = ''
        for n_line in range(n_start, len(lines)):
            line = lines[n_line]
            
            try:
                tok0 = line.split(' ',1)[1].split(':',1)[0].strip()
                tok1 = line.split(' ',1)[1].split(':',1)[1].strip()
                for var in variable_list: 
                    if tok0 == var:
                        if isNumber(tok1):
                            exec('global '+ tok0 + '; ' + tok0 + ' = ' + tok1)
                        else:
                            exec('global '+ tok0 + '; ' + tok0 + ' = ' + "tok1")
                            
            except:
                pass
                           
#==============================================================================
#             if line.startswith('# Title'):
#                 title = line.split(':', 1)[1]
#                 Title = title
#             
#             if line.startswith('# function'):
#                 function = line.split(':', 1)[1]
#                 
#             if line.startswith('# Npoints'):
#                 n_points = int(line.split(':')[1])
#            
#             if line.startswith('# Baseline'):
#                 baseline = float(line.split(':')[1])
#                 
#==============================================================================
            if line.startswith('# t_min'):
            # line will look like # t_min, t_max = :  (12, 50)\n
            # split on :, take second element, strip space, ( ) and \n) gives 12, 50
            # split on , gives list ['1','2']
               t_min, t_max = line.split(':')[1].strip(' ()\n').split(',')
               t_min = float(t_min)
               t_max = float(t_max)
                                               
            # break out of loop when we get to Begin Data
            if line.startswith('# Begin data'):
                pass
                break
        
        temp = temp_label.split(',')[0].strip() + ' K'        
        mass_factor ={'H2':1.0, 'HD':1.22474, 'D2':1.41421}            
        ion_tof = float(ion_tof_label.split(',')[0]) * mass_factor[molecule]
        print('ion_tof =', ion_tof)

        #------------------------------------------------------------------------------------------
        # get the time signal and fit arrays plot
        #------------------------------------------------------------------------------------------
        time   = []
        sig    = []
        fit    = []
        
        n_start = n_line + 1
        for n_line in range(n_start, n_start + n_points ):
            line = lines[n_line]
            if float(line.split()[0]) < ion_tof:
                continue
            time.append  (float(line.split()[0]))
            sig.append   (float(line.split()[1]))
            fit.append   (float(line.split()[2]))
        
        time = np.array(time)
        sig  = np.array(sig)            
        fit  = np.array(fit)

        #------------------------------------------------------------------------------------------
        # get the cutoff array over the range wher fit > fit.max()
        # use the same range of lines used for the time, signal, and fit arrays
        #     i.e. do not reset n_start
        # note: didn't use n_line as index variable so that n_line remains 
        #       the last data line
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
        
        cutoff_max = 0.8 * sig.max() + 0.2 * fit.max()
        cutoff = np.array(cutoff) * cutoff_max
        time_cutoff = np.array(time_cutoff)
        
        #------------------------------------------------------------------------------------------
        #  create the plot label based on plot type
        #------------------------------------------------------------------------------------------
        sig_file = os.path.basename(signal_file)
        
        if function == 'cal':
            label = 'Calibration, '
            label += avg_type_label                     + '\n'
            label += molecule + str(states) + ' ' + temp+ '\n'
            label += date                               + '\n'
            label += '\n'
            label += 'ecutm   ' + ecutm_label           + '\n'
            label += 'ecuts   ' + ecuts_label           + '\n'
            label += 'ffr     ' + ffr_label             + '\n'                
            label += 'ion_tof ' + ion_tof_label         + '\n'
                       
        if function == 'erf':
            label = 'erf, ' 
            label += avg_type_label + '\n'
            label += molecule + str(states) + ' ' + temp+ '\n'
            label += date                               + '\n'
            label += '\n'
            label += 'E0      ' + e0_label              + '\n'
            label += 'W       ' + w_label               + '\n'
            try:
                label += 'tcutc   ' + tcutc_label       + '\n'
                label += 'tcutw   ' + tcutw_label       + '\n'
            except:
                pass
            
            try:
                label += 'ecutm   ' + ecutm_label       + '\n'
                label += 'ecuts   ' + ecuts_label       + '\n'
            except:
                pass
                

        
        sig_max = sig.max()
        sig_min = sig.min()
    
        fig = plt.figure(figsize = (6,6), dpi = 200)
        fig = plt.figure()
        
        fig.suptitle(Title, fontsize=14, y=1.01)
        plt.title(sig_file, fontsize = 10 )
        
        ax = plt.subplot(111)
        ax.set_xlabel('TOF [' + unicodedata.lookup('micro sign') + 's]', fontsize=16)
        ax.set_ylabel('Signal', fontsize=16)
    
        ax.annotate(label, xy = [0.55, 0.95, ], xycoords = 'axes fraction', 
                    va = 'top', family='monospace', )
        ax.annotate('Cutoff', xy = (4., cutoff_max * 1.03), 
                    xycoords ='data', va='bottom')
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
        
#figtext(.5,.85,'Lorem ipsum dolor sit amet, consectetur adipiscing elit',fontsize=10,ha='center')
        
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
        
    plot_file_path = 'G:\\Spyder\\Au111-Fitting\\'
    plot_file_name = plot_file_path + 'Fit0005_H2_v0j1_erf.fit_out'
    
    plot_fit(plot_file_name)