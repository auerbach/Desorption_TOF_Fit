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
    multi_page_pdf = PdfPages(plot_filename + '.pdf')
    
    # read the fit_out file    
    with open(filename,'r') as fit_out:
        lines = fit_out.readlines()
     
    #--------------------------------------------------------------------------
    # get the number of plots to make
    # note: n_line will be at the Number of plots line which is the first line
    #       of the 'data for plots' section of the fit_out file
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
                     'tcutc_label', 'tcutw_label',
                     'temp_label', 
                     'Title',
                     'w_label',                                                              
                     'y_scale_label']
        
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
                if tok0 in variable_list: 
                    if isNumber(tok1):
                        exec('global '+ tok0 + '; ' + tok0 + ' = ' + tok1)
                    else:
                        exec('global '+ tok0 + '; ' + tok0 + ' = ' + "tok1")
            except:
                pass
            

            #----------------------------------------------------------------------
            #  process # t_min, t_max : (4.98, 25.00)           
            #  
            #----------------------------------------------------------------------                           
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
        cutoff = []
        
        n_start = n_line + 1
        for n_line in range(n_start, n_start + n_points ):
            line = lines[n_line]
            if float(line.split()[0]) < ion_tof:
                continue
            time.append   (float(line.split()[0]))
            sig.append    (float(line.split()[1]))
            fit.append    (float(line.split()[2]))
            cutoff.append (float(line.split()[3]))
        
        time   = np.array(time)
        sig    = np.array(sig)            
        fit    = np.array(fit)
        cutoff = np.array(cutoff)
        
        cutoff_max = 0.8 * sig.max() + 0.2 * fit.max()
        cutoff = np.array(cutoff) * cutoff_max
        
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
#==============================================================================
#             # -------------------------------------------------------------------------------------
#             # remove the , from E0 and W label and adjust alignment
#             # -------------------------------------------------------------------------------------
#             e0, e0_pm = e0_label.split(',')
#             w,  w_pm  = w_label.split(',')
#             new_e0_label = e0 + '  ' + e0_pm.strip()            
#             new_w_label  = w  + ' '  + w_pm.strip()
#==============================================================================

            label = 'Erf fit, ' 
            label += avg_type_label + '\n\n'
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

        # ---------------------------------------------------------------------
        # Set the size of the figure
        # ---------------------------------------------------------------------    
        fig = plt.figure(figsize = (8,6), dpi = 200)
        #fig = plt.figure()
        
        # ---------------------------------------------------------------------         
        # make extra room on top for title and suptitle
        # top sets the location of the top axis and suptitle
        # y = 1.01 in the axis title moves it up a little so _ looks OK
        # ---------------------------------------------------------------------
        fig.subplots_adjust(top=0.9)
        
        fig.suptitle(Title, fontsize=14)
        plt.title(sig_file, fontsize = 10, y=1.01 )
        
        # ---------------------------------------------------------------------
        # set up the axis labels
        # ---------------------------------------------------------------------
        ax = plt.subplot(111)
        ax.set_xlabel('TOF [' + unicodedata.lookup('micro sign') + 's]', fontsize=16)
        ax.set_ylabel('Signal', fontsize=16)

        # ---------------------------------------------------------------------
        # add the label to the plot
        # annotate the Cutoff line
        # add bars at t_min and t_max use for the fit
        # ---------------------------------------------------------------------    
        #ax.annotate(label, xy = [0.65, 0.95, ], xycoords = 'axes fraction', 
        #            va = 'top', family='monospace' )
        
        # these are matplotlib.patch.Patch properties
        props = dict(facecolor='lightyellow', alpha=.8)
        
        ax.text(0.6, 0.95, label,  family = 'monospace', fontsize=11, bbox = props,
                transform=ax.transAxes, va='top')
        
        ax.annotate('Cutoff', xy = (time[0], cutoff_max * 1.01), 
                    xycoords ='data', va='bottom')
        ax.annotate('', xycoords = 'data', xy = (t_min , sig_min), xytext =(t_min, sig_max*.3), 
                    arrowprops=dict(linewidth = 2.5, linestyle = '-', arrowstyle = '-')) 
        ax.annotate('', xycoords = 'data', xy = (t_max , sig_min), xytext =(t_max, sig_max*.3), 
                    arrowprops=dict(linewidth = 2.5, linestyle = '-', arrowstyle = '-'))   

        # ---------------------------------------------------------------------
        # set the min and max times for the plot
        # the criteria used here are possiblyt to simple but they seem to work
        # calibration plots will have max time of 60 microseconds
        # fit plots will have a max time of 30 microseconds
        # ---------------------------------------------------------------------        
        t_min_plot = 2        
        if t_max <= 20:
            t_max_plot = 30
        else: 
            t_max_plot = 60
            
        # ---------------------------------------------------------------------        
        # plot the signal and fit
        # ---------------------------------------------------------------------        
        plt.xlim((t_min_plot, t_max_plot))
        plt.ylim((sig_min, sig_max * 1.05))
        plt.plot(time, sig, 'b.')
        plt.plot(time, fit, 'r', linewidth = 2)
        
        # ---------------------------------------------------------------------        
        # plot the cutoff up to 0.6  of plot range to avoid obscuring the plot label
        # ---------------------------------------------------------------------        
        index = np.where(time > (t_max_plot - t_min_plot) * 0.6)[0][0]
        index = len(time)
        plt.plot(time[0:index], cutoff[0:index], 'g', linestyle = '--', linewidth = 2)
                
        # ---------------------------------------------------------------------        
        # save the plot as a multipage .pdf file and display the plot in console
        # ---------------------------------------------------------------------        
        multi_page_pdf.savefig()
        plt.show()
    
    multi_page_pdf.close()    
    
    
    
if __name__ == '__main__':
        
    plot_file_path = 'C:\\Users\\dja\\Desktop\\All DJA local\\Permeation Experiment\\Example Fits\\'
    plot_file_name = plot_file_path + 'Fit0036_D2_v1j2_erf.fit_out'
    plot_file_name = plot_file_path + 'Fit0001_D2_v0j2_ERF.fit_out'
    
    plot_fit(plot_file_name)