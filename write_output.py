#!/usr/bin/env python
"""
write_output(glbl, data, pathToFits)

write the .fit_out file with results of fit and data for plotting

"""

import lmfit
import numpy as np
import scipy as sp
import unicodedata
import openpyxl as pyxl
import compute_tof
from Cutoff import cutoff_function


def write_fit_out(glbl, data, path, control_filename, fit_number, fit_result, PlotDataSets):
    state_string = 'v' + str(data.states[0][0]) + 'j' + str(data.states[0][1])
    result_filename = 'Fit' + fit_number + '_' + data.molecules[0] + \
                      '_' + state_string + '_' + glbl.functions[0] + '.fit_out'

    parms = fit_result.params

    with open(path + result_filename, 'w') as result_file:
        result_file.write('#\n')
        result_file.write('# Fit {} Results\n'.format(fit_number))
        result_file.write('#\n')
        result_file.write('#' + 60 * '-' + '\n')
        result_file.write('# Control file: {}\n'.format(control_filename))
        result_file.write('#' + 60 * '-' + '\n')

        # --------------------------------------------------------------------------
        # write the control file to the result file
        # --------------------------------------------------------------------------
        with open(control_filename, 'r') as cmd:
            cmd_lines = cmd.readlines()

        for line in cmd_lines:
            result_file.write('# ' + line)

        result_file.write('\n#' + 60 * '-' + '\n')
        result_file.write('# End Control File\n')
        result_file.write('#' + 60 * '-' + '\n')
        result_file.write('\n')

        # --------------------------------------------------------------------------
        # write the angular averaging parameters
        # --------------------------------------------------------------------------
        result_file.write('#' + 60 * '-' + '\n')
        result_file.write('# Angular Averaging Parameters\n')
        result_file.write('#' + 60 * '-' + '\n')
        result_file.write('# points_source    : ' + str(glbl.points_source) + '\n')
        result_file.write('# points_detector  : ' + str(glbl.points_detector) + '\n')
        result_file.write('# z_source         : ' + str(glbl.z_source) + '\n')
        result_file.write('# r_source         : ' + str(glbl.r_source) + '\n')
        result_file.write('# z_aperture       : ' + str(glbl.z_aperture) + '\n')
        result_file.write('# r_aperture       : ' + str(glbl.r_aperture) + '\n')
        result_file.write('# z_final          : ' + str(glbl.z_final) + '\n')
        result_file.write('# r_final          : ' + str(glbl.r_final) + '\n')
        # result_file.write('# it mse\n')

        for i in range(0, len(glbl.angles_list), 10):
            result_file.write('# angles_list : ' + str(glbl.angles_list[i: i + 10]) + '\n')
        result_file.write('#' + 60 * '-' + '\n')
        result_file.write('# End Angular Averaging Parameter\n')
        result_file.write('#' + 60 * '-' + '\n')
        result_file.write('\n')

        # --------------------------------------------------------------------------
        # write the fit report to the result file
        # --------------------------------------------------------------------------
        result_file.write('#' + 60 * '-' + '\n')
        result_file.write('# Fit Report\n')
        result_file.write('#' + 60 * '-' + '\n')

        report = lmfit.fit_report(fit_result).split('\n')
        for line in report:
            result_file.write('# ' + line + '\n')
        result_file.write('#' + 60 * '-' + '\n')
        result_file.write('# End Fit Report\n')
        result_file.write('#' + 60 * '-' + '\n')
        result_file.write('\n')

        # --------------------------------------------------------------------------
        # the number of datasets to plot
        # --------------------------------------------------------------------------
        result_file.write('#' + 60 * '-' + '\n')
        result_file.write('# Labels and data for plots\n')
        result_file.write('#' + 60 * '-' + '\n')
        result_file.write('# Number of plots : ' + str(len(PlotDataSets)) + '\n')
        result_file.write('#' + 60 * '-' + '\n')
        result_file.write('\n')

        # --------------------------------------------------------------------------
        # for each data set, write plot title, labels, and data
        # --------------------------------------------------------------------------
        for i in range(len(PlotDataSets)):
            n_dataset = i + 1

            Nmin = data.fit_index_ranges[i][0]
            Nmax = data.fit_index_ranges[i][1]
            Tmin = data.datasets[i][0][Nmin]
            Tmax = data.datasets[i][0][Nmax]
            ProbCurveType = glbl.functions[i]
            averaging_type = glbl.averaging_types[i]
            # DataFile = signal_filenames[i]
            cutoff_type = glbl.cutoff_types[i]

            result_file.write('#' + 60 * '-' + '\n')
            result_file.write('# Plot ' + str(n_dataset) + '\n')
            result_file.write('#' + 60 * '-' + '\n')

            # --------------------------------------------------------------------------
            # write the plot title
            # --------------------------------------------------------------------------
            result_file.write('# Title    : Fit {:04d}_{}'.format(int(fit_number), n_dataset) + '\n')

            # --------------------------------------------------------------------------
            # write the plot label
            # --------------------------------------------------------------------------
            # result_file.write('# Label    : ' + data.original_signal_names[i] + '\n')
            result_file.write('# Label    : \n')
            # result_file.write('# Label    : ' + '----------------------\n')

            if averaging_type == 'none':
                avg_type_label = 'No Angular Averaging'
            elif averaging_type == 'point_detector':
                avg_type_label = 'Point Detector'
            elif averaging_type == 'line_detector':
                avg_type_label = 'Line Detector'

            # result_file.write('# Label    : Averaging: ' + avg_type + '\n')
            result_file.write('# Label    : function:  ' + ProbCurveType + '\n')
            result_file.write('# Label    : ' + avg_type_label + '\n')
            result_file.write('# Label    : \n')
            final_params = fit_result.params

            parm_for_plot_label = ['e0', 'w', 'tcutc', 'tcutw', 'ecutm', 'ecuts', 'ffr']
            # if glbl.functions[i].lower().startswith('cal'):
            #     parm_for_plot_label.append('FFR')
            for p in parm_for_plot_label:
                p_ = p + '_' + str(i + 1)
                try:
                    result_file.write('# Label    : {:6s}{:>7.3f}'.format(p, final_params[p_].value))

                    if parms[p + '_1'].glbl and i > 0:
                        result_file.write(' (global)\n')
                    if not final_params[p_].vary:
                        result_file.write(' (fixed)\n')
                    else:
                        stderr = '{:6.3f}'.format(final_params[p_].stderr)
                        result_file.write(' (' + unicodedata.lookup('plus-minus sign')
                                          + str(stderr) + ')\n')
                except:
                    pass

            # --------------------------------------------------------------------------
            # Get the point where time > 2E-6 sec and ignore earlier points to avoid noise
            # --------------------------------------------------------------------------
            i_start = np.where(PlotDataSets[i][0] > 2.0E-6)[0][0]

            # --------------------------------------------------------------------------
            # Get the plot data and convert time to microseconds
            # --------------------------------------------------------------------------
            Time = PlotDataSets[i][0][i_start:]
            Signal = PlotDataSets[i][1][i_start:]
            Fit = compute_tof.TOF(Time, n_dataset, fit_result.params, data, glbl, averaging_type, glbl.angles_list,
                                  ProbCurveType, cutoff_type, data.mass_molecules)

            mass_factor = np.sqrt(data.mass_molecules[i] / glbl.massH2)
            ion_tof = fit_result.params['ion_tof_%i' % n_dataset].value * mass_factor

            # Compute cutoff for plotting.
            #   Create a new time array which is corrected for ion time of flight
            #   Replace time=0 by time = 1e-8 to remove singularity in velocity and
            #   energy if t=0
            Time2 = Time - ion_tof * 1E-6
            Time2 = np.where(Time2 != 0, Time2, np.repeat(0.01E-6, len(Time)))
            Cutoff = cutoff_function(fit_result.params, data, glbl,
                                            n_dataset, Time2, cutoff_type)
            # Convert time to microseconds for plotting
            Time = Time * 1E6  # convert to microseconds for plotting

            # --------------------------------------------------------------------------
            # write the number of data lines and range of lines fit
            # --------------------------------------------------------------------------
            result_file.write('# Npoints    : ' + str(len(Time)) + '\n')
            result_file.write('# n_min, n_max : ' + str(Nmin) + ',' + str(Nmax) + '\n')
            result_file.write('# t_min, t_max : ' + str(Tmin * 1E6) + ',' + str(Tmax * 1E6) + '\n')
            result_file.write('# Baseline   : ' +
                              str(fit_result.params['baseline_' + str(n_dataset)].value) + '\n')

            result_file.write('#' + 68 * '-' + '\n')

            # --------------------------------------------------------------------------
            # write the plot data
            # --------------------------------------------------------------------------
            result_file.write('# time            Signal                Fit                 Cutoff\n')
            result_file.write('#' + 68 * '-' + '\n')
            result_file.write('# Begin data\n')

            for j in range(len(Time)):
                result_file.write('{:6.2f} {:20.5e} {:20.5e} {:20.5e}\n'
                                  .format(Time[j], Signal[j], Fit[j], Cutoff[j]))

            result_file.write('# End data\n')
            result_file.write('#' + 68 * '-' + '\n')
            result_file.write('# End Plot ' + str(n_dataset) + '\n')

    # end comment out for profiling

    # print("--- %s seconds ---" % (time.time() - start_time))

    return result_filename


# ------------------------------------------------------------------------------
# Enter fit results in results.xlsx
# ------------------------------------------------------------------------------
def write_results_excel_file(glbl, data, path, fit_number, fit_result):
    xls_filename = path + 'fit_results.xlsx'
    wb = pyxl.load_workbook(xls_filename)
    ws = wb.active

    i_found = None
    next_row = None
    fit_name = 'Fit' + fit_number

    for i in range(1, ws.max_row + 1):
        test = ws.cell(row=i, column=1).value
        try:
            if ws.cell(row=i, column=1).value.startswith(fit_name):
                i_found = i
        except:
            pass

    if i_found:
        print('An entry for fit', fit_name, 'already exists')

        while True:
            ans = input('Overwrite (O)  Append (A)  or Skip (S): ').upper()
            if ans.startswith('O'):
                next_row = i_found
                break
            elif ans.startswith('A'):
                next_row = ws.max_row + 1
                break
            elif ans.startswith('S'):
                next_row = None
                break
            else:
                print('Please enter "O", "A", or "S" : ')

    else:
        next_row = ws.max_row + 2

    if (next_row):
        for n in range(len(glbl.signalFiles)):
            ws.cell(row=next_row, column=1).value = fit_name + '_' + str(n + 1)
            ws.cell(row=next_row, column=2).value = data.molecules[n]  # molecule
            ws.cell(row=next_row, column=3).value = data.states[n][0]
            ws.cell(row=next_row, column=4).value = data.states[n][1]
            ws.cell(row=next_row, column=5).value = glbl.functions[n]
            ws.cell(row=next_row, column=6).value = glbl.averaging_types[n]
            if glbl.functions[0].lower() == 'erf':
                ws.cell(row=next_row, column=7).value = fit_result.params['e0_' + str(n + 1)].value
                ws.cell(row=next_row, column=8).value = fit_result.params['e0_' + str(n + 1)].stderr
                ws.cell(row=next_row, column=9).value = fit_result.params['w_' + str(n + 1)].value
                ws.cell(row=next_row, column=10).value = fit_result.params['w_' + str(n + 1)].stderr
            ws.cell(row=next_row, column=11).value = glbl.Tmins[0]
            ws.cell(row=next_row, column=12).value = glbl.Tmaxs[0]
            ws.cell(row=next_row, column=13).value = fit_result.params['ffr_' + str(n + 1)].value
            ws.cell(row=next_row, column=14).value = fit_result.params['ecutm_' + str(n + 1)].value
            ws.cell(row=next_row, column=15).value = fit_result.params['ecuts_' + str(n + 1)].value
            if n == 0:
                ws.cell(row=next_row, column=16).value = glbl.comment_xlsx
            next_row += 1
    wb.save(xls_filename)
