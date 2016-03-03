#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob
# import string
import time
import re

# data_path = 'd:\\users\\dja\\desktop\\permeation\\data\\2016.02.10\\Permeation Data\\Data'
# data_path = 'd:\\users\\dja\\desktop\\permeation\\data\\2016.01.18\\Permeation Data\\Reference Line'
# data_path = 'D:\\Users\\dja\\Git\\Desorption_TOF_Fit\\Data\\2016.02.10'
# data_path = 'd:\\users\\dja\desktop\\Permeation\\All-DJA\\2016.02.26\\Dataset\\D2v0'
# data_path = 'd:\\users\\dja\desktop\\Permeation\\All-DJA\\2016.02.26\\Dataset\\H2v1'
# data_path = 'd:\\users\\dja\desktop\\Permeation\\All-DJA\\2016.02.26\\Dataset\\Knudsen'
# data_path = 'd:\\users\\dja\desktop\\Permeation\\All-DJA\\2016.02.26\\Dataset\\Reference Lines'
data_path = 'd:\\users\\dja\desktop\\Permeation\\All-DJA\\2016.03.02'

# find all files with extension .dat in data_path
os.chdir(data_path)
for file in glob.glob('*.dat'):
    default_temperature = None
    print('processing file: ',file)
    
    # get the filename without path or extension
    fn_with_path = os.path.splitext(file)[0]
    fn = os.path.basename(fn_with_path)

    
    time_modified = time.ctime(os.path.getmtime(file))
    
    # find the type of the file
    #   if filename starts with numbers it is a raw data file with header
    #   if not, it is a corrected file
    #   if filename contains 'Kn' it is a Knudsen data file
    if fn[0:3].isdigit():
        data_col = 2
        
#==============================================================================
#         if 'KN' in fn.upper():
#             file_type = 'Knudsen v1'
#         else:
#             file_type = 'Permeation v1'
#==============================================================================
    
    else:
        file_format = 'Processed_data v1'
        # processed data files do not have a temperature in the header: assume 1061
        temperature = 1061
        data_col = 3
        print('file format is processed data:  assuming temperature = 1061K')
            
#==============================================================================
#     # Get v and J from the file name
#     print('processing file', file)
#     v = fn[fn.find('v')+1]
#     ij = fn.find('J')
#     
#     if file[ij + 2].isdigit():
#         j = fn[ij + 1 : ij + 3]
#     else:
#         j = fn[ij + 1]
#         
#     # get molecule from the filename
#     if 'D2' in fn.upper(): 
#         molecule = 'D2'
#     elif 'H2' in fn.upper():
#         molecule = 'H2'
#     else:
#         print('opps no molecule name in filename ', fn)
#==============================================================================
    
    infile  = open(fn_with_path + '.dat'   , 'r')        
    outfile = open(fn_with_path + '.datv2' , 'w')
    lines = infile.readlines()

    for n_line, line in enumerate(lines):
        if line.upper() in ['H2,HD,D2']:
            molecule = line.upper()
        
        if line.lower().startswith('v'):
            v = int(line[1:])

        if line.lower().startswith('j'):
            j = int(line[1:])
        
        if line.startswith('Knnew'):
            t_knudsen = float(re.search(r'(\d+)', line).group(1)) + 273.
            pass

        if line.lower()[0:2] in ['ag', 'au', 'cu']:
            t_xtal = float(re.search(r'(\d+)', line).group(1))
            
        if line.lower()[0:2] in ['h2', 'hd', 'd2']:
            molecule = line.upper()[0:2]
            
        if line.lower().startswith('time'):
            data_row = n_line +2            # +2 becuase lines[0] is line 1 of the file
            header_rows = n_line + 1
        
    try:
        if t_xtal > 500:
            temperature = t_xtal
        else:
            temperature = t_knudsen
    except:
        if not default_temperature:
            default_temperature = input('No temperatures found. Please input temperature: ')
        temperature = float(default_temperature)
        
    data_start = None
    data_lines = []        
    
    data_lines.append('#' + 47*'-' + '\n')
    data_lines.append('# Original Header\n')
    data_lines.append('#' + 47*'-' + '\n')
    for i in range(data_row -1 ):
        data_lines.append('# ' + lines[i])
    data_lines.append('\n')
    data_lines.append('#' + 47*'-' + '\n')
    data_lines.append('# Begin Data\n')
    data_lines.append('#' + 47*'-' + '\n')
    
    header_length = len(data_lines)
    
                            
    for i in range(data_row -1, len(lines)):
        
        # find the start of valid data        
        #   some data sets have 'time - -' rather than numbers for several lines
        #   test if we can convert second column to a float; if not ignore line        
        try: 
            y1 = float(lines[i].split()[1])
            break
        except:
            pass
        
    data_start = i

    
    for i in range(data_start, len(lines)):
        x1 = float(lines[i].split()[0])
        y1 = float(lines[i].split()[1])
        
        if data_col == 2:
            data_lines.append('{0:6.2f} {1:20.5e} \n'.format(x1, y1))                
        
        if data_col == 3:
            y2 = float(lines[i].split()[2])        
            data_lines.append('{0:6.2f} {1:20.5e} {2:20.5e}\n'.format(x1, y1, y2))    
    
    data_row = header_length + 13
    # write the .datv2 header 
    outfile.write('#' + 47*'-' +                             '\n')
    outfile.write('# data_format   :  ' + '2.1'            + '\n')        
    outfile.write('# original_file :  ' + file             + '\n')
    outfile.write('# date          :  ' + time_modified    + '\n')   
    outfile.write('# molecule      :  ' + molecule         + '\n')
    outfile.write('# temperature   :  ' + str(temperature) + '\n')
    outfile.write('# v             :  ' + str(v)           + '\n')
    outfile.write('# J             :  ' + str(j)           + '\n')
    outfile.write('# data_col      :  ' + str(data_col)    + '\n')
    outfile.write('# data_row      :  ' + str(data_row)    + '\n')
    outfile.write('#' + 47*'-' +                             '\n')
    outfile.write('\n')

    # write the original header and valid data    
    for line in data_lines:
        outfile.write(line)
    
    infile.close()
    outfile.close()