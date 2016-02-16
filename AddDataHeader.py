#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob
# import string
import time

data_path = 'd:\\users\\dja\\desktop\\permeation\\data\\2016.02.10\\Permeation Data\\Data'
data_path = 'd:\\users\\dja\\desktop\\permeation\\data\\2016.01.18\\Permeation Data\\Reference Line'
data_path = 'D:\\Users\\dja\\Git\\Desorption_TOF_Fit\\Data\\ReferenceLine'
os.chdir(data_path)
#for file in os.listdir(data_path):
for file in glob.glob('*.dat'):
    # fn, ft = file.split('.')[0], file.split('.')[1]
    fn = os.path.splitext(file)[0] 
    
    if fn[0:3].isdigit():
        orig_data_row = 28
        header_rows   = 12 + 3  # 3 for text block 'original header'
        data_row      = orig_data_row + header_rows + 4        
        temperature = 1063
        data_col = 2
    else:
        orig_data_row = 4
        header_rows   = 12 + 3  # 3 for text block 'original header'
        data_row      = orig_data_row + header_rows + 4        
        temperature = 1063
        data_col = 3
            
        
    
    if not os.path.exists(data_path + '\\' + fn + '.datv2') or True:
        print('processing file', file)
        v = fn[fn.find('v')+1]
        ij = fn.find('J')
        
        if file[ij + 2].isdigit():
            j = fn[ij + 1 : ij + 3]
        else:
            j = fn[ij + 1]

        if 'D2' in fn: 
            molecule = 'D2'
        elif 'H2' in fn:
            molecule = 'H2'
        else:
            print('opps no molecule name in filename ', fn)
        
        time_modified = time.ctime(os.path.getmtime(file))
        
        infile  = open(fn + '.dat' , 'r')        
        outfile = open(fn + '.datv2' , 'w')
        lines = infile.readlines()

        data_start = None
        data_lines = []        
        
        data_lines.append('#' + 47*'-' + '\n')
        data_lines.append('# Original Header\n')
        data_lines.append('#' + 47*'-' + '\n')
        for i in range(0,orig_data_row -1 ):
            data_lines.append('# ' + lines[i])

        data_lines.append('\n')
        data_lines.append('#' + 47*'-' + '\n')
        data_lines.append('# Begin Data\n')
        data_lines.append('#' + 47*'-' + '\n')
                                
        for i in range(orig_data_row -1, len(lines)):
            
            x1 = float(lines[i].split()[0])
            try: 
                y1 = float(lines[i].split()[1])
                #y2 = float(lines[i].split()[2])  
                if not data_start:
                    data_start = i + header_rows + 1
                    
                #data_lines.append('{0:6.2f} {1:20.5e} {2:20.5e}\n'.format(x1, y1, y2))    
                data_lines.append('{0:6.2f} {1:20.5e}\n'.format(x1, y1))
            except:
                y1 = lines[i].split()[1]
                #y2 = lines[i].split()[2]
                # following line writes the lines of data with no y values
                # data_lines.append('{0:6.2f} {1:20s} {2:20s}\n'.format(x1, y1, y2))
                data_lines.append('{0:6.2f} {1:20s}\n'.format(x1, y1))
            
        data_start = data_row  # since we don't write the lines with no y values
        outfile.write('#' + 47*'-' +                             '\n')
        outfile.write('# data_format   :  ' + '2.1'            + '\n')        
        outfile.write('# original_file :  ' + file             + '\n')
        outfile.write('# date          :  ' + time_modified    + '\n')   
        outfile.write('# molecule      :  ' + molecule         + '\n')
        outfile.write('# temperature   :  ' + str(temperature) + '\n')
        outfile.write('# v             :  ' + str(v)           + '\n')
        outfile.write('# J             :  ' + str(j)           + '\n')
        outfile.write('# data_col      :  ' + str(data_col)     + '\n')
        outfile.write('# data_row      :  ' + str(data_start)  + '\n')
        outfile.write('#' + 47*'-' +                             '\n')
        outfile.write('\n')
        
        for line in data_lines:
            outfile.write(line)
        
        infile.close()
        outfile.close()