#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Module to parse the control file for FitDesorbtionTOF

"""

__author__ = 'Daniel Auerbach'  
__email__  = 'daniel@djauerbach.com'

from Parameters2 import Parameters2
import GlobalVariables as glbl

# Global Variables
backgroundFiles = []
backgroundFile = None

errors = []    

function = None
functions =  []

globalParms=[]
# parms = Parameters2()
parms = Parameters2()


signalFile = None
signalFiles = []

# names of parameters and constants
parmList = ['A', 'E0', 'W', 
            'TCutC', 'TCutW', 'Temp', 'FFR', 'IonTOF', 'Yscale', 'Baseline']

constList =['testConst', 'test1', 'test2', 'test3',
            'Label',
            'ThetaStep', 'ZDiffWall', 'RDiffWall', 
            'ZRef', 'ZAperture', 'RAperture', 
            'ZSource', 'RSource', 'ZLaser', 'LLaser', 'ZFinal', 'RFinal', 
            'NPointsDetector', 'NPointsSource', 'GridType',
            'DataLine', 'MassLine', 'TemperatureLine', 'VibStateLine', 'RotStateLine',
            'Tmin', 'Tmax']
            
            
#==============================================================================
# open file for reading; abort if file not found
#==============================================================================  
def openForReading(filename):
    try:
        file=open(filename, 'r')
        return file
    except FileNotFoundError:
        print('\n ***')
        print(' *** File ', filename, ' not found')
        print(' *** Ending program')
        print(' *** \n')
        raise SystemExit
        
    except:
        quit


#==============================================================================
# break input line into tokens
#=========================================================================
def getTokens(line):

    global tokens       
    # split the line on ',' 
    tokens = line.split(',')
    for j in range(len(tokens)):
    # replace nule tokens (i.e. '') with None
        tokens[j] = tokens[j].strip()
        if tokens[j] == '': tokens[j]=None
        # remove comments and strip whitespace
        try:
            ii = tokens[j].index('#')
            tokens[j] = tokens[j][0:ii].strip()
        except:
            pass
    return tokens


#==============================================================================
# check if cmd is an add Parameter item 
#==============================================================================
def isParameterCmd(tokens):
    return any(tokens[0] in s for s in parmList)


#==============================================================================
# check if cmd is an add Global Variable item 
#==============================================================================
def isGlobalConst(tokens):         
    matching = [tokens[0] == s for s in constList]   
    return any(matching)


#==============================================================================
# add Parameter to the Parameters class instance 
#==============================================================================
def addParameter(runNumber, tokens, line, lineNumber):
    
    global backgroundFiles, backgroundFile, errors, function, functions, \
           globalParms, parms, signalFile, signalFiles
    
    # check if number of tokesn is correct    
    #   there should be name, and 2, 3 or 4 parameters
    #   i.e. 3 to 5 tokens on the line
    if len(tokens) < 3 or len(tokens) > 5:
        print('\n *** wrong number of items on line #', lineNumber)
        print(' *** items found =', len(tokens), 
              '. Should be 3, 4, or 5')
        print(' *** line =', line)
        print(' *** items =', tokens)
        errors.append(('wrong number of items', lineNumber))    
        return
    
    # set defaults
    vary1 = False
    min1 = None
    max1 = None
    
    # set vary and global        
    if tokens[2] != '0': 
            vary1 = True
    if tokens[2] == '2':
        globalParms.append(tokens[0])
        
    # set bounds
    try: 
        min1=float(tokens[3])
    except: 
        pass
    try: 
        max1 = float(tokens[4])
    except: 
        pass

    parms.add(tokens[0] + '_' + str(runNumber),
              value = float(tokens[1]),
              vary = vary1,
              min = min1,
              max = max1
              )
              
              
#==============================================================================
# add const
#==============================================================================             
def addGlobalConst(tokens, line, lineNumber):
    if len(tokens) != 2:
        error = 'wrong number of items for const on line #'
        print('\n ***' + error, lineNumber)
        print(' *** items found =', len(tokens), 
              '. Should be 2')
        print(' *** line =', line)
        print(' *** items =', tokens)
        errors.append((error, lineNumber))    
        return
    exec('glbl.' + tokens[0] + '=' + tokens[1])
    


#==============================================================================
# add Parameter to the Parameters class instance 
#==============================================================================
def addGlobalParameters(runNumber):
    
    global backgroundFiles, backgroundFile, errors, function, functions, \
           globalParms, parms, signalFile, signalFiles
    
    for parmName in globalParms:
        oldparm = parms[parmName + '_1']
        parms.add(parmName + '_' + str(runNumber),
                  value = oldparm.value,
                  vary = oldparm.vary,
                  min = oldparm.min,
                  max = oldparm.max,
                  expr = parmName + '_1')
                            

def processEndSection(runNumber, lineNumber):
    global backgroundFiles, backgroundFile, errors, function, functions, \
           globalParms, parms, signalFile, signalFiles
    
    # check if signalFile specified
    if signalFile == None:
        print('\n *** Error: No signal file name for')
        print(' *** data section ending at line ', lineNumber)
        errors.append(('no signal file',lineNumber))
    
    # check function type is specified.  If not specified in this 
    # section of command file - use previously specified function type.
    if function == None:
        if len(functions) ==0:
            print('\n *** Error: No function specified for')
            print(' *** data section ending at line ', lineNumber)
            errors.append(('no function',lineNumber))
        else:
            print('\n *** Info: using previously specified function')
            print(' *** data section ending at line ', lineNumber)
            # runNumber starts at 1; list indexing starts at 0
            functions.append(functions[runNumber - 2])
        
    # if this is not the first section, add global parameters
    if runNumber > 1:
        addGlobalParameters(runNumber)
        
    # save filenames and reset to None for entry to next section
    signalFiles.append(signalFile)
    backgroundFiles.append(backgroundFile)
    signalFile = None
    backgroundFile = None
    function = None
    
    # increment run number for use in next section
    runNumber += 1
    return

                            
#==============================================================================
# parse command file
#==============================================================================
def parseCmdFile(filename):
    
    global backgroundFiles, backgroundFile, errors, function, functions, \
           globalParms, parms, signalFile, signalFiles

    # note initialization in global scope
    runNumber = 1
    
    # open and read command file
    cmdFile=openForReading(filename)
    lines = cmdFile.readlines()
    
    # pasre file, line by line
    for i in range(len(lines)):
        line=lines[i].strip()
        lineNumber = i+1     
        
        # skip blank lines and comments
        if line.startswith('#') or len(line) == 0:
            continue
        
        # break line into tokens        
        tokens = getTokens(line)
        
        # check if line specifies a parameter
        if isParameterCmd(tokens):
            addParameter(runNumber, tokens, line, lineNumber)
        
        # check if line specifies a constant
        elif isGlobalConst(tokens):
            addGlobalConst(tokens, line, lineNumber)
            
        # check if line specifies a function form
        elif tokens[0].upper() == 'FUNCTION':
            function = tokens[1]
            functions.append(function)
        
        # check if line specifies a signal or background file name
        elif tokens[0].upper() == 'SIGNAL':
            signalFile = tokens[1]
            continue
        elif tokens[0].upper() == 'BACKGROUND':
            backgroundFile = tokens[1]        
            continue
    
        # check if line indicates end of dataset section
        elif tokens[0].upper().startswith('END'):
            processEndSection(runNumber, lineNumber)
            
        else:
            print('\n *** unknown command on line', i+1)
            print(' *** line =', line)
            errors.append(('unknown command',i+1))
    
    return parms, functions, signalFiles, backgroundFiles, errors

#==============================================================================
# test parseCmdFile if executed as main
#==============================================================================
if __name__ == '__main__':         
    
    filename = 'Fits\\fit001.tof_in'
    const_filename = 'testSetVariables.dat'
    
    # parms, functions, signalFiles, backgroundFiles, errors = parseCmdFile(const_filename)
    # print('test1, test2, test3 =', glbl.test1, glbl.test2, glbl.test3)
    
    parms, functions, signalFiles, backgroundFiles, errors = parseCmdFile(filename)
    
    print('\nFuncions: ', functions)    
    print('\nsignalFiles\n', signalFiles)
    print('\nbackgroundFiles\n', backgroundFiles)
    print()
    parms.pretty_print()
    
    if len(errors) > 0:
        print('Errors\n', errors)
        raise SystemExit ('Error in command file' )
    
    print('Tmin, Tmax =', glbl.Tmin, glbl.Tmax)
    print(' all done')
