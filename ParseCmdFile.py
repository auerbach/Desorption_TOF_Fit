#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Module to parse the control file for FitDesorbtionTOF

"""

__author__ = 'Daniel Auerbach'  
__email__  = 'daniel@djauerbach.com'

from glob import glob
from Parameters2 import Parameters2
# import ast
import GlobalVariables as glbl


# names of parameters constants and lists
parmList = ['A', 'E0', 'W', 
            'TCutC', 'TCutW', 'ECutM', 'ECutS',
            'Temp', 'FFR', 'IonTOF', 'Yscale', 'Baseline']

constList =['ThetaStep', 'ZDiffWall', 'RDiffWall', 
            'ZRef', 'ZAperture', 'RAperture', 
            'ZSource', 'RSource', 'ZLaser', 'LLaser', 'ZFinal', 'RFinal',
            'NPointsDetector', 'NPointsSource', 'GridType']
            
listList = ['DataColLine', 'DataRowLine', 'MoleculeLine', 'TemperatureLine', 
            'VibStateLine', 'RotStateLine',
            'Tmin', 'Tmax', 'Function', 'AveragingType', 'cutoff_type']

            

#==============================================================================
# initialize
#==============================================================================
def initialize():
    
    # variables in module global scope

    global errors,  globalParms, parms
    
    glbl.backgroundFile = None
    glbl.backgroundFiles = []   
    errors = []      
    glbl.Function = None
    glbl.Functions =  []  
    globalParms=[]
    parms = Parameters2()
    glbl.signalFile = None
    glbl.signalFiles = []
    glbl.Tmax = None    
    glbl.Tmin = None

    glbl.DataLines = []
    glbl.MassLines = []
    glbl.TemperatureLines = []
    glbl.VibStateLines = []
    glbl.RotStateLines = []
    glbl.Tmaxs = []
    glbl.Tmins = []
    

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

    # split the line on ',' 
    tokens = line.split(',')
    for j in range(len(tokens)):
    # replace nule tokens (i.e. '') with None
        tokens[j] = tokens[j].strip()
        if tokens[j] == '': tokens[j] = None
        
        try:
            #test if token has an embedded '#'
            ii = tokens[j].index('#')
            # if so, remove rest of token and strip white space
            tokens[j] = tokens[j][0:ii].strip()
            # ignore rest of tokens
            for k in range(j + 1, len(tokens)):
                tokens.pop(-1)
            break 
        except:
            pass
    return tokens


#==============================================================================
# check if cmd is an add Parameter item 
#==============================================================================
def isParameter(tokens):
    return any(tokens[0] in s for s in parmList)


#==============================================================================
# check if cmd is an add Global Variable item 
#==============================================================================
def isGlobalConst(tokens):         
    matching = [tokens[0] == s for s in constList]   
    return any(matching)


#==============================================================================
# check if cmd is a List item 
#==============================================================================
def isListItem(tokens):
    matching = [tokens[0].upper() == s.upper() for s in listList]   
    #print('isListItem: tokens=',tokens, 'any=', any(matching)) #, 'matching=', matching)
    return any(matching)

#==============================================================================
# check if token is a number 
#==============================================================================
def is_number(s):
    try:
        float(s)
        return True
    except:
        return False

#==============================================================================
# add Parameter to the Parameters class instance 
#==============================================================================
def addParameter(runNumber, tokens, line, lineNumber):
    
    global errors,  globalParms, parms
    
    # check if number of tokesn is correct    
    #   there should be name, and 2, 3 or 4 parameters
    #   i.e. 3 to 5 tokens on the line
    if len(tokens) < 3 or len(tokens) > 5:
        print('\n *** wrong number of items on line ', lineNumber)
        print(' *** items found =', len(tokens), 
              '. Should be 3, 4, or 5')
        print(' *** line =', line)
        print(' *** items =', tokens)
        errors.append(('wrong number of items', 'line ' + str(lineNumber)))   
        return
    
    # set defaults
    vary1 = False
    glbl1 = False
    min1 = None
    max1 = None
    
    # set vary and global        
    if tokens[2] != '0': 
            vary1 = True
    if tokens[2] == '2':
        glbl1 = True
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
              glbl = glbl1,
              min = min1,
              max = max1
              )
    # print('addParameter runNumber, tokens =', runNumber, tokens)
              
              
#==============================================================================
# add const
#==============================================================================             
def addGlobalConst(tokens, line, lineNumber):
    if len(tokens) != 2:
        error = 'wrong number of items for const on line ' + str(lineNumber)
        print('\n ***' + error, lineNumber)
        print(' *** items found =', len(tokens), 
              '. Should be 2')
        print(' *** line =', line)
        print(' *** items =', tokens)
        errors.append(error)    
        return
    exec('glbl.' + tokens[0] + '=' + tokens[1])
    


#==============================================================================
# add list item
#==============================================================================             
def addListItem(tokens, line, lineNumber):
    
    global errors,  globalParms, parms
           
    if len(tokens) != 2:
        error = 'wrong number of arguments items on line ' + str(lineNumber)
        errors.append(error)
        
        print('\n ***' + error)
        print(' *** items found =', len(tokens), 
              '. Should be 2')
        print(' *** line =', line)
        print(' *** items =', tokens)
        

    # item is a string, need to add quotes to exec, e.g. ERF become 'ERF'
    if is_number(tokens[1]):
        exec('glbl.' + tokens[0] + '=' + tokens[1])
    else:
        exec('glbl.' + tokens[0] + '=' + "tokens[1]")
        
#==============================================================================
#     # append the item to list
#         cmd = 'glbl.' + tokens[0] + 's.append(glbl.' + tokens[0] +')'
#         exec(cmd)
#         # exec('glbl.' + tokens[0] + 's.append(glbl.' + tokens[0])
#         pass
#==============================================================================

#==============================================================================
# add Parameter to the Parameters class instance 
#==============================================================================
def addGlobalParameters(runNumber):
    
    global errors,  globalParms, parms
    
    for parmName in globalParms:
        oldparm = parms[parmName + '_1']
        parms.add(parmName + '_' + str(runNumber),
                  value = oldparm.value,
                  vary = oldparm.vary,
                  min = oldparm.min,
                  max = oldparm.max,
                  expr = parmName + '_1')

                            
#==============================================================================
# process end of section
#==============================================================================
def processEndSection(runNumber, lineNumber):
    global errors,  globalParms, parms
    
    optionalList = ['Tmin', 'Tmax']
   
    # if this is not the first section, add global parameters 
    #   and duplicate old list items if needed
    if runNumber > 1:
        addGlobalParameters(runNumber)
        
        # if item was not supplied, use value from previous run if available
        for item in listList:
            if not eval('glbl.' + item):
                try:
                    exec('glbl.' + item + ' = glbl.' + item + 's[' + str(runNumber -2) + ']')
                except:
                    pass
            
        
    
        
    # add item to items list and check if required items have been supplied        
    for item in listList:
        exec('glbl.' + item + 's.append(glbl.' + item + ')')
        if not eval('glbl.' + item) and not item in optionalList:
            print('\n *** Error: No ', item, ' for')
            print(' *** data section ending at line ', lineNumber)
            errors.append('no ' + str(item) + ' for section ending at line ' +str(lineNumber))
            
    # check if required parameters have been supplied
    required_parms = ['Yscale', 'Baseline', 'IonTOF', 'FFR', 'Temp']

    if glbl.cutoff_type.lower() == 'exp':
        required_parms.append('ECutM')
        required_parms.append('ECutS')
    
    elif glbl.cutoff_type.lower() == 'tanh':
        required_parms.append('TCutC')
        required_parms.append('TCutW')
        
    if glbl.Function.lower() == 'erf':
        required_parms.append('E0')
        required_parms.append('W')
        
    
    for p in required_parms:
        try:
            parms[p + '_' + str(runNumber) ].value
            continue
        except:
            if runNumber > 1:
                oldparm = parms[p + '_' + str(runNumber - 1)]
                parms.add(p + '_' + str(runNumber),
                  value = oldparm.value,
                  vary = oldparm.vary,
                  min = oldparm.min,
                  max = oldparm.max,
                  expr = oldparm.expr)
            else:
                print('\n *** Error: No parameter', p, ' for')
                print(' *** data section ending at line ', lineNumber)
                errors.append('no parameter' + p + ' for section ending at line ' +str(lineNumber))
                
        
    # save filenames and reset to None for entry to next section
    glbl.signalFiles.append(glbl.signalFile)
    glbl.backgroundFiles.append(glbl.backgroundFile)     
    glbl.signalFile = None
    glbl.backgroundFile = None
    
    glbl.Function = None
    glbl.AveragingType = None
    glbl.cutoff_type = None
    glbl.Tmin = None
    glbl.Tmax = None



          
                            
#==============================================================================
# parse command file
#==============================================================================
def parseCmdFile(filename):
    
    global errors,  globalParms, parms
    
    initialize()    
    runNumber = 1
        
    # open and read command file
    cmdFile=openForReading(filename)
    lines = cmdFile.readlines()
    cmdFile.close()
    
    # parse command file, line by line
    for i in range(len(lines)):
        line=lines[i].strip()
        lineNumber = i+1
        
        # skip blank lines and comments
        if line.startswith('#') or len(line) == 0:
            continue
        
        # break line into tokens        
        tokens = getTokens(line)
        
        # check if line specifies a parameter
        if isParameter(tokens):
            addParameter(runNumber, tokens, line, lineNumber)
        
        # check if line specifies a constant
        elif isGlobalConst(tokens):
            addGlobalConst(tokens, line, lineNumber)
            
        elif isListItem(tokens):
            addListItem(tokens, line, lineNumber)
            pass
        
        # check if line specifies a signal file name
        # note this could be consolidatae with listList items
        elif tokens[0].upper() == 'SIGNAL':
            #use following code to input signal file name without quotes
            filelist = []
            for file in glob(tokens[1]):
                filelist.append(file)
            if len(filelist) > 1:
                print('***** error nonunique match to file pattern \n') 
                print('*****', tokens[1])
                print('***** filelist =', filelist)
                errors.append(' nonunique match for file patern' + tokens[1] + 
                                ' on line ' + str(lineNumber))
            glbl.signalFile = file
            
            # use the following line to input signal file name with quotes
            #exec('glbl.signalFile =' + tokens[1])
            continue
        
        # check if line specifies a background file name
        # note this could be consolidatae with listList items
        elif tokens[0].upper() == 'BACKGROUND':
            # use following code to input background file name without quotes
            filelist = []
            for file in glob(tokens[1]):
                filelist.append(file)
            if len(filelist) > 1:
                print('***** error nonunique match to file pattern \n') 
                print('*****', tokens[1])
                print('***** filelist =', filelist)
                errors.append(' nonunique match for file patern' + tokens[1] + 
                                ' on line ' + str(lineNumber))
            glbl.backgroundFile = file
            # use the following line to file name with quotes
            # exec('glbl.backgroundFile =' + tokens[1])
            continue
    
        # check if line indicates end of dataset section
        elif tokens[0].upper().startswith('END'):
            processEndSection(runNumber, lineNumber)
            runNumber += 1
            
        else:
            print('\n *** unknown command on line', i+1)
            print(' *** command =', tokens[0])
            errors.append('unknown command on line ' + str(i+1) + ' command = ' + tokens[0])
    
    return parms, glbl.Functions, glbl.signalFiles, glbl.backgroundFiles, errors


#--------------------------------------------------------------------------------------------------
# print list items
#--------------------------------------------------------------------------------------------------
def print_list_items():
    for item in listList:
        print('{:16s} = {}'.format(item + 's', eval('glbl.' + item + 's')))
        
        
        
        
    
#==============================================================================
# test parseCmdFile if executed as main
#==============================================================================
if __name__ == '__main__':         
    
    
    global errors,  globalParms, parms
    
    print()
    print(' Test ParseCmdFile')
    print()
    print('=====================')
    print('    Initial State    ')
    print('=====================')
    print_list_items()
    print()    

#==============================================================================
#     # parse the constants file    
#     const_filename = 'testSetVariables.dat'   
#     parms, functions, signalFiles, backgroundFiles, errors = parseCmdFile(const_filename)
#==============================================================================

    filename = 'Fits\\fit021.tof_in'



    
    parms, functions, signalFiles, backgroundFiles, errors = parseCmdFile(filename)
    
    print()
    print('=========================')    
    print('       Final State       ')
    print('=========================')
    print()
    print('functions:            ', functions)
    print('glbl.Functions:       ', glbl.Functions)    
    print('signalFiles:          ', signalFiles)
    print('glbl.signalFiles:     ', glbl.signalFiles)
    print('backgroundFiles:      ', backgroundFiles)
    print('glbl.backgroundFiles: ', glbl.backgroundFiles)
    print('Tmins:                ', glbl.Tmins)
    print('Tmaxs:                ', glbl.Tmaxs)
    print()
    
    
    print_list_items()    
    print()        
    
    is_global_fit = any([parms[p].glbl for p in parms])
    print('{:16s} = {}'.format('Global Fit', is_global_fit))
    print()
    
    parms.pretty_print()
    
    if len(errors) > 0:
        print('*****')
        print('***** Errors in command file found')
        for error in errors:
            print('***** ', error)
        print('*****')
        
        raise SystemExit ('Error in command file')
    
    